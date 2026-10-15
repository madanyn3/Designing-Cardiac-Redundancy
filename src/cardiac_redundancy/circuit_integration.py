"""Generalized valve-count integrator for M3 circuit fixtures.

Derived from the verified four-valve implementation; historical implementation
is retained unchanged. Valve isolation applies only at segment boundaries.
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .mechanics import ModelDomainError
from .source import EDGES, VALVES


@dataclass
class Trajectory:
    """Sampled result, not the solver's internal step history.

    time: length-N sample times [s].
    states: shape (number of volumes + 4,N), physical/source mL then mL/s.
    events: ordered valve action records with localized times [s].
    nfev: RHS evaluation count reported by solve_ivp, summed over segments;
        not a count of every pressure/root solve or numerical-Jacobian call.
    segments: count of positive-duration dense solution pieces.

    Dense interpolants are used during assembly but are not retained here.
    """
    time: np.ndarray
    states: np.ndarray
    events: list
    nfev: int
    segments: int
    termination: dict | None = None


def integrate_circuit(model, duration, *, initial_state=None, start=0., active=True,
                     method="Radau", rtol=None, max_step=None, sample_step=.001):
    """Fixed native-activation run. Caller must split any scheduled intervention.

    Only a localized closing-valve flow is reset to zero. Volumes are never reset.
    Returned samples are interpolated from dense solver segments, not integration
    steps; output spacing is independent of the maximum solver step.

    Parameters:
        model: SourceModel with a compatible frozen baseline configuration.
        duration, start: interval length and absolute start time [s]. Driver
            phase follows absolute time; restarting at zero also resets phase.
        initial_state: optional model-sized state; copied before use.
        active: use the native periodic driver if True, passive walls if False.
        method: solve_ivp method; default Radau is an implicit stiff-capable method.
        rtol, max_step: optional relative tolerance / maximum step [s] overrides;
            absolute tolerances always come from the model configuration.
        sample_step: upper bound on uniform output spacing [s]. linspace includes
            both endpoints; actual spacing can be smaller for a noninteger ratio.

    Returns a Trajectory. Invalid states, failed integration or switching stalls
    raise exceptions rather than silently returning a successful shortened run.
    A regional domain crossing returns the trajectory up to the boundary with
    an explicit termination record; callers must check it before measuring.
    At an initial zero margin the regional run stops conservatively at time zero.
    Initial valve modes are inferred, not accepted as restart inputs. Dedicated
    simultaneous/near-coincident event and restart validation remains pending.
    """
    if not np.isfinite(duration) or duration <= 0 or not np.isfinite(sample_step) or sample_step <= 0:
        raise ValueError("Positive finite duration and sample step required")
    settings = model.settings
    valve_edges = model.valve_edges
    valve_names = model.valve_names
    enabled = np.asarray(model.valve_enabled, dtype=bool)
    nq = len(valve_names)
    nv = getattr(model, "volume_count", 6)
    termination = None
    y = model.initial_state.copy() if initial_state is None else np.array(initial_state, dtype=float, copy=True)
    if np.any(y[nv:] < 0):
        raise ModelDomainError("Initial healthy-valve flow cannot be negative")
    if hasattr(model, "domain_margins"):
        margins=model.domain_margins(y)
        if min(margins.values()) <= 0:
            node=min(margins,key=margins.get)
            return Trajectory(np.array([start]),y[:,None],[],0,0,
                              {"kind":"model_domain_limit","node":node,"time_s":start})
    p, _ = model.pressures(start, y, active=active)
    # Existing forward flow keeps a valve open even under adverse pressure:
    # inertance must decelerate it first. A zero-flow valve opens for P_up>P_down.
    modes = np.array([y[nv+i] > 0 or p[a] > p[b] for i,(a,b) in enumerate(valve_edges)])
    modes &= enabled
    if np.any(y[nv:][~enabled] != 0):
        raise ValueError("Cannot isolate a valve with nonzero inertial flow")
    t, end = float(start), float(start + duration)
    dense_segments, events = [], []
    nfev = 0
    no_progress = 0
    while t < end:
        # Freeze mode semantics for this ODE call; switch only after it returns.
        mode_snapshot = modes.copy()
        event_functions = []
        for i, (a,b) in enumerate(valve_edges):
            if not enabled[i]:
                def event(time, state): return 1.0
                event.direction = 0
            elif modes[i]:
                # Close only at a downward flow crossing, not merely because
                # the pressure gradient has become adverse.
                def event(time, state, index=i): return state[nv+index]
                event.direction = -1
            else:
                # A closed valve waits for an upward pressure-gradient crossing.
                # Default arguments bind this loop's endpoints/index immediately;
                # otherwise Python closures would all refer to the final valve.
                def event(time, state, donor=a, recipient=b):
                    pressures, _ = model.pressures(time, state, active=active)
                    return pressures[donor] - pressures[recipient]
                event.direction = 1
            event.terminal = True  # Stop at the first event and rebuild modes.
            event_functions.append(event)
        domain_names = []
        if hasattr(model, "domain_margins"):
            domain_names=list(model.domain_margins(y))
            for node in domain_names:
                def boundary(time,state,node=node): return model.domain_margins(state)[node]
                boundary.direction=-1
                boundary.terminal=True
                event_functions.append(boundary)
        # Per-coordinate absolute tolerances account for volume/flow units.
        # dense_output supplies interpolants for later sampling; sample_step
        # does not constrain solver steps or improve event localization by itself.
        sol = solve_ivp(lambda time,state: model.rhs(time,state,mode_snapshot,active=active),
                        (t,end), y, method=method, dense_output=True, events=event_functions,
                        rtol=rtol if rtol is not None else settings["relative_tolerance"],
                        atol=[settings["absolute_volume_tolerance_mL"]]*nv + [settings["absolute_flow_tolerance_mL_s"]]*nq,
                        max_step=max_step if max_step is not None else settings["maximum_step_s"])
        nfev += sol.nfev
        if not sol.success:
            raise RuntimeError(f"Integrator failed at {sol.t[-1]}: {sol.message}")
        next_t = float(sol.t[-1])
        y = sol.y[:,-1].copy()
        if next_t > t:
            dense_segments.append((t,next_t,sol.sol))
            no_progress = 0
        else:
            no_progress += 1
            # Fail explicitly if a root repeatedly retriggers at the same time.
            # This guard diagnoses stalls; it is not a proof against chattering.
            if no_progress > 8:
                raise RuntimeError(f"Valve events stalled at {t}; modes={modes}")
        t = next_t
        if sol.status == 0:  # Reached the requested end, rather than an event.
            break
        # Process events returned by this terminal solve. Do not assume that
        # every mathematically coincident event is necessarily returned together.
        triggered = [i for i, times in enumerate(sol.t_events) if len(times)]
        if not triggered:
            raise RuntimeError("Terminal stop without a valve event")
        domain_hits=[i for i in triggered if i>=nq]
        if domain_hits:
            termination={"kind":"model_domain_limit","node":domain_names[domain_hits[0]-nq],"time_s":t}
            end=t
            break
        for i in triggered:
            closing = bool(modes[i])
            if closing:
                # Hard-coded localization guard [mL/s], matching the current
                # nominal flow atol. It is not a pressure/injury threshold.
                if abs(y[nv+i]) > 1e-6:
                    raise RuntimeError("Localized valve closure exceeds flow tolerance")
                y[nv+i] = 0.0  # Enforce the located Q=0 event, not volume clipping.
            modes[i] = not modes[i]
            events.append({"time_s": t, "valve": valve_names[i], "action": "close" if closing else "open"})
        if len(events) > 10000:
            raise RuntimeError("Excess valve switching; possible chattering")
    # Build one uniform output grid from all positive-duration segments. At a
    # shared endpoint a later segment overwrites the preceding segment's sample.
    count = int(np.ceil((end-start)/sample_step)) + 1
    times = np.linspace(start,end,count)
    states = np.empty((len(y),count))
    assigned = np.zeros(count,dtype=bool)
    for left,right,dense in dense_segments:
        mask = (times >= left) & (times <= right)
        # A short event segment can contain no requested output samples. SciPy's
        # OdeSolution cannot evaluate an empty array; skip sampling, not physics.
        # Its integration and valve event have already been retained above.
        if not np.any(mask):
            continue
        states[:,mask] = dense(times[mask])
        assigned |= mask
    if not np.all(assigned):
        raise RuntimeError("Dense output has unassigned samples")
    states[:,-1] = y  # Preserve the final post-event state if end is an event.
    return Trajectory(times, states, events, nfev, len(dense_segments), termination)
