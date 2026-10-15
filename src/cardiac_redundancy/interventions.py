"""Scheduled native-arrest orchestration, with one-sided event observations.

The existing integrator solves each constant-activation interval. This layer
ends that interval exactly at arrest and carries every continuous coordinate
into the passive interval. It does not erase inertial flow or reset pressure.
Pressure is recomputed from the new constitutive law and can jump.
"""
from dataclasses import dataclass
import math
import numpy as np

from .integration import integrate_source, Trajectory


@dataclass
class ArrestRun:
    # Keep segments separate so downstream pressure/oxygen processing cannot
    # accidentally interpolate across a discontinuity in activation.
    segments: list[tuple[bool, Trajectory]]
    interventions: list[dict]
    termination: dict | None


def observe(model, time, state, active):
    """Record the same state under one side's pressure and flow laws."""
    p, meta = model.pressures(time, state, active=active)
    # Valve flow is already a continuous state. The regional flow accessor
    # returns it directly; modes affect its derivative, not this observation.
    q = model.flows(time, state, [True]*4, active=active, pressures=p)
    return dict(state=np.asarray(state).tolist(), pressures_mmHg=p,
                flows_mL_s=q, activation=meta['activation'])


def integrate_arrest(model, *, initial_state, start, arrest_time, end, **solver):
    """Run F0: active until arrest_time, passive afterwards, without support.

    All times are absolute seconds. Arrest may coincide with start; a positive
    post-arrest interval is required. A prior domain exit suppresses arrest.
    The returned event contains left/right observations at identical states.
    """
    if not all(math.isfinite(v) for v in (start, arrest_time, end)) or not start <= arrest_time < end:
        raise ValueError('Require finite start <= arrest_time < end')
    y = np.array(initial_state, dtype=float, copy=True)
    segments = []
    if arrest_time > start:
        pre = integrate_source(model, arrest_time-start, initial_state=y,
                               start=start, active=True, **solver)
        segments.append((True, pre))
        if pre.termination is not None:
            return ArrestRun(segments, [], pre.termination)
        y = pre.states[:, -1].copy()
    # Respect zero-time domain violations before executing an intervention.
    if hasattr(model, 'domain_margins'):
        margins = model.domain_margins(y)
        if min(margins.values()) <= 0:
            stop = {'kind':'model_domain_limit', 'node':min(margins,key=margins.get),
                    'time_s':arrest_time}
            return ArrestRun(segments, [], stop)
    event = {'kind':'native_arrest', 'time_s':arrest_time,
             'before':observe(model, arrest_time, y, True),
             'after':observe(model, arrest_time, y, False)}
    post = integrate_source(model, end-arrest_time, initial_state=y,
                            start=arrest_time, active=False, **solver)
    segments.append((False, post))
    return ArrestRun(segments, [event], post.termination)
