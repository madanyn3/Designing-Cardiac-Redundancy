"""Prescribed, finite metabolic transitions; no suppression mechanism implied.

All demand vectors use mL O2 STPD/min. Times and durations are seconds.
Hydraulic flow and valve state are deliberately absent from this interface.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class DemandCommand:
    time_s: float
    target_mL_min: object
    duration_s: object


class FiniteDemandSchedule:
    """Continuous linear ramps with independent per-organ durations.

    A new command interrupts any active ramp, starting from its current value.
    Returning demand to baseline is not synonymous with functional recovery.
    Positive durations are mandatory, including for a return command.
    """
    def __init__(self, baseline_mL_min, commands=()):
        self.baseline = np.array(baseline_mL_min, dtype=float, copy=True)
        if self.baseline.ndim != 1 or not self.baseline.size or np.any(~np.isfinite(self.baseline)) or np.any(self.baseline < 0):
            raise ValueError('Finite nonnegative baseline vector required')
        self.baseline.setflags(write=False)
        self.events = []
        previous = -np.inf
        for command in commands:
            time = float(command.time_s)
            target = np.array(command.target_mL_min, dtype=float, copy=True)
            duration = np.broadcast_to(np.asarray(command.duration_s, dtype=float), self.baseline.shape).copy()
            if not np.isfinite(time) or time < 0 or time <= previous:
                raise ValueError('Command times must be finite, nonnegative and strictly increasing')
            if target.shape != self.baseline.shape or np.any(~np.isfinite(target)) or np.any(target < 0):
                raise ValueError('Invalid target demand vector')
            if np.any(~np.isfinite(duration)) or np.any(duration <= 0):
                raise ValueError('Finite positive ramp durations required')
            start = self.at(time)
            for array in (start, target, duration): array.setflags(write=False)
            self.events.append((time, start, target, duration))
            previous = time
        self.events = tuple(self.events)

    def at(self, time_s):
        """Demand at one finite time, retaining the actual value at interruption."""
        if not np.isfinite(time_s): raise ValueError('Finite time required')
        for time, start, target, duration in reversed(self.events):
            if time_s >= time:
                fraction = np.minimum((time_s-time)/duration, 1.)
                return start+(target-start)*fraction
        return self.baseline.copy()

    def breakpoints(self, start_s, end_s):
        """Boundaries for integration, including interrupted/individual ramp ends."""
        if not np.isfinite(start_s) or not np.isfinite(end_s) or end_s < start_s:
            raise ValueError('Finite ordered interval required')
        points = {float(start_s), float(end_s)}
        for time, _, _, duration in self.events:
            for t in (time, *(time+duration)):
                if start_s < t < end_s: points.add(float(t))
        return sorted(points)

    def integral_mL(self, start_s, end_s):
        """Exact integral of piecewise-linear demand, including seconds/minute.

        Adding unused old-ramp end points after interruption is harmless: the
        actual schedule remains linear between all resulting boundaries.
        """
        points = self.breakpoints(start_s, end_s)
        total = np.zeros_like(self.baseline)
        for a, b in zip(points[:-1], points[1:]):
            total += (b-a)*(self.at(a)+self.at(b))/120
        return total


def scenario_vector(organs, budget, scenario, *, auxiliary_multiplier=1.):
    """Carry M1 parent reductions onto fixed organ splits; keep BIO independent.

    The 5%-residual scenario applies to all native tissue pools, including
    respiratory muscles: it is an extreme hypothetical comparator. Prescribed
    ventilation in that comparator is NOT evidence those muscles sustain it.
    """
    if not np.isfinite(auxiliary_multiplier) or auxiliary_multiplier < 0:
        raise ValueError('Invalid independent auxiliary multiplier')
    if scenario != 'hypothetical_5_percent' and scenario not in budget['scenarios']:
        raise ValueError('Unknown demand scenario')
    values = []
    for organ in organs:
        if organ.parent_pool == 'auxiliary_heart':
            multiplier = auxiliary_multiplier
        elif scenario == 'hypothetical_5_percent':
            multiplier = .05
        else:
            base = budget['baseline_mL_min'][organ.parent_pool]
            if base <= 0: raise ValueError('Positive reference parent demand required')
            multiplier = budget['scenarios'][scenario]['organ_mL_min'][organ.parent_pool]/base
        values.append(organ.demand_mL_min*multiplier)
    return np.asarray(values)
