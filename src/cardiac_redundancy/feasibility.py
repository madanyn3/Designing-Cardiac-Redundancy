"""T8.1 endpoint evaluator: one declared verdict per evaluated design.

Every screen here is reused from constraints.py, so an M8 verdict is the same
arithmetic the earlier milestones reported, assembled in one place and applied
identically to every candidate.

Three protocol rules are enforced rather than left to the caller.

The startup allowance is a fixed interval measured from arrest, declared before
the search. A design that detects or activates late therefore pays for the
extra exposure; it does not receive a longer grace period because its own
release came later.

A pressure verdict uses the full declared triple (mean, maximum and minimum
over complete fixed windows), not the mean alone.

Satisfaction is always "no violation of the declared screens through the
tested horizon". It is never support duration, survival, recovery, or a result
outside the sampled window. Numerical failure and model-domain exit stay
separate from physiological violation, so a design is never credited with a
pass it obtained by failing to integrate.
"""
from dataclasses import dataclass, field
import math
import numpy as np
from .constraints import (combined_status, first_sampled_violation, oxygen_screen,
                          positive_exposure, pressure_windows, series)

#: Declared joint pressure screen from configs/auxiliary-heart-design.yaml:
#: MAP >= 65, SBP >= 90, DBP >= 50 over complete 1 s windows, no violations
#: allowed. A transfer from post-ROSC references, not a validated minimum.
PRESSURE_SCREEN_MMHG = dict(mean=65., maximum=90., minimum=50.)


@dataclass(frozen=True)
class FeasibilityPolicy:
    """Frozen evaluation policy. Freeze under T8.1 before any main run."""
    thresholds_mmHg: dict = field(default_factory=lambda: dict(PRESSURE_SCREEN_MMHG))
    window_s: float = 1.
    startup_allowance_s: float = 2.
    oxygen_tolerance_mL_s: float = 1e-8
    horizon_s: float = 10.
    unresolved_requirements: tuple = ()

    def __post_init__(self):
        if set(self.thresholds_mmHg) != {'mean', 'maximum', 'minimum'}:
            raise ValueError('Pressure screen requires mean, maximum and minimum')
        if any(not math.isfinite(v) for v in self.thresholds_mmHg.values()):
            raise ValueError('Finite pressure thresholds required')
        if not self.thresholds_mmHg['minimum'] <= self.thresholds_mmHg['mean'] <= self.thresholds_mmHg['maximum']:
            raise ValueError('Ordered DBP <= MAP <= SBP screen required')
        for name in ('window_s', 'oxygen_tolerance_mL_s', 'horizon_s'):
            v = getattr(self, name)
            if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or v <= 0:
                raise ValueError(f'Positive finite {name} required')
        a = self.startup_allowance_s
        if isinstance(a, bool) or not isinstance(a, (int, float)) or not math.isfinite(a) or a < 0:
            raise ValueError('Finite nonnegative startup allowance required')
        if abs(a/self.window_s-round(a/self.window_s)) > 1e-9:
            raise ValueError('Startup allowance must be a whole number of windows')
        if a >= self.horizon_s:
            raise ValueError('Startup allowance must leave a sustained interval')


def numerical_status(termination, paired_errors, paired_tolerance):
    """Separate a domain exit and a failed refinement check from physiology."""
    if termination is not None:
        return 'model_domain_limit'
    values = [float(e) for e in paired_errors]
    if any(not math.isfinite(e) for e in values) or (values and max(values) > paired_tolerance):
        return 'numerical_failure'
    return 'valid'


def pressure_endpoint(policy, time_s, ao_mmHg, release_s=None):
    """Startup exposure and sustained screen, reported separately.

    Only complete windows beginning at or after the declared allowance bind
    the verdict. Startup windows are reported with their shortfall exposure so
    the pre-activation interval is never silently discarded.
    """
    if not isinstance(policy, FeasibilityPolicy):
        raise ValueError('Frozen feasibility policy required')
    t, p = series(time_s, ao_mmHg)
    full = pressure_windows(t, p, policy.thresholds_mmHg, policy.window_s)
    start = float(t[0])
    startup, sustained = [], []
    for record in full['windows']:
        (startup if record['start_s'] < start+policy.startup_allowance_s-1e-9 else sustained).append(record)
    violating = [r for r in sustained if 'violated' in r['status'].values()]
    status = ('violated' if violating else
              'satisfied' if sustained and full['incomplete_tail_s'] is None else 'unassessed')
    shortfall = {k: float(sum(max(policy.thresholds_mmHg[k]-r['metrics_mmHg'][k], 0.)*policy.window_s
                              for r in startup)) for k in policy.thresholds_mmHg}
    return dict(status=status, thresholds_mmHg=dict(policy.thresholds_mmHg),
                window_s=policy.window_s, startup_allowance_s=policy.startup_allowance_s,
                release_s=None if release_s is None else float(release_s),
                startup_windows=len(startup), sustained_windows=len(sustained),
                startup_shortfall_mmHg_s=shortfall,
                first_violating_sustained_window=violating[0] if violating else None,
                sustained_exposure=({k: dict(
                    violating_window_duration_s=float(sum(policy.window_s for r in sustained if r['status'][k] == 'violated')),
                    shortfall_mmHg_s=float(sum(max(policy.thresholds_mmHg[k]-r['metrics_mmHg'][k], 0.)*policy.window_s
                                               for r in sustained))) for k in policy.thresholds_mmHg}),
                incomplete_tail_s=full['incomplete_tail_s'],
                interpretation='central AO fixed-window screens with a declared startup allowance; '
                               'not regional pressure adequacy and not demonstrated viability')


def oxygen_endpoint(policy, time_s, groups, baseline_mL_min):
    """Per-group demand-satisfaction screens, including auxiliary self-supply.

    `groups` maps a role name ('active', 'preservation', 'pump') to the demand,
    uptake, delivered and cumulative-deficit series for that group.
    """
    if not isinstance(policy, FeasibilityPolicy):
        raise ValueError('Frozen feasibility policy required')
    if not groups:
        raise ValueError('At least one organ group required')
    out = {}
    for role, s in groups.items():
        screen = oxygen_screen(time_s, s['demand'], s['uptake'], s['delivered'],
                               s['cumulative_deficit'], baseline_mL_min,
                               policy.oxygen_tolerance_mL_s)
        t, _ = series(time_s, s['demand'])
        start = float(t[0])
        # The startup and sustained intervals are both closed at the allowance
        # boundary, so neither integral loses the partial interval beside it.
        mask = t >= start+policy.startup_allowance_s
        early = t <= start+policy.startup_allowance_s+1e-12
        unmet = np.asarray(s['demand'], dtype=float)-np.asarray(s['uptake'], dtype=float)
        screen['sustained'] = dict(
            status='violated' if np.any(unmet[mask] > policy.oxygen_tolerance_mL_s) else 'satisfied',
            first_violation=first_sampled_violation(t[mask], unmet[mask] > policy.oxygen_tolerance_mL_s),
            excess_over_numerical_tolerance=positive_exposure(t[mask], unmet[mask]-policy.oxygen_tolerance_mL_s))
        screen['startup_unmet_mL'] = float(np.trapezoid(np.maximum(unmet[early], 0.), t[early])) if early.sum() > 1 else 0.
        out[role] = screen
    return out


def verdict(policy, numerical, pressure, oxygen, reported=None):
    """Assemble the single declared status, keeping every qualification visible."""
    if not isinstance(policy, FeasibilityPolicy):
        raise ValueError('Frozen feasibility policy required')
    statuses = [o['sustained']['status'] for o in oxygen.values()]
    combined = combined_status(numerical, pressure['status'], statuses,
                               list(policy.unresolved_requirements))
    binding = []
    if pressure['status'] == 'violated': binding.append('pressure')
    binding += [f'oxygen:{role}' for role, o in oxygen.items() if o['sustained']['status'] == 'violated']
    return dict(status=combined['status'], numerical_status=numerical,
                binding_constraints=binding, horizon_s=policy.horizon_s,
                pressure=pressure, oxygen=oxygen, reported=dict(reported or {}),
                biological_recovery=combined['biological_recovery'],
                unresolved_requirements=combined['unresolved_requirements'],
                meaning='no violation of the declared screens through the tested horizon; '
                        'not support duration, survival or recovery',
                interpretation=combined['interpretation'])
