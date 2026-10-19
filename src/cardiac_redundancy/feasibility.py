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


#: Requirements the model does not evaluate, from configs/m4-constraint-policy.json.
#: Carried on every policy so no run can report an unqualified pass while they
#: stand: combined_status downgrades a clean numeric result to 'indeterminate'.
#: Pass an empty tuple only to inspect the screens in isolation, never to
#: report a verdict.
UNRESOLVED_REQUIREMENTS = (
    'regional organ pressure minima',
    'whole-brain regional distribution and function',
    'recovery-linked oxygen-deficit thresholds for every tissue',
    'endogenous suppression mechanism and reversibility',
    'respiratory drive, airway, ventilation-demand coupling and CO2 clearance',
    'tissue oxygen accessibility and microvascular exchange',
    'glycogen/glucose stock and delivery limits',
    'blood compatibility, endothelial injury and reperfusion injury',
    'biological auxiliary self-demand calibration where applicable',
)


#: Declared evaluation criteria. The strict triple is the M2 research screen.
#: The MAP-primary criterion enforces only the mean, because adult post-arrest
#: haemodynamic guidance is expressed as a minimum MAP (AHA 2025 Part 11: MAP
#: at least 65 mmHg, with four randomized trials finding no benefit from higher
#: targets) and states no adult systolic target, while SBP 90 and DBP 50 derive
#: from one perioperative observational cohort with radial measurement and a
#: five-minute exposure definition. Both criteria are always reported; neither
#: replaces the other, and selecting one after seeing which passes is
#: forbidden by the protocol.
CRITERIA = dict(strict_triple=('mean', 'maximum', 'minimum'), map_primary=('mean',))


@dataclass(frozen=True)
class FeasibilityPolicy:
    """Frozen evaluation policy. Freeze under T8.1 before any main run."""
    thresholds_mmHg: dict = field(default_factory=lambda: dict(PRESSURE_SCREEN_MMHG))
    window_s: float = 1.
    startup_allowance_s: float = 2.
    oxygen_tolerance_mL_s: float = 1e-8
    horizon_s: float = 10.
    unresolved_requirements: tuple = UNRESOLVED_REQUIREMENTS
    #: 'fixed' uses the declared 1 s window throughout. 'cycle' follows
    #: docs/m2-pressure-requirements.md: the integral mean, maximum and minimum
    #: over each auxiliary activation cycle once the pump is running, with the
    #: fixed window retained for the nonpulsatile pre-release interval.
    window_mode: str = 'fixed'
    #: Which criterion in CRITERIA decides `status`. The others are still
    #: computed and reported under `status_by_criterion`.
    criterion: str = 'strict_triple'
    #: Spread the last five sustained cycle means must fall within for the
    #: transient to count as settled at the horizon.
    plateau_tolerance_mmHg: float = .05

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
        if self.window_mode not in ('fixed', 'cycle'):
            raise ValueError('Window mode must be fixed or cycle')
        if self.criterion not in CRITERIA:
            raise ValueError('Criterion must name an entry in CRITERIA')
        if (isinstance(self.plateau_tolerance_mmHg, bool)
                or not isinstance(self.plateau_tolerance_mmHg, (int, float))
                or not math.isfinite(self.plateau_tolerance_mmHg)
                or self.plateau_tolerance_mmHg <= 0):
            raise ValueError('Positive finite plateau tolerance required')


def numerical_status(termination, paired_errors, paired_tolerance):
    """Separate a domain exit and a failed refinement check from physiology."""
    if termination is not None:
        return 'model_domain_limit'
    values = [float(e) for e in paired_errors]
    if any(not math.isfinite(e) for e in values) or (values and max(values) > paired_tolerance):
        return 'numerical_failure'
    return 'valid'


def _windowed(policy, t, p, release_s, period_s):
    """Window records and any incomplete tail, honouring the window mode.

    In cycle mode the pre-release interval is nonpulsatile, because the native
    heart is arrested and the auxiliary pump has not been released, so the
    declared fixed window is kept there and cycle-aligned windows of one
    activation period begin at release. This follows the M2 evaluation
    contract rather than imposing one grid on two different regimes.
    """
    if policy.window_mode == 'fixed' or release_s is None:
        full = pressure_windows(t, p, policy.thresholds_mmHg, policy.window_s)
        return full['windows'], full['incomplete_tail_s']
    if period_s is None or not math.isfinite(period_s) or period_s <= 0:
        raise ValueError('Cycle mode requires the auxiliary activation period')
    windows = []
    tail = None
    pre = t <= release_s+1e-12
    if pre.sum() > 1 and release_s-float(t[0]) >= policy.window_s:
        early = pressure_windows(t[pre], p[pre], policy.thresholds_mmHg, policy.window_s)
        windows += early['windows']
    post = t >= release_s-1e-12
    if post.sum() > 1:
        late = pressure_windows(t[post], p[post], policy.thresholds_mmHg, period_s)
        windows += late['windows']
        tail = late['incomplete_tail_s']
    return windows, tail


def pressure_endpoint(policy, time_s, ao_mmHg, release_s=None, period_s=None):
    """Startup exposure and sustained screen, reported separately.

    Only complete windows beginning at or after the declared allowance bind
    the verdict. Startup windows are reported with their shortfall exposure so
    the pre-activation interval is never silently discarded.

    Every criterion in CRITERIA is evaluated and reported; `status` follows the
    policy's declared one.
    """
    if not isinstance(policy, FeasibilityPolicy):
        raise ValueError('Frozen feasibility policy required')
    t, p = series(time_s, ao_mmHg)
    records, tail = _windowed(policy, t, p, release_s, period_s)
    start = float(t[0])
    startup, sustained = [], []
    for record in records:
        (startup if record['start_s'] < start+policy.startup_allowance_s-1e-9 else sustained).append(record)
    enforced = CRITERIA[policy.criterion]

    # A leftover fraction of a window at the end of a completed trajectory is a
    # grid remainder, not an unobserved horizon: cycle-aligned windows begin at
    # release, so a partial cycle is left over whenever the horizon is not an
    # exact multiple of the period. A run that stopped early is caught by the
    # numerical status instead, so voiding these would discard physiology.
    # Anything larger than one window is still unassessed.
    widths = [r['end_s']-r['start_s'] for r in records] or [policy.window_s]
    remainder = 0. if tail is None else float(tail[1]-tail[0])
    # Only cycle mode earns this allowance. In fixed mode the grid starts at the
    # first sample, so a remainder means the trajectory itself ended mid-window,
    # which the M4 incomplete-horizon policy requires to stay unassessed.
    covered = tail is None or (policy.window_mode == 'cycle' and remainder <= max(widths)+1e-9)

    def decide(metrics):
        bad = [r for r in sustained if any(r['status'][k] == 'violated' for k in metrics)]
        return ('violated' if bad else
                'satisfied' if sustained and covered else 'unassessed'), bad

    by_criterion = {name: decide(metrics)[0] for name, metrics in CRITERIA.items()}
    status, violating = decide(enforced)
    span = lambda r: r['end_s']-r['start_s']
    shortfall = {k: float(sum(max(policy.thresholds_mmHg[k]-r['metrics_mmHg'][k], 0.)*span(r)
                              for r in startup)) for k in policy.thresholds_mmHg}
    # Whether the transient has flattened by the horizon. Cycle-mean pressure
    # declines toward a steady state after activation, so a horizon that ends
    # mid-transient credits a value the design does not maintain. This reports
    # the spread of the final cycles rather than assuming settlement.
    tail_means = [r['metrics_mmHg']['mean'] for r in sustained[-5:]]
    spread = float(max(tail_means)-min(tail_means)) if len(tail_means) >= 2 else None
    plateau = None if spread is None else bool(spread <= policy.plateau_tolerance_mmHg)
    return dict(status=status, status_by_criterion=by_criterion, criterion=policy.criterion,
                plateau_reached=plateau, plateau_spread_mmHg=spread,
                plateau_tolerance_mmHg=policy.plateau_tolerance_mmHg,
                final_cycle_means_mmHg=[float(v) for v in tail_means],
                enforced_metrics=list(enforced),
                thresholds_mmHg=dict(policy.thresholds_mmHg),
                window_mode=policy.window_mode, cycle_period_s=period_s,
                window_s=policy.window_s, startup_allowance_s=policy.startup_allowance_s,
                release_s=None if release_s is None else float(release_s),
                startup_windows=len(startup), sustained_windows=len(sustained),
                # A window counts once however many metrics it violates, so a
                # per-metric count is never mistaken for the overall total.
                violating_sustained_windows=len(violating),
                first_satisfied_sustained_window_s=next(
                    (r['start_s'] for r in sustained if 'violated' not in r['status'].values()), None),
                windows=[dict(r, interval='startup') for r in startup]
                        +[dict(r, interval='sustained') for r in sustained],
                startup_shortfall_mmHg_s=shortfall,
                first_violating_sustained_window=violating[0] if violating else None,
                sustained_exposure=({k: dict(
                    violating_window_duration_s=float(sum(span(r) for r in sustained if r['status'][k] == 'violated')),
                    shortfall_mmHg_s=float(sum(max(policy.thresholds_mmHg[k]-r['metrics_mmHg'][k], 0.)*span(r)
                                               for r in sustained))) for k in policy.thresholds_mmHg}),
                incomplete_tail_s=tail, unassessed_remainder_s=remainder,
                remainder_within_one_window=covered,
                interpretation='central AO fixed-window screens with a declared startup allowance; '
                               'not regional pressure adequacy and not demonstrated viability')


def pressure_margin(endpoint):
    """Signed mmHg margin of the worst sustained window against the screen.

    The minimum over sustained windows and over the three metrics of
    (value - threshold). Nonnegative exactly when the screen is satisfied, so
    a root-find on this quantity locates the same boundary the boolean verdict
    reports, in a quarter of the evaluations a bisection on the boolean needs.

    It is a continuous but only piecewise-smooth function: the window and
    metric attaining the minimum switch along the boundary, so a safeguarded
    bracketing method is required and a raw secant or Newton step is not.
    Returns None when no complete sustained window exists to judge.
    """
    sustained = [w for w in endpoint.get('windows', ()) if w['interval'] == 'sustained']
    if not sustained:
        return None
    thresholds = endpoint['thresholds_mmHg']
    enforced = endpoint.get('enforced_metrics') or list(thresholds)
    return min(w['metrics_mmHg'][k]-thresholds[k] for w in sustained for k in enforced)


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
