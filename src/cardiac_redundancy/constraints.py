"""Research screens, exposure proxies and explicitly unassessed requirements.

No screen here predicts injury or recovery. Sampled crossings are brackets,
not solver-localized events. All time inputs are seconds.
"""
import numpy as np


def series(time, values):
    t = np.asarray(time, dtype=float); v = np.asarray(values, dtype=float)
    if t.ndim != 1 or len(t) < 2 or v.shape != t.shape or np.any(~np.isfinite(t)) or np.any(~np.isfinite(v)) or np.any(np.diff(t) <= 0):
        raise ValueError('Finite aligned series with strictly increasing time required')
    return t, v


def first_sampled_violation(time, failed):
    """First observed crossing bracket, including an initial violation."""
    idx = np.flatnonzero(failed)
    if not len(idx): return dict(observed=False, horizon_s=float(time[-1]), interpretation='not observed within sampled horizon')
    k = int(idx[0])
    return dict(observed=True, initial=bool(k == 0), bracket_s=[float(time[max(0, k-1)]), float(time[k])],
                interpretation='sample bracket; not exact onset or injury time')


def positive_exposure(time, margin):
    """Integral of positive piecewise-linear margin; split zero crossings.

    Used for unmet rate or pressure shortfall. Duration is exposure on the
    interpolant, not a validated tissue tolerance or unsampled event measure.
    """
    t, v = series(time, margin); duration = area = 0.
    for a, b, x, y in zip(t[:-1], t[1:], v[:-1], v[1:]):
        dt = b-a
        if x >= 0 and y >= 0:
            area += dt*(x+y)/2
            if x > 0 or y > 0: duration += dt
        elif x > 0 or y > 0:
            positive = max(x, y); fraction = positive/abs(y-x)
            duration += dt*fraction; area += .5*positive*dt*fraction
    return dict(duration_s=float(duration), integral=float(area))


def pressure_windows(time, pressure, thresholds, window_s=1.):
    t, p = series(time, pressure)
    if not np.isfinite(window_s) or window_s <= 0 or set(thresholds) != {'mean', 'maximum', 'minimum'} or any(not np.isfinite(x) for x in thresholds.values()):
        raise ValueError('Invalid pressure policy')
    count = int(np.floor((t[-1]-t[0])/window_s+1e-10)); records = []
    for i in range(count):
        a = t[0]+i*window_s; b = a+window_s
        if b > t[-1]+1e-9: break
        grid = np.r_[a, t[(t > a) & (t < b)], b]; values = np.interp(grid, t, p)
        metrics = dict(mean=float(np.trapezoid(values, grid)/window_s), maximum=float(values.max()), minimum=float(values.min()))
        records.append(dict(start_s=float(a), end_s=float(b), metrics_mmHg=metrics,
            window_pressure_range_mmHg=metrics['maximum']-metrics['minimum'],
            status={k:'satisfied' if v >= thresholds[k] else 'violated' for k, v in metrics.items()}))
    tail_start = t[0]+len(records)*window_s
    coverage_complete = abs(t[-1]-tail_start) <= 1e-9
    violated = [r for r in records if 'violated' in r['status'].values()]
    status = 'violated' if violated else 'satisfied' if records and coverage_complete else 'unassessed'
    exposures = {k:dict(violating_window_duration_s=float(sum(window_s for r in records if r['status'][k]=='violated')),
        shortfall_mmHg_s=float(sum(max(thresholds[k]-r['metrics_mmHg'][k], 0)*window_s for r in records))) for k in thresholds}
    return dict(status=status, window_s=window_s, thresholds_mmHg=thresholds, windows=records,
        first_violating_window=violated[0] if violated else None, exposure=exposures,
        incomplete_tail_s=None if coverage_complete else [float(tail_start),float(t[-1])],
        interpretation='central AO fixed-window research screens; not regional pressure adequacy')


def oxygen_screen(time, demand, uptake, delivered, cumulative_deficit, baseline_mL_min, tolerance_mL_s, ledger_tolerance_mL=1e-7, blood_extraction=None):
    t, d = series(time, demand)
    _, u = series(time, uptake); _, supply = series(time, delivered); _, debt = series(time, cumulative_deficit)
    _, extracted = series(time, uptake if blood_extraction is None else blood_extraction)
    if not np.isfinite(tolerance_mL_s) or tolerance_mL_s < 0 or not np.isfinite(baseline_mL_min) or baseline_mL_min <= 0 or not np.isfinite(ledger_tolerance_mL) or ledger_tolerance_mL < 0:
        raise ValueError('Invalid numerical tolerance or baseline')
    if np.any(d < 0) or np.any(u < 0) or np.any(supply < 0) or np.any(u > d+1e-10) or np.any(extracted < 0) or np.any(extracted > supply+1e-10) or np.any(np.diff(debt) < -ledger_tolerance_mL) or np.any(debt < -ledger_tolerance_mL):
        raise ValueError('Invalid uptake/demand ledger')
    unmet = d-u
    ratio = np.divide(extracted, supply, out=np.full_like(u, np.nan), where=supply > 0)
    finite = ratio[np.isfinite(ratio)]
    return dict(status='violated' if np.any(unmet > tolerance_mL_s) else 'satisfied',
        tolerance_mL_s=float(tolerance_mL_s), peak_unmet_mL_s=float(unmet.max()),
        accumulated_unmet_mL=float(debt[-1]), baseline_demand_equivalent_min=float(debt[-1]/baseline_mL_min),
        ledger_tolerance_mL=float(ledger_tolerance_mL),maximum_ledger_backward_step_mL=float(max(0.,-np.min(np.diff(debt)))),
        first_violation=first_sampled_violation(t, unmet > tolerance_mL_s),
        excess_over_numerical_tolerance=positive_exposure(t, unmet-tolerance_mL_s),
        extraction_fraction_range=None if not len(finite) else [float(finite.min()),float(finite.max())],
        extraction_undefined_samples=int(np.sum(supply == 0)),
        interpretation='demand satisfaction under assumed extraction; accumulated deficit is not injury dose')


def combined_status(numerical_status, pressure_status, oxygen_statuses, unresolved_requirements):
    """A known violation and missing requirements must both remain visible."""
    allowed = {'satisfied', 'violated', 'unassessed'}
    if pressure_status not in allowed or any(s not in allowed for s in oxygen_statuses): raise ValueError('Invalid screen status')
    if numerical_status not in {'valid', 'numerical_failure', 'model_domain_limit'}: raise ValueError('Invalid numerical status')
    statuses = [pressure_status, *oxygen_statuses]
    if numerical_status != 'valid': overall = numerical_status
    elif 'violated' in statuses: overall = 'violated_research_screen'
    elif unresolved_requirements or 'unassessed' in statuses: overall = 'indeterminate'
    else: overall = 'satisfied_within_declared_model'
    return dict(status=overall, numerical_status=numerical_status,
        biological_recovery='unassessed', unresolved_requirements=list(unresolved_requirements),
        interpretation='no clinical safety or recovery conclusion')
