"""T8.3 latency relation: how long a configuration takes to become adequate.

Activation latency is not a design control a search may minimize; it is a
capability of whatever detection and recruitment path exists. What the model
can establish is the relation between that latency, the design, and the time
the circulation needs to reach the declared screens.

Reported per case: the first sustained window that satisfies the full triple,
the settling time measured from release, and the shortfall exposure accrued
before adequacy. A case that never satisfies the screen within the horizon is
reported as such, not as an infinite settling time.

The fixed two-second startup allowance of the frozen policy is one declared
choice of tolerable exposure. The relation here is what a different allowance
would have to be justified against, per condition, since the three pressure
metrics are not equally supported.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json, sys
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import optimize_m8_design as opt
from cardiac_redundancy.feasibility import FeasibilityPolicy, pressure_endpoint

PROTOCOL = json.loads((ROOT / 'configs/m8-protocol.json').read_text(encoding='utf-8'))


def declared_policy(horizon_s):
    """The frozen evaluation policy, so latency is judged like everything else."""
    declared = dict(PROTOCOL['evaluation_policy'])
    declared['unresolved_requirements'] = tuple(declared['unresolved_requirements'])
    declared['horizon_s'] = horizon_s
    return FeasibilityPolicy(**declared)

OUT = 'reports/m8-optimize'
RELEASES = [.5, 1., 1.5, 2., 2.5, 3.]


def per_metric_adequacy(endpoint, release_s):
    """First window at or after release in which each metric is satisfied.

    Separating the metrics matters because they are not equally supported: the
    MAP 65 screen is an AHA post-ROSC recommendation, while SBP 90 and DBP 50
    are transferred from one observational surgical cohort, and AHA withdrew
    its systolic target. A design adequate on MAP but not SBP is a different
    statement from one inadequate on all three.
    """
    windows = endpoint['windows']
    out = {}
    for metric in endpoint['thresholds_mmHg']:
        first = next((w['start_s'] for w in windows
                      if w['start_s'] >= release_s-1e-9 and w['status'][metric] == 'satisfied'), None)
        out[metric] = dict(first_satisfied_window_s=first,
                           settling_from_release_s=None if first is None else first-release_s)
    enforced = endpoint['enforced_metrics']
    joint = next((w['start_s'] for w in windows
                  if w['start_s'] >= release_s-1e-9
                  and all(w['status'][k] == 'satisfied' for k in enforced)), None)
    out['joint'] = dict(first_satisfied_window_s=joint,
                        settling_from_release_s=None if joint is None else joint-release_s)
    return out


def one(job):
    scale, period, coronary, kappa, release = job
    d = opt.design(scale, period, coronary)
    record = opt.hydraulics(d, kappa, release_s=release)
    pressure = record['pressure']
    adequacy = per_metric_adequacy(pressure, release)
    return dict(scale=scale, rate_bpm=60/period, coronary_conductance=coronary,
                compression_kappa=kappa, release_s=release,
                pressure_status=pressure['status'],
                numerical_status=record['numerical_status'],
                violating_sustained_windows=pressure.get('violating_sustained_windows'),
                startup_shortfall_mmHg_s=pressure['startup_shortfall_mmHg_s'],
                sustained_shortfall_mmHg_s={k: v['shortfall_mmHg_s']
                                            for k, v in pressure['sustained_exposure'].items()},
                status_by_criterion=pressure['status_by_criterion'],
                plateau_reached=pressure.get('plateau_reached'),
                horizon_s=record['duration_s'],
                MAP_mmHg=record['runs'][-1]['final_MAP_mmHg'],
                auxiliary_demand_mL_min=opt.objective(record), adequacy=adequacy)


def main(workers=6):
    # The v4 boundary designs at the settled horizon, plus a larger comparator,
    # plus both designs protocol v5 reports. The first two hold only +1.33 and
    # +0.69 mmHg of pressure margin, below the declared 3.0 mmHg minimum, so the
    # relation has to be re-established on designs that clear it: a slope read
    # off a design sitting on the threshold describes the threshold crossing
    # rather than the design's latency tolerance. They are kept for comparison.
    designs = [(.475, 1., 1., .5), (.387, .75, 1., .5), (.65, 1., 1., .5),
               (.536, 1.043, 1.748, .5), (.749, 1.304, .953, .5)]
    jobs = [(s, p, c, k, r) for (s, p, c, k) in designs for r in RELEASES]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(one, jobs))
    result = dict(scope='activation-latency relation at fixed designs; prescribed release stands '
                        'for detection plus recruitment',
                  releases_s=RELEASES, designs=designs, rows=rows,
                  interpretation='settling times are window-resolution brackets on the declared '
                                 'screens, not tissue tolerances')
    (ROOT / OUT / 'latency.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('%6s %5s %6s  %-9s %8s %8s %8s %8s' % (
        'scale', 'rate', 'release', 'pressure', 'joint_s', 'MAP_s', 'SBP_s', 'DBP_s'), flush=True)
    for r in sorted(rows, key=lambda x: (x['scale'], -x['rate_bpm'], x['release_s'])):
        g = lambda k: (r['adequacy'][k]['settling_from_release_s'])
        fmt = lambda v: '    none' if v is None else '%8.1f' % v
        print('%6.3f %5.0f %6.1f  %-9s %s %s %s %s' % (
            r['scale'], r['rate_bpm'], r['release_s'], r['pressure_status'],
            fmt(g('joint')), fmt(g('mean')), fmt(g('maximum')), fmt(g('minimum'))), flush=True)
    return result


if __name__ == '__main__':
    main()
