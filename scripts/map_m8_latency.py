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
from cardiac_redundancy.constraints import pressure_windows
from cardiac_redundancy.feasibility import PRESSURE_SCREEN_MMHG

OUT = 'reports/m8-optimize'
RELEASES = [.5, 1., 1.5, 2., 2.5, 3.]


def per_metric_adequacy(time_s, ao_mmHg, release_s):
    """First window at or after release in which each metric is satisfied.

    Separating the metrics matters because they are not equally supported: the
    MAP 65 screen is an AHA post-ROSC recommendation, while SBP 90 and DBP 50
    are transferred from one observational surgical cohort, and AHA withdrew
    its systolic target. A design adequate on MAP but not SBP is a different
    statement from one inadequate on all three.
    """
    screen = pressure_windows(time_s, ao_mmHg, PRESSURE_SCREEN_MMHG, 1.)
    out = {}
    for metric in PRESSURE_SCREEN_MMHG:
        first = next((w['start_s'] for w in screen['windows']
                      if w['start_s'] >= release_s-1e-9 and w['status'][metric] == 'satisfied'), None)
        out[metric] = dict(first_satisfied_window_s=first,
                           settling_from_release_s=None if first is None else first-release_s)
    joint = next((w['start_s'] for w in screen['windows']
                  if w['start_s'] >= release_s-1e-9
                  and 'violated' not in w['status'].values()), None)
    out['joint'] = dict(first_satisfied_window_s=joint,
                        settling_from_release_s=None if joint is None else joint-release_s)
    return out


def one(job):
    scale, period, coronary, kappa, release = job
    d = opt.design(scale, period, coronary)
    record = opt.hydraulics(d, kappa, release_s=release)
    name = opt.name_of(d, kappa, release)
    a = np.load(ROOT / OUT / name / 'hydraulic-refined.npz')
    t = a['time_s']; ao = a['pressures_mmHg'][:, 2]
    adequacy = per_metric_adequacy(t, ao, release)
    pressure = record['pressure']
    return dict(scale=scale, rate_bpm=60/period, coronary_conductance=coronary,
                compression_kappa=kappa, release_s=release,
                pressure_status=pressure['status'],
                numerical_status=record['numerical_status'],
                violating_sustained_windows=pressure.get('violating_sustained_windows'),
                startup_shortfall_mmHg_s=pressure['startup_shortfall_mmHg_s'],
                sustained_shortfall_mmHg_s={k: v['shortfall_mmHg_s']
                                            for k, v in pressure['sustained_exposure'].items()},
                MAP_mmHg=record['runs'][-1]['final_MAP_mmHg'],
                auxiliary_demand_mL_min=opt.objective(record), adequacy=adequacy)


def main(workers=6):
    designs = [(.651, 1., .969, .5), (.8, 1., .969, .5), (.651, .75, .969, .5)]
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
