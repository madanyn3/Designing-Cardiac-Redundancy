"""Re-map Arm A's feasible envelope across the valley, not just its cheapest corner.

T8.4 tested the four cheapest designs, all at cavity scale 0.536 to 0.566, and
concluded that none survives compression kappa 0. A cross-check then found that
a smaller design, scale 0.443, survives kappa 0 in Arm A. If that holds, the
earlier conclusion describes **which corner of the valley was sampled**, not a
property of the arm, because a smaller cavity carries a smaller coronary bed and
so loses less flow when compression is removed.

This takes a spread of feasible designs across the valley rather than its
cheapest few, and tests each against both compression alternatives. The result
is a boundary in cavity scale rather than a single pass or fail.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json, sys, time
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import optimize_m8_design as opt

OUT = 'reports/m8-envelope'
KAPPAS = (0., .5, 1.)
#: A spread across the valley: the cheapest design, the small design the
#: cross-check found feasible, and points between and above them.
DESIGNS = ((.443, 60/67.8, 1.591), (.496, 1., 1.203), (.536, 60/53.4, 1.508),
           (.566, 60/49.2, 1.027), (.650, 1., 1.203))


def one(job):
    scale, period, coronary, kappa = job
    row = dict(scale=scale, period_s=period, rate_bpm=60/period,
               coronary_conductance=coronary, kappa=kappa)
    try:
        status = opt.pair_status(scale, period, coronary, kappa)
    except Exception as exc:
        row.update(feasible=False, outcome='evaluation_failed', cause=repr(exc)[:140])
        return row
    row.update({k: status.get(k) for k in ('feasible', 'pressure', 'self_supply',
                                           'MAP_mmHg', 'auxiliary_demand_mL_min',
                                           'sustained_unmet_mL')})
    return row


def main(workers=6):
    start = time.perf_counter()
    jobs = [(s, p, c, k) for (s, p, c) in DESIGNS for k in KAPPAS]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(one, jobs))
    by_kappa = {k: sorted([r for r in rows if r['kappa'] == k], key=lambda r: r['scale'])
                for k in KAPPAS}
    survivors = {k: [r['scale'] for r in v if r.get('feasible')] for k, v in by_kappa.items()}
    result = dict(scope='Arm A envelope re-mapped across the valley at each compression case',
                  designs=[dict(scale=s, period_s=p, coronary_conductance=c) for s, p, c in DESIGNS],
                  kappas=list(KAPPAS), rows=rows, survivors=survivors,
                  wall_s=time.perf_counter()-start,
                  interpretation='boundaries of the sampled designs, not of the possible')
    dest = ROOT / OUT; dest.mkdir(parents=True, exist_ok=True)
    (dest / 'envelope.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n',
                                        encoding='utf-8')
    print('envelope re-map: %d designs x %d kappa, %.1f min' % (
        len(DESIGNS), len(KAPPAS), result['wall_s']/60), flush=True)
    print('  %6s %8s %8s %8s' % ('scale', 'kappa 0', 'kappa 0.5', 'kappa 1'), flush=True)
    for s, p, c in DESIGNS:
        cells = []
        for k in KAPPAS:
            r = next(x for x in rows if x['scale'] == s and x['kappa'] == k)
            cells.append('pass' if r.get('feasible') else
                         ('press' if r.get('pressure') != 'satisfied' else 'self'))
        print('  %6.3f %8s %8s %8s' % (s, *cells), flush=True)
    print('  (press = pressure screen failed, self = auxiliary self-supply failed)', flush=True)
    return result


if __name__ == '__main__':
    main()
