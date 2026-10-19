"""T8.0 verification: a matched-budget search beside the structural argument.

The boundary-tracing search rests on three measured structural claims -- the
objective rises monotonically in cavity scale and rate, conductance is nearly
free energetically, and the constraints nearly separate. A structural argument
cannot verify itself, so the protocol requires an independent search at a
matched evaluation budget.

Two baselines run against the located optimum:

  * constrained differential evolution over the same bounds, which assumes no
    structure at all;
  * a Sobol-style space-filling sample at the same budget, which is the honest
    floor: if a search cannot beat uniform sampling, it has shown nothing.

Both are given the same number of evaluations and both report the best feasible
design found. Infeasible and failed evaluations are recorded with their cause
and are never counted as passes. A design is only compared on the declared
objective, auxiliary oxygen demand, and only when it satisfies both screens.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json, sys, time
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import optimize_m8_design as opt
from cardiac_redundancy.design import BOUNDS, DesignPoint

OUT = 'reports/m8-verification'
PROTOCOL = json.loads((ROOT / 'configs/m8-protocol.json').read_text(encoding='utf-8'))
#: Searched jointly: cavity scale, period, coronary conductance. The remaining
#: design controls stay pinned, as in the boundary trace, so the comparison is
#: like for like rather than a different problem.
AXES = ('scale', 'period_s', 'coronary_conductance')


def _bounds(centre, width):
    """A neighbourhood of the located optimum, clipped to the declared bounds."""
    out = {}
    for axis, key in (('scale', 'left_scale'), ('period_s', 'period_s'),
                      ('coronary_conductance', 'coronary_conductance')):
        lo, hi = BOUNDS[key][0]
        c = centre[axis]
        out[axis] = (max(lo, c-width[axis]), min(hi, c+width[axis]))
    return out


def evaluate_point(args):
    """One candidate: both screens, or a recorded reason it is not a result."""
    scale, period, coronary, kappa, tag = args
    d = opt.design(round(scale, 3), round(period, 3), round(coronary, 3))
    row = dict(tag=tag, scale=d.left_scale, period_s=d.period_s, rate_bpm=d.rate_bpm,
               coronary_conductance=d.coronary_conductance)
    try:
        status = opt.pair_status(d.left_scale, d.period_s, d.coronary_conductance, kappa)
    except Exception as exc:
        row.update(feasible=False, outcome='evaluation_failed', cause=repr(exc)[:160])
        return row
    row.update(feasible=bool(status.get('feasible')), pressure=status.get('pressure'),
               self_supply=status.get('self_supply'),
               auxiliary_demand_mL_min=status.get('auxiliary_demand_mL_min'),
               MAP_mmHg=status.get('MAP_mmHg'),
               outcome='feasible' if status.get('feasible') else 'infeasible')
    return row


def sample(bounds, budget, seed):
    """Space-filling sample: a scrambled Sobol sequence over the box."""
    from scipy.stats import qmc
    engine = qmc.Sobol(d=len(AXES), scramble=True, seed=seed)
    raw = engine.random(budget)
    lo = np.array([bounds[a][0] for a in AXES]); hi = np.array([bounds[a][1] for a in AXES])
    return lo+raw*(hi-lo)


def run(centre, kappa=.5, budget=32, seed=8101, width=None, workers=8):
    width = width or dict(scale=.08, period_s=.25, coronary_conductance=.5)
    bounds = _bounds(centre, width)
    points = sample(bounds, budget, seed)
    jobs = [(p[0], p[1], p[2], kappa, 'space_filling') for p in points]
    start = time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(evaluate_point, jobs))
    feasible = [r for r in rows if r['feasible']]
    best = min(feasible, key=lambda r: r['auxiliary_demand_mL_min']) if feasible else None
    result = dict(
        scope='matched-budget space-filling verification of the boundary-traced optimum',
        protocol_version=PROTOCOL['version'], horizon_s=PROTOCOL['evaluation_policy']['horizon_s'],
        compression_kappa=kappa, centre=centre, bounds={k: list(v) for k, v in bounds.items()},
        budget=budget, seed=seed, wall_s=time.perf_counter()-start,
        evaluated=len(rows), feasible=len(feasible),
        failures=[r for r in rows if r['outcome'] == 'evaluation_failed'],
        best=best, rows=rows,
        claim='best found within the tested neighbourhood and budget; not a global optimum')
    if best and centre.get('auxiliary_demand_mL_min'):
        result['beats_traced_optimum'] = bool(
            best['auxiliary_demand_mL_min'] < centre['auxiliary_demand_mL_min'])
        result['margin_mL_O2_min'] = centre['auxiliary_demand_mL_min']-best['auxiliary_demand_mL_min']
    dest = ROOT / OUT; dest.mkdir(parents=True, exist_ok=True)
    (dest / f'space-filling-kappa{kappa:g}.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    # usage: verify_m8_search.py <scale> <period> <coronary> <demand> [budget]
    a = sys.argv[1:]
    centre = dict(scale=float(a[0]), period_s=float(a[1]), coronary_conductance=float(a[2]),
                  auxiliary_demand_mL_min=float(a[3]))
    r = run(centre, budget=int(a[4]) if len(a) > 4 else 32)
    print('space-filling: %d evaluated, %d feasible, %.1f min' % (
        r['evaluated'], r['feasible'], r['wall_s']/60), flush=True)
    if r['best']:
        b = r['best']
        print('  best found: scale %.3f, %.0f/min, conductance %.3f -> %.3f mL O2/min' % (
            b['scale'], b['rate_bpm'], b['coronary_conductance'], b['auxiliary_demand_mL_min']), flush=True)
        print('  traced optimum %.3f -> %s by %.3f' % (
            centre['auxiliary_demand_mL_min'],
            'BEATEN' if r.get('beats_traced_optimum') else 'not beaten',
            abs(r.get('margin_mL_O2_min', 0.))), flush=True)
    else:
        print('  no feasible design found in the neighbourhood at this budget', flush=True)
    for f in r['failures']:
        print('  failed: scale %.3f %s' % (f['scale'], f['cause'][:90]), flush=True)
