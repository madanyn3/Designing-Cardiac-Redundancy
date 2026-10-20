"""T8.0 joint optimization over cavity scale, rate and coronary conductance.

Boundary tracing located its optimum on a grid of three discrete rates, and a
matched-budget sample beat it by finding 49 beats/min between those grid lines.
The rate axis is continuous and the objective is shallow along it, so this
searches all three design axes jointly instead.

    minimize    auxiliary oxygen demand (mL O2/min)
    over        cavity scale, activation period, coronary conductance
    subject to  the declared pressure screen and auxiliary self-supply

Two stages, both respecting the constraints rather than penalising them:

  * a scrambled Sobol sample over the declared box, which makes no structural
    assumption and is embarrassingly parallel;
  * a compass (pattern) search from the best feasible point, which the protocol
    names as a fallback method for expensive black-box problems with hidden
    constraints. Each iteration evaluates two steps per axis in parallel and
    halves the step when none improves.

Evaluation is staged: a candidate that fails the pressure screen costs only the
hydraulic solve, because self-supply is never evaluated on a design that has
already failed. Infeasible and failed points are recorded with their cause and
are never counted as passes.

Nothing here claims a global optimum. The result is the best feasible design
found within the declared box and the stated evaluation budget.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json, sys, time
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import optimize_m8_design as opt
from cardiac_redundancy.design import BOUNDS

OUT = 'reports/m8-joint'
PROTOCOL = json.loads((ROOT / 'configs/m8-protocol.json').read_text(encoding='utf-8'))
AXES = ('scale', 'period_s', 'coronary_conductance')
#: Search box. Wider than the region the traced search covered, especially on
#: rate, where the sample found the optimum outside the traced grid. Clipped to
#: the declared parameter bounds; a wider box costs only failed pressure solves.
BOX = dict(scale=(.30, .80), period_s=(.80, 1.50), coronary_conductance=(.70, 1.90))
STEP0 = dict(scale=.04, period_s=.08, coronary_conductance=.12)


def clip(point):
    out = {}
    for axis, key in (('scale', 'left_scale'), ('period_s', 'period_s'),
                      ('coronary_conductance', 'coronary_conductance')):
        lo, hi = BOUNDS[key][0]
        out[axis] = float(min(max(point[axis], max(lo, BOX[axis][0])), min(hi, BOX[axis][1])))
    return out


def evaluate(args):
    point, kappa, tag = args
    p = clip(point)
    d = opt.design(round(p['scale'], 3), round(p['period_s'], 3),
                   round(p['coronary_conductance'], 3))
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
               sustained_unmet_mL=status.get('sustained_unmet_mL'),
               outcome='feasible' if status.get('feasible') else 'infeasible')
    return row


def sobol(n, seed):
    from scipy.stats import qmc
    raw = qmc.Sobol(d=3, scramble=True, seed=seed).random(n)
    lo = np.array([BOX[a][0] for a in AXES]); hi = np.array([BOX[a][1] for a in AXES])
    scaled = lo+raw*(hi-lo)
    return [dict(zip(AXES, row)) for row in scaled]


def best_of(rows):
    feasible = [r for r in rows if r.get('feasible')]
    return min(feasible, key=lambda r: r['auxiliary_demand_mL_min']) if feasible else None


def compass(centre, kappa, pool, iterations, step, history):
    """Pattern search: step each axis both ways, shrink when none improves."""
    current = dict(centre); trace = []
    best = None
    for it in range(iterations):
        candidates = []
        for axis in AXES:
            for sign in (+1, -1):
                point = dict(current); point[axis] = current[axis]+sign*step[axis]
                candidates.append((point, kappa, f'compass{it}:{axis}{sign:+d}'))
        rows = list(pool.map(evaluate, candidates))
        history.extend(rows)
        contender = best_of(rows)
        incumbent = best or best_of([r for r in history if r['tag'] == 'incumbent']) or None
        reference = best['auxiliary_demand_mL_min'] if best else centre['auxiliary_demand_mL_min']
        if contender and contender['auxiliary_demand_mL_min'] < reference:
            best = contender
            current = {a: contender[{'scale': 'scale', 'period_s': 'period_s',
                                     'coronary_conductance': 'coronary_conductance'}[a]] for a in AXES}
            trace.append(dict(iteration=it, moved=True, step=dict(step),
                              demand=contender['auxiliary_demand_mL_min']))
        else:
            step = {a: step[a]/2 for a in AXES}
            trace.append(dict(iteration=it, moved=False, step=dict(step)))
    return best, trace


def main(kappa=.5, budget=48, seed=8101, iterations=4, workers=8):
    start = time.perf_counter()
    history = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        points = sobol(budget, seed)
        rows = list(pool.map(evaluate, [(p, kappa, 'sobol') for p in points]))
        history.extend(rows)
        seed_point = best_of(rows)
        refined, trace = (None, [])
        if seed_point:
            centre = dict(scale=seed_point['scale'], period_s=seed_point['period_s'],
                          coronary_conductance=seed_point['coronary_conductance'],
                          auxiliary_demand_mL_min=seed_point['auxiliary_demand_mL_min'])
            refined, trace = compass(centre, kappa, pool, iterations, dict(STEP0), history)
    overall = best_of(history)
    result = dict(
        scope='joint search over cavity scale, activation period and coronary conductance',
        protocol_version=PROTOCOL['version'], horizon_s=PROTOCOL['evaluation_policy']['horizon_s'],
        compression_kappa=kappa, box={k: list(v) for k, v in BOX.items()},
        sobol_budget=budget, seed=seed, compass_iterations=iterations,
        evaluations=len(history), feasible=sum(1 for r in history if r.get('feasible')),
        failures=[r for r in history if r.get('outcome') == 'evaluation_failed'],
        sobol_best=seed_point, compass_best=refined, compass_trace=trace, best=overall,
        rows=history, wall_s=time.perf_counter()-start,
        claim='best feasible design found within the declared box and budget; not a global optimum')
    dest = ROOT / OUT; dest.mkdir(parents=True, exist_ok=True)
    (dest / f'joint-kappa{kappa:g}.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('joint search: %d evaluations, %d feasible, %.1f min' % (
        result['evaluations'], result['feasible'], result['wall_s']/60), flush=True)
    for label, row in (('sobol best', seed_point), ('compass best', refined), ('overall', overall)):
        if row:
            print('  %-13s scale %.3f, %5.1f/min, conductance %.3f -> %.3f mL O2/min  MAP %.2f' % (
                label, row['scale'], row['rate_bpm'], row['coronary_conductance'],
                row['auxiliary_demand_mL_min'], row['MAP_mmHg']), flush=True)
        else:
            print('  %-13s none' % label, flush=True)
    for t in trace:
        print('  compass %d: %s step %s' % (t['iteration'], 'moved' if t['moved'] else 'shrank',
              {k: round(v, 4) for k, v in t['step'].items()}), flush=True)
    return result


if __name__ == '__main__':
    main(kappa=float(sys.argv[1]) if len(sys.argv) > 1 else .5,
         budget=int(sys.argv[2]) if len(sys.argv) > 2 else 48,
         iterations=int(sys.argv[3]) if len(sys.argv) > 3 else 4)
