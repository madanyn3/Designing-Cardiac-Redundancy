"""Arm B joint optimization over cavity scale, rate and coronary conductance.

The same search Arm A received, run on the circulation that represents visceral
adipose and the metabolism that may oxidise fat. The objective and constraints
are unchanged:

    minimize    auxiliary oxygen demand (mL O2/min)
    over        cavity scale, activation period, coronary conductance
    subject to  the declared pressure screen and auxiliary self-supply

The adipose policy is **not** searched. It is a categorical support-policy
choice over existing vasculature, so the search runs once per policy and the
two are compared, rather than letting an optimizer pick the one that flatters
the result.

Evaluation is staged as in Arm A: a design failing the pressure screen costs
only the hydraulic solve, because self-supply is never evaluated on a design
that has already failed. The cache key carries the policy as well as the
design, model form and horizon, since a preserved and a perfused depot are
different trajectories for the same design.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json, sys, time
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import optimize_m8_design as opt
import optimize_m8_joint as joint
import run_m8_hydraulics as hyd
import run_m8b_metabolism as m8b
from cardiac_redundancy.adipose_policy import AdiposePolicy
from cardiac_redundancy.design import ScenarioPoint

OUT = 'reports/m8b-joint'
PROTOCOL = json.loads((ROOT / 'configs/m8-protocol.json').read_text(encoding='utf-8'))
AXES = joint.AXES
BOX = joint.BOX
STEP0 = joint.STEP0


def name_of(design, kappa, policy, duration_s):
    """Identity including the policy: two policies are two trajectories."""
    base = opt.name_of(design, kappa, .5, duration_s)
    return base+('-perfused' if policy.perfused else '-preserved')


def pair_status(scale, period, coronary, kappa, policy):
    """Both Arm B constraints at one design, or why it is not a result."""
    design = opt.design(round(scale, 3), round(period, 3), round(coronary, 3))
    horizon = PROTOCOL['evaluation_policy']['horizon_s']
    name = name_of(design, kappa, policy, horizon)
    folder = ROOT / OUT / name
    row = dict(scale=design.left_scale, period_s=design.period_s, rate_bpm=design.rate_bpm,
               coronary_conductance=design.coronary_conductance, perfused=policy.perfused)
    if (folder / 'hydraulics.json').is_file():
        record = json.loads((folder / 'hydraulics.json').read_text(encoding='utf-8'))
    else:
        record = hyd.evaluate(design, name, OUT, duration_s=horizon, policy='fixed',
                              kappa=kappa, release_s=.5, network=hyd.ARM_B_NETWORK,
                              additional_mandatory=policy.additional_mandatory)
    row.update(pressure=record['pressure']['status'], numerical=record['numerical_status'],
               MAP_mmHg=record['runs'][-1]['final_MAP_mmHg'])
    if record['pressure']['status'] != 'satisfied' or record['numerical_status'] != 'valid':
        row.update(feasible=False, outcome='pressure_violated')
        return row
    cached = folder / 'm8b-coupled.json'
    if cached.is_file():
        coupled = json.loads(cached.read_text(encoding='utf-8'))
    else:
        coupled = m8b.evaluate(design, ScenarioPoint(), policy, name, Path(OUT) / name, out=OUT)
    if coupled.get('status') != 'verified':
        row.update(feasible=False, outcome=coupled.get('status'))
        return row
    groups = coupled['groups']
    row.update(outcome='evaluated',
               feasible=all(g['sustained_status'] == 'satisfied' for g in groups.values()),
               self_supply=groups['pump']['sustained_status'],
               sustained_unmet_mL=groups['pump']['unmet_sustained_mL'],
               auxiliary_demand_mL_min=sum(coupled['auxiliary_demand_mL_min'].values()),
               adipose_flow_mL_min=coupled['adipose_flow_mL_min'],
               fat_supply_mL_O2_eq_min=coupled['fat_supply']['supply_mL_O2_eq_min'],
               fat_binding=coupled['fat_supply']['binding'],
               covered_fraction=coupled['glucose_shares']['covered_fraction'])
    return row


def evaluate(args):
    point, kappa, perfused, tag = args
    p = joint.clip(point)
    try:
        row = pair_status(p['scale'], p['period_s'], p['coronary_conductance'],
                          kappa, AdiposePolicy(perfused=perfused))
    except Exception as exc:
        return dict(tag=tag, feasible=False, outcome='evaluation_failed',
                    cause=repr(exc)[:160], **p)
    row['tag'] = tag
    return row


def main(perfused=True, kappa=.5, budget=32, seed=8201, iterations=3, workers=8):
    start = time.perf_counter()
    history = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        points = joint.sobol(budget, seed)
        rows = list(pool.map(evaluate, [(p, kappa, perfused, 'sobol') for p in points]))
        history.extend(rows)
        seed_point = joint.best_of(rows)
        refined, trace = None, []
        if seed_point:
            centre = {a: seed_point[a] for a in AXES}
            centre['auxiliary_demand_mL_min'] = seed_point['auxiliary_demand_mL_min']
            step = dict(STEP0)
            best = None
            for it in range(iterations):
                candidates = []
                for axis in AXES:
                    for sign in (+1, -1):
                        q = dict(centre); q[axis] = centre[axis]+sign*step[axis]
                        candidates.append((q, kappa, perfused, 'compass%d:%s%+d' % (it, axis, sign)))
                got = list(pool.map(evaluate, candidates))
                history.extend(got)
                contender = joint.best_of(got)
                reference = (best or seed_point)['auxiliary_demand_mL_min']
                if contender and contender['auxiliary_demand_mL_min'] < reference:
                    best = contender
                    centre = {a: contender[a] for a in AXES}
                    trace.append(dict(iteration=it, moved=True,
                                      demand=contender['auxiliary_demand_mL_min']))
                else:
                    step = {a: step[a]/2 for a in AXES}
                    trace.append(dict(iteration=it, moved=False))
            refined = best
    overall = joint.best_of(history)
    result = dict(arm='B', scope='joint search on the Arm B circulation and fat-capable metabolism',
                  adipose_policy=AdiposePolicy(perfused=perfused).describe(),
                  protocol_version=PROTOCOL['version'],
                  horizon_s=PROTOCOL['evaluation_policy']['horizon_s'], compression_kappa=kappa,
                  box={k: list(v) for k, v in BOX.items()}, sobol_budget=budget, seed=seed,
                  evaluations=len(history),
                  feasible=sum(1 for r in history if r.get('feasible')),
                  failures=[r for r in history if r.get('outcome') == 'evaluation_failed'],
                  sobol_best=seed_point, compass_best=refined, compass_trace=trace,
                  best=overall, rows=history, wall_s=time.perf_counter()-start,
                  claim='best feasible design found within the declared box and budget')
    dest = ROOT / OUT; dest.mkdir(parents=True, exist_ok=True)
    tag = 'perfused' if perfused else 'preserved'
    (dest / ('joint-%s-kappa%g.json' % (tag, kappa))).write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('Arm B joint (%s): %d evaluations, %d feasible, %.1f min' % (
        tag, result['evaluations'], result['feasible'], result['wall_s']/60), flush=True)
    if overall:
        print('  best: scale %.3f, %.1f/min, conductance %.3f -> %.3f mL O2/min  MAP %.2f' % (
            overall['scale'], overall['rate_bpm'], overall['coronary_conductance'],
            overall['auxiliary_demand_mL_min'], overall['MAP_mmHg']), flush=True)
        print('       depot %.1f mL/min, fat %.1f mL O2-eq/min (%s binds), covers %.1f%%' % (
            overall['adipose_flow_mL_min'], overall['fat_supply_mL_O2_eq_min'],
            overall['fat_binding'], 100*overall['covered_fraction']), flush=True)
    else:
        print('  no feasible design found', flush=True)
    return result


if __name__ == '__main__':
    main(perfused=(len(sys.argv) < 2 or sys.argv[1].lower() in ('1', 'true', 'perfused')),
         budget=int(sys.argv[2]) if len(sys.argv) > 2 else 32,
         iterations=int(sys.argv[3]) if len(sys.argv) > 3 else 3)
