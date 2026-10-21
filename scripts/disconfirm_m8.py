"""T8.4 disconfirmation: look for where the result fails and where rankings invert.

A design chosen as cheapest under one set of assumptions is not thereby the
best design. This takes the leading feasible designs and re-evaluates them
under the declared stresses, then asks two questions the optimization cannot
ask of itself:

  * **Does feasibility survive?** Each design is tested against the adverse and
    favourable physiology triples and against the compression-kappa cases.
  * **Does the ranking survive?** If the order by objective inverts under a
    stress, the cheapest design at reference is not the one to report, and
    selecting on reference alone was the wrong criterion.

Scenario changes replay over the saved circulation, because no scenario factor
enters the hydraulic model. Kappa changes do not: they alter coronary flow, so
those cases are solved from the circulation up.

Boundaries found here are boundaries of the sampled domain. A design that fails
at kappa 1 is reported as failing at kappa 1, not as failing in general, and
nothing is claimed outside the declared ranges.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json, sys, time
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import optimize_m8_design as opt
import run_m8_metabolism as m8
from cardiac_redundancy.design import DesignPoint, ScenarioPoint

OUT = 'reports/m8-disconfirm'
PROTOCOL = json.loads((ROOT / 'configs/m8-protocol.json').read_text(encoding='utf-8'))
#: Stresses applied to every leading design. Scenario cases replay over the
#: reference circulation; kappa cases are solved from the circulation up.
STRESSES = (('reference', dict(scenario='reference', kappa=.5)),
            ('adverse', dict(scenario='adverse', kappa=.5)),
            ('favourable', dict(scenario='favourable', kappa=.5)),
            ('kappa0', dict(scenario='reference', kappa=0.)),
            ('kappa1', dict(scenario='reference', kappa=1.)))


def leading_designs(limit=4):
    """The best feasible designs found, pooled across the searches."""
    rows = []
    for path in ('reports/m8-joint/joint-kappa0.5.json',
                 'reports/m8-verification/space-filling-kappa0.5.json'):
        f = ROOT / path
        if f.is_file():
            rows += [r for r in json.loads(f.read_text(encoding='utf-8'))['rows']
                     if r.get('feasible')]
    rows.sort(key=lambda r: r['auxiliary_demand_mL_min'])
    seen, out = set(), []
    for r in rows:
        key = (round(r['scale'], 3), round(r['period_s'], 3), round(r['coronary_conductance'], 3))
        if key in seen:
            continue
        seen.add(key); out.append(r)
        if len(out) >= limit:
            break
    return out


def one(job):
    row, label, spec = job
    design = opt.design(row['scale'], row['period_s'], row['coronary_conductance'])
    scenario = ScenarioPoint(**PROTOCOL['scenarios'][spec['scenario']])
    kappa = spec['kappa']
    tag = 'scale%g-period%g-coron%g-%s' % (design.left_scale, design.period_s,
                                           design.coronary_conductance, label)
    record = dict(design=vars(design), stress=label, scenario=spec['scenario'], kappa=kappa,
                  reference_demand_mL_O2_min=row['auxiliary_demand_mL_min'])
    try:
        status = opt.pair_status(design.left_scale, design.period_s,
                                 design.coronary_conductance, kappa)
    except Exception as exc:
        record.update(outcome='evaluation_failed', cause=repr(exc)[:160], feasible=False)
        return record
    if not status.get('pressure') == 'satisfied':
        record.update(outcome='pressure_violated', feasible=False,
                      pressure=status.get('pressure'), MAP_mmHg=status.get('MAP_mmHg'))
        return record
    if spec['scenario'] == 'reference':
        record.update(outcome='evaluated', feasible=bool(status.get('feasible')),
                      pressure=status.get('pressure'), self_supply=status.get('self_supply'),
                      auxiliary_demand_mL_min=status.get('auxiliary_demand_mL_min'),
                      sustained_unmet_mL=status.get('sustained_unmet_mL'),
                      MAP_mmHg=status.get('MAP_mmHg'))
        return record
    # Scenario stresses replay over the circulation already solved above.
    folder = opt.name_of(design, kappa, .5, PROTOCOL['evaluation_policy']['horizon_s'])
    try:
        coupled = m8.evaluate(design, scenario, tag, Path(opt.OUT) / folder, out=OUT)
    except Exception as exc:
        record.update(outcome='coupled_failed', cause=repr(exc)[:160], feasible=False)
        return record
    if coupled.get('status') != 'verified':
        record.update(outcome=coupled.get('status'), feasible=False)
        return record
    groups = {k: dict(status=g['sustained_status'], unmet=g['unmet_sustained_mL'])
              for k, g in coupled['groups'].items()}
    record.update(outcome='evaluated',
                  feasible=all(g['status'] == 'satisfied' for g in groups.values()),
                  groups=groups, MAP_mmHg=status.get('MAP_mmHg'),
                  auxiliary_demand_mL_min=sum(coupled['auxiliary_demand_mL_min'].values()),
                  violated=[k for k, g in groups.items() if g['status'] != 'satisfied'])
    return record


def rankings(rows):
    """Order designs by objective within each stress, and find inversions."""
    order = {}
    for label, _ in STRESSES:
        here = [r for r in rows if r['stress'] == label and r.get('auxiliary_demand_mL_min')]
        here.sort(key=lambda r: r['auxiliary_demand_mL_min'])
        order[label] = [r['design']['left_scale'] for r in here]
    reference = order.get('reference', [])
    reversals = []
    for label, seq in order.items():
        if label == 'reference':
            continue
        common = [s for s in reference if s in seq]
        if [s for s in seq if s in common] != common:
            reversals.append(dict(stress=label, reference_order=common,
                                  stressed_order=[s for s in seq if s in common]))
    return order, reversals


def main(limit=4, workers=5):
    start = time.perf_counter()
    designs = leading_designs(limit)
    jobs = [(row, label, spec) for row in designs for label, spec in STRESSES]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(one, jobs))
    order, reversals = rankings(rows)
    survived = {label: [r for r in rows if r['stress'] == label and r.get('feasible')]
                for label, _ in STRESSES}
    result = dict(scope='T8.4 disconfirmation of the leading designs',
                  protocol_version=PROTOCOL['version'], designs=designs,
                  stresses=[label for label, _ in STRESSES], rows=rows,
                  survival={k: len(v) for k, v in survived.items()},
                  tested=len(designs), rankings=order, ranking_reversals=reversals,
                  wall_s=time.perf_counter()-start,
                  interpretation='fractions of tested cases, never probabilities; boundaries '
                                 'are boundaries of the sampled domain, not of the possible')
    dest = ROOT / OUT; dest.mkdir(parents=True, exist_ok=True)
    (dest / 'disconfirmation.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('disconfirmation: %d designs x %d stresses, %.1f min' % (
        len(designs), len(STRESSES), result['wall_s']/60), flush=True)
    for label, _ in STRESSES:
        print('  %-12s %d of %d designs survive' % (label, len(survived[label]), len(designs)), flush=True)
    for r in rows:
        if not r.get('feasible'):
            print('    scale %.3f %-11s %s %s' % (r['design']['left_scale'], r['stress'],
                  r.get('outcome'), ','.join(r.get('violated', []))), flush=True)
    if reversals:
        print('  RANKING REVERSALS:', flush=True)
        for rv in reversals:
            print('    under %s: %s -> %s' % (rv['stress'], rv['reference_order'],
                                              rv['stressed_order']), flush=True)
    else:
        print('  no ranking reversal among surviving designs', flush=True)
    return result


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
