"""T8.4 for Arm B: stress the fat-capable designs, including their own policy.

Arm A's five stresses apply unchanged, and Arm B adds one it alone can have.
The adipose policy is a capability claim about the support circulation, not a
property of the subject, so a design selected while the depot is perfused must
be asked what happens when that perfusion is not achieved. If feasibility
depends on it, the claim rests on the policy and must say so.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json, sys, time
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import optimize_m8b_joint as jb
from cardiac_redundancy.adipose_policy import AdiposePolicy

OUT = 'reports/m8b-disconfirm'
PROTOCOL = json.loads((ROOT / 'configs/m8-protocol.json').read_text(encoding='utf-8'))
STRESSES = (('reference', dict(kappa=.5, perfused=True)),
            ('kappa0', dict(kappa=0., perfused=True)),
            ('kappa1', dict(kappa=1., perfused=True)),
            ('depot_suppressed', dict(kappa=.5, perfused=False)))


def leading(limit=3):
    path = ROOT / 'reports/m8b-joint/joint-perfused-kappa0.5.json'
    if not path.is_file():
        raise SystemExit('Arm B joint search has not finished')
    rows = [r for r in json.loads(path.read_text(encoding='utf-8'))['rows'] if r.get('feasible')]
    rows.sort(key=lambda r: r['auxiliary_demand_mL_min'])
    seen, out = set(), []
    for r in rows:
        key = (round(r['scale'], 3), round(r['period_s'], 3), round(r['coronary_conductance'], 3))
        if key not in seen:
            seen.add(key); out.append(r)
        if len(out) >= limit:
            break
    return out


def one(job):
    row, label, spec = job
    record = dict(stress=label, scale=row['scale'], period_s=row['period_s'],
                  coronary_conductance=row['coronary_conductance'],
                  reference_demand_mL_O2_min=row['auxiliary_demand_mL_min'], **spec)
    try:
        status = jb.pair_status(row['scale'], row['period_s'], row['coronary_conductance'],
                                spec['kappa'], AdiposePolicy(perfused=spec['perfused']))
    except Exception as exc:
        record.update(feasible=False, outcome='evaluation_failed', cause=repr(exc)[:160])
        return record
    record.update({k: status.get(k) for k in
                   ('feasible', 'outcome', 'pressure', 'self_supply', 'MAP_mmHg',
                    'auxiliary_demand_mL_min', 'adipose_flow_mL_min',
                    'fat_supply_mL_O2_eq_min', 'covered_fraction', 'sustained_unmet_mL')})
    return record


def main(limit=3, workers=6):
    start = time.perf_counter()
    designs = leading(limit)
    jobs = [(row, label, spec) for row in designs for label, spec in STRESSES]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(one, jobs))
    survival = {label: sum(1 for r in rows if r['stress'] == label and r.get('feasible'))
                for label, _ in STRESSES}
    result = dict(arm='B', scope='T8.4 disconfirmation of the leading Arm B designs',
                  protocol_version=PROTOCOL['version'], tested=len(designs),
                  stresses=[label for label, _ in STRESSES], survival=survival, rows=rows,
                  wall_s=time.perf_counter()-start,
                  interpretation='fractions of tested cases, never probabilities')
    dest = ROOT / OUT; dest.mkdir(parents=True, exist_ok=True)
    (dest / 'disconfirmation.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('Arm B disconfirmation: %d designs x %d stresses, %.1f min' % (
        len(designs), len(STRESSES), result['wall_s']/60), flush=True)
    for label, _ in STRESSES:
        print('  %-18s %d of %d survive' % (label, survival[label], len(designs)), flush=True)
    for r in rows:
        if not r.get('feasible'):
            print('    scale %.3f %-18s %s %s' % (r['scale'], r['stress'], r.get('outcome'),
                  '' if r.get('self_supply') is None else 'self_supply=%s' % r['self_supply']),
                  flush=True)
    return result


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
