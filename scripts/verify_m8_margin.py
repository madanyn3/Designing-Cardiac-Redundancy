"""Verify the margin-constrained candidates on auxiliary self-supply.

The margin frontier was computed from pressure-feasible designs only, so each
candidate still has to pass the coupled stage. Protocol v5 requires a declared
minimum pressure margin, and this reports the cheapest design that holds it
*and* feeds its own myocardium.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import glob, json, sys, time
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import optimize_m8_design as opt
import run_m8_metabolism as m8
from cardiac_redundancy.design import DesignPoint, ScenarioPoint
from cardiac_redundancy.feasibility import pressure_margin

OUT = 'reports/m8-margin'
PROTOCOL = json.loads((ROOT / 'configs/m8-protocol.json').read_text(encoding='utf-8'))
FLOOR = PROTOCOL['minimum_pressure_margin_mmHg']


def candidates(limit=14):
    """Cheapest pressure-feasible designs holding the declared margin."""
    rows = []
    for path in glob.glob(str(ROOT / 'reports/m8-optimize/*h40*/hydraulics.json')):
        d = json.loads(Path(path).read_text(encoding='utf-8'))
        if d.get('compression_kappa') != .5 or d['pressure']['status'] != 'satisfied':
            continue
        margin = pressure_margin(d['pressure'])
        if margin is None or margin < FLOOR:
            continue
        design = DesignPoint(**d['design'])
        try:
            _, sides = m8.auxiliary_demand(design, ScenarioPoint(), d)
        except Exception:
            continue
        rows.append(dict(scale=design.left_scale, period_s=design.period_s,
                         rate_bpm=design.rate_bpm,
                         coronary_conductance=design.coronary_conductance,
                         margin_mmHg=margin,
                         demand=sum(v['total_mL_min'] for v in sides.values())))
    rows.sort(key=lambda r: r['demand'])
    seen, out = set(), []
    for r in rows:
        key = (round(r['scale'], 3), round(r['period_s'], 3), round(r['coronary_conductance'], 3))
        if key in seen:
            continue
        seen.add(key); out.append(r)
        if len(out) >= limit:
            break
    return out


def one(row):
    try:
        status = opt.pair_status(row['scale'], row['period_s'],
                                 row['coronary_conductance'], .5)
    except Exception as exc:
        return dict(row, feasible=False, outcome='evaluation_failed', cause=repr(exc)[:140])
    return dict(row, **{k: status.get(k) for k in
                        ('feasible', 'pressure', 'self_supply', 'MAP_mmHg',
                         'auxiliary_demand_mL_min', 'sustained_unmet_mL')})


def main(limit=14, workers=7):
    start = time.perf_counter()
    rows = candidates(limit)
    if not rows:
        raise SystemExit('no design holds the declared margin')
    with ProcessPoolExecutor(max_workers=workers) as pool:
        out = list(pool.map(one, rows))
    feasible = [r for r in out if r.get('feasible')]
    best = min(feasible, key=lambda r: r['auxiliary_demand_mL_min']) if feasible else None
    result = dict(scope='self-supply verification of the margin-constrained candidates',
                  minimum_margin_mmHg=FLOOR, protocol_version=PROTOCOL['version'],
                  tested=len(out), feasible=len(feasible), rows=out, best=best,
                  wall_s=time.perf_counter()-start)
    dest = ROOT / OUT; dest.mkdir(parents=True, exist_ok=True)
    (dest / 'margin-candidates.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('margin >= %.1f mmHg: %d candidates, %d feasible, %.1f min' % (
        FLOOR, len(out), len(feasible), result['wall_s']/60), flush=True)
    print('  %7s %7s %8s %8s %9s %s' % ('scale', 'rate', 'coron', 'margin', 'demand', 'self-supply'), flush=True)
    for r in sorted(out, key=lambda r: r['demand']):
        print('  %7.3f %7.1f %8.3f %+8.2f %9.3f %s' % (
            r['scale'], r['rate_bpm'], r['coronary_conductance'], r['margin_mmHg'],
            r['demand'], r.get('self_supply') or r.get('outcome')), flush=True)
    return result


if __name__ == '__main__':
    main()
