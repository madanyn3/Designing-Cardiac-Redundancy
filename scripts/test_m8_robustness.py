"""T8.4 scenario robustness for designs that passed at reference physiology.

No scenario factor enters the hydraulic model, so a change of physiology
replays over the saved trajectory and costs one coupled solve rather than a
full evaluation. The design is held fixed: re-optimizing it per scenario would
answer a different question and hide fragility.

Protocol v2 judges a candidate on its worst declared scenario, so a design is
robust only if every scenario here satisfies the screens.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json, sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import run_m8_metabolism as met
from cardiac_redundancy.design import DesignPoint, ScenarioPoint

OUT = 'reports/m8-optimize'
PROTOCOL = json.loads((ROOT / 'configs/m8-protocol.json').read_text(encoding='utf-8'))


def one(job):
    folder, scenario_name = job
    source = ROOT / OUT / folder
    record = json.loads((source / 'hydraulics.json').read_text(encoding='utf-8'))
    design = DesignPoint(**record['design'])
    scenario = ScenarioPoint(**PROTOCOL['scenarios'][scenario_name])
    name = f'{folder}-{scenario_name}'
    cached = ROOT / OUT / name / 'm8-coupled.json'
    if cached.is_file():
        result = json.loads(cached.read_text(encoding='utf-8'))
    else:
        try:
            result = met.evaluate(design, scenario, name, Path(OUT) / folder, out=OUT)
        except Exception as exc:
            return dict(folder=folder, scenario=scenario_name, status='error', error=repr(exc)[:200])
    if result.get('status') != 'verified':
        return dict(folder=folder, scenario=scenario_name, status=result.get('status'))
    groups = {role: dict(sustained_status=g['sustained_status'],
                         unmet_sustained_mL=g['unmet_sustained_mL']) for role, g in result['groups'].items()}
    return dict(folder=folder, scenario=scenario_name, status='verified',
                scale=design.left_scale, rate_bpm=design.rate_bpm,
                coronary_conductance=design.coronary_conductance,
                compression_kappa=record.get('compression_kappa'),
                auxiliary_demand_mL_min=sum(result['auxiliary_demand_mL_min'].values()),
                coupled_mass_ratio=result['coupled_mass_ratio'], groups=groups,
                feasible=all(g['sustained_status'] == 'satisfied' for g in groups.values()))


def main(scenarios=('adverse', 'favourable'), workers=4):
    phase = json.loads((ROOT / OUT / 'phase-b-kappa0.5.json').read_text(encoding='utf-8'))
    folders = []
    for pair in phase['feasible_pairs']:
        d = DesignPoint(left_scale=pair['scale'], right_scale=pair['scale'],
                        septal_scale=pair['scale'], pericardial_scale=pair['scale'],
                        period_s=pair['period_s'], coronary_conductance=pair['coronary_conductance'])
        folders.append('scale%g-period%g-coron%g-kappa%g'
                       % (d.left_scale, d.period_s, d.coronary_conductance, pair['compression_kappa']))
    jobs = [(f, s) for f in folders for s in scenarios]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(one, jobs))
    survivors = [r for r in rows if r.get('feasible')]
    result = dict(scope='scenario robustness of reference-feasible designs; design held fixed',
                  scenarios=list(scenarios), rows=rows,
                  survivors=[dict(scale=r['scale'], rate_bpm=r['rate_bpm'], scenario=r['scenario'],
                                  auxiliary_demand_mL_min=r['auxiliary_demand_mL_min']) for r in survivors],
                  interpretation='fractions of tested scenarios, never probabilities')
    (ROOT / OUT / 'robustness.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    for r in rows:
        if r.get('status') != 'verified':
            print('  %-46s %-11s %s' % (r['folder'], r['scenario'], r.get('status')), flush=True); continue
        bad = [k for k, g in r['groups'].items() if g['sustained_status'] != 'satisfied']
        print('  scale %.3f @ %3.0f/min  %-11s demand %7.3f  %s' % (
            r['scale'], r['rate_bpm'], r['scenario'], r['auxiliary_demand_mL_min'],
            'FEASIBLE' if r['feasible'] else 'violated: '+','.join(bad)), flush=True)
    print('survivors: %d of %d' % (len(survivors), len(rows)), flush=True)
    return result


if __name__ == '__main__':
    main(tuple(sys.argv[1:]) or ('adverse', 'favourable'))
