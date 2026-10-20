"""T8.2 core matrix: one design held fixed across the declared operating cases.

The design is never re-tuned per case. Re-optimizing for each adverse
circumstance would answer a different question and hide fragility, so the
matrix measures how one configuration behaves when circumstances change.

Cases come from the frozen protocol's operating_cases. Hydraulic flags
(activation policy, suppression override, native restart and reopening, sensor
loss, reduced filling) are applied in the circulation; metabolic flags (apnea,
demand restoration) are applied in the replay. A case whose hydraulics exit the
model domain is recorded with its cause and is not counted as a result.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json, sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import run_m8_hydraulics as hyd
import run_m8_metabolism as met
from cardiac_redundancy.design import DesignPoint, ScenarioPoint
from cardiac_redundancy.duration import supported_duration

OUT = 'reports/m8-matrix'
PROTOCOL = json.loads((ROOT / 'configs/m8-protocol.json').read_text(encoding='utf-8'))
#: The T8.0 joint optimum at the settled 40 s horizon under protocol v4.
OPTIMUM = dict(left_scale=.496, right_scale=.496, septal_scale=.496, pericardial_scale=.496,
               period_s=1., coronary_conductance=1.203)
HYDRAULIC_FLAGS = ('preservation_override', 'deplete_PU', 'restart', 'reopen', 'sensor_loss',
                   'inlet_multiplier', 'rearrest', 'withdraw', 'withdrawal_s')
METABOLIC_FLAGS = ('apnea', 'demand_restore')


def one(job):
    case_name, spec, kappa = job
    design = DesignPoint(**OPTIMUM)
    policy = spec['policy']
    case = {k: v for k, v in spec.items() if k in HYDRAULIC_FLAGS}
    name = f'{case_name}-kappa{kappa:g}'
    folder = ROOT / OUT / name
    if (folder / 'hydraulics.json').is_file():
        record = json.loads((folder / 'hydraulics.json').read_text(encoding='utf-8'))
    else:
        record = hyd.evaluate(design, name, OUT, policy=policy, kappa=kappa,
                              release_s=spec.get('release_s', .5), case=case)
    row = dict(case=case_name, meaning=spec['meaning'], policy=policy,
               compression_kappa=kappa, numerical_status=record['numerical_status'],
               pressure_status=record['pressure']['status'],
               status_by_criterion=record['pressure']['status_by_criterion'],
               release_s=record['runs'][-1]['release_s'],
               MAP_mmHg=record['runs'][-1]['final_MAP_mmHg'],
               violating_sustained_windows=record['pressure']['violating_sustained_windows'],
               sustained_windows=record['pressure']['sustained_windows'])
    if record['numerical_status'] != 'valid':
        row.update(outcome='model_domain_limit', note='not counted as a result')
        return row
    coupled_name = name
    cached = ROOT / OUT / coupled_name / 'm8-coupled.json'
    if cached.is_file():
        coupled = json.loads(cached.read_text(encoding='utf-8'))
    else:
        extra = {k: v for k, v in spec.items() if k in METABOLIC_FLAGS}
        coupled = met.evaluate(design, ScenarioPoint(), coupled_name, Path(OUT) / name,
                               out=OUT, extra_case=extra)
    if coupled.get('status') != 'verified':
        row.update(outcome=coupled.get('status'), note='coupled stage unusable')
        return row
    total = sum(coupled['groups'][k]['requested_mL'] for k in coupled['groups'])*60/record['duration_s']
    duration = supported_duration(total)
    row.update(outcome='evaluated',
               groups={k: dict(sustained_status=g['sustained_status'],
                               unmet_sustained_mL=g['unmet_sustained_mL'])
                       for k, g in coupled['groups'].items()},
               auxiliary_demand_mL_min=sum(coupled['auxiliary_demand_mL_min'].values()),
               total_demand_mL_O2_min=total,
               supported_duration_min=duration['minutes'],
               duration_limited_by=duration['limited_by'],
               all_groups_satisfied=all(g['sustained_status'] == 'satisfied'
                                        for g in coupled['groups'].values()))
    return row


def main(kappa=.5, workers=5):
    jobs = [(name, spec, kappa) for name, spec in PROTOCOL['operating_cases'].items()]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(one, jobs))
    result = dict(scope='T8.2 core matrix: the T8.0 optimum held fixed across declared '
                        'operating cases, reference physiology',
                  design=OPTIMUM, compression_kappa=kappa,
                  protocol_version=PROTOCOL['version'], rows=rows,
                  interpretation='fractions of tested cases, never probabilities; a pass means '
                                 'no modeled violation observed through the horizon')
    dest = ROOT / OUT; dest.mkdir(parents=True, exist_ok=True)
    (dest / f'matrix-kappa{kappa:g}.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('%-24s %-10s %-11s %8s %9s %s' % ('case', 'pressure', 'organs', 'MAP', 'minutes', 'release'), flush=True)
    for r in rows:
        if r.get('outcome') != 'evaluated':
            print('%-24s %-10s %s' % (r['case'], r.get('pressure_status', '-'), r.get('outcome')), flush=True)
            continue
        bad = [k for k, g in r['groups'].items() if g['sustained_status'] != 'satisfied']
        print('%-24s %-10s %-11s %8.2f %9s %s' % (
            r['case'], r['pressure_status'], 'satisfied' if not bad else ','.join(bad),
            r['MAP_mmHg'],
            'none' if r['supported_duration_min'] is None else '%.1f' % r['supported_duration_min'],
            '-' if r['release_s'] is None else '%.2f s' % r['release_s']), flush=True)
    return result


if __name__ == '__main__':
    main(float(sys.argv[1]) if len(sys.argv) > 1 else .5)
