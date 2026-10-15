"""Static M2 gate checks and materialized hydraulic fixture; no ODE integration."""
from pathlib import Path
from copy import deepcopy
import hashlib
import itertools
import json
import math
import sys

ROOT = Path(__file__).resolve().parents[1]
FILES = ['configs/baseline.yaml', 'configs/project-extension.yaml',
         'configs/auxiliary-heart-design.yaml', 'configs/m2-frozen-overlay.yaml',
         'docs/model-specification.md', 'docs/m2-frozen-specification.md',
         'docs/verification-plan.md', 'docs/m2-circuits.md',
         'docs/m2-pressure-requirements.md']


def read(path):
    return json.loads((ROOT/path).read_text(encoding='utf-8'))


def waterfall(pa, pv, pc, r):
    if r <= 0:
        raise ValueError('Resistance must be positive')
    return (max(pa-max(pv, pc), 0)-max(pv-max(pa, pc), 0))/r


def audit():
    base = read('configs/baseline.yaml')
    parent = read('configs/project-extension.yaml')
    overlay = read('configs/m2-frozen-overlay.yaml')
    design = read('configs/auxiliary-heart-design.yaml')
    edges = deepcopy(parent['branches'])
    brain = next(e for e in edges if e['name']=='brain')
    pc = overlay['brain']['ICP_mmHg'] + overlay['brain']['tone_offset_mmHg']
    rbrain = (100-max(2, pc))/brain['target_flow_mL_s']
    brain['R_mmHg_s_mL'] = rbrain
    brain['law'] = 'signed_waterfall; Pc=ICP+tone_offset'
    assert math.isclose(waterfall(100, 2, pc, rbrain), brain['target_flow_mL_s'])
    nchecks = 0
    for pa, pv, pc_test in itertools.product([-10, 0, 10, 20, 80, 100], repeat=3):
        q = waterfall(pa, pv, pc_test, 2)
        assert math.isclose(q, -waterfall(pv, pa, pc_test, 2), abs_tol=1e-12)
        assert (pa-pv)*q >= -1e-12
        if min(pa, pv) >= pc_test:
            assert math.isclose(q, (pa-pv)/2, abs_tol=1e-12)
        if max(pa, pv) <= pc_test:
            assert q == 0
        nchecks += 1
    native = next(e for e in edges if e['name']=='heart')
    kidney = next(e for e in edges if e['name']=='kidneys')
    edges = [e for e in edges if e['name'] not in ['heart', 'kidneys']]
    for name, w in zip(overlay['kidneys']['replace_pair_with'], overlay['kidneys']['baseline_weights']):
        e = deepcopy(kidney)
        e.update(name=name, R_mmHg_s_mL=kidney['R_mmHg_s_mL']/w,
                 baseline_fraction=kidney['baseline_fraction']*w,
                 target_flow_mL_s=kidney['target_flow_mL_s']*w)
        edges.append(e)
    assert math.isclose(sum(1/e['R_mmHg_s_mL'] for e in edges if e['name'].startswith('kidney_')),
                        1/kidney['R_mmHg_s_mL'])
    alpha = overlay['auxiliary_fixture']['alpha']
    weights = [overlay['coronary']['left_conductance_fraction'],
               overlay['coronary']['right_conductance_fraction']]
    assert math.isclose(sum(weights), 1)
    for side, w in zip(['L', 'R'], weights):
        for owner, scale in [('native', 1), ('auxiliary', alpha)]:
            edges.append({'name': f'{owner}_coronary_{side}', 'from':'AO', 'to':'VC',
                          'R_mmHg_s_mL':native['R_mmHg_s_mL']/w/scale,
                          'law':'signed_waterfall_compression_sensitivity',
                          'classification':'A; benchmark conductance, not fixed achieved flow'})
    assert math.isclose(sum(1/e['R_mmHg_s_mL'] for e in edges if e['name'].startswith('native_coronary')),
                        1/native['R_mmHg_s_mL'])
    volumes = deepcopy(parent['initial_physical_volumes_mL'])
    for chamber, donor in overlay['auxiliary_fixture']['donors'].items():
        v = alpha*base['initial_volumes_mL'][chamber[1:]]
        volumes[chamber] = v
        volumes[donor] -= v
        assert volumes[donor] > parent['unstressed_offsets_mL'][donor]
    assert math.isclose(sum(volumes.values()), parent['total_blood_volume_mL'])
    assert all(v > 0 for v in volumes.values())
    assert math.isclose(overlay['auxiliary_fixture']['driver_dimensionless_width'],
                        base['driver']['B_per_s2']*overlay['auxiliary_fixture']['period_s']**2)
    graph = [(e['from'], e['to']) for e in edges]
    graph += [('PU','LV'), ('LV','AO'), ('VC','RV'), ('RV','PA'), ('PA','PU'),
              ('PU','aLV'), ('aLV','AO'), ('VC','aRV'), ('aRV','PA')]
    nodes = set(sum(([a,b] for a,b in graph), []))
    net = dict.fromkeys(nodes, 0.)
    for i, (a,b) in enumerate(graph):
        q = (-1)**i*(i+1)/3
        net[a] -= q
        net[b] += q
    assert abs(sum(net.values())) < 1e-10
    for origin in nodes:
        seen={origin}
        while True:
            new=seen|{b for a,b in graph if a in seen}
            if new==seen: break
            seen=new
        assert seen==nodes  # Post-trigger all-open topology, not instantaneous valve states.
    pressure=design['pressure']
    assert pressure['proven_minimum_SBP_mmHg'] is None
    assert pressure['proven_minimum_DBP_mmHg'] is None
    sbp, dbp=pressure['initial_waveform_guess_SBP_DBP_mmHg']
    assert dbp <= pressure['mean_AO_nominal_mmHg'] <= sbp
    fixture={'status':'materialized pre-warmup biological hydraulic fixture; not simulated',
             'initial_physical_volumes_mL':volumes, 'regional_edges':edges,
             'overlay':'configs/m2-frozen-overlay.yaml',
             'note':'Auxiliary coronary beds included; abstract fixtures omit auxiliary chambers and their coronary beds. Controls apply at trigger.'}
    (ROOT/'data/m2-hydraulic-fixture.json').write_text(json.dumps(fixture,indent=2)+'\n',encoding='utf-8')
    return {'status':'passed static freeze checks; no ODE or physiology validation',
            'waterfall_pressure_combinations_checked':nchecks,
            'brain_R_mmHg_s_mL':rbrain, 'brain_critical_pressure_A_mmHg':pc,
            'total_physical_blood_mL':sum(volumes.values()),
            'auxiliary_initial_volumes_mL':{k:volumes[k] for k in ['aLV','aRV']},
            'expanded_regional_edge_count':len(edges),
            'checks':['waterfall antisymmetry, resistor limit, closure and dissipativity',
                      'cerebral anchor','kidney/coronary conductance preservation',
                      'auxiliary blood allocation','normalized driver width',
                      'post-trigger graph paths and incidence','pressure policy'],
            'not_done':['ODE','periodic warm-up','organ adequacy','recovery','pump sizing']}


if __name__=='__main__':
    result=audit()
    manifest_path=ROOT/'data/m2-freeze-manifest.json'
    digests={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES}
    if '--freeze' in sys.argv:
        manifest_path.write_text(json.dumps({'version':'M2 v1.0','date':'2026-10-03','sha256':digests},indent=2)+'\n',encoding='utf-8')
    else:
        expected=json.loads(manifest_path.read_text(encoding='utf-8'))['sha256']
        assert digests==expected, 'Frozen specification changed: review and explicitly refreeze'
    result['snapshot_files_verified']=len(digests)
    (ROOT/'data/m2-freeze-audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
