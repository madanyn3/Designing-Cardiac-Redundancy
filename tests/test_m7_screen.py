"""Check design isolation before interpreting endpoint effects."""
import importlib.util
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import run_m7_metabolic_screen as screen

def test_oat_design_changes_exactly_one_factor():
    cases=screen.cases();assert len(cases)==13
    nominal=cases[0]['knobs']
    assert len({c['name'] for c in cases})==13
    for c in cases[1:]:
        assert sum(c['knobs'][k]!=v for k,v in nominal.items())==1

def test_output_isolation_and_fixed_hydraulics():
    assert screen.OUT==ROOT/'reports/m7-metabolic-screen'
    segments,paths=screen.hydraulics({'name':'a_different_output_case'})
    assert paths==['reports/m6-coupled/observed_arrest/hydraulic-refined.npz']
    assert segments[0]['time'][0]==0
    assert abs(segments[0]['time'][-1]-10)<1e-10


def test_completed_nominal_reproduces_frozen_m6():
    import json
    current=screen.OUT/'nominal/result.json'
    assert current.exists(), 'Run the paired screen before this result regression'
    new=json.loads(current.read_text(encoding='utf-8'))
    old=screen.read('reports/m6-coupled/observed_arrest/result.json')
    assert new['status']=='verified'
    assert abs(new['records'][-1]['unmet_oxygen_mL']-old['records'][-1]['unmet_oxygen_mL'])<1e-7
    for name in ('brain','coronary_left','coronary_right','kidney_protected','hepatic_outflow','respiratory_muscle','lung_supported'):
        assert new['profile'][name]==old['profile'][name]


def test_stock_endpoints_preserve_release_accounting():
    import json
    for name in ('glycogen_mmol_L_liver_low','glycogen_mmol_L_liver_high'):
        r=json.loads((screen.OUT/name/'result.json').read_text(encoding='utf-8'))
        assert r['status']=='verified'
        p=r['parameters'];initial=p['liver_volume_L']*p['glycogen_mmol_L_liver']*p['accessible_fraction']
        final=r['records'][-1]
        assert abs(initial-final['glycogen_remaining_mmol']-final['released_mmol'])<1e-8
