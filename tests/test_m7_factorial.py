import itertools,json
from pathlib import Path
import numpy as np
import pytest
from cardiac_redundancy.factorial_analysis import decompose,interpolate
from cardiac_redundancy.sensitivity_circulation import SensitivityCirculation
from cardiac_redundancy.controlled_circulation import ControlledCirculation
ROOT=Path(__file__).resolve().parents[1]

def model(cls=SensitivityCirculation,**kwargs):
    return cls(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml',arterial_fraction=.2,inlet_resistance_fraction=.1,**kwargs)

def test_known_main_and_interaction_effects():
    x=np.array(list(itertools.product((-1,1),repeat=3)))
    y=4+2*x[:,0]-3*x[:,1]+5*x[:,0]*x[:,2]
    d=decompose(x,y,['a','b','c'])
    assert d['terms']['a']['effect']==4
    assert d['terms']['b']['effect']==-6
    assert d['terms']['a * c']['effect']==10
    assert d['reconstruction_error']<1e-12
    assert sum(t['variance_share'] for t in d['terms'].values())==pytest.approx(1)
    assert interpolate(d,dict(a=0,b=0,c=0),dict(a=(-1,1),b=(-1,1),c=(-1,1)))==4

def test_incomplete_or_duplicate_design_rejected():
    with pytest.raises(ValueError):decompose([[1],[1]],[2,3],['a'])

def test_constant_response_has_no_variance_ranking():
    assert decompose([[-1],[1]],[2,2],['a'])['total_shares']['a'] is None

def test_coronary_multiplier_is_local_and_conservative():
    base=model(ControlledCirculation);m=model(coronary_conductance=2)
    state=base.initial_state.copy();p=base.pressures(0,state)[0]
    a=base.flows(0,state,[True]*8,pressures=p);b=m.flows(0,state,[True]*8,pressures=p)
    for name in a:assert b[name]==pytest.approx(a[name]*(2 if name.startswith('aux_coronary_') else 1))
    delta=np.array([b[n]-a[n] for n in m.flow_names]);dv=m.incidence@delta
    assert abs(dv.sum())<1e-12
    assert np.count_nonzero(abs(dv)>1e-12)<=2

def test_identity_multiplier_reproduces_reference():
    a=model(ControlledCirculation);b=model()
    assert a.flows(0,a.initial_state,[True]*8)==b.flows(0,b.initial_state,[True]*8)

@pytest.mark.parametrize('value',[0,-1,float('nan'),float('inf'),True])
def test_invalid_coronary_multiplier(value):
    with pytest.raises(ValueError):model(coronary_conductance=value)

def test_protocol_complete_and_nonoverlapping_blocks():
    p=json.loads((ROOT/'configs/m7-factorial.json').read_text())
    assert len(p['cases'])==17
    for block,fs in p['factors'].items():
        rows=[c for c in p['cases'] if c['block']==block and c['signs'] is not None]
        assert set(tuple(c['signs']) for c in rows)==set(itertools.product((-1,1),repeat=3))
        for c in rows:
            for f,sign in zip(fs,c['signs']):assert c['knobs'][f['id']]==f['low' if sign==-1 else 'high']
