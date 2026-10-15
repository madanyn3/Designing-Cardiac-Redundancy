import json
from pathlib import Path
import numpy as np
import pytest
from cardiac_redundancy.selective import SelectiveModel
from cardiac_redundancy.whole_network import WholeNetworkModel,validate_roles
ROOT=Path(__file__).resolve().parents[1]
ARGS=(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml')
KW=dict(arterial_fraction=.2,inlet_resistance_fraction=.1)

def test_all_open_exact_regression():
    a=SelectiveModel(*ARGS,**KW);b=WholeNetworkModel(*ARGS,**KW)
    for active in (True,False):
        for t in (0,.2,.375):
            np.testing.assert_array_equal(a.rhs(t,a.initial_state,[True]*4,active=active),b.rhs(t,b.initial_state,[True]*4,active=active))

def test_independent_renal_control_and_native_coronary_closure():
    m=WholeNetworkModel(*ARGS,**KW);y=m.initial_state.copy();old=m.flows(.375,y,[True]*4)
    m.set_regional_controls(time=0,factors={'kidney_other':0,'coronary_left':.25})
    new=m.flows(.375,y,[True]*4)
    assert new['kidney_other']==0
    assert new['kidney_protected']==old['kidney_protected']
    assert new['coronary_left']==pytest.approx(.25*old['coronary_left'])
    assert new['brain']==old['brain']
    assert sum(m.rhs(.375,y,[True]*4)[:16])==pytest.approx(0,abs=1e-10)

def test_controlled_hepatic_balance_and_reversal():
    m=WholeNetworkModel(*ARGS,**KW);y=m.initial_state
    gut=['stomach_oesophagus','small_intestine','large_intestine','pancreas','spleen']
    m.set_regional_controls(time=0,factors=dict.fromkeys(gut,0))
    p,_=m.pressures(0,y);q=m.flows(0,y,[False]*4,pressures=p)
    assert all(q[n]==0 for n in gut)
    assert q['hepatic_outflow']==pytest.approx(q['liver_arterial'])
    p['AO']=1;p['VC']=10
    branches=m.base.hepatic_branches;f=m.regional_at(0)
    p['H']=sum(p[b.recipient if b.donor=='H' else b.donor]*f[b.name]/b.resistance for b in branches)/sum(f[b.name]/b.resistance for b in branches)
    q=m.flows(0,y,[False]*4,pressures=p)
    assert q['liver_arterial']<0 and q['hepatic_outflow']==pytest.approx(q['liver_arterial'])
    m.set_regional_controls(time=0,factors={'liver_arterial':0})
    assert m.flows(0,y,[False]*4)['hepatic_outflow']==pytest.approx(0,abs=1e-12)
    with pytest.raises(ValueError):m.set_regional_controls(time=0,factors={'hepatic_outflow':0})

def test_ramps_roles_and_unknown_controls():
    m=WholeNetworkModel(*ARGS,**KW)
    m.set_regional_controls(time=1,factors={'kidney_other':0},ramp_s=2)
    assert m.regional_at(2)['kidney_other']==.5
    for t in (1,2,3):assert sum(m.rhs(t,m.initial_state,[True]*4)[:16])==pytest.approx(0,abs=1e-10)
    with pytest.raises(ValueError):m.set_regional_controls(time=0,factors={'support_left':0})
    scenario=json.loads((ROOT/'configs/m3-whole-network-protocol.json').read_text())
    assert validate_roles(m,scenario)
    scenario['preservation_supply_edges'].remove('coronary_left')
    with pytest.raises(ValueError):validate_roles(m,scenario)
