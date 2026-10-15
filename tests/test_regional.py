from pathlib import Path
import numpy as np
import pytest

from cardiac_redundancy.regional import RegionalModel
from cardiac_redundancy.integration import integrate_source

ROOT=Path(__file__).resolve().parents[1]


@pytest.fixture
def model():
    return RegionalModel(ROOT/'configs/baseline.yaml',ROOT/'configs/project-extension.yaml',ROOT/'configs/m2-frozen-overlay.yaml')


@pytest.fixture
def amended():
    return RegionalModel(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml')


def test_evidence_allocation_and_named_complements(amended):
    test_regional_anchors_and_split_conductance(amended)
    arterial=[b for b in amended.project['branches'] if b['from']=='AO']
    assert sum(b['target_flow_mL_s'] for b in arterial)==pytest.approx(5600/60)
    assert 'gut' not in amended.branch_by_name
    assert sum(amended.project['remaining_flow_components_mL_min'].values())==pytest.approx(
        next(b['target_flow_mL_s']*60 for b in arterial if b['name']=='remaining'))
    assert sum(next(b['target_flow_mL_s']*60 for b in arterial if b['name']==n)
               for n in amended.project['portal_branches'])==pytest.approx(1064)


@pytest.mark.parametrize('reverse',[False,True])
def test_expanded_portal_balance_and_conservation(amended,reverse):
    y=amended.initial_state.copy()
    if reverse: y[2]=151.; y[3]+=132.3381
    q=amended.flows(0,y,[True]*4)
    assert sum(q[n] for n in amended.project['portal_branches'])+q['liver_arterial']==pytest.approx(q['hepatic_outflow'],abs=1e-11)
    assert amended.rhs(0,y,[True]*4)[:11].sum()==pytest.approx(0,abs=1e-11)
    if reverse: assert all(q[n]<0 for n in amended.project['portal_branches'])


def test_regional_anchors_and_split_conductance(model):
    p={'AO':100.,'VC':2.,'PU':4.,'PA':20.,'H':12.}
    p.update({n:5. for n in model.volume_names[6:]})
    flows=model.flows(0,model.initial_state,[False]*4,pressures=p)
    for original in model.project['branches']:
        name=original['name']
        if name=='heart': q=flows['coronary_left']+flows['coronary_right']
        elif name=='kidneys': q=flows['kidney_protected']+flows['kidney_other']
        else: q=flows[name]
        assert q==pytest.approx(original['target_flow_mL_s'])
    assert flows['kidney_protected']==pytest.approx(flows['kidney_other'])
    assert flows['coronary_left']/flows['coronary_right']==pytest.approx(7/3)


@pytest.mark.parametrize('reverse',[False,True])
def test_hepatic_junction_and_regional_conservation(model,reverse):
    y=model.initial_state.copy()
    if reverse:
        # Raise venous pressure above AO using a test-state redistribution.
        # Algebraic balance remains valid for reversed splanchnic flow.
        y[2]=151.; y[3]+=132.3381
    before=y.copy()
    p,_=model.pressures(0,y)
    q=model.flows(0,y,[True]*4,pressures=p)
    assert q['gut']+q['liver_arterial']==pytest.approx(q['hepatic_outflow'],abs=1e-12)
    if reverse: assert q['gut']<0 and q['hepatic_outflow']<0
    dy=model.rhs(0,y,[True]*4)
    assert dy[:11].sum()==pytest.approx(0,abs=1e-11)
    np.testing.assert_array_equal(y,before)
    np.testing.assert_array_equal(model.incidence.sum(axis=0),0)


def test_regional_cycle_conserves_all_pools_and_localizes_boundary(model):
    tr=integrate_source(model,.75)
    assert tr.termination is None
    assert tr.states.shape[0]==15
    assert np.max(abs(tr.states[:11].sum(axis=0)-5250))/5250<1e-8
    assert min(min(model.domain_margins(y).values()) for y in tr.states.T)>0
    assert {e['valve'] for e in tr.events}=={'MT','AV','TC','PV'}


def test_domain_crossing_is_reported_at_boundary_without_clipping(model):
    # A conservative test-only transfer drains VC to AO at 100 mL/s while
    # native valve derivatives are zero. It has an exact event time.
    model.initial_state[11:]=0
    model.initial_state[3]=model.inventory.unstressed_volumes[3]+1.
    def rhs(t,y,modes,active=True):
        dy=np.zeros(15); dy[3]=-100.; dy[2]=100.; return dy
    model.rhs=rhs
    tr=integrate_source(model,.02,active=False,max_step=.001)
    assert tr.termination['kind']=='model_domain_limit'
    assert tr.termination['node']=='VC'
    assert tr.time[-1]==pytest.approx(.01,abs=1e-10)
    assert tr.states[3,-1]==pytest.approx(model.inventory.unstressed_volumes[3],abs=1e-9)
    assert tr.states[:11,-1].sum()==pytest.approx(model.initial_state[:11].sum(),abs=1e-9)
