"""Circuit incidence, similarity, isolation and generalized-event regression."""
import numpy as np
import pytest
from cardiac_redundancy.circuit import CircuitModel
from cardiac_redundancy.circuit_integration import integrate_circuit
from cardiac_redundancy.integration import integrate_source
from cardiac_redundancy.mechanics import septal_root
from test_backup_pump import ARGS, KW


def test_general_integrator_exact_four_valve_regression():
    m=CircuitModel(*ARGS,**KW)
    m.native.set_controls(time=0,ramp_s=.1,pump_on=True,right_head=20,left_head=100)
    a=integrate_source(m.native,.12,active=False,sample_step=.01)
    b=integrate_circuit(m,.12,active=False,sample_step=.01)
    np.testing.assert_array_equal(a.states,b.states)
    assert a.events==b.events and a.termination==b.termination


@pytest.mark.parametrize('architecture',['B1','E1','E2','BIO'])
def test_topology_conservation_and_disabled_support(architecture):
    m=CircuitModel(*ARGS,architecture=architecture,**KW)
    assert sum(m.initial_state[:m.volume_count])==pytest.approx(5250)
    assert np.all(m.incidence.sum(axis=0)==0)
    q=m.flows(0,m.initial_state,[True]*len(m.valve_names))
    assert sum(m.rhs(0,m.initial_state,[True]*len(m.valve_names))[:m.volume_count])==pytest.approx(0,abs=1e-10)
    assert all(q[n]==0 for n in q if n.startswith('support_') or n=='oxygenator')
    if architecture=='E1':assert 'support_left' not in m.flow_names
    if architecture=='E2':assert ('support_left','PU','AO') in m.edges
    if architecture in ('E1','E2'):assert ('oxygenator','VC','AO') in m.edges


def test_auxiliary_similarity_allocation_isolation_and_release():
    m=CircuitModel(*ARGS,architecture='BIO',**KW)
    source=m.native.base.source
    np.testing.assert_allclose(m.initial_state[16:18],m.alpha*source.initial_state[:2])
    for donor,idx in (('PU',16),('VC',17)):
        i=m.volume_names.index(donor)
        assert m.initial_state[i]+m.initial_state[idx]==pytest.approx(m.native.initial_state[i])
    tr=integrate_circuit(m,.05,sample_step=.01)
    assert tr.termination is None
    np.testing.assert_array_equal(tr.states[16:18],np.repeat(m.initial_state[16:18,None],tr.states.shape[1],axis=1))
    assert np.all(tr.states[22:]==0)
    y=tr.states[:,-1].copy();m.release_auxiliary(.05)
    p,meta=m.pressures(.425,y,active=False)
    assert meta['aux_activation']==pytest.approx(.375)
    # Normalized cavity coordinates give the frozen E/lambda/V scaling exactly.
    lv,rv=y[16:18]/m.alpha
    v,_=septal_root(lv,rv,meta['aux_activation'],*source.walls)
    assert p['aLV']-p['aPERI']==pytest.approx(source.walls[0].pressure(lv-v,meta['aux_activation']))
    np.testing.assert_array_equal(y,tr.states[:,-1])
    with pytest.raises(ValueError):m.release_auxiliary(.1)
    m.withdrawal=.425
    assert m.pressures(1.425,y,False)[1]['aux_activation']==0


def test_initial_domain_exit_is_explicit_and_preserves_state():
    m=CircuitModel(*ARGS,**KW);y=m.initial_state.copy()
    i=m.volume_names.index('PU');delta=y[i]-m.native.U[i];y[i]-=delta;y[3]+=delta
    tr=integrate_circuit(m,1,initial_state=y,active=False)
    assert tr.termination['node']=='PU' and tr.termination['kind']=='model_domain_limit'
    np.testing.assert_array_equal(tr.states[:,0],y)


def test_native_restart_preserves_forward_inertial_states():
    m=CircuitModel(*ARGS,**KW)
    a=integrate_circuit(m,.05,active=False,sample_step=.01)
    b=integrate_circuit(m,.05,start=.05,initial_state=a.states[:,-1],active=True,sample_step=.01)
    np.testing.assert_array_equal(a.states[:,-1],b.states[:,0])
    assert b.termination is None
