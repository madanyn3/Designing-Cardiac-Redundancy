from dataclasses import replace
import numpy as np
import pytest
from scipy.integrate import solve_ivp
from cardiac_redundancy.oxygen import OxygenParameters
from cardiac_redundancy.organ_oxygen import OrganDemand,OrganOxygenTransport
from cardiac_redundancy.respiratory import RespiratoryParameters,RespiratoryOxygenTransport
from cardiac_redundancy.preservation import PreservationChoice,apply_policy


def model(store=0.,oxygenator=False,**kwargs):
    edges=[('pulmonary_exchange','PA','PU'),('left','PU','AO'),('resp','AO','VC'),('right','VC','PA')]
    if oxygenator:edges.append(('oxygenator','VC','AO'))
    return RespiratoryOxygenTransport(['PA','PU','AO','VC'],edges,
        [OrganDemand('resp','resp',60,.75,'respiratory_muscle')],store_seconds=store,**kwargs)


def test_zero_store_maintained_regresses_to_old_model():
    m=model();old=OrganOxygenTransport(m.names,m.edges,m.organs,m.parameters)
    y=m.initial_state([200]*4);q=[10]*4
    r,d=m.evaluate([200]*4,y,q)
    np.testing.assert_allclose(r[:m.base_size],old.rhs([200]*4,y[:m.base_size],q,4690/60),atol=1e-12)
    assert m.balance(r)==pytest.approx(0,abs=1e-12)


@pytest.mark.parametrize('store',[0,5,15])
def test_finite_store_and_apnea_conserve_oxygen(store):
    m=model(store);y=m.initial_state([200]*4)
    sol=solve_ivp(lambda t,x:m.rhs([200]*4,x,[10]*4,mode='apnea'),(0,30),y,method='Radau',rtol=1e-8,atol=1e-9,max_step=.2)
    assert sol.success
    np.testing.assert_allclose(m.balance(sol.y.T),m.balance(y),atol=1e-7)
    assert sol.y[-1].min()>=-1e-7 and sol.y[-1].max()<=store+1e-7
    assert sol.y[5,-1]==0 # no external native gas source during apnea


def test_trapped_store_consumption_is_finite_and_not_deleted():
    m=model(5);y=m.initial_state([200]*4)
    _,d=m.evaluate([200]*4,y,[0]*4,mode='apnea')
    assert d['blood_extracted'][0]==0 and d['consumed'][0]==1
    sol=solve_ivp(lambda t,x:m.rhs([200]*4,x,[0]*4,mode='apnea'),(0,30),y,rtol=1e-8,atol=1e-9,max_step=.1)
    assert sol.success and sol.y[6,-1]<=5+1e-7
    assert sol.y[7,-1]>24.9
    np.testing.assert_allclose(m.balance(sol.y.T),m.balance(y),atol=1e-7)


def test_feedback_depends_on_reference_need_not_suppressed_demand():
    m=model();y=m.initial_state([200]*4)
    _,d=m.evaluate([200]*4,y,[10]*4,mode='feedback',demand_mL_min=[3])
    assert d['respiratory_availability']==pytest.approx(.05)
    assert d['ventilation_target_mL_s_BTPS']==pytest.approx(4690/60*.05)
    r,d=m.evaluate([200]*4,y,[10]*4,mode='feedback',drive_fraction=0)
    assert d['ventilation_target_mL_s_BTPS']==0 and r[m.base_size]<0


@pytest.mark.parametrize('sweep',[0,100])
def test_oxygenator_finite_gas_and_external_ledger(sweep):
    m=model(oxygenator=True,parameters=replace(OxygenParameters(),exchange_fraction=0),
            respiratory=replace(RespiratoryParameters(),oxygenator_gas_volume_mL_STPD=10))
    y=m.initial_state([200]*4);q=[10,10,20,10,10]
    sol=solve_ivp(lambda t,x:m.rhs([200]*4,x,q,mode='apnea',sweep_mL_s_STPD=sweep),
                  (0,30),y,method='Radau',rtol=1e-8,atol=1e-9,max_step=.1)
    assert sol.success
    np.testing.assert_allclose(m.balance(sol.y.T),m.balance(y),atol=1e-7)
    b=m.base_size
    assert sol.y[b+1].min()>=-1e-7
    assert sol.y[b+1].max()<=m.respiratory.oxygenator_gas_volume_mL_STPD+1e-7
    if not sweep:assert sol.y[b+2,-1]==0 and sol.y[b+1,-1]<y[b+1]
    else:assert sol.y[b+2,-1]>0


def test_reverse_oxygenator_transfer_and_disabled_exchange():
    m=model(oxygenator=True);y=m.initial_state([200]*4)
    rate,d=m.evaluate([200]*4,y,[-10,-10,-20,-10,-10],mode='apnea')
    assert d['oxygenator_transfer_mL_s']>0
    assert m.balance(rate)==pytest.approx(0,abs=1e-12)
    m=model(oxygenator=True,respiratory=replace(RespiratoryParameters(),oxygenator_exchange_fraction=0))
    _,d=m.evaluate([200]*4,m.initial_state([200]*4),[10,10,20,10,10])
    assert d['oxygenator_transfer_mL_s']==0


def test_outcome_exclusion_does_not_erase_demand_or_apply_controls():
    m=model();r=apply_policy(m.organs,{'resp':PreservationChoice(1,0,False)})
    assert r['demand_mL_min']==[60] and r['outcome_required']==[False]
    assert r['requested_conductance']=={'resp':0}
    with pytest.raises(ValueError):apply_policy(m.organs,{})
    with pytest.raises(ValueError):RespiratoryParameters(response_time_s=0)
    with pytest.raises(ValueError):m.rhs([200]*4,m.initial_state([200]*4),[10]*4,sweep_mL_s_STPD=1)


def test_time_varying_demand_with_isolated_store_preserves_both_ledgers():
    m=model(5);y=m.initial_state([200]*4)
    # D(t)=1+t/10 mL/s integrates analytically to 15 mL over ten seconds.
    sol=solve_ivp(lambda t,x:m.rhs([200]*4,x,[0]*4,mode='apnea',
        demand_mL_min=[60*(1+t/10)]),(0,10),y,rtol=1e-9,atol=1e-11,max_step=.02)
    assert sol.success
    assert sol.y[6,-1]+sol.y[7,-1]==pytest.approx(15,abs=1e-8)
    assert sol.y[6,-1]<=5+1e-8
    np.testing.assert_allclose(m.balance(sol.y.T),m.balance(y),atol=1e-8)
