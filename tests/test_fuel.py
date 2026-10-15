from dataclasses import replace
import numpy as np
import pytest
from scipy.integrate import solve_ivp

from cardiac_redundancy.fuel import (FuelParameters,GlucoseTransport,FuelOxygenTransport,
    O2_ML_STPD_PER_MMOL_GLUCOSE,GLUCOSE_MG_PER_MMOL,GLUCOSYL_MG_PER_MMOL)
from cardiac_redundancy.organ_oxygen import OrganDemand
from cardiac_redundancy.respiratory import RespiratoryOxygenTransport
from cardiac_redundancy.transport import marker_rhs


def fixture(parameters=FuelParameters(),extraction=1.):
    names=('PA','PU','AO','VC')
    edges=(('pulmonary_exchange','PA','PU'),('gut','AO','H'),
           ('arterial','AO','H'),('hepatic_outflow','H','VC'))
    organs=(OrganDemand('gut','gut',30,.75,'gut'),
            OrganDemand('liver','hepatic_outflow',60,.75,'liver'))
    g=GlucoseTransport(names,edges,organs,parameters,extraction)
    y=g.initial_state([100]*4);y[:4]=[1,1,1,.5]
    return g,y


def test_stock_and_blood_units_do_not_confuse_polymer_and_free_glucose():
    p=FuelParameters();assert p.glycogen_mmol==pytest.approx(329.76)
    assert p.glycogen_mmol*GLUCOSE_MG_PER_MMOL/1000==pytest.approx(59.40824256)
    assert p.glycogen_mmol*GLUCOSYL_MG_PER_MMOL/1000==pytest.approx(53.46761616)
    assert 5250*p.blood_concentration_mmol_mL==pytest.approx(23.64864864865)
    assert O2_ML_STPD_PER_MMOL_GLUCOSE==pytest.approx(134.484)


def test_portal_consumption_precedes_liver_with_separate_export():
    g,y=fixture();r,d=g.evaluate([100]*4,y,[0,5,5,10],[.02,.06])
    np.testing.assert_allclose(d['delivered_mmol_s'],[.05,.08])
    np.testing.assert_allclose(d['consumed_mmol_s'],[.02,.06])
    assert d['release_mmol_s']==pytest.approx(.0035)
    np.testing.assert_allclose(r[:4],[0,0,-.1,.0235])
    assert g.balance(r)==pytest.approx(0,abs=1e-15)


@pytest.mark.parametrize('flows',[[0,0,0,0],[0,-5,-5,-10]])
def test_zero_or_reverse_hepatic_flow_blocks_export_but_not_demand(flows):
    g,y=fixture();r,d=g.evaluate([100]*4,y,flows,[.02,.06])
    assert d['release_mmol_s']==0
    assert np.sum(d['unmet_mmol_s'])>0
    assert g.balance(r)==pytest.approx(0,abs=1e-15)


def test_zero_demand_and_zero_release_regresses_to_passive_transport():
    g,y=fixture(replace(FuelParameters(),release_cap_mmol_min=0))
    q=[3,5,5,10];r,_=g.evaluate([100]*4,y,q,[0,0])
    np.testing.assert_allclose(r[:4],marker_rhs([100]*4,y[:4],g.edges,q,g.names))
    np.testing.assert_array_equal(r[4:],0)


def test_empty_store_functional_loss_and_flow_limit_are_distinct():
    g,y=fixture();q=[0,5,5,10]
    assert g.evaluate([100]*4,y,q,[.1,.1],liver_function_fraction=0)[1]['release_mmol_s']==0
    y[4]=0
    assert g.evaluate([100]*4,y,q,[.1,.1])[1]['release_mmol_s']==0
    y[4]=100
    assert g.evaluate([100]*4,y,[0,.1,.1,.2],[.1,.1])[1]['release_mmol_s']==pytest.approx(.0002)


def test_nonoxidative_sink_consumes_trapped_blood_without_energy_credit():
    g,y=fixture(replace(FuelParameters(),nonoxidative_mmol_min=.12))
    r,_=g.evaluate([100]*4,y,[0]*4,[0,0])
    assert r[-2]==pytest.approx(.002)
    assert r[:4].sum()==pytest.approx(-.002)
    assert g.balance(r)==pytest.approx(0,abs=1e-15)


def test_coupling_cannot_burn_oxygen_when_glucose_is_empty():
    g,y=fixture(replace(FuelParameters(),release_cap_mmol_min=0))
    o=RespiratoryOxygenTransport(g.names,g.edges,g.organs,store_seconds=5)
    model=FuelOxygenTransport(o,g);z=model.initial_state([100]*4,po2_mmHg=100)
    z[model.oxygen_size:model.oxygen_size+4]=0
    r,d=model.evaluate([100]*4,z,[0,5,5,10])
    np.testing.assert_array_equal(d['consumed'],0)
    np.testing.assert_allclose(d['unmet'],[.5,1])
    assert model.balances(r)==pytest.approx((0,0),abs=1e-13)


def test_oxygen_limitation_returns_unused_glucose_and_retains_deficits():
    g,_=fixture();o=RespiratoryOxygenTransport(g.names,g.edges,g.organs)
    model=FuelOxygenTransport(o,g);z=model.initial_state([100]*4,po2_mmHg=0,alveolar_fraction=0)
    r,d=model.evaluate([100]*4,z,[0,5,5,10],mode='apnea')
    np.testing.assert_array_equal(d['consumed'],0)
    np.testing.assert_array_equal(d['glucose']['consumed_mmol_s'],0)
    np.testing.assert_allclose(d['unmet'],[.5,1])
    assert model.balances(r)==pytest.approx((0,0),abs=1e-13)


def test_serial_reservations_can_be_withdrawn_at_lower_oxygen_use():
    g,y=fixture();q=[0,5,5,10];d=[.08,.1]
    _,reservation=g.evaluate([100]*4,y,q,d)
    actual=reservation['consumed_mmol_s']*np.array([.1,.7])
    r,diag=g.evaluate([100]*4,y,q,d,actual_mmol_s=actual)
    np.testing.assert_allclose(diag['consumed_mmol_s'],actual)
    assert g.balance(r)==pytest.approx(0,abs=1e-15)


def test_joint_dynamic_conservation_and_original_demand_identity():
    g,_=fixture(replace(FuelParameters(),glycogen_mmol_L_liver=.001,initial_plasma_glucose_mmol_L=.001))
    o=RespiratoryOxygenTransport(g.names,g.edges,g.organs,store_seconds=5)
    model=FuelOxygenTransport(o,g);z=model.initial_state([100]*4,po2_mmHg=100)
    # Flow-balanced reduced circuit avoids changing the prescribed blood volumes.
    q=[0,0,0,0]
    sol=solve_ivp(lambda t,y:model.rhs([100]*4,y,q),(0,10),z,method='Radau',rtol=1e-9,atol=1e-11)
    assert sol.success
    for balances,initial in zip(model.balances(sol.y.T),model.balances(z)):
        np.testing.assert_allclose(balances,initial,atol=1e-9)
    n=g.n;m=g.m
    np.testing.assert_allclose(sol.y[n+2:n+2+m,-1]+sol.y[n+2+m:n+2+2*m,-1],[5,10],atol=1e-9)


def test_finite_stock_depletion_and_refinement_in_closed_circulating_fixture():
    organs=(OrganDemand('t','t',1,1,'t'),)
    g=GlucoseTransport(('A','V'),(('pump','V','A'),('t','A','H'),('hepatic_outflow','H','V')),
        organs,replace(FuelParameters(),liver_volume_L=1,glycogen_mmol_L_liver=.01,
                       initial_plasma_glucose_mmol_L=.01,release_cap_mmol_min=6),1.)
    y=g.initial_state([100,100]);time=np.linspace(0,60,121);rows=[]
    for tol in [1e-7,1e-10]:
        sol=solve_ivp(lambda t,z:g.evaluate([100,100],z,[10,10,10],[.01])[0],(0,60),y,
                      method='Radau',rtol=tol,atol=tol*.01,max_step=.25,t_eval=time)
        assert sol.success;rows.append(sol.y)
        assert sol.y[:3].min()>-1e-9
        np.testing.assert_allclose(g.balance(sol.y.T),g.balance(y),atol=1e-10)
        assert sol.y[2,-1]<1e-10
        assert sol.y[3,-1]+sol.y[4,-1]==pytest.approx(.6,abs=1e-9)
    np.testing.assert_allclose(rows[0],rows[1],atol=1e-7,rtol=1e-6)


def test_invalid_inputs_and_partial_fuel_scope_are_visible():
    with pytest.raises(ValueError):FuelParameters(release_cap_mmol_min=-1)
    with pytest.raises(ValueError):FuelParameters(accessible_fraction=1.1)
    g,y=fixture()
    with pytest.raises(ValueError):g.evaluate([100]*4,y,[0,1,1,3],[0,0])
    with pytest.raises(ValueError):g.evaluate([100]*4,y,[0,1,1,2],[0,0],actual_mmol_s=[1,0])
    o=RespiratoryOxygenTransport(g.names,g.edges,g.organs)
    model=FuelOxygenTransport(o,g,.5)
    assert model.evaluate([100]*4,model.initial_state([100]*4),[0,1,1,2])[1]['other_substrates_unassessed']


def test_numerical_extension_does_not_clip_negative_trial_state():
    g,y=fixture();y[4]=-.1;original=y.copy()
    with pytest.raises(ValueError):g.evaluate([100]*4,y,[0,5,5,10],[.1,.1])
    r,d=g.evaluate([100]*4,y,[0,5,5,10],[.1,.1],numerical_trial=True)
    np.testing.assert_array_equal(y,original)
    assert d['release_mmol_s']==0
    assert g.balance(r)==pytest.approx(0,abs=1e-15)
