import numpy as np
import pytest
from cardiac_redundancy.bilateral_support import LungAllocation,BilateralSupport,BilateralRespiration
from cardiac_redundancy.respiratory import RespiratoryOxygenTransport
from test_backup_pump import ARGS,KW

def test_equal_open_lungs_reproduce_single_pool_without_extra_store():
    allocation=LungAllocation(preserved_perfusion=1,preserved_nutrition=1,preserved_ventilation=1)
    a=BilateralRespiration(('PA','PU'),(('pulmonary_exchange','PA','PU'),('pulmonary_preserved','PA','PU')),(),allocation=allocation)
    b=RespiratoryOxygenTransport(('PA','PU'),(('pulmonary_exchange','PA','PU'),),())
    v=np.array([100.,100.]);ya=a.initial_state(v);yb=b.initial_state(v)
    da,_=a.evaluate(v,ya,[2.,2.]);db,_=b.evaluate(v,yb,[4.])
    np.testing.assert_allclose(da[:a.base_size],db[:b.base_size])
    assert a.balance(ya)==pytest.approx(b.balance(yb))
    assert a.balance(da)==pytest.approx(0,abs=1e-12)

def test_unventilated_lung_uses_its_own_gas_and_conserves_oxygen():
    a=BilateralRespiration(('PA','PU'),(('pulmonary_exchange','PA','PU'),('pulmonary_preserved','PA','PU')),())
    v=np.array([100.,100.]);y=a.initial_state(v)
    cap=a.parameters.alveolar_volume_mL_STPD
    y[2]=.1*cap;y[a.base_size+1]=.02*cap
    d,z=a.evaluate(v,y,[3.,.3])
    assert z['lung_ventilation_mL_s'][1]==0
    assert z['lung_fractions'][0]==pytest.approx(.16)
    assert z['lung_fractions'][1]==pytest.approx(.04)
    assert a.balance(d)==pytest.approx(0,abs=1e-12)
    assert d[a.base_size+1]==pytest.approx(-z['lung_transfer_mL_s'][1])

@pytest.mark.parametrize('device',['BIO','mechanical'])
def test_regional_paths_conserve_blood_and_protect_mandatory_organs(device):
    m=BilateralSupport(*ARGS,device=device,**KW);m.apply_support_controls(0)
    f=m.native.regional_at(0)
    assert f['coronary_left']==f['coronary_right']==f['kidney_protected']==1
    assert f['kidney_other']==.25 and f['pulmonary_exchange']==pytest.approx(.55)
    q=m.flows(0,m.initial_state,[True]*8,False)
    assert q['pulmonary_exchange']/q['pulmonary_preserved']==pytest.approx(10)
    assert q['lung_supported']/q['lung_preserved']==pytest.approx(4)
    assert abs(sum(m.rhs(0,m.initial_state,[True]*8,False)[:18]))<1e-8
    assert ('aux_coronary_left' in m.flow_names)==(device=='BIO')

@pytest.mark.parametrize('kw',[{'supported_share':0},{'preserved_perfusion':-1},{'supported_perfusion':0,'preserved_perfusion':0}])
def test_invalid_allocation(kw):
    with pytest.raises(ValueError):LungAllocation(**kw)
