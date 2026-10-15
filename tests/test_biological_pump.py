import numpy as np
import pytest
from cardiac_redundancy.biological_pump import BiologicalParameters, BiologicalPumpModel
from cardiac_redundancy.circuit import CircuitModel
from cardiac_redundancy.circuit_integration import integrate_circuit
from test_backup_pump import ARGS,KW

@pytest.mark.parametrize('kwargs',[{'left_scale':0},{'right_scale':1.1},{'period_s':float('nan')},{'period_s':True},{'left_contractility':2.1},{'septal_scale':.09}])
def test_bounds(kwargs):
    with pytest.raises(ValueError):BiologicalParameters(**kwargs)

def test_similarity_and_native_regression():
    old=CircuitModel(*ARGS,architecture='BIO',**KW);new=BiologicalPumpModel(*ARGS,**KW)
    np.testing.assert_array_equal(old.initial_state,new.initial_state)
    old.release_auxiliary(0);new.release_auxiliary(0)
    for t in (0,.2,.375,1.125):
        a,_=old.pressures(t,old.initial_state);b,_=new.pressures(t,new.initial_state)
        for key in a:assert a[key]==pytest.approx(b[key],abs=1e-7)
    a=integrate_circuit(old,.12,active=False,sample_step=.01)
    b=integrate_circuit(new,.12,active=False,sample_step=.01)
    np.testing.assert_allclose(a.states,b.states,atol=2e-6,rtol=1e-7)

def test_independent_allocation_and_conservation():
    p=BiologicalParameters(left_scale=.3,right_scale=.2,left_contractility=1.5)
    m=BiologicalPumpModel(*ARGS,biological=p,**KW)
    assert sum(m.initial_state[:18])==pytest.approx(5250)
    np.testing.assert_allclose(m.initial_state[16:18],np.array([.3,.2])*m.native.base.source.initial_state[:2])
    m.release_auxiliary(0)
    assert abs(sum(m.rhs(.375,m.initial_state,[True]*8)[:18]))<1e-9
    assert abs(m.pressures(.375,m.initial_state)[1]['aux_septal_residual_mmHg'])<1e-8
