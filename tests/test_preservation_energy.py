import numpy as np
import pytest
from cardiac_redundancy.preservation_energy import native_fraction,integrated_native_fraction,oxygen_energy

def test_ramp_integral_and_plateau():
    assert native_fraction(0,.3,5)==1
    assert native_fraction(5,.3,5)==pytest.approx(.3)
    assert integrated_native_fraction(8,.3,5)==pytest.approx(4.15)
    x=np.linspace(0,8,8001)
    assert np.trapezoid([native_fraction(t,.3,5) for t in x],x)==pytest.approx(4.15)

def test_step_and_no_reduction():
    assert native_fraction(0,.3,0)==.3
    assert integrated_native_fraction(8,.3,0)==pytest.approx(2.4)
    assert integrated_native_fraction(8,1,5)==8

def test_energy_is_linear_and_not_free_credit():
    np.testing.assert_allclose(oxygen_energy([0,1,2]),[0,20.2924,40.5848])
    with pytest.raises(ValueError):oxygen_energy(-1)

@pytest.mark.parametrize('args',[(-1,.3,5),(1,0,5),(1,1.1,5),(1,.3,-1)])
def test_invalid_transition(args):
    with pytest.raises(ValueError):native_fraction(*args)
