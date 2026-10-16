import numpy as np
import pytest
from cardiac_redundancy.structural_closures import flow_extraction, transition_fraction


def test_extraction_reference_zero_and_high_flow_limits():
    assert np.isclose(flow_extraction(5, 5, .75), .75)
    assert flow_extraction(0, 5, .75) == 1
    assert flow_extraction(0, 5, 0) == 0
    ps = -5*np.log(.25)
    assert np.isclose(1e9*flow_extraction(1e9, 5, .75), ps)


def test_reverse_flow_uses_magnitude_and_uptake_is_bounded():
    q = np.array([-100, -1, 0, 1, 100])
    e = flow_extraction(q, 5, .75)
    np.testing.assert_allclose(e, e[::-1])
    assert np.all((e >= 0) & (e <= 1))
    assert np.all(np.abs(q)*e <= np.abs(q))


def test_transition_integrals_and_limits():
    t = np.linspace(0, 10, 10001)
    linear = transition_fraction(t, .25, 5, 'linear')
    exponential = transition_fraction(t, .25, 5, 'exponential')
    assert np.isclose(np.trapezoid(linear, t), 4.375)
    assert np.isclose(np.trapezoid(exponential, t), .25*10+.75*5*(1-np.exp(-2)))
    assert linear[0] == exponential[0] == 1
    assert linear[-1] == .25
    np.testing.assert_allclose(transition_fraction(t, .25, 0, 'linear'), .25)


def test_invalid_parameters_rejected():
    with pytest.raises(ValueError):
        flow_extraction(1, 0, .5)
    with pytest.raises(ValueError):
        flow_extraction(1, 1, 1)
    with pytest.raises(ValueError):
        transition_fraction(1, .25, 5, 'unknown')
