from dataclasses import replace

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from cardiac_redundancy.oxygen import (
    ConservativeOxygenTransport, OxygenParameters, oxygen_content,
    oxygen_pressure, saturation,
)


def fixture(parameters=OxygenParameters()):
    return ConservativeOxygenTransport(
        ['PA', 'PU'], [('pulmonary_exchange', 'PA', 'PU'), ('return', 'PU', 'PA')], parameters)


def test_content_units_half_saturation_and_inverse():
    p = OxygenParameters()
    assert saturation(26.8) == pytest.approx(.5)
    assert oxygen_content(26.8) == pytest.approx((1.34*15*.5+.0031*26.8)/100)
    assert oxygen_content(0) == 0
    pressures = np.array([0, .1, 1, 5, 20, 26.8, 40, 100, 500])
    contents = oxygen_content(pressures)
    assert np.all(np.diff(contents) > 0)
    np.testing.assert_allclose([oxygen_pressure(c) for c in contents], pressures, atol=1e-9)
    assert oxygen_content(100, replace(p, hemoglobin_g_dL=0)) == pytest.approx(.0031)


def test_gas_units_and_no_flow_ventilation_analytic_solution():
    m = fixture(); p = m.parameters; volume = np.array([100., 100.])
    assert p.alveolar_volume_mL_STPD == pytest.approx(2500*713/760*273.15/310.15)
    y = m.initial_state(volume)
    sol = solve_ivp(lambda t, x: m.rhs(volume, x, [0, 0], 80), (0, 20), y, rtol=1e-10, atol=1e-11)
    assert sol.success
    expected = p.alveolar_volume_mL_STPD*(.21+(.14-.21)*np.exp(-80*20/2500))
    assert sol.y[-2, -1] == pytest.approx(expected, rel=1e-9)
    np.testing.assert_allclose(sol.y[:2, -1], y[:2])
    np.testing.assert_allclose(m.balance(sol.y.T), m.balance(y), atol=1e-10)


@pytest.mark.parametrize('flow', [20., -20., 0.])
def test_closed_apnea_conserves_oxygen_and_equilibrates(flow):
    m = fixture(); volume = np.array([100., 100.]); y = m.initial_state(volume)
    rate = m.rhs(volume, y, [flow, flow])
    assert sum(rate[:-1]) == pytest.approx(0, abs=1e-12)
    sol = solve_ivp(lambda t, x: m.rhs(volume, x, [flow, flow]), (0, 30), y, rtol=1e-9, atol=1e-10)
    assert sol.success and np.min(sol.y[:-1]) >= 0
    np.testing.assert_allclose(m.balance(sol.y.T), m.balance(y), atol=1e-9)
    if flow:
        assert sol.y[-2, -1] < y[-2]  # blood gains only what finite lung loses
    else:
        np.testing.assert_array_equal(sol.y[:, -1], y)


def test_exchange_disabled_and_deoxygenated_gas():
    m = fixture(replace(OxygenParameters(), exchange_fraction=0))
    y = m.initial_state([100, 100])
    np.testing.assert_allclose(m.rhs([100, 100], y, [10, 10]), 0, atol=1e-12)
    m = fixture(); y = m.initial_state([100, 100], alveolar_fraction=0)
    r = m.rhs([100, 100], y, [10, 10])
    assert r[-2] > 0 and r[1] < 0  # oxygen can leave blood for empty gas
    assert sum(r[:-1]) == pytest.approx(0, abs=1e-12)


def test_hepatic_mixing_uses_actual_inflows():
    m = ConservativeOxygenTransport(['PA', 'PU', 'A', 'B', 'C'],
        [('pulmonary_exchange', 'PA', 'PU'), ('a', 'A', 'H'), ('b', 'B', 'H'), ('out', 'H', 'C')])
    y = m.initial_state([100]*5); y[:5] = [0, 0, 20, 10, 0]
    r = m.rhs([100]*5, y, [0, 2, 1, 3])
    np.testing.assert_allclose(r[:5], [0, 0, -.4, -.1, .5])
    assert m.balance(r) == pytest.approx(0)


def test_invalid_domain_is_reported_not_clipped():
    m = fixture(); y = m.initial_state([100, 100])
    with pytest.raises(ValueError): m.rhs([100, 100], y, [1])
    with pytest.raises(ValueError): m.rhs([0, 100], y, [1, 1])
    y[-2] = -1
    with pytest.raises(ValueError): m.rhs([100, 100], y, [1, 1])
    with pytest.raises(ValueError): OxygenParameters(exchange_fraction=1.1)
    with pytest.raises(ValueError): oxygen_content(-1)
