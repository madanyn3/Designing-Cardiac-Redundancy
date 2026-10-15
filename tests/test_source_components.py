from pathlib import Path
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from cardiac_redundancy.mechanics import (ModelDomainError, Wall, driver,
    inertial_derivative, resistive_flow, septal_root, vascular_pressure)
from cardiac_redundancy.source import SourceModel, INCIDENCE

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def model():
    return SourceModel(ROOT / "configs/baseline.yaml")


@pytest.mark.parametrize("e", [0, .25, .5, 1])
def test_symmetric_ventricles_have_zero_septal_displacement(e):
    wall = Wall(2, 0, .1, .03, 0)
    septum = Wall(10, 0, 1, .1, 0)
    root, error = septal_root(80, 80, e, wall, wall, septum, hint=7)
    assert abs(root) < 1e-11
    assert abs(error) < 1e-8


def test_fully_active_septum_has_analytic_linear_solution(model):
    left, right, septum = model.walls
    lv, rv = model.initial_state[:2]
    exact = (septum.elastance * septum.dead_volume + left.elastance * (lv-left.dead_volume)
             - right.elastance * (rv-right.dead_volume)) / (septum.elastance+left.elastance+right.elastance)
    root, error = septal_root(lv, rv, 1, left, right, septum)
    assert root == pytest.approx(exact, abs=1e-12)
    assert abs(error) < 1e-8


def test_passive_reference_and_negative_wall_coordinate():
    wall = Wall(2, 0, 1, .1, 0)
    assert wall.pressure(0, 0) == 0
    assert wall.pressure(-1, 0) < 0  # Geometric coordinate is not cavity blood volume.


@pytest.mark.parametrize("e", [-.1, 1.1, float("nan")])
def test_invalid_activation_rejected(e):
    with pytest.raises(ModelDomainError):
        Wall(2, 0, 1, .1, 0).pressure(1, e)


def test_overflow_is_explicit_failure():
    with pytest.raises(ModelDomainError):
        Wall(2, 0, 1, .1, 0).pressure(1e9, 0)


def test_driver_periodicity_and_peak():
    assert driver(.375, .75, 1, 80, .375) == 1
    for t in [-.2, .01, .2, .5]:
        assert driver(t, .75, 1, 80, .375) == pytest.approx(driver(t+750, .75, 1, 80, .375), abs=1e-12)


def test_external_pressure_and_offset_translation():
    assert vascular_pressure(100, .2, -4, 0) == 16
    assert vascular_pressure(600, .2, -4, 500) == 16
    assert resistive_flow(5, 10, 2) == -2.5


def test_rc_relaxation_matches_analytic_solution():
    r, c, initial = 2., 3., 12.
    def rhs(t, y):
        return [-resistive_flow(vascular_pressure(y[0], 1/c), 0, r)]
    times = np.linspace(0, 12, 101)
    sol = solve_ivp(rhs, (0, 12), [initial], t_eval=times, rtol=1e-10, atol=1e-12)
    assert sol.success
    np.testing.assert_allclose(sol.y[0], initial*np.exp(-times/(r*c)), rtol=1e-8)


def test_inertial_forward_flow_decays_before_closure():
    # Adverse pressure does not erase positive inertial flow instantaneously.
    gradient, resistance, inertance, q0 = -2., 1., .5, 4.
    exact_closure = inertance/resistance * math.log(3)
    def closure(t, y): return y[0]
    closure.terminal, closure.direction = True, -1
    sol = solve_ivp(lambda t,y: [inertial_derivative(gradient,y[0],resistance,inertance)],
                    (0, 1), [q0], events=closure, rtol=1e-10, atol=1e-12)
    assert sol.t_events[0][0] == pytest.approx(exact_closure, abs=1e-8)
    assert sol.y[0, 1] > 0


@pytest.mark.parametrize("active", [False, True])
def test_source_rhs_conserves_volume_without_mutation(model, active):
    y = model.initial_state.copy()
    before = y.copy()
    dy = model.rhs(.375, y, [True, False, True, False], active=active)
    assert abs(sum(dy[:6])) < 1e-10
    np.testing.assert_array_equal(y, before)
    assert dy[7] == dy[9] == 0


def test_incidence_conserves_arbitrary_signed_transfers():
    np.testing.assert_array_equal(INCIDENCE.sum(axis=0), np.zeros(6))
    assert abs(sum(INCIDENCE @ np.array([1, -2, 3, -4, 5, -6.]))) < 1e-12


def test_arrest_removes_activation_but_retains_passive_pressure(model):
    y = model.initial_state.copy()
    normal, _ = model.pressures(.375, y)
    arrested, meta = model.pressures(.375, y, active=False)
    assert meta["activation"] == 0
    assert abs(meta["septal_residual_mmHg"]) < 1e-8
    assert arrested["LV"] != 0
    assert arrested["LV"] < normal["LV"]
    assert normal["AO"] == arrested["AO"]


def test_source_volume_and_closed_flow_contract(model):
    assert sum(model.initial_state[:6]) == pytest.approx(1500, abs=1e-9)
    y = model.initial_state.copy()
    y[6:] = 0
    expected = model.rhs(0, y, [False]*4)
    perturbed = y.copy()
    perturbed[6:] = 1e-10
    np.testing.assert_array_equal(model.rhs(0, perturbed, [False]*4), expected)
    bad = model.initial_state.copy()
    bad[0] = 0
    with pytest.raises(ModelDomainError):
        model.pressures(0, bad)
