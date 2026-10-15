import json
from pathlib import Path
import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp
from cardiac_redundancy.demand_schedule import DemandCommand, FiniteDemandSchedule, scenario_vector
from cardiac_redundancy.organ_oxygen import OrganDemand, OrganOxygenTransport


@pytest.mark.parametrize('duration', [10, 60, 300])
def test_finite_induction_hold_return_and_exact_integral(duration):
    s = FiniteDemandSchedule([60, 120], [DemandCommand(5, [6, 12], duration),
        DemandCommand(10+duration, [60, 120], duration)])
    np.testing.assert_allclose(s.at(5), [60, 120])
    np.testing.assert_allclose(s.at(5+duration/2), [33, 66])
    np.testing.assert_allclose(s.at(5+duration), [6, 12])
    end = 10+2*duration
    np.testing.assert_allclose(s.at(end), [60, 120])
    expected = (np.array([60, 120])*5+np.array([6, 12])*5+duration*np.array([66, 132]))/60
    np.testing.assert_allclose(s.integral_mL(0, end), expected)
    for knot in s.breakpoints(0, end):
        np.testing.assert_allclose(s.at(knot-1e-8), s.at(knot+1e-8), atol=1e-6)


def test_interrupted_induction_returns_from_actual_demand_without_jump():
    s = FiniteDemandSchedule([100], [DemandCommand(0, [0], 100), DemandCommand(20, [100], 10)])
    assert s.at(20)[0] == 80
    assert s.at(25)[0] == 90
    assert s.at(30)[0] == 100
    assert s.integral_mL(0, 40)[0] == pytest.approx((1800+900+1000)/60)


def test_independent_organ_durations_and_quadrature():
    s = FiniteDemandSchedule([100, 60], [DemandCommand(3, [20, 30], [10, 60]),
        DemandCommand(25, [100, 60], [60, 10])])
    np.testing.assert_allclose(s.at(13), [20, 55])
    knots = s.breakpoints(0, 100)
    expected = [quad(lambda t: s.at(t)[i]/60, 0, 100, points=knots)[0] for i in range(2)]
    np.testing.assert_allclose(s.integral_mL(0, 100), expected, atol=1e-10)


def test_m1_targets_and_independent_working_auxiliary_heart():
    budget = json.loads((Path(__file__).resolve().parents[1]/'data/m1-demand-budget.json').read_text())
    organs = [OrganDemand(p, p, v, .75, p) for p, v in budget['baseline_mL_min'].items()]
    organs.append(OrganDemand('aux', 'aux', 10, .75, 'auxiliary_heart'))
    for scenario in budget['scenarios']:
        target = scenario_vector(organs, budget, scenario)
        assert target[:-1].sum() == pytest.approx(budget['scenarios'][scenario]['whole_body_mL_min'])
        assert target[-1] == 10
    target = scenario_vector(organs, budget, 'hypothetical_5_percent')
    assert target[:-1].sum() == pytest.approx(sum(budget['baseline_mL_min'].values())*.05)
    assert target[-1] == 10


def test_full_long_transition_coupled_to_finite_oxygen_and_deficit_ledgers():
    m = OrganOxygenTransport(['PA', 'PU'], [('pulmonary_exchange', 'PA', 'PU'), ('tissue', 'PU', 'PA')],
        [OrganDemand('tissue', 'tissue', 120, .5, 'test')])
    s = FiniteDemandSchedule([120], [DemandCommand(5, [30], 300), DemandCommand(310, [120], 300)])
    y = m.initial_state([200, 200]); invariant = m.balance(y)
    for a, b in zip(s.breakpoints(0, 620)[:-1], s.breakpoints(0, 620)[1:]):
        sol = solve_ivp(lambda t, x: m.rhs([200, 200], x, [10, 10], 80, s.at(t)),
                        (a, b), y, method='Radau', rtol=1e-8, atol=1e-9, max_step=2)
        assert sol.success and np.min(sol.y[:3]) >= 0
        np.testing.assert_allclose(m.balance(sol.y.T), invariant, atol=1e-7)
        y = sol.y[:, -1]
    assert y[-2]+y[-1] == pytest.approx(s.integral_mL(0, 620)[0], abs=1e-6)


def test_no_suppression_and_invalid_commands():
    s = FiniteDemandSchedule([60])
    assert s.integral_mL(0, 100)[0] == 100
    with pytest.raises(ValueError): FiniteDemandSchedule([60], [DemandCommand(0, [0], 0)])
    with pytest.raises(ValueError): FiniteDemandSchedule([60], [DemandCommand(0, [-1], 10)])
    with pytest.raises(ValueError): FiniteDemandSchedule([60], [DemandCommand(1, [0], 10), DemandCommand(1, [60], 10)])
    with pytest.raises(ValueError): s.at(float('nan'))
    with pytest.raises(ValueError): s.integral_mL(10, 0)
