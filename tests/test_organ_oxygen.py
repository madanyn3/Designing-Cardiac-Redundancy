import json
from pathlib import Path
import numpy as np
import pytest
from scipy.integrate import solve_ivp

from cardiac_redundancy.organ_oxygen import OrganDemand, OrganOxygenTransport, map_demands
from cardiac_redundancy.oxygen import ConservativeOxygenTransport
from cardiac_redundancy.whole_network import WholeNetworkModel

ROOT = Path(__file__).resolve().parents[1]


def hepatic():
    names = ['PA', 'PU', 'AO', 'VC']
    edges = [('pulmonary_exchange', 'PA', 'PU'), ('gut', 'AO', 'H'),
             ('arterial', 'AO', 'H'), ('liver', 'H', 'VC')]
    organs = [OrganDemand('gut', 'gut', 30, 1, 'gut'), OrganDemand('liver', 'liver', 60, 1, 'liver')]
    m = OrganOxygenTransport(names, edges, organs)
    y = m.initial_state([100]*4)
    y[:4] = [20, 20, 20, 10]  # AO .2, VC .1 mL O2/mL blood
    return m, y


def test_portal_uptake_precedes_mixing_and_liver():
    m, y = hepatic()
    r, d = m.evaluate([100]*4, y, [0, 5, 5, 10])
    assert d['hepatic_content'] == pytest.approx(.15)
    np.testing.assert_allclose(d['delivered'], [1, 1.5])
    np.testing.assert_allclose(d['consumed'], [.5, 1])
    np.testing.assert_allclose(r[:4], [0, 0, -2, .5])
    assert m.balance(r) == pytest.approx(0)


def test_reversed_liver_uptake_precedes_reversed_portal_uptake():
    m, y = hepatic()
    r, d = m.evaluate([100]*4, y, [0, -5, -5, -10])
    # VC supplies 1 mL/s, all consumed by liver before reaching H/gut.
    np.testing.assert_allclose(d['consumed'], [0, 1])
    np.testing.assert_allclose(d['unmet'], [.5, 0])
    assert d['hepatic_content'] == 0
    assert m.balance(r) == pytest.approx(0)


def test_closed_branch_demand_persists_and_supply_limited_uptake():
    m, y = hepatic()
    r, d = m.evaluate([100]*4, y, [0, 0, 1, 1])
    np.testing.assert_allclose(d['consumed'], [0, .2])
    np.testing.assert_allclose(d['unmet'], [.5, .8])
    assert np.all(d['consumed'] <= d['delivered'])
    assert m.balance(r) == pytest.approx(0)


def test_zero_uptake_regresses_to_foundation_with_deficit_ledger():
    m, y = hepatic(); base = ConservativeOxygenTransport(m.names, m.edges)
    r, d = m.evaluate([100]*4, y, [3, 2, 3, 5], 50, demand_mL_min=[0, 0])
    np.testing.assert_allclose(r[:6], base.rhs([100]*4, y[:6], [3, 2, 3, 5], 50), atol=1e-14)
    np.testing.assert_array_equal(r[6:], 0)


def test_dynamic_consumption_and_deficit_ledgers_are_independent():
    names = ['PA', 'PU']; edges = [('pulmonary_exchange', 'PA', 'PU'), ('tissue', 'PU', 'PA')]
    m = OrganOxygenTransport(names, edges, [OrganDemand('tissue', 'tissue', 120, .5, 'test')])
    initial = m.initial_state([200, 200])
    sol = solve_ivp(lambda t, y: m.rhs([200, 200], y, [10, 10]), (0, 10), initial,
                    method='Radau', rtol=1e-9, atol=1e-10)
    assert sol.success and np.min(sol.y[:-2]) >= 0
    np.testing.assert_allclose(m.balance(sol.y.T), m.balance(initial), atol=1e-8)
    assert sol.y[-2, -1]+sol.y[-1, -1] == pytest.approx(20, abs=1e-8)
    assert sol.y[-1, -1] > 0


def project():
    return WholeNetworkModel(ROOT/'configs/baseline.yaml', ROOT/'configs/m3-evidence-regional.yaml',
        ROOT/'configs/m2-frozen-overlay.yaml', arterial_fraction=.2, inlet_resistance_fraction=.1)


@pytest.mark.parametrize('scenario', ['unsuppressed', 'modest_reduction', 'stronger_reduction'])
def test_mapping_preserves_every_m1_pool_and_all_tissue_routes(scenario):
    m = project()
    budget = json.loads((ROOT/'data/m1-demand-budget.json').read_text())['scenarios'][scenario]['organ_mL_min']
    mapping = json.loads((ROOT/'configs/m4-organ-mapping.json').read_text())['mapping']
    organs = map_demands(budget, mapping, m.flow_names, .75)
    for pool, demand in budget.items():
        assert sum(o.demand_mL_min for o in organs if o.parent_pool == pool) == pytest.approx(demand)
    supplies = {b.name for b in m.base.branches if b.donor == 'AO'}
    expected = {n+'_exchange' if n.startswith('limb_') else 'hepatic_outflow' if n == 'liver_arterial' else n for n in supplies}
    assert {o.edge for o in organs} == expected
    assert all(not o.edge.endswith('_return') for o in organs)
    assert sum(o.demand_mL_min for o in organs) == pytest.approx(sum(budget.values()))


def test_biological_pump_cannot_silently_have_zero_demand():
    m = project(); budget = json.loads((ROOT/'data/m1-demand-budget.json').read_text())['baseline_mL_min']
    mapping = json.loads((ROOT/'configs/m4-organ-mapping.json').read_text())['mapping']
    flows = (*m.flow_names, 'aux_coronary_left', 'aux_coronary_right')
    with pytest.raises(ValueError): map_demands(budget, mapping, flows, .75)
    organs = map_demands(budget, mapping, flows, .75, auxiliary_total_mL_min=10)
    assert sum(o.demand_mL_min for o in organs) == pytest.approx(sum(budget.values())+10)
    with pytest.raises(ValueError): map_demands(budget, mapping, m.flow_names, .75, auxiliary_total_mL_min=10)


def test_invalid_mapping_or_demand_rejected():
    with pytest.raises(ValueError): OrganDemand('a', 'a', -1, .5, 'test')
    with pytest.raises(ValueError): OrganDemand('a', 'a', 1, 1.01, 'test')
    with pytest.raises(ValueError): map_demands({'a': 1}, {'a': {'x': .5}}, ['x'], .5)
    with pytest.raises(ValueError): map_demands({'a': 1}, {}, [], .5)
    m, y = hepatic()
    with pytest.raises(ValueError): m.rhs([100]*4, y, [0, 1, 1, 3])
    with pytest.raises(ValueError): m.rhs([100]*4, y, [0, 1, 1, 2], demand_mL_min=[-1, 1])
