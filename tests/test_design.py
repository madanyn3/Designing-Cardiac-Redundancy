import json
from pathlib import Path
import pytest
from cardiac_redundancy.design import (BOUNDS, DesignPoint, ScenarioPoint, contractility,
    demand_assumptions, mass_ratio, registry)
from cardiac_redundancy.biological_pump import BiologicalParameters

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = json.loads((ROOT/'configs/m5-myocardial-demand.json').read_text(encoding='utf-8'))['reference']


def test_columns_partition_the_declared_bounds():
    r = registry()
    assert set(r['design_controls']) | set(r['scenario_factors']) == set(BOUNDS)
    assert not set(r['design_controls']) & set(r['scenario_factors'])


def test_reference_point_reproduces_the_verified_fixture():
    d = DesignPoint()
    assert d.biological() == BiologicalParameters()
    assert d.rate_bpm == pytest.approx(80.)
    assert mass_ratio(d, ScenarioPoint()) == pytest.approx(.25)
    a = demand_assumptions(d, ScenarioPoint(), REFERENCE)
    assert a.mass_ratio == pytest.approx(REFERENCE['mass_ratio'])
    assert a.work_cost_multiplier == pytest.approx(REFERENCE['work_cost_multiplier'])
    assert (a.maintenance_fraction, a.activation_fraction, a.work_fraction) == (
        REFERENCE['maintenance_fraction'], REFERENCE['activation_fraction'], REFERENCE['work_fraction'])


def test_larger_cavity_is_charged_more_mass_and_faster_rate_more_activation():
    small, large = DesignPoint(), DesignPoint(left_scale=.6, right_scale=.6, septal_scale=.6, pericardial_scale=.6)
    assert mass_ratio(large, ScenarioPoint()) > mass_ratio(small, ScenarioPoint())
    assert mass_ratio(large, ScenarioPoint()) == pytest.approx(.6)
    assert DesignPoint(period_s=.5).rate_bpm > small.rate_bpm


def test_mass_coupling_is_a_scenario_not_a_design_control():
    assert 'mass_coupling' in registry()['scenario_factors']
    assert 'mass_coupling' not in registry()['design_controls']
    assert mass_ratio(DesignPoint(), ScenarioPoint(mass_coupling=1.6)) == pytest.approx(.4)


@pytest.mark.parametrize('field,value', [('left_scale', 1.01), ('left_scale', .09),
    ('period_s', .39), ('period_s', 1.51), ('coronary_conductance', 2.01),
    ('preservation_conductance', .51), ('left_contractility', 2.01)])
def test_out_of_bound_design_controls_are_rejected(field, value):
    with pytest.raises(ValueError):
        DesignPoint(**{field: value})


@pytest.mark.parametrize('field,value', [('extraction_limit', .8), ('mass_coupling', .7),
    ('preservation_fraction', .51), ('glycogen_mmol_L_liver', 500.),
    ('release_cap_mmol_min', .5), ('liver_function_fraction', 1.2)])
def test_out_of_bound_scenario_factors_are_rejected(field, value):
    with pytest.raises(ValueError):
        ScenarioPoint(**{field: value})


@pytest.mark.parametrize('bad', [True, float('nan'), float('inf'), '0.3'])
def test_nonnumeric_and_nonfinite_values_are_rejected(bad):
    with pytest.raises(ValueError):
        DesignPoint(left_scale=bad)


def test_per_side_contractility_is_not_pooled():
    d = DesignPoint(left_contractility=1.4, right_contractility=.8)
    assert contractility(d, 'left') == pytest.approx(1.4)
    assert contractility(d, 'right') == pytest.approx(.8)
    with pytest.raises(ValueError):
        contractility(d, 'septal')


def test_demand_assumptions_rejects_unvalidated_inputs():
    with pytest.raises(ValueError):
        mass_ratio(dict(left_scale=.25), ScenarioPoint())
    with pytest.raises(ValueError):
        demand_assumptions(DesignPoint(), dict(mass_coupling=1.), REFERENCE)
