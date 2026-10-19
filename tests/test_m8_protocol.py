"""The frozen protocol file must keep agreeing with the validated code."""
import json
from dataclasses import asdict
from pathlib import Path
import pytest
from cardiac_redundancy.design import BOUNDS, DesignPoint, ScenarioPoint, registry
from cardiac_redundancy.feasibility import FeasibilityPolicy, PRESSURE_SCREEN_MMHG

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = json.loads((ROOT/'configs/m8-protocol.json').read_text(encoding='utf-8'))


def test_registry_and_bounds_match_the_validators():
    assert PROTOCOL['registry'] == registry()
    assert PROTOCOL['declared_bounds'] == {k: list(v[0]) for k, v in BOUNDS.items()}


def test_declared_policy_is_constructible_and_states_the_v3_choices():
    declared = dict(PROTOCOL['evaluation_policy'])
    declared['unresolved_requirements'] = tuple(declared['unresolved_requirements'])
    policy = FeasibilityPolicy(**declared)
    assert policy.thresholds_mmHg == PRESSURE_SCREEN_MMHG
    # v3 deliberately departs from the stricter dataclass defaults, and says so.
    assert policy.window_mode == 'cycle' and policy.criterion == 'map_primary'
    assert FeasibilityPolicy().window_mode == 'fixed'
    assert FeasibilityPolicy().criterion == 'strict_triple'


def test_v3_amendments_record_the_evidence_and_the_withdrawal():
    changes = PROTOCOL['amendments_v3']
    assert len(changes) == 3
    criterion = next(c for c in changes if 'map_primary' in c['change'])
    assert any('10.1161/CIR.0000000000001375' in s for s in criterion['sources'])
    assert 'withdrawn' in criterion['withdrawn_claim']
    assert any('cycle-aligned' in c['change'] for c in changes)
    assert PROTOCOL['primary_endpoint'].startswith('supported duration')


def test_reference_points_are_the_verified_fixture():
    assert PROTOCOL['reference_design'] == asdict(DesignPoint())
    assert PROTOCOL['reference_scenario'] == asdict(ScenarioPoint())


@pytest.mark.parametrize('name', ['reference', 'adverse', 'favourable'])
def test_every_declared_scenario_is_constructible(name):
    point = ScenarioPoint(**PROTOCOL['scenarios'][name])
    for field, value in PROTOCOL['scenarios'][name].items():
        lo, hi = BOUNDS[field][0]
        assert lo <= value <= hi
    assert isinstance(point, ScenarioPoint)


def test_scenario_factors_are_never_design_controls():
    for name, entries in PROTOCOL['scenarios'].items():
        assert not set(entries) & set(PROTOCOL['registry']['design_controls']), name


def test_operating_cases_declare_a_policy_and_a_meaning():
    assert set(PROTOCOL['operating_cases']) >= {
        'no_support', 'support_no_suppression', 'support_with_suppression',
        'ventilation_degraded', 'delayed_activation', 'partial_recovery', 'shared_failure'}
    for name, case in PROTOCOL['operating_cases'].items():
        assert case['policy'] in ('none', 'fixed', 'observed'), name
        assert case['meaning'], name


def test_budget_is_derived_from_a_measured_cost():
    b = PROTOCOL['budget']
    assert b['measured_single_evaluation_s'] > 0 and b['measured_hydraulic_s'] > 0
    for stage in ('stage1', 'stage2'):
        expected = b[stage]['evaluations']*(b['measured_hydraulic_s'] if stage == 'stage1'
                                            else b['measured_single_evaluation_s'])/b['workers']/3600
        assert b[stage]['wall_hours'] == pytest.approx(expected, abs=.01)


def test_numerical_settings_match_the_verified_runners():
    h = PROTOCOL['numerical']['hydraulic']
    assert [(p['rtol'], p['max_step_s']) for p in h['pairs']] == [(1e-7, .01), (1e-9, .005)]
    m = PROTOCOL['numerical']['metabolic']
    assert [(p['rtol'], p['max_step_s']) for p in m['pairs']] == [(1e-9, .005), (1e-11, .0025)]


def test_code_version_and_inputs_are_recorded():
    assert PROTOCOL['code']['commit'] and len(PROTOCOL['code']['commit']) == 40
    assert len(PROTOCOL['inputs']) >= 18
    assert all(len(d) == 64 for d in PROTOCOL['inputs'].values())


def test_inherited_limitations_and_forbidden_moves_stay_declared():
    text = ' '.join(PROTOCOL['inherited_limitations'])
    assert 'L-M4-F01' in text and 'R-M3-01' in text and 'supply-limited' in text
    assert any('efficiency' in f for f in PROTOCOL['objective']['forbidden'])
    assert 'never as indefinite support' in PROTOCOL['objective']['censoring']
