import pytest
from dataclasses import replace
from cardiac_redundancy.fat import (LipolysisParameters, O2_ML_STPD_PER_MMOL_FFA,
                                    fat_supply, fuel_split)
from cardiac_redundancy.duration import SubstrateBudget, supported_duration

#: Measured at the T8.0 optimum over the 60 s coupled run.
TOTAL, BRAIN = 108.62, 28.866
#: Visceral adipose reference flow, 1.3333 mL/s, at the two candidate roles.
FLOW_PRESERVED, FLOW_FULL = 1.3333*60*.25, 1.3333*60


def test_palmitate_stoichiometry():
    assert O2_ML_STPD_PER_MMOL_FFA == pytest.approx(23*22.414)
    assert O2_ML_STPD_PER_MMOL_FFA == pytest.approx(515.5, abs=.2)


def test_perfusion_binds_at_the_preserved_role():
    s = fat_supply(FLOW_PRESERVED)
    assert s['binding'] == 'adipose perfusion'
    assert s['supply_mmol_min'] == pytest.approx(FLOW_PRESERVED/1000.)
    assert s['supply_mL_O2_eq_min'] == pytest.approx(10.3, abs=.2)


def test_perfusion_still_binds_at_full_perfusion():
    """Even unsuppressed, transport is below lipolytic capacity."""
    s = fat_supply(FLOW_FULL)
    assert s['binding'] == 'adipose perfusion'
    assert s['perfusion_transport_mmol_min'] < s['lipolytic_capacity_mmol_min']
    assert s['supply_mL_O2_eq_min'] == pytest.approx(41.2, abs=.3)


def test_lipolytic_capacity_binds_only_at_implausibly_high_flow():
    p = LipolysisParameters()
    crossover = p.lipolytic_capacity_mmol_min*1000./p.arteriovenous_difference_mmol_L
    assert crossover == pytest.approx(189., abs=1.)       # mL/min
    assert fat_supply(crossover*1.2)['binding'] == 'lipolytic capacity'


def test_the_brain_share_is_a_floor_on_glucose_demand():
    plenty = dict(supply_mL_O2_eq_min=1e3)
    split = fuel_split(TOTAL, BRAIN, plenty)
    assert split['glucose_demand_mL_O2_min'] == pytest.approx(BRAIN)
    assert split['glucose_fraction'] == pytest.approx(BRAIN/TOTAL)
    assert split['glucose_fraction'] == pytest.approx(split['floor_fraction'])
    assert split['non_brain_covered_fraction'] == pytest.approx(1.)


def test_no_fat_reproduces_the_arm_a_assumption():
    split = fuel_split(TOTAL, BRAIN, dict(supply_mL_O2_eq_min=0.))
    assert split['glucose_fraction'] == pytest.approx(1.)
    assert split['glucose_demand_mL_O2_min'] == pytest.approx(TOTAL)


def test_the_two_adipose_roles_give_different_glucose_fractions():
    preserved = fuel_split(TOTAL, BRAIN, fat_supply(FLOW_PRESERVED))
    full = fuel_split(TOTAL, BRAIN, fat_supply(FLOW_FULL))
    assert preserved['non_brain_covered_fraction'] == pytest.approx(.13, abs=.02)
    assert full['non_brain_covered_fraction'] == pytest.approx(.52, abs=.02)
    assert full['glucose_fraction'] < preserved['glucose_fraction'] < 1.


def test_duration_responds_to_the_adipose_role():
    """The point of Arm B: perfusing fat buys supported duration."""
    budget = SubstrateBudget()
    def minutes(flow):
        split = fuel_split(TOTAL, BRAIN, fat_supply(flow))
        b = replace(budget, glucose_fraction=split['glucose_fraction'])
        return supported_duration(TOTAL, b)['minutes']
    arm_a = supported_duration(TOTAL, budget)['minutes']
    assert arm_a == pytest.approx(33.8, abs=.3)
    assert minutes(FLOW_PRESERVED) > arm_a
    assert minutes(FLOW_FULL) > minutes(FLOW_PRESERVED)


def test_arteriovenous_difference_scales_supply_linearly():
    p = LipolysisParameters()
    base = fat_supply(FLOW_FULL, p)['supply_mmol_min']
    doubled = fat_supply(FLOW_FULL, replace(p, arteriovenous_difference_mmol_L=2.))['supply_mmol_min']
    assert doubled == pytest.approx(2*base)


@pytest.mark.parametrize('bad', [-1., float('nan'), True, '80'])
def test_invalid_flow_rejected(bad):
    with pytest.raises(ValueError):
        fat_supply(bad)


@pytest.mark.parametrize('kwargs', [dict(body_mass_kg=0.), dict(accessible_fraction=1.4),
                                    dict(whole_body_capacity_umol_kg_min=-1.)])
def test_invalid_parameters_rejected(kwargs):
    with pytest.raises(ValueError):
        LipolysisParameters(**kwargs)


def test_brain_demand_cannot_exceed_total():
    with pytest.raises(ValueError):
        fuel_split(50., 80., dict(supply_mL_O2_eq_min=0.))


ORGANS = ['brain', 'aux_coronary_left', 'remaining', 'respiratory_muscle']
DEMANDS = [28.866, 7.190, 8.071, 5.000]


def test_zero_supply_reproduces_arm_a_exactly():
    from cardiac_redundancy.fat import glucose_shares
    g = glucose_shares(ORGANS, DEMANDS, dict(supply_mL_O2_eq_min=0.))
    assert g['shares'] == [1., 1., 1., 1.]
    assert g['covered_fraction'] == pytest.approx(0.)


def test_the_obligate_organ_never_leaves_glucose():
    from cardiac_redundancy.fat import glucose_shares
    for supply in (0., 10., 1e4):
        g = glucose_shares(ORGANS, DEMANDS, dict(supply_mL_O2_eq_min=supply))
        assert g['shares'][ORGANS.index('brain')] == pytest.approx(1.)
        assert g['obligate_organs'] == ['brain']


def test_abundant_supply_takes_non_obligate_organs_off_glucose():
    from cardiac_redundancy.fat import glucose_shares
    g = glucose_shares(ORGANS, DEMANDS, dict(supply_mL_O2_eq_min=1e4))
    assert g['covered_fraction'] == pytest.approx(1.)
    assert g['shares'][ORGANS.index('aux_coronary_left')] == pytest.approx(0.)
    assert g['met_by_fat_mL_O2_eq_min'] == pytest.approx(sum(DEMANDS)-28.866)


def test_partial_supply_is_shared_proportionally():
    from cardiac_redundancy.fat import glucose_shares
    non_obligate = sum(DEMANDS)-28.866
    g = glucose_shares(ORGANS, DEMANDS, dict(supply_mL_O2_eq_min=non_obligate/2))
    assert g['covered_fraction'] == pytest.approx(.5)
    for name in ('aux_coronary_left', 'remaining', 'respiratory_muscle'):
        assert g['shares'][ORGANS.index(name)] == pytest.approx(.5)


def test_the_pump_group_benefits_which_is_the_arm_b_prediction():
    """Auxiliary self-supply binds Arm A; myocardium prefers fatty acid."""
    from cardiac_redundancy.fat import fat_supply, glucose_shares
    perfused = glucose_shares(ORGANS, DEMANDS, fat_supply(53.6))
    suppressed = glucose_shares(ORGANS, DEMANDS, fat_supply(53.6*.25))
    i = ORGANS.index('aux_coronary_left')
    assert perfused['shares'][i] < suppressed['shares'][i] < 1.


def test_shares_are_valid_for_fuel_oxygen_transport():
    from cardiac_redundancy.fat import fat_supply, glucose_shares
    g = glucose_shares(ORGANS, DEMANDS, fat_supply(53.6))
    assert all(0. <= s <= 1. for s in g['shares'])
    assert len(g['shares']) == len(ORGANS)


def test_mismatched_lengths_and_bad_values_rejected():
    from cardiac_redundancy.fat import glucose_shares
    with pytest.raises(ValueError):
        glucose_shares(ORGANS, DEMANDS[:2], dict(supply_mL_O2_eq_min=0.))
    with pytest.raises(ValueError):
        glucose_shares(ORGANS, [-1.]+DEMANDS[1:], dict(supply_mL_O2_eq_min=0.))
    with pytest.raises(ValueError):
        glucose_shares(ORGANS, DEMANDS, dict(supply_mL_O2_eq_min=float('nan')))
