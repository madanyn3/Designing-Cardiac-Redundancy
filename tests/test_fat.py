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
