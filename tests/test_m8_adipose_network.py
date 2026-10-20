"""The Arm B network must differ from Arm A by exactly one declared change.

Disaggregation, not addition: visceral adipose comes out of the residual
branch, so total cardiac output, the baseline fraction sum and the total oxygen
demand are all unchanged. Adding a branch instead would inflate the closed
M1/M3 accounting, and the difference would be invisible in any single result.
"""
import json
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
FROZEN = json.loads((ROOT/'configs/m3-evidence-regional.yaml').read_text(encoding='utf-8'))
ARM_B = json.loads((ROOT/'configs/m8-adipose-regional.yaml').read_text(encoding='utf-8'))
FROZEN_MAP = json.loads((ROOT/'configs/m4-organ-mapping.json').read_text(encoding='utf-8'))
ARM_B_MAP = json.loads((ROOT/'configs/m8-adipose-organ-mapping.json').read_text(encoding='utf-8'))
BUDGET = json.loads((ROOT/'data/m1-demand-budget.json').read_text(encoding='utf-8'))['baseline_mL_min']


def ao(net):
    return [b for b in net['branches'] if b['from'] == 'AO']


def test_the_frozen_network_is_untouched():
    """Arm A must stay reproducible; Arm B is a separate file."""
    assert 'visceral_adipose' not in {b['name'] for b in FROZEN['branches']}
    assert len(FROZEN['branches']) == 23


def test_exactly_one_branch_is_added_and_it_is_a_split():
    names = {b['name'] for b in ARM_B['branches']}
    assert names-{b['name'] for b in FROZEN['branches']} == {'visceral_adipose'}
    assert len(ARM_B['branches']) == len(FROZEN['branches'])+1


def test_total_flow_and_fraction_sum_are_preserved():
    assert sum(b['target_flow_mL_s'] for b in ao(ARM_B)) == pytest.approx(
        sum(b['target_flow_mL_s'] for b in ao(FROZEN)), abs=1e-9)
    assert sum(b['baseline_fraction'] for b in ao(ARM_B)) == pytest.approx(1., abs=1e-12)


def test_only_the_residual_branch_changed():
    frozen = {b['name']: b for b in FROZEN['branches']}
    for b in ARM_B['branches']:
        if b['name'] in ('visceral_adipose', 'remaining'):
            continue
        assert b == frozen[b['name']], b['name']


def test_the_residual_branch_lost_exactly_what_adipose_gained():
    f = {b['name']: b for b in FROZEN['branches']}
    a = {b['name']: b for b in ARM_B['branches']}
    moved = f['remaining']['target_flow_mL_s']-a['remaining']['target_flow_mL_s']
    assert moved == pytest.approx(a['visceral_adipose']['target_flow_mL_s'], abs=1e-12)
    moved_fraction = f['remaining']['baseline_fraction']-a['remaining']['baseline_fraction']
    assert moved_fraction == pytest.approx(a['visceral_adipose']['baseline_fraction'], abs=1e-15)


def test_resistance_construction_survives_for_branches_sharing_endpoints():
    """R = (P_from - P_to)/Q; liver_arterial drains to the hepatic node and has
    its own driving pressure, so only AO-to-VC branches share this one."""
    drive = {b['R_mmHg_s_mL']*b['target_flow_mL_s']
             for b in ARM_B['branches'] if b['from'] == 'AO' and b['to'] == 'VC'}
    assert max(drive)-min(drive) < 1e-6
    reference = json.loads((ROOT/'configs/m3-evidence-regional.yaml').read_text(
        encoding='utf-8'))['regional_reference']
    assert next(iter(drive)) == pytest.approx(
        reference['P_AO_mmHg']-reference['P_VC_mmHg'], abs=1e-6)


def test_adipose_flow_follows_the_declared_mass_and_perfusion():
    d = ARM_B['arm_b_disaggregation']
    expected = d['mass_kg']*10.*d['flow_mL_per_100g_min']/60.
    assert d['visceral_flow_mL_s'] == pytest.approx(expected, rel=1e-12)
    assert 3. <= d['flow_mL_per_100g_min'] <= 5.          # the reported fasting range


def test_demand_comes_from_the_projects_own_adipose_ratio():
    """Not an imported specific metabolic rate: the ratio is the project's."""
    d = ARM_B['arm_b_disaggregation']
    limb = next(b for b in FROZEN['branches'] if b['name'] == 'limb_adipose')
    assert d['demand_per_flow_mL_O2_min_per_mL_s'] == pytest.approx(
        BUDGET['limb_adipose']/limb['target_flow_mL_s'], rel=1e-12)
    assert d['visceral_demand_mL_O2_min'] == pytest.approx(
        d['demand_per_flow_mL_O2_min_per_mL_s']*d['visceral_flow_mL_s'], rel=1e-12)


def test_mapping_weights_still_sum_to_one_and_only_residual_changed():
    for entry, weights in ARM_B_MAP['mapping'].items():
        assert sum(weights.values()) == pytest.approx(1., abs=1e-12), entry
        if entry != 'remaining':
            assert weights == FROZEN_MAP['mapping'][entry], entry
    arm_b = ARM_B_MAP['mapping']['remaining']
    assert 'visceral_adipose' in arm_b
    moved = FROZEN_MAP['mapping']['remaining']['remaining']-arm_b['remaining']
    assert moved == pytest.approx(arm_b['visceral_adipose'], abs=1e-12)


def test_total_oxygen_demand_is_unchanged_by_the_split():
    entry = BUDGET['remaining']
    before = sum(entry*w for w in FROZEN_MAP['mapping']['remaining'].values())
    after = sum(entry*w for w in ARM_B_MAP['mapping']['remaining'].values())
    assert after == pytest.approx(before, rel=1e-12)
