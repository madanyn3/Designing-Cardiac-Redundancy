"""The seams between Arm A and Arm B must agree before any expensive run.

Arm A's runners were written around one regional network, one organ mapping and
one preservation profile. Arm B varies all three, and four integration defects
appeared in quick succession because each seam was discovered by crashing a
forty-minute run rather than by a check costing milliseconds:

  * `flow_names` read in a scope that did not have the model;
  * the coupled cache keyed without the horizon;
  * the network path hardcoded while the mapping was swapped;
  * the depot present in the mapping but absent from the profile.

These tests build both arms' models and assert the three artefacts are mutually
consistent, so a mismatch fails in seconds and names itself.
"""
import json
from pathlib import Path
import pytest
from cardiac_redundancy.adipose_policy import ADIPOSE_BRANCH, AdiposePolicy
from cardiac_redundancy.bilateral_support import BilateralSupport, LungAllocation
from cardiac_redundancy.organ_oxygen import map_demands

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT/'configs/m5-final-comparison.json').read_text(encoding='utf-8'))
BUDGET = json.loads((ROOT/'data/m1-demand-budget.json').read_text(encoding='utf-8'))['baseline_mL_min']
PROFILE = json.loads((ROOT/'configs/m5-organ-preservation.json').read_text(
    encoding='utf-8'))['cases'][1]['profile']

ARMS = {
    'A': ('configs/m3-evidence-regional.yaml', 'configs/m4-organ-mapping.json'),
    'B': ('configs/m8-adipose-regional.yaml', 'configs/m8-adipose-organ-mapping.json'),
}


def build(network):
    return BilateralSupport(ROOT/'configs/baseline.yaml', ROOT/network,
                            ROOT/'configs/m2-frozen-overlay.yaml', device='BIO',
                            allocation=LungAllocation(**CFG['lung_allocation']),
                            arterial_fraction=.2, inlet_resistance_fraction=.1)


def organs_for(arm):
    network, mapping_path = ARMS[arm]
    model = build(network)
    mapping = json.loads((ROOT/mapping_path).read_text(encoding='utf-8'))['mapping']
    share = CFG['lung_allocation']['supported_share']
    mapping = dict(mapping, lungs={'lung_supported': share, 'lung_preserved': 1-share})
    return model, map_demands(BUDGET, mapping, model.flow_names, .75, 10)


@pytest.mark.parametrize('arm', sorted(ARMS))
def test_every_mapped_route_exists_in_the_hydraulic_model(arm):
    """The defect this guards: Arm B's mapping named a flow Arm A's model
    lacked, which surfaced only part way into a coupled replay."""
    model, organs = organs_for(arm)
    assert organs, arm
    for organ in organs:
        assert organ.name in model.flow_names or organ.name.startswith('lung'), organ.name


@pytest.mark.parametrize('arm', sorted(ARMS))
def test_every_organ_has_a_preservation_profile_entry(arm):
    """Arm B must supply the depot's entry; Arm A must need no additions."""
    _, organs = organs_for(arm)
    entries = dict(PROFILE)
    if arm == 'B':
        entries[ADIPOSE_BRANCH] = AdiposePolicy(False).organ_entry(.25)
    missing = [o.name for o in organs
               if o.name not in entries and not o.name.startswith('lung')]
    assert missing == [], (arm, missing)


def test_arm_a_needs_no_depot_entry():
    _, organs = organs_for('A')
    assert ADIPOSE_BRANCH not in {o.name for o in organs}


def test_arm_b_adds_exactly_the_depot():
    _, a = organs_for('A')
    _, b = organs_for('B')
    assert {o.name for o in b}-{o.name for o in a} == {ADIPOSE_BRANCH}


def test_total_organ_demand_is_identical_across_arms():
    """Disaggregation moves demand between organs; it must not create any."""
    _, a = organs_for('A')
    _, b = organs_for('B')
    assert sum(o.demand_mL_min for o in b) == pytest.approx(
        sum(o.demand_mL_min for o in a), rel=1e-9)


def test_glucose_shares_line_up_with_the_organ_order():
    """The share array is positional; a mismatch would silently mis-assign fuel."""
    from cardiac_redundancy.fat import fat_supply, glucose_shares
    _, organs = organs_for('B')
    names = [o.name for o in organs]
    shares = glucose_shares(names, [o.demand_mL_min for o in organs], fat_supply(53.6))
    assert shares['organ_names'] == names
    assert len(shares['shares']) == len(organs)
    assert shares['shares'][names.index('brain')] == pytest.approx(1.)


def test_the_runners_do_not_hardcode_the_regional_network():
    """Arm B passes its own network through the case; the default is frozen."""
    import inspect, sys
    sys.path.insert(0, str(ROOT/'scripts'))
    import run_m7_global_metabolism as met
    source = inspect.getsource(met.network)
    assert "case.get('network'" in source
    assert 'configs/m3-evidence-regional.yaml' in source     # the default survives
