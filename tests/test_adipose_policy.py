"""The adipose policy must move circulation and metabolism together."""
import json
from pathlib import Path
import pytest
from cardiac_redundancy.adipose_policy import (ADIPOSE_BRANCH, AdiposePolicy, apply,
                                               organ_entries)
from cardiac_redundancy.bilateral_support import BilateralSupport, LungAllocation

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT/'configs/m5-final-comparison.json').read_text(encoding='utf-8'))
BASE = json.loads((ROOT/'configs/m5-organ-preservation.json').read_text(
    encoding='utf-8'))['cases'][1]['profile']


def model(network='configs/m8-adipose-regional.yaml'):
    return BilateralSupport(ROOT/'configs/baseline.yaml', ROOT/network,
                            ROOT/'configs/m2-frozen-overlay.yaml', device='BIO',
                            allocation=LungAllocation(**CFG['lung_allocation']),
                            arterial_fraction=.2, inlet_resistance_fraction=.1)


def test_policy_requires_an_explicit_boolean():
    with pytest.raises(ValueError):
        AdiposePolicy(perfused=1)
    with pytest.raises(ValueError):
        AdiposePolicy(perfused='yes')


def test_suppressed_policy_leaves_the_depot_with_preserved_territories():
    factors = apply(model(), AdiposePolicy(False), preservation_conductance=.25)
    assert factors[ADIPOSE_BRANCH] == pytest.approx(.25)
    assert factors['remaining'] == pytest.approx(.25)
    assert factors['brain'] == pytest.approx(1.)


def test_perfused_policy_opens_the_depot_fully_and_nothing_else():
    suppressed = apply(model(), AdiposePolicy(False), preservation_conductance=.25)
    perfused = apply(model(), AdiposePolicy(True), preservation_conductance=.25)
    assert perfused[ADIPOSE_BRANCH] == pytest.approx(1.)
    changed = {k for k in perfused if perfused[k] != suppressed[k]}
    assert changed == {ADIPOSE_BRANCH}


def test_the_depot_resistance_is_unchanged_by_either_policy():
    """Policy opens a fraction of an existing bed; it does not redesign it."""
    net = json.loads((ROOT/'configs/m8-adipose-regional.yaml').read_text(encoding='utf-8'))
    before = next(b['R_mmHg_s_mL'] for b in net['branches'] if b['name'] == ADIPOSE_BRANCH)
    for perfused in (False, True):
        m = model()
        apply(m, AdiposePolicy(perfused))
        assert m.native.base.branch_by_name[ADIPOSE_BRANCH].resistance == pytest.approx(before)


def test_conductance_factor_matches_the_hydraulic_result():
    for perfused, expected in ((False, .25), (True, 1.)):
        policy = AdiposePolicy(perfused)
        assert policy.conductance_factor(.25) == pytest.approx(expected)
        factors = apply(model(), policy, preservation_conductance=.25)
        assert factors[ADIPOSE_BRANCH] == pytest.approx(policy.conductance_factor(.25))


def test_metabolic_entry_follows_the_same_flag():
    suppressed = organ_entries(BASE, AdiposePolicy(False), .25)[ADIPOSE_BRANCH]
    perfused = organ_entries(BASE, AdiposePolicy(True), .25)[ADIPOSE_BRANCH]
    assert suppressed == dict(role='preservation', fraction=.25, ramp_s=0.)
    assert perfused == dict(role='active', fraction=1., ramp_s=0.)


def test_entries_are_copied_not_mutated():
    out = organ_entries(BASE, AdiposePolicy(True), .25)
    assert ADIPOSE_BRANCH not in BASE
    assert out['brain'] == BASE['brain']
    out['brain']['fraction'] = 999.
    assert BASE['brain']['fraction'] != 999.


def test_both_sides_agree_for_every_policy():
    """The failure this guards: a depot perfused hydraulically but still
    carrying a suppressed demand, which no single result would reveal."""
    for perfused in (False, True):
        policy = AdiposePolicy(perfused)
        factors = apply(model(), policy, preservation_conductance=.25)
        entry = organ_entries(BASE, policy, .25)[ADIPOSE_BRANCH]
        hydraulically_full = factors[ADIPOSE_BRANCH] == pytest.approx(1.)
        metabolically_full = entry['role'] == 'active' and entry['fraction'] == 1.
        assert hydraulically_full == metabolically_full


def test_unvalidated_policy_refused():
    with pytest.raises(ValueError):
        apply(model(), dict(perfused=True))
    with pytest.raises(ValueError):
        organ_entries(BASE, 'perfused', .25)


def test_arm_a_network_rejects_the_depot_as_mandatory():
    """Arm A has no such branch, so the policy cannot be applied there."""
    with pytest.raises(ValueError):
        apply(model('configs/m3-evidence-regional.yaml'), AdiposePolicy(True))
