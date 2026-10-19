import json
from dataclasses import replace
from pathlib import Path
import pytest
from cardiac_redundancy.duration import SubstrateBudget, duration_sensitivity, supported_duration
from cardiac_redundancy.fuel import O2_ML_STPD_PER_MMOL_GLUCOSE

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = json.loads((ROOT/'configs/m4-fuel-reference.json').read_text(encoding='utf-8'))['parameters']


def test_budget_defaults_match_the_declared_fuel_reference():
    b = SubstrateBudget()
    assert b.release_cap_mmol_min == REFERENCE['release_cap_mmol_min']
    assert b.plasma_glucose_mmol_L == REFERENCE['initial_plasma_glucose_mmol_L']
    assert b.plasma_to_whole_blood_ratio == REFERENCE['plasma_to_whole_blood_ratio']
    assert b.glycogen_mmol_L_liver == REFERENCE['glycogen_mmol_L_liver']
    assert b.liver_volume_L == REFERENCE['liver_volume_L']


def test_pool_and_store_arithmetic():
    b = SubstrateBudget()
    assert b.circulating_pool_mmol == pytest.approx(5./1.11*5.25, rel=1e-12)
    assert b.glycogen_total_mmol == pytest.approx(229.*1.44, rel=1e-12)


def test_reference_optimum_duration():
    """The measured 109.5 mL O2/min demand of the T8.0 optimum."""
    r = supported_duration(109.5)
    assert r['glucose_demand_mmol_min'] == pytest.approx(109.5/O2_ML_STPD_PER_MMOL_GLUCOSE)
    assert r['net_drain_mmol_min'] == pytest.approx(109.5/O2_ML_STPD_PER_MMOL_GLUCOSE-.21)
    assert r['minutes'] == pytest.approx(39.1, abs=.2)
    assert 'circulating glucose pool' in r['limited_by']


def test_the_store_is_not_the_limit():
    r = supported_duration(109.5)
    assert r['store_limited_minutes'] > 1000        # over a day at the export cap
    assert r['minutes'] < r['store_limited_minutes']/30


def test_export_covering_demand_is_not_substrate_limited():
    covered = .21*O2_ML_STPD_PER_MMOL_GLUCOSE
    r = supported_duration(covered*.9)
    assert r['minutes'] is None and 'not substrate-limited' in r['limited_by']
    assert r['margin_mmol_min'] > 0
    assert supported_duration(covered*1.01)['minutes'] is not None


def test_a_larger_pump_shortens_duration():
    """Only the auxiliary term differs: 1.233, 1.885 and 3.110 mL over 10 s."""
    organs = (11.159+5.203)*6
    minutes = [supported_duration(organs+pump*6)['minutes'] for pump in (1.233, 1.885, 3.110)]
    assert minutes == sorted(minutes, reverse=True)
    assert minutes[0] == pytest.approx(41.1, abs=.3)
    assert minutes[2] == pytest.approx(35.9, abs=.3)


def test_higher_export_cap_lengthens_duration():
    b = SubstrateBudget()
    assert (supported_duration(109.5, replace(b, release_cap_mmol_min=.4))['minutes']
            > supported_duration(109.5, b)['minutes'])
    assert (supported_duration(109.5, replace(b, release_cap_mmol_min=.1))['minutes']
            < supported_duration(109.5, b)['minutes'])


def test_demand_is_the_stronger_lever_than_hepatic_export():
    """Measured elasticities: demand about -1.2, export about +0.36.

    Demand dominates because the glucose draw (0.814 mmol/min at the optimum)
    is roughly four times the export cap (0.21), so the pool is emptied mostly
    by demand rather than refilled by export. Anything that lowers total
    demand buys more duration than an equal fractional gain in export.
    """
    s = duration_sensitivity(109.5)
    assert s['demand_elasticity'] < 0 < s['export_elasticity']
    assert abs(s['demand_elasticity']) > abs(s['export_elasticity'])
    assert s['demand_elasticity'] == pytest.approx(-1.19, abs=.05)
    assert s['export_elasticity'] == pytest.approx(.36, abs=.05)


def test_accessible_fraction_scales_the_pool_only():
    half = supported_duration(109.5, SubstrateBudget(accessible_fraction=.5))
    full = supported_duration(109.5)
    assert half['minutes'] == pytest.approx(full['minutes']/2, rel=1e-9)


@pytest.mark.parametrize('bad', [-1., float('nan'), float('inf'), True, '109'])
def test_invalid_demand_rejected(bad):
    with pytest.raises(ValueError):
        supported_duration(bad)


@pytest.mark.parametrize('kwargs', [dict(blood_volume_L=0.), dict(plasma_to_whole_blood_ratio=0.),
                                    dict(accessible_fraction=1.5), dict(release_cap_mmol_min=-.1),
                                    dict(liver_volume_L=float('nan'))])
def test_invalid_budgets_rejected(kwargs):
    with pytest.raises(ValueError):
        SubstrateBudget(**kwargs)


def test_unvalidated_budget_refused():
    with pytest.raises(ValueError):
        supported_duration(109.5, dict(release_cap_mmol_min=.21))
