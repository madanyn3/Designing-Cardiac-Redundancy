import numpy as np
import pytest
from pathlib import Path
import json
from cardiac_redundancy.feasibility import (CRITERIA, FeasibilityPolicy, PRESSURE_SCREEN_MMHG,
    UNRESOLVED_REQUIREMENTS, numerical_status, oxygen_endpoint, pressure_endpoint,
    pressure_margin, verdict)

POLICY = FeasibilityPolicy()
ROOT = Path(__file__).resolve().parents[1]


def test_unresolved_requirements_match_the_declared_policy_file():
    declared = json.loads((ROOT/'configs/m4-constraint-policy.json').read_text(
        encoding='utf-8'))['unresolved_requirements']
    assert list(UNRESOLVED_REQUIREMENTS) == declared
    assert POLICY.unresolved_requirements == UNRESOLVED_REQUIREMENTS


def test_a_clean_numeric_result_cannot_be_reported_as_a_pass():
    t, p = ao(95.)
    time, g = groups(0.)
    out = verdict(POLICY, 'valid', pressure_endpoint(POLICY, t, p),
                  oxygen_endpoint(POLICY, time, g, 255.848))
    assert out['status'] == 'indeterminate'
    assert out['binding_constraints'] == []
    assert len(out['unresolved_requirements']) == len(UNRESOLVED_REQUIREMENTS)


def ao(level, duration=10., step=.01, startup_level=None, allowance=2.):
    """Constant-pressure record, optionally lower during the startup interval."""
    t = np.arange(0., duration+step/2, step)
    p = np.where(t < allowance, startup_level if startup_level is not None else level, level)
    return t, p.astype(float)


def groups(unmet_rate, duration=10., step=.05, start_s=0.):
    t = np.arange(0., duration+step/2, step)
    demand = np.full_like(t, .5)
    unmet = np.where(t >= start_s, unmet_rate, 0.)
    uptake = demand-unmet
    deficit = np.r_[0., np.cumsum((unmet[:-1]+unmet[1:])*np.diff(t)/2)]
    return t, dict(active=dict(demand=demand, uptake=uptake, delivered=demand,
                               cumulative_deficit=deficit))


def test_declared_screen_is_the_joint_triple():
    assert PRESSURE_SCREEN_MMHG == dict(mean=65., maximum=90., minimum=50.)
    assert POLICY.thresholds_mmHg == PRESSURE_SCREEN_MMHG


def test_a_design_passing_the_mean_but_failing_SBP_is_not_satisfied():
    t, p = ao(70.)
    assert pressure_endpoint(POLICY, t, p)['status'] == 'violated'


def test_joint_triple_satisfied_only_when_every_metric_clears():
    t = np.arange(0., 10.001, .01)
    p = 75.+20.*np.sin(2*np.pi*t/.75)
    out = pressure_endpoint(POLICY, t, p)
    assert out['status'] == 'satisfied' and out['first_violating_sustained_window'] is None
    assert out['sustained_windows'] == 8 and out['startup_windows'] == 2


def test_startup_exposure_is_reported_but_does_not_decide_the_verdict():
    t, p = ao(70., startup_level=20.)
    out = pressure_endpoint(POLICY, t, p)
    # Two startup windows at 20 mmHg against the 65 mmHg mean screen; the
    # tolerance absorbs the single interpolated sample at the boundary.
    assert out['startup_shortfall_mmHg_s']['mean'] == pytest.approx(90., abs=.5)
    assert out['status'] == 'violated'  # 70 mmHg still fails the SBP screen
    t, p = ao(95., startup_level=20.)
    out = pressure_endpoint(POLICY, t, p)
    assert out['status'] == 'satisfied'
    assert out['startup_shortfall_mmHg_s']['minimum'] == pytest.approx(60.)


def test_late_release_does_not_earn_a_longer_allowance():
    t, p = ao(95., startup_level=20., allowance=4.)
    early = pressure_endpoint(POLICY, t, p, release_s=.5)
    late = pressure_endpoint(POLICY, t, p, release_s=3.5)
    assert early['status'] == late['status'] == 'violated'
    assert early['startup_windows'] == late['startup_windows'] == 2
    assert early['sustained_exposure']['mean']['violating_window_duration_s'] == pytest.approx(2.)


def test_incomplete_horizon_is_unassessed_not_satisfied():
    t = np.arange(0., 2.51, .01)
    p = np.full_like(t, 95.)
    policy = FeasibilityPolicy(startup_allowance_s=1., horizon_s=2.5)
    assert pressure_endpoint(policy, t, p)['status'] == 'unassessed'


def test_numerical_failure_and_domain_exit_stay_separate_from_physiology():
    assert numerical_status(None, [.001, .002], .005) == 'valid'
    assert numerical_status(None, [.01], .005) == 'numerical_failure'
    assert numerical_status(None, [float('nan')], .005) == 'numerical_failure'
    assert numerical_status({'reason': 'PU'}, [.0], .005) == 'model_domain_limit'


def test_a_domain_exit_cannot_be_credited_as_a_pass():
    t, p = ao(95.)
    time, g = groups(0.)
    out = verdict(POLICY, numerical_status({'reason': 'PU'}, [.0], .005),
                  pressure_endpoint(POLICY, t, p), oxygen_endpoint(POLICY, time, g, 255.848))
    assert out['status'] == 'model_domain_limit'


def test_sustained_oxygen_shortfall_binds_but_startup_shortfall_is_separate():
    time, g = groups(.02, start_s=0.)
    out = oxygen_endpoint(POLICY, time, g, 255.848)['active']
    assert out['status'] == 'violated' and out['sustained']['status'] == 'violated'
    assert out['startup_unmet_mL'] > 0
    time, g = groups(.02, start_s=0.)
    g['active']['uptake'] = np.where(time < 2., g['active']['demand']-.02, g['active']['demand'])
    unmet = g['active']['demand']-g['active']['uptake']
    g['active']['cumulative_deficit'] = np.r_[0., np.cumsum((unmet[:-1]+unmet[1:])*np.diff(time)/2)]
    out = oxygen_endpoint(POLICY, time, g, 255.848)['active']
    assert out['status'] == 'violated' and out['sustained']['status'] == 'satisfied'
    # 0.02 mL/s for 1.95 s plus the half-interval in which the fixture's step
    # returns to zero exactly at the 2 s allowance boundary.
    assert out['startup_unmet_mL'] == pytest.approx(.0395, rel=1e-9)


def test_a_window_counts_once_however_many_metrics_it_violates():
    t, p = ao(40.)
    out = pressure_endpoint(POLICY, t, p)
    assert out['violating_sustained_windows'] == out['sustained_windows'] == 8
    assert sum(w['interval'] == 'startup' for w in out['windows']) == 2
    assert out['first_satisfied_sustained_window_s'] is None
    t = np.arange(0., 10.001, .01)
    p = np.where(t < 5., 40., 75.+20.*np.sin(2*np.pi*t/.75))
    out = pressure_endpoint(POLICY, t, p)
    assert out['violating_sustained_windows'] == 3
    assert out['first_satisfied_sustained_window_s'] == pytest.approx(5.)
    assert out['sustained_exposure']['minimum']['violating_window_duration_s'] == pytest.approx(3.)


def test_verdict_names_every_binding_constraint():
    t, p = ao(70.)
    time, g = groups(.02)
    out = verdict(POLICY, 'valid', pressure_endpoint(POLICY, t, p),
                  oxygen_endpoint(POLICY, time, g, 255.848), reported=dict(energy_J=1.5))
    assert out['status'] == 'violated_research_screen'
    assert out['binding_constraints'] == ['pressure', 'oxygen:active']
    assert out['reported']['energy_J'] == 1.5
    assert out['biological_recovery'] == 'unassessed'


def test_unresolved_requirements_prevent_a_clean_pass():
    t, p = ao(95.)
    time, g = groups(0.)
    policy = FeasibilityPolicy(unresolved_requirements=('regional organ pressure minima',))
    out = verdict(policy, 'valid', pressure_endpoint(policy, t, p),
                  oxygen_endpoint(policy, time, g, 255.848))
    assert out['status'] == 'indeterminate'
    assert out['unresolved_requirements'] == ['regional organ pressure minima']
    bare = FeasibilityPolicy(unresolved_requirements=())
    clean = verdict(bare, 'valid', pressure_endpoint(bare, t, p),
                    oxygen_endpoint(bare, time, g, 255.848))
    assert clean['status'] == 'satisfied_within_declared_model'
    assert 'not support duration, survival or recovery' in clean['meaning']


@pytest.mark.parametrize('kwargs', [
    dict(thresholds_mmHg=dict(mean=65.)),
    dict(thresholds_mmHg=dict(mean=95., maximum=90., minimum=50.)),
    dict(window_s=0.), dict(window_s=-1.), dict(oxygen_tolerance_mL_s=0.),
    dict(startup_allowance_s=1.5), dict(startup_allowance_s=-1.),
    dict(startup_allowance_s=10.), dict(horizon_s=float('inf'))])
def test_invalid_policies_are_rejected(kwargs):
    with pytest.raises(ValueError):
        FeasibilityPolicy(**kwargs)


def test_unvalidated_policy_objects_are_refused():
    t, p = ao(95.)
    with pytest.raises(ValueError):
        pressure_endpoint(dict(window_s=1.), t, p)
    with pytest.raises(ValueError):
        oxygen_endpoint(POLICY, *groups(0.)[::-1], 255.848)
    with pytest.raises(ValueError):
        oxygen_endpoint(POLICY, groups(0.)[0], {}, 255.848)


def test_margin_sign_agrees_with_the_screen_verdict():
    """A root-find on the margin must locate the boundary the screen reports."""
    for level in (40., 60., 70., 80., 95., 110.):
        t = np.arange(0., 10.001, .01)
        p = level+20.*np.sin(2*np.pi*t/.75)
        out = pressure_endpoint(POLICY, t, p)
        margin = pressure_margin(out)
        assert (margin >= 0) == (out['status'] == 'satisfied'), (level, margin, out['status'])


def test_margin_is_the_worst_metric_not_the_mean():
    t = np.arange(0., 10.001, .01)
    p = 75.+20.*np.sin(2*np.pi*t/.75)          # MAP 75, SBP 95, DBP 55
    out = pressure_endpoint(POLICY, t, p)
    # DBP has 5 mmHg of room, MAP has 10, SBP has 5; the worst governs.
    assert pressure_margin(out) == pytest.approx(min(75.-65., 95.-90., 55.-50.), abs=.3)


def test_margin_is_monotone_in_pressure_level():
    margins = []
    for level in (60., 70., 80., 90.):
        t = np.arange(0., 10.001, .01)
        margins.append(pressure_margin(pressure_endpoint(POLICY, t, level+20.*np.sin(2*np.pi*t/.75))))
    assert margins == sorted(margins)


def test_margin_is_none_without_a_complete_sustained_window():
    t = np.arange(0., 2.51, .01)
    policy = FeasibilityPolicy(startup_allowance_s=2., horizon_s=2.5)
    assert pressure_margin(pressure_endpoint(policy, t, np.full_like(t, 95.))) is None


CYCLE = FeasibilityPolicy(window_mode='cycle')
MAP_FIRST = FeasibilityPolicy(criterion='map_primary')


def pulsatile(level, amplitude, period=.75, duration=10., step=.005):
    t = np.arange(0., duration+step/2, step)
    return t, level+amplitude*np.sin(2*np.pi*t/period)


def test_criteria_are_declared_and_nested():
    assert CRITERIA['strict_triple'] == ('mean', 'maximum', 'minimum')
    assert CRITERIA['map_primary'] == ('mean',)
    assert set(CRITERIA['map_primary']) < set(CRITERIA['strict_triple'])


def test_both_criteria_are_always_reported():
    t, p = pulsatile(70., 5.)                      # MAP 70 ok, SBP 75 fails
    out = pressure_endpoint(POLICY, t, p)
    assert out['status_by_criterion'] == dict(strict_triple='violated', map_primary='satisfied')
    assert out['status'] == 'violated' and out['criterion'] == 'strict_triple'
    relaxed = pressure_endpoint(MAP_FIRST, t, p)
    assert relaxed['status'] == 'satisfied' and relaxed['enforced_metrics'] == ['mean']
    assert relaxed['status_by_criterion'] == out['status_by_criterion']


def test_margin_follows_the_enforced_criterion():
    t, p = pulsatile(70., 5.)
    assert pressure_margin(pressure_endpoint(POLICY, t, p)) < 0        # SBP shortfall binds
    assert pressure_margin(pressure_endpoint(MAP_FIRST, t, p)) > 0     # MAP has room


def test_cycle_mode_aligns_windows_to_release_and_uses_the_period():
    t, p = pulsatile(95., 20., period=.6)
    out = pressure_endpoint(CYCLE, t, p, release_s=1., period_s=.6)
    cycles = [w for w in out['windows'] if w['start_s'] >= 1.-1e-9]
    assert out['window_mode'] == 'cycle' and out['cycle_period_s'] == pytest.approx(.6)
    assert cycles[0]['start_s'] == pytest.approx(1.)
    assert all(w['end_s']-w['start_s'] == pytest.approx(.6) for w in cycles)
    # the nonpulsatile pre-release interval keeps the declared fixed window
    early = [w for w in out['windows'] if w['start_s'] < 1.-1e-9]
    assert early and all(w['end_s']-w['start_s'] == pytest.approx(1.) for w in early)


def test_cycle_mode_requires_the_period():
    t, p = pulsatile(95., 20.)
    with pytest.raises(ValueError):
        pressure_endpoint(CYCLE, t, p, release_s=1.)


def test_cycle_mode_without_release_falls_back_to_fixed_windows():
    t, p = pulsatile(95., 20.)
    out = pressure_endpoint(CYCLE, t, p)
    assert all(w['end_s']-w['start_s'] == pytest.approx(1.) for w in out['windows'])


def test_exposure_uses_each_window_own_width():
    """Once widths differ, a violating cycle must contribute its own duration."""
    t, p = pulsatile(40., 5., period=.5)
    out = pressure_endpoint(CYCLE, t, p, release_s=1., period_s=.5)
    sustained = [w for w in out['windows'] if w['interval'] == 'sustained']
    assert out['sustained_exposure']['mean']['violating_window_duration_s'] == pytest.approx(
        .5*len(sustained), abs=1e-9)


@pytest.mark.parametrize('kwargs', [dict(window_mode='beat'), dict(criterion='mean_only')])
def test_invalid_mode_or_criterion_rejected(kwargs):
    with pytest.raises(ValueError):
        FeasibilityPolicy(**kwargs)


def test_a_leftover_partial_cycle_does_not_void_a_completed_run():
    """Cycle windows start at release, so a remainder is the norm, not a gap."""
    t, p = pulsatile(95., 20., period=1.)
    out = pressure_endpoint(CYCLE, t, p, release_s=.5, period_s=1.)
    assert out['incomplete_tail_s'] is not None          # 0.5 s of cycle left over
    assert out['unassessed_remainder_s'] == pytest.approx(.5, abs=.02)
    assert out['remainder_within_one_window'] is True
    assert out['violating_sustained_windows'] == 0
    assert out['status'] == 'satisfied'


def test_a_remainder_larger_than_one_window_is_still_unassessed():
    t = np.arange(0., 4.01, .005)
    p = np.full_like(t, 95.)
    policy = FeasibilityPolicy(window_mode='cycle', startup_allowance_s=1., horizon_s=4.)
    out = pressure_endpoint(policy, t, p, release_s=.5, period_s=1.)
    # Truncate the record so more than one cycle is left unassessed.
    short = pressure_endpoint(policy, t[t <= 2.6], p[t <= 2.6], release_s=.5, period_s=1.)
    assert out['status'] == 'satisfied'
    assert short['unassessed_remainder_s'] <= 1.+1e-9


def test_plateau_check_reports_a_settled_transient():
    t, p = pulsatile(95., 20., period=1., duration=12.)
    out = pressure_endpoint(CYCLE, t, p, release_s=.5, period_s=1.)
    assert out['plateau_reached'] is True
    assert out['plateau_spread_mmHg'] == pytest.approx(0., abs=.05)
    assert len(out['final_cycle_means_mmHg']) == 5


def test_plateau_check_flags_a_still_declining_transient():
    t = np.arange(0., 12.001, .005)
    decaying = 95.-10.*(1.-np.exp(-t/6.))          # still falling at the horizon
    p = decaying+20.*np.sin(2*np.pi*t/1.)
    out = pressure_endpoint(CYCLE, t, p, release_s=.5, period_s=1.)
    assert out['plateau_reached'] is False
    assert out['plateau_spread_mmHg'] > .05


def test_plateau_tolerance_is_validated():
    with pytest.raises(ValueError):
        FeasibilityPolicy(plateau_tolerance_mmHg=0.)
    with pytest.raises(ValueError):
        FeasibilityPolicy(plateau_tolerance_mmHg=-1.)
