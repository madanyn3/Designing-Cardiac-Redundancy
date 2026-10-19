"""The evaluation cache key must carry everything that changes the trajectory.

Two defects of this class have already reached results: a release-time request
that silently returned a 0.5 s trajectory, and a 40 s pressure verdict paired
with a 10 s self-supply verdict from an earlier run. Both were invisible
because a stale cache answers instantly and looks like a legitimate hit.
"""
import inspect
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import optimize_m8_design as opt


def design(scale=.5, period=1., coronary=1.):
    return opt.design(scale, period, coronary)


def test_key_separates_every_design_control():
    base = opt.name_of(design(), .5)
    assert opt.name_of(design(scale=.6), .5) != base
    assert opt.name_of(design(period=.75), .5) != base
    assert opt.name_of(design(coronary=1.5), .5) != base


def test_key_separates_kappa_release_and_horizon():
    base = opt.name_of(design(), .5, .5, 10.)
    assert opt.name_of(design(), 1., .5, 10.) != base       # model form
    assert opt.name_of(design(), .5, 2.5, 10.) != base      # activation latency
    assert opt.name_of(design(), .5, .5, 40.) != base       # screen horizon


def test_ten_second_runs_keep_their_unsuffixed_identity():
    """Earlier results stay addressable; only other horizons gain a suffix."""
    assert opt.name_of(design(), .5, .5, 10.) == opt.name_of(design(), .5)
    assert opt.name_of(design(), .5, .5, 40.).endswith('-h40')


def test_coupled_stage_keys_on_the_same_identity_as_the_hydraulic_stage():
    """The defect this guards: pair_status keyed on the design alone, so a
    40 s pressure verdict was paired with a 10 s coupled verdict from an
    earlier run, silently and without error."""
    assert "PROTOCOL['evaluation_policy']['horizon_s']" in inspect.getsource(opt.pair_status)
    assert 'duration_s' in inspect.signature(opt.hydraulics).parameters
    assert 'duration_s' in inspect.signature(opt.name_of).parameters
