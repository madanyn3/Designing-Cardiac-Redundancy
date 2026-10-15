"""Analytic support characteristics, integration regression and conservation."""
import json
from pathlib import Path
import numpy as np
import pytest
from cardiac_redundancy.backup_pump import BackupPumpModel, PumpParameters, PumpStage
from cardiac_redundancy.whole_network import WholeNetworkModel
from cardiac_redundancy.integration import integrate_source

ROOT = Path(__file__).resolve().parents[1]
ARGS = (ROOT/'configs/baseline.yaml', ROOT/'configs/m3-evidence-regional.yaml', ROOT/'configs/m2-frozen-overlay.yaml')
KW = dict(arterial_fraction=.2, inlet_resistance_fraction=.1)


def test_default_exact_regression_with_support_and_regional_ramps():
    old = WholeNetworkModel(*ARGS, **KW)
    new = BackupPumpModel(*ARGS, **KW, pump_parameters=PumpParameters.from_json(ROOT/'configs/m3-backup-pump.json'))
    for m in (old, new):
        m.set_controls(time=0, ramp_s=.1, pump_on=True, right_head=20, left_head=100)
        m.set_regional_controls(time=0, ramp_s=.1, factors={'kidney_other':.25, 'small_intestine':0})
    for active in (True, False):
        for t in (0, .05, .1):
            np.testing.assert_array_equal(old.rhs(t, old.initial_state, [True]*4, active), new.rhs(t, new.initial_state, [True]*4, active))
    a = integrate_source(old, .12, active=False, sample_step=.01)
    b = integrate_source(new, .12, active=False, sample_step=.01)
    assert a.termination is None and b.termination is None
    np.testing.assert_array_equal(a.states, b.states)
    assert a.events == b.events


def test_independent_stage_limits_and_passive_head_behavior():
    m = BackupPumpModel(*ARGS, **KW, pump_parameters=PumpParameters(PumpStage(.2, 30), PumpStage(.5, 80)))
    y = m.initial_state
    p, _ = m.pressures(0, y)
    p.update(VC=5., PA=20., PU=5., AO=90.)
    m.set_controls(time=0, pump_on=True, right_head=20, left_head=100)
    q = m.flows(0, y, [False]*4, pressures=p)
    assert q['support_right'] == 25  # (5+20-20)/.2
    assert q['support_left'] == 30   # (5+100-90)/.5
    m.set_controls(time=0, right_head=100, left_head=200)
    q = m.flows(0, y, [False]*4, pressures=p)
    assert q['support_right'] == 30 and q['support_left'] == 80
    m.set_controls(time=0, right_head=0, left_head=0)
    assert m.flows(0, y, [False]*4, pressures=p)['support_left'] == 0
    # Zero head is not isolation: favorable pressure can drive passive flow.
    p['PU'] = 100
    assert m.flows(0, y, [False]*4, pressures=p)['support_left'] == 20
    m.set_controls(time=0, pump_on=False)
    assert m.flows(0, y, [False]*4, pressures=p)['support_left'] == 0


def test_one_stage_isolation_and_integrated_blood_conservation():
    m = BackupPumpModel(*ARGS, **KW, pump_parameters=PumpParameters(PumpStage(.2, 0), PumpStage(.3, 40)))
    m.set_controls(time=0, pump_on=True, right_head=20, left_head=100)
    tr = integrate_source(m, .12, active=False, sample_step=.01)
    assert tr.termination is None
    assert np.max(abs(tr.states[:m.volume_count].sum(axis=0)-5250))/5250 < 1e-8
    for t, y in zip(tr.time, tr.states.T):
        q = m.flows(t, y, [True]*4, active=False)
        assert q['support_right'] == 0 and 0 <= q['support_left'] <= 40


@pytest.mark.parametrize('r,c', [(0,1), (-1,1), (.1,-1), (float('nan'),1), (.1,float('inf')), (True,1), (.1,'120')])
def test_reject_invalid_parameters(r,c):
    with pytest.raises(ValueError):
        PumpStage(r,c)


def test_reject_incomplete_or_misspelled_configuration(tmp_path):
    p = tmp_path/'pump.json'
    for data in ({'right':{}}, {'right':{}, 'left':{}}, {'right':{'capacity':1}, 'left':{}}, []):
        p.write_text(json.dumps(data))
        with pytest.raises(ValueError):
            PumpParameters.from_json(p)
