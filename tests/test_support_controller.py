from dataclasses import replace
import numpy as np
import pytest
from cardiac_redundancy.support_controller import ControllerParameters, SupportController


def samples(controller, start, stop, pressure=50., flow=0., **kwargs):
    return [controller.update(float(t), pressure, flow, **kwargs)
            for t in np.round(np.arange(start, stop+.0001, .01), 8)]


def test_confirmation_release_and_ramp_are_distinct():
    c = SupportController()
    u = samples(c, 0, 1.5)
    assert c.events[1]['time_s'] == .25
    assert c.release_s == .5
    assert u[50] == 0 and u[100] == pytest.approx(.5)
    assert c.state == 'SUPPORT' and u[-1] == 1
    assert all(0 <= v <= 1 for v in u)


def test_short_pressure_dip_cancels_and_native_flow_blocks_false_trigger():
    c = SupportController()
    samples(c, 0, .1)
    c.update(.11, 66, 0)
    assert c.state == 'MONITOR' and c.release_s is None
    samples(c, .12, .5, pressure=40, flow=3)
    assert c.state == 'MONITOR' and c.command == 0


def test_recovery_requires_adequacy_and_aborts_continuously():
    c = SupportController()
    samples(c, 0, 1.5)
    samples(c, 1.51, 4, pressure=80, flow=3)
    assert c.state == 'SUPPORT'
    samples(c, 4.01, 6.51, pressure=80, flow=3, adequate_for_recovery=True)
    assert c.state == 'RECOVER' and c.command == pytest.approx(.75)
    before = c.command
    c.update(6.52, 40, 0)
    assert c.state == 'ACTIVATE' and abs(c.command-before) <= .00500001
    samples(c, 6.53, 7.52)
    assert c.state == 'SUPPORT' and c.command == 1


def test_complete_withdrawal_returns_to_monitor_and_can_retrigger():
    c = SupportController()
    samples(c, 0, 1.5)
    samples(c, 1.51, 5.51, pressure=80, flow=3, adequate_for_recovery=True)
    assert c.state == 'MONITOR' and c.command == 0
    c.update(5.52, 40, 0)
    assert c.state == 'CONFIRM'


@pytest.mark.parametrize('kind', ['sensor', 'gap', 'domain'])
def test_fault_is_latched_and_holds_command(kind):
    c = SupportController()
    samples(c, 0, 1)
    before = c.command
    t = 1.2 if kind == 'gap' else 1.01
    c.update(t, None if kind == 'sensor' else 50, 0, model_valid=kind != 'domain')
    assert c.state == 'FAULT' and c.command == before
    c.update(t+.01, 90, 5, adequate_for_recovery=True)
    assert c.state == 'FAULT' and c.command == before


def test_repeatability_and_no_time_reversal():
    a, b = SupportController(), SupportController()
    assert samples(a, 0, 2) == samples(b, 0, 2)
    assert a.events == b.events
    with pytest.raises(ValueError): a.update(2, 50, 0)


@pytest.mark.parametrize('change', [dict(ramp_s=0), dict(confirmation_s=-1),
    dict(reset_MAP_mmHg=50), dict(recovery_native_L_min=.1), dict(signal_s=float('nan'))])
def test_invalid_configuration(change):
    with pytest.raises(ValueError): replace(ControllerParameters(), **change)
