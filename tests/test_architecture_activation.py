import numpy as np
import pytest
from cardiac_redundancy.architecture_energy import external_input
from cardiac_redundancy.activation_timing import ActivationTiming
from cardiac_redundancy.biological_pump import BiologicalPumpModel
from cardiac_redundancy.circuit_integration import integrate_circuit
from test_backup_pump import ARGS,KW

def test_stage_returns_do_not_cancel_other_stage_input():
    r=external_input([0,2],[[2,-1],[2,-1]],.5,.25)
    assert r['positive_port_J']==4 and r['returned_port_J']==2
    assert r['signed_port_J']==2 and r['external_input_J']==8.5

@pytest.mark.parametrize('eta',[0,-1,1.1,float('nan')])
def test_invalid_efficiency(eta):
    with pytest.raises(ValueError):external_input([0,1],[[1,1],[1,1]],eta)

def test_delays_and_ramp_are_separate():
    s=ActivationTiming()
    assert s.latency_s==pytest.approx(.5)
    assert s.amplitude(.49)==0 and s.amplitude(.5)==0
    assert s.amplitude(1)==pytest.approx(.5) and s.amplitude(1.5)==1

@pytest.mark.parametrize('kwargs',[{'detection_s':-1},{'ramp_s':0},{'signal_s':float('nan')},{'recruitment_s':True}])
def test_invalid_timing(kwargs):
    with pytest.raises(ValueError):ActivationTiming(**kwargs)

def test_release_does_not_teleport_volume_and_uses_ramp():
    m=BiologicalPumpModel(*ARGS,**KW);s=ActivationTiming(.01,.01,.01,.2)
    pre=integrate_circuit(m,.03,active=False,sample_step=.01)
    assert pre.termination is None and not any(m.valve_enabled[4:])
    y=pre.states[:,-1].copy();before=y.copy()
    with pytest.raises(ValueError):s.release(m,0,.02)
    s.release(m,0,.03)
    np.testing.assert_array_equal(y,before)
    assert m.trigger==.03 and m.auxiliary_ramp_s==.2
    assert m.pressures(.03,y,False)[1]['aux_activation']==0
    post=integrate_circuit(m,.05,start=.03,initial_state=y,active=False,sample_step=.01)
    assert post.termination is None
    assert np.max(abs(post.states[:18].sum(axis=0)-5250))<1e-5
    with pytest.raises(ValueError):s.release(m,0,.03)
