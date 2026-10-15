import numpy as np
import pytest
from cardiac_redundancy.controlled_circulation import Envelope, ControlledCirculation
from cardiac_redundancy.bilateral_support import BilateralSupport
from cardiac_redundancy.circuit_integration import integrate_circuit
from test_backup_pump import ARGS, KW


def test_envelope_history_abort_and_invalid_events():
    e=Envelope();e.set(1,1,2)
    assert e.value(0)==0 and e.value(2)==.5
    e.set(2,0,1)
    assert e.value(1.5)==.25 and e.value(2)==.5 and e.value(3)==0
    with pytest.raises(ValueError):e.set(1,1,1)
    with pytest.raises(ValueError):e.set(4,2,1)


def test_full_envelopes_reproduce_native_and_auxiliary_pressures():
    a=ControlledCirculation(*ARGS,**KW);b=BilateralSupport(*ARGS,**KW)
    a.release_auxiliary(0);b.release_auxiliary(0)
    a.aux_envelope.set(0,1,1);a.native_envelope.set(0,1,0)
    for t in (.1,.4,1.1):
        pa,_=a.pressures(t,a.initial_state);pb,_=b.pressures(t,b.initial_state,True)
        for key in pa:assert pa[key]==pytest.approx(pb[key],abs=1e-8)


def test_native_recovery_preserves_states_and_conserves_blood():
    a=ControlledCirculation(*ARGS,**KW);a.native_envelope.set(0,1,.2)
    y=a.initial_state.copy()
    tr=integrate_circuit(a,.25,sample_step=.01)
    assert tr.termination is None
    np.testing.assert_array_equal(y,tr.states[:,0])
    assert np.max(abs(tr.states[:18].sum(axis=0)-y[:18].sum()))<1e-5
