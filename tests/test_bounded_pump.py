from pathlib import Path
import numpy as np
import pytest
from cardiac_redundancy.bounded_pump import PreloadGate,bounded_stage_flow,PreloadPumpModel
from cardiac_redundancy.backup_pump import BackupPumpModel,PumpParameters,PumpStage
from cardiac_redundancy.selective import finite_head_flow

ROOT=Path(__file__).resolve().parents[1]
ARGS=(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml')


def test_preload_closure_cannot_be_overcome_by_large_head():
    assert bounded_stage_flow(-2,80,10000)==0
    assert bounded_stage_flow(0,80,10000)==0
    assert bounded_stage_flow(2,80,10000)==60


def test_filled_gate_matches_existing_law_and_preserves_afterload_dependence():
    for pin in [4,8]:
        for pout in [10,80,120]:
            for head in [0,20,100]:
                assert bounded_stage_flow(pin,pout,head)==finite_head_flow(pin,pout,head,.1,120)
    assert bounded_stage_flow(4,120,100)==0
    assert bounded_stage_flow(4,100,100)==40


def test_external_pressure_and_isolation_are_explicit():
    assert PreloadGate(external_pressure_mmHg=-4).fraction(-2)==.5
    assert bounded_stage_flow(8,0,0)>0
    assert bounded_stage_flow(8,0,0,enabled=False)==0
    with pytest.raises(ValueError):PreloadGate(full_above_mmHg=0)
    with pytest.raises(ValueError):PreloadGate(closed_below_mmHg=True)
    with pytest.raises(ValueError):bounded_stage_flow(1,2,float('nan'))


def test_adapter_fully_open_regression_against_existing_hydraulics():
    base=BackupPumpModel(*ARGS,arterial_fraction=.2,inlet_resistance_fraction=.1)
    gate=PreloadGate(-100,-99)
    model=PreloadPumpModel(*ARGS,arterial_fraction=.2,inlet_resistance_fraction=.1,right_gate=gate,left_gate=gate)
    for x in (base,model):x.set_controls(time=0,ramp_s=1,pump_on=True,right_head=20,left_head=100)
    np.testing.assert_allclose(model.rhs(2,model.initial_state,[False]*4,active=False),
                               base.rhs(2,base.initial_state,[False]*4,active=False),rtol=0,atol=0)


def test_independent_stage_limits_preserve_global_volume_balance():
    model=PreloadPumpModel(*ARGS,arterial_fraction=.2,inlet_resistance_fraction=.1,
        pump_parameters=PumpParameters(PumpStage(.1,20),PumpStage(.1,10)),
        right_gate=PreloadGate(100,101),left_gate=PreloadGate(-100,-99))
    model.set_controls(time=0,ramp_s=1,pump_on=True,right_head=100,left_head=200)
    q=model.flows(2,model.initial_state,[False]*4,active=False)
    assert q['support_right']==0
    assert q['support_left']==10
    derivative=model.rhs(2,model.initial_state,[False]*4,active=False)
    assert derivative[:len(model.volume_names)].sum()==pytest.approx(0,abs=1e-10)


def test_pulmonary_inlet_gate_uses_transmural_filling_convention():
    model=PreloadPumpModel(*ARGS,arterial_fraction=.2,inlet_resistance_fraction=.1,
        left_gate=PreloadGate(external_pressure_mmHg=-4))
    y=model.initial_state.copy();j=model.volume_names.index('PU');v=model.volume_names.index('VC')
    target=model.U[j]+2*model.C['PU'];transfer=y[j]-target;y[j]=target;y[v]+=transfer
    p,_=model.pressures(0,y,active=False)
    assert p['PU']==pytest.approx(-2)
    assert model.left_gate.fraction(p['PU'])==pytest.approx(.5)
    assert PreloadGate().fraction(p['PU'])==0
    assert y[:len(model.volume_names)].sum()==pytest.approx(5250)


def test_equal_support_flows_do_not_alone_balance_pulmonary_inventory():
    model=PreloadPumpModel(*ARGS,arterial_fraction=.2,inlet_resistance_fraction=.1)
    q=dict.fromkeys(model.flow_names,0.)
    q.update(support_right=20.,support_left=20.,PV=3.,MT=1.,lungs=.5)
    derivative=model.incidence@np.array([q[n] for n in model.flow_names])
    assert derivative[model.volume_names.index('PA')]+derivative[model.volume_names.index('PU')]==pytest.approx(2.5)
