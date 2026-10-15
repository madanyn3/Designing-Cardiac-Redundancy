import numpy as np
import pytest
from cardiac_redundancy.biological_energy import chamber_power,biological_power
from cardiac_redundancy.biological_pump import BiologicalPumpModel
from cardiac_redundancy.pump_energy import MMHG_ML_TO_J as K
from test_backup_pump import ARGS,KW

def test_pressure_driven_steady_valve_balance():
    r=chamber_power(10,5,2,0,5,3,1,1,2,2)
    assert r['valve_storage_rate_W']==pytest.approx(0)
    assert r['chamber_W']==pytest.approx(-10*K)
    assert r['port_W']==pytest.approx(-44*K)
    assert r['valve_loss_W']==pytest.approx(34*K)
    assert abs(r['balance_residual_W'])<1e-12

def test_inertial_energy_and_pericardial_partition():
    r=chamber_power(8,70,90,4,2,10,.1,.2,.01,.02)
    assert r['valve_stored_J']==pytest.approx(.5*K*(.01*4+.02*100))
    assert r['transmural_W']+r['pericardial_exchange_W']==pytest.approx(r['chamber_W'])
    assert abs(r['balance_residual_W'])<1e-12

def test_zero_flow_not_zero_metabolic_claim():
    r=chamber_power(5,100,80,0,0,0,1,1,1,1)
    assert all(v==0 for v in r.values())

@pytest.mark.parametrize('bad',[float('nan'),float('inf')])
def test_invalid_pressure(bad):
    with pytest.raises(ValueError):chamber_power(bad,1,2,0,0,0,1,1,1,1)

def test_live_model_balance():
    m=BiologicalPumpModel(*ARGS,**KW);m.release_auxiliary(0)
    y=m.initial_state.copy();y[22:]=[1,2,3,4]
    for r in biological_power(m,.375,y).values():assert abs(r['balance_residual_W'])<1e-12
