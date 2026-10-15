import numpy as np
import pytest
from cardiac_redundancy.pump_energy import stage_power, integrate_power, pump_power, MMHG_ML_TO_J as K
from test_backup_pump import ARGS, KW
from cardiac_redundancy.backup_pump import BackupPumpModel


@pytest.mark.parametrize('head,cap,enabled', [(100,120,True),(100,10,True),(0,120,True),(100,0,True),(100,120,False)])
def test_power_balance_including_limiter_and_isolation(head,cap,enabled):
    r=stage_power(100,90,head,.5,cap,enabled)
    assert r['head_W']==pytest.approx(r['port_W']+r['resistance_loss_W']+r['limiter_loss_W'],abs=1e-14)
    assert r['resistance_loss_W']>=0 and r['limiter_loss_W']>=0
    if head==0:
        assert r['port_W']<0 and r['head_W']==0  # no energy-recovery credit


def test_known_power_and_constant_work():
    # 30 mL/s through 90 mmHg: net 0.359970446 W; 10 mmHg resistor drop.
    r=stage_power(5,95,100,1/3,120)
    assert r['flow_mL_s']==30
    assert r['port_W']==pytest.approx(K*2700)
    assert r['limiter_loss_W']==pytest.approx(0,abs=1e-14)
    work=integrate_power([0,600],[r,r])
    assert work['port_J'][-1]==pytest.approx(K*2700*600)


def test_linear_power_quadrature_and_invalid_grid():
    a=stage_power(0,90,100,.5,120)
    b={k:2*v for k,v in a.items()}
    assert integrate_power([2,5],[a,b])['port_J'][-1]==pytest.approx(4.5*a['port_W'])
    for t in ([0,0],[1,0],[0,float('nan')],[0]):
        with pytest.raises(ValueError):integrate_power(t,[a]*len(t))


def test_model_power_matches_actual_support_flow_without_state_mutation():
    m=BackupPumpModel(*ARGS,**KW)
    m.set_controls(time=0,ramp_s=1,pump_on=True,right_head=20,left_head=100)
    y=m.initial_state.copy()
    for t in (0,.5,1):
        powers=pump_power(m,t,y,active=False)
        flows=m.flows(t,y,[False]*4,active=False)
        for side in ('right','left'):
            assert powers[side]['flow_mL_s']==flows['support_'+side]
    np.testing.assert_array_equal(y,m.initial_state)
