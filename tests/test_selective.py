from pathlib import Path
import numpy as np
import pytest
from scipy.integrate import solve_ivp
from cardiac_redundancy.selective import SelectiveModel,finite_head_flow,transport_rhs

ROOT=Path(__file__).resolve().parents[1]
def model(f=.2):
    return SelectiveModel(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',
                         ROOT/'configs/m2-frozen-overlay.yaml',arterial_fraction=f,inlet_resistance_fraction=.1)


def test_split_inventory_and_reference_series_resistance():
    for f in (.1,.2,.3):
        m=model(f); assert sum(m.initial_state[:16])==pytest.approx(5250)
        for pool,a in zip(m.pool_names,m.arterial_names):
            i=m.volume_names.index(pool); j=m.volume_names.index(a)
            assert m.U[i]+m.U[j]==pytest.approx(m.base.inventory.unstressed_volumes[i])
            assert m.initial_state[i]+m.initial_state[j]==pytest.approx(m.base.initial_state[i])
        p,_=m.pressures(0,m.initial_state);p['AO']=100
        q=m.flows(0,m.initial_state,[False]*4,pressures=p)
        for n in m.beds:
            target=95/m.base.branch_by_name[n].resistance
            assert q[n]==pytest.approx(target)
            assert q[n+'_exchange']==pytest.approx(target)


def test_gates_leakage_trapping_and_conservation():
    m=model();y=m.initial_state.copy()
    m.set_controls(time=0,inlet=0,outlet=0)
    q=m.flows(.2,y,[True]*4);dy=m.rhs(.2,y,[True]*4)
    for bed,pool,a in zip(m.beds,m.pool_names,m.arterial_names):
        assert q[bed]==q[bed+'_return']==0
        assert dy[m.volume_names.index(pool)]+dy[m.volume_names.index(a)]==pytest.approx(0,abs=1e-12)
    m.set_controls(time=0,inlet=0,outlet=1,collateral=.02,pump_on=True,right_head=20,left_head=100)
    q=m.flows(.2,y,[True]*4)
    assert any(q[n+'_collateral']!=0 for n in m.beds)
    for active in (True,False):assert sum(m.rhs(.2,y,[True]*4,active=active)[:16])==pytest.approx(0,abs=1e-10)
    np.testing.assert_array_equal(y,m.initial_state)


def test_ramp_and_pump_limits():
    m=model();m.set_controls(time=2,ramp_s=2,inlet=0)
    assert m.control_at(2)['inlet']==1 and m.control_at(3)['inlet']==.5 and m.control_at(4)['inlet']==0
    assert finite_head_flow(2,100,20,.1,120)==0
    assert finite_head_flow(2,100,100,.1,120)==20
    assert finite_head_flow(2,100,200,.1,120)==120
    with pytest.raises(ValueError):m.set_controls(time=4,outlet=-.1)


def test_passive_marker_washout_matches_analytic_and_reverses_donor():
    # Two equal reservoirs with equal counter-circulating flow: conserved
    # marker tends to uniform concentration with time constant V/(2Q).
    def rhs(t,a):return transport_rhs([10.,10.],a,[(0,1),(1,0)],[1.,1.])
    tr=solve_ivp(rhs,(0,10),[10.,0.],rtol=1e-10,atol=1e-11)
    assert tr.y[0,-1]==pytest.approx(5+5*np.exp(-2),abs=1e-8)
    assert np.max(abs(tr.y.sum(axis=0)-10))<1e-10
    np.testing.assert_array_equal(transport_rhs([10,10],[10,0],[(1,0)],[-1]),[-1,1])
    np.testing.assert_array_equal(transport_rhs([10,10],[10,0],[(0,1)],[0]),[0,0])
