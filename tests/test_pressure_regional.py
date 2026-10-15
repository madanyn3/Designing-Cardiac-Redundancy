from pathlib import Path
import numpy as np
import pytest
from scipy.integrate import solve_ivp

from cardiac_redundancy.pressure_regional import PressureRegionalModel, waterfall
from cardiac_redundancy.regional import RegionalModel

ROOT=Path(__file__).resolve().parents[1]
PATHS=(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml')


def test_waterfall_limits_symmetry_and_passivity():
    assert waterfall(100,2,20,8)==10
    assert waterfall(10,2,20,8)==0
    assert waterfall(100,2,-5,8)==(100-2)/8
    for a in (-5,0,2,20,100):
        for v in (-5,0,2,20,100):
            for c in (0,10,20,40):
                q=waterfall(a,v,c,8)
                assert q==-waterfall(v,a,c,8)
                assert (a-v)*q>=0


@pytest.mark.parametrize('method',['Radau','BDF'])
def test_threshold_crossings_integrated_against_analytic_volume(method):
    # Prescribed inlet ramp -1..3 crosses Pc=1 at t=1. Exact forward
    # integral is 1/2 mL; reversing endpoints must reverse the transfer.
    for direction in (-1,1):
        def rhs(t,y):
            q=direction*waterfall(2*t-1,0,1,2)
            return [-q,q]
        tr=solve_ivp(rhs,(0,2),[10.,10.],method=method,rtol=1e-9,atol=1e-10,max_step=.01)
        assert tr.success
        assert tr.y[1,-1]==pytest.approx(10+direction*.5,abs=1e-7)
        assert np.max(abs(tr.y.sum(axis=0)-20))<1e-10


def test_revised_anchor_and_no_recalibration_in_sensitivities():
    nominal=PressureRegionalModel(*PATHS)
    assert nominal.brain_resistance==pytest.approx(80/11.2)
    assert waterfall(100,2,20,nominal.brain_resistance)==pytest.approx(11.2)
    for icp in (5,10,15):
        for tone in (0,10,20):
            m=PressureRegionalModel(*PATHS,closure={'ICP_mmHg':icp,'tone_offset_mmHg':tone})
            assert m.brain_resistance==nominal.brain_resistance
    for factor in (.5,1,1.5):
        m=PressureRegionalModel(*PATHS,closure={'brain_resistance_multiplier':factor})
        assert m.brain_resistance==pytest.approx(factor*nominal.brain_resistance)


def test_disabled_closures_reproduce_comparator_and_conserve():
    ref=RegionalModel(*PATHS)
    m=PressureRegionalModel(*PATHS,closure={'brain_enabled':False,'coronary_enabled':False})
    for t in (0,.2,.375,.6):
        np.testing.assert_allclose(m.rhs(t,m.initial_state,[True]*4),ref.rhs(t,ref.initial_state,[True]*4),rtol=0,atol=0)
    m=PressureRegionalModel(*PATHS)
    for active in (True,False):
        for t in (0,.2,.375,.6):
            assert sum(m.rhs(t,m.initial_state,[True]*4,active=active)[:11])==pytest.approx(0,abs=1e-10)


def test_coronary_compression_and_passive_arrest_pressure():
    m=PressureRegionalModel(*PATHS)
    p={'AO':100.,'VC':2.,'PERI':5.,'LV':105.,'RV':25.}
    c=m.closing_pressures(p)
    assert c['coronary_left']==55 and c['coronary_right']==15
    zero=PressureRegionalModel(*PATHS,closure={'kappa':0})
    assert zero.closing_pressures(p)['coronary_left']==5
    # At kappa=0 pericardial compression can remain; it is not in general
    # identical to disabling coronary closure altogether.
    passive,meta=m.pressures(.375,m.initial_state,active=False)
    assert meta['activation']==0
    assert all(np.isfinite(list(m.closing_pressures(passive).values())))
    for k in (0,.5,1):
        c=PressureRegionalModel(*PATHS,closure={'kappa':k}).closing_pressures(p)
        assert waterfall(100,2,c['coronary_left'],30)<=98/30


@pytest.mark.parametrize('args',[(1,0,0,0),(float('nan'),0,0,1),(1,0,0,-1)])
def test_invalid_waterfall(args):
    with pytest.raises(ValueError): waterfall(*args)
