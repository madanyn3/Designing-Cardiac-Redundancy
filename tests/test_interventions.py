from pathlib import Path
import numpy as np
import pytest

from cardiac_redundancy.pressure_regional import PressureRegionalModel
from cardiac_redundancy.integration import integrate_source
from cardiac_redundancy.interventions import integrate_arrest

ROOT=Path(__file__).resolve().parents[1]


def model():
    return PressureRegionalModel(ROOT/'configs/baseline.yaml',
        ROOT/'configs/m3-evidence-regional.yaml', ROOT/'configs/m2-frozen-overlay.yaml')


def test_arrest_preserves_state_and_removes_activation():
    m=model()
    run=integrate_arrest(m,initial_state=m.initial_state,start=0,
                         arrest_time=.375,end=.45)
    event=run.interventions[0]
    assert event['before']['activation']>0.99
    assert event['after']['activation']==0
    np.testing.assert_array_equal(event['before']['state'],event['after']['state'])
    np.testing.assert_allclose(run.segments[0][1].states[:,-1],
                               run.segments[1][1].states[:,0],rtol=0,atol=1e-12)
    assert event['before']['pressures_mmHg']['LV']!=event['after']['pressures_mmHg']['LV']
    for active,tr in run.segments:
        assert np.max(abs(tr.states[:11].sum(axis=0)-5250))<1e-6
        assert tr.termination is None


def test_unmodified_restart_matches_uninterrupted():
    m=model()
    whole=integrate_source(m,.6,sample_step=.001)
    a=integrate_source(m,.3,sample_step=.001)
    b=integrate_source(m,.3,start=.3,initial_state=a.states[:,-1],sample_step=.001)
    scale=np.maximum(1,np.max(abs(whole.states[:,300:]),axis=1))
    assert np.max(abs(b.states-whole.states[:,300:])/scale[:,None])<5e-5


def test_domain_exit_prevents_future_arrest():
    m=model(); y=m.initial_state.copy()
    y[3]=m.inventory.unstressed_volumes[3]
    for arrest in (0,.1):
        run=integrate_arrest(m,initial_state=y,start=0,arrest_time=arrest,end=.2)
        assert run.termination['kind']=='model_domain_limit'
        assert run.interventions==[]


@pytest.mark.parametrize('start,arrest,end',[(1,0,2),(0,1,1),(0,float('nan'),2)])
def test_invalid_event_schedule(start,arrest,end):
    m=model()
    with pytest.raises(ValueError):
        integrate_arrest(m,initial_state=m.initial_state,start=start,arrest_time=arrest,end=end)
