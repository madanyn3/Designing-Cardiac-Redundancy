import json
from pathlib import Path

import numpy as np
import pytest

from cardiac_redundancy.experiments import apply_knobs, periodic_error, waveform_summary


@pytest.fixture
def baseline():
    return json.loads((Path(__file__).resolve().parents[1]/"configs/baseline.yaml").read_text())


def test_knobs_copy_and_preserve_inventory(baseline):
    before=json.dumps(baseline,sort_keys=True)
    changed=apply_knobs(baseline,{"initial_VC_to_PU_mL":50,"systemic_resistance_factor":2})
    assert sum(changed["initial_volumes_mL"].values())==pytest.approx(1500)
    assert changed["source_parameters"]["heart_parameters"]["R_sys"]==2*baseline["source_parameters"]["heart_parameters"]["R_sys"]
    assert json.dumps(baseline,sort_keys=True)==before


def test_rate_preserves_normalized_activation_and_compliance_is_inverse(baseline):
    changed=apply_knobs(baseline,{"heart_rate_bpm":60,"aortic_compliance_factor":2})
    d=changed["driver"]
    assert d["period_s"]==1
    assert d["C_s"]/d["period_s"]==pytest.approx(.5)
    assert d["B_per_s2"]*d["period_s"]**2==pytest.approx(45)
    assert changed["source_parameters"]["aorta"]["E_es_ao"]==baseline["source_parameters"]["aorta"]["E_es_ao"]/2


def test_thoracic_pressure_updates_actual_three_offsets(baseline):
    c=apply_knobs(baseline,{"thoracic_pressure_mmHg":0})
    assert all(c["pressure_external_mmHg"][n]==0 for n in ["PA","PU","pericardium"])


def test_thoracic_shift_enters_chambers_and_pulmonary_pressures_once(baseline,tmp_path):
    from cardiac_redundancy.source import SourceModel
    original=tmp_path/"original.json"
    modified=tmp_path/"modified.json"
    original.write_text(json.dumps(baseline))
    modified.write_text(json.dumps(apply_knobs(baseline,{"thoracic_pressure_mmHg":0})))
    a,b=SourceModel(original),SourceModel(modified)
    p,_=a.pressures(.375,a.initial_state)
    q,_=b.pressures(.375,b.initial_state)
    for node in ['LV','RV','PA','PU']: assert q[node]-p[node]==pytest.approx(4)
    for node in ['AO','VC']: assert q[node]==p[node]


@pytest.mark.parametrize("knobs",[{"bad_knob":1},{"heart_rate_bpm":0},{"systemic_resistance_factor":-1},
                                   {"activation_amplitude":1.1},{"initial_VC_to_PU_mL":1e6},{"lv_elastance_factor":float('nan')}])
def test_invalid_knob_cases_rejected(baseline,knobs):
    with pytest.raises(ValueError): apply_knobs(baseline,knobs)


def test_periodicity_checks_internal_phases_not_only_endpoints():
    previous=np.zeros((10,101))
    current=previous.copy()
    current[7,50]=.01
    assert periodic_error(current,previous)==.01
    assert periodic_error(previous,previous)==0


def test_pressure_mean_constant_and_sinusoid():
    t=np.linspace(0,1,1001)
    assert waveform_summary(t,np.full_like(t,80))=={"mean":80,"minimum":80,"maximum":80}
    s=waveform_summary(t,80+20*np.sin(2*np.pi*t))
    assert s["mean"]==pytest.approx(80,abs=1e-10)
    assert s["minimum"]==pytest.approx(60)
    assert s["maximum"]==pytest.approx(100)


def test_nonuniform_time_mean_uses_weights():
    t=np.array([0,.1,1.])
    assert waveform_summary(t,2*t)["mean"]==pytest.approx(1)


@pytest.mark.parametrize("limit,converges",[(4,True),(3,False)])
def test_warmup_requires_three_consecutive_comparisons(monkeypatch,limit,converges):
    from types import SimpleNamespace
    import cardiac_redundancy.experiments as exp
    initial=np.r_[np.full(6,100.),np.zeros(4)]
    model=SimpleNamespace(config={"driver":{"period_s":1}},initial_state=initial,
        settings={"warmup_max_cycles":limit,"periodic_state_relative_tolerance":1e-5,
                  "mass_drift_relative_limit":1e-8,"absolute_flow_tolerance_mL_s":1e-6})
    def steady(model,duration,initial_state,start,sample_step,**solver):
        return SimpleNamespace(states=np.repeat(initial_state[:,None],101,axis=1))
    monkeypatch.setattr(exp,"integrate_source",steady)
    protocol={"warmup_phase_intervals":100,"consecutive_periodic_comparisons":3}
    if converges:
        state,time,history=exp.warm_up(model,protocol,{})
        assert len(history)==4 and time==4
        assert history[0]["phase_error"] is None
        assert history[-1]["passing_streak"]==3
    else:
        with pytest.raises(RuntimeError,match="did not converge"):
            exp.warm_up(model,protocol,{})
