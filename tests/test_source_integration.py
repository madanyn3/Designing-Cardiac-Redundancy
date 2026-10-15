from pathlib import Path

import numpy as np

from cardiac_redundancy.integration import integrate_source
from cardiac_redundancy.source import SourceModel


def test_one_cycle_switches_all_valves_and_preserves_blood():
    model = SourceModel(Path(__file__).resolve().parents[1] / "configs/baseline.yaml")
    before = model.initial_state.copy()
    result = integrate_source(model, .75)
    assert abs(result.time[-1] - .75) < 1e-12
    assert np.max(np.abs(result.states[:6].sum(axis=0) - 1500)) / 1500 < 1e-8
    assert result.states[6:].min() >= -1e-6
    assert result.states[:6].min() > 0
    assert {e["valve"] for e in result.events} == {"MT","AV","TC","PV"}
    for valve in {e["valve"] for e in result.events}:
        actions = [e["action"] for e in result.events if e["valve"] == valve]
        assert all(a != b for a,b in zip(actions, actions[1:]))
    np.testing.assert_array_equal(model.initial_state, before)
    np.testing.assert_allclose(result.states[:,0], before, atol=1e-12, rtol=0)


def test_sparse_output_does_not_drop_short_valve_segments():
    model = SourceModel(Path(__file__).resolve().parents[1] / "configs/baseline.yaml")
    # Only endpoints are sampled. Intermediate valve segments still integrate
    # and contribute events, even when no output time falls inside them.
    result = integrate_source(model, .75, sample_step=.75)
    assert result.states.shape == (10,2)
    assert len(result.events) >= 4
    assert abs(result.states[:6,-1].sum()-1500)/1500 < 1e-8
