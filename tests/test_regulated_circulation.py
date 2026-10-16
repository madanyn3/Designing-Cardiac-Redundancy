import numpy as np
import pytest
from cardiac_redundancy.regulated_circulation import regulated_flow, RegulatedCirculation
from pathlib import Path


def test_regulation_limits_and_passivity():
    assert regulated_flow(0, 10, 4) == 0
    assert regulated_flow(-2, 10, 4) == -2
    assert regulated_flow(1, 10, 4) == 4
    assert regulated_flow(5, 10, 4) == 10
    assert regulated_flow(30, 10, 4) == 15
    with pytest.raises(ValueError):
        regulated_flow(1, 10, .5)


def test_regulation_conserves_volume_and_changes_only_brain_flow():
    root = Path(__file__).resolve().parents[1]
    m = RegulatedCirculation(root/'configs/baseline.yaml', root/'configs/m3-evidence-regional.yaml', root/'configs/m2-frozen-overlay.yaml', arterial_fraction=.2, inlet_resistance_fraction=.1, regulation_gain=4.)
    y = np.load(root/'reports/m5-joint-preservation/reference_g25/fine-refined.npz')['states'][:,0]
    q = m.flows(0, y, [True]*8)
    m.regulation_gain = None
    base = m.flows(0, y, [True]*8)
    for name in q:
        if name != 'brain':
            assert q[name] == base[name]
    dv = m.incidence@np.array([q[n] for n in m.flow_names])
    assert abs(dv.sum()) < 1e-10
    assert abs(dv[-1]) < 1e-9
