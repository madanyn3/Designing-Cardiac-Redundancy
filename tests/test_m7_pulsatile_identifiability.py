"""Checks for the declared observation/parameter and SVD conventions."""
import importlib.util
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT/'scripts'/f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_log_knobs_scale_all_active_walls_only():
    module = load('run_m7_pulsatile_identifiability')
    knobs = module.knobs(np.log([1.1, .9, 1.2, .8]))
    assert np.isclose(knobs['aortic_compliance_factor'], 1.2)
    assert all(np.isclose(knobs[k], .8) for k in ['lv_elastance_factor', 'rv_elastance_factor', 'septal_elastance_factor'])
    assert len(knobs) == 6


def test_two_observations_have_two_unseen_directions():
    module = load('report_m7_pulsatile_identifiability')
    j = np.array([[1., 2, 0, 0], [0, 0, 3, 4]])
    result = module.svd_audit(j)
    assert result['nullity'] == 2
    assert result['condition_number'] is None
    null_vectors = np.asarray(result['right_singular_vectors'])[2:]
    np.testing.assert_allclose(j@null_vectors.T, 0, atol=1e-14)


def test_full_rank_reports_conditioning():
    module = load('report_m7_pulsatile_identifiability')
    result = module.svd_audit(np.diag([4., 3, 2, .01]))
    assert result['nullity'] == 0
    assert result['condition_number'] == 400
