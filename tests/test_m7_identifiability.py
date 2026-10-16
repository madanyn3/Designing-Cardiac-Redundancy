"""Analytical checks of observation rank and nonunique synthetic recovery."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

spec = importlib.util.spec_from_file_location('identifiability', Path(__file__).resolve().parents[1]/'scripts/report_m7_identifiability.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_repeated_total_flow_does_not_separate_parallel_paths():
    a = module.observation_matrix([70, 80, 90, 100], 5)
    assert np.linalg.matrix_rank(a) == 1
    np.testing.assert_allclose(a @ np.array([1, -1, 0, 0, 0]), 0)


def test_regional_measurements_recover_unique_conductances():
    a = module.observation_matrix([70, 80, 90], 5, (0, 1, 2, 3))
    truth = np.array([.1, .3, .2, .05, .15])
    assert np.linalg.matrix_rank(a) == 5
    np.testing.assert_allclose(np.linalg.lstsq(a, a @ truth, rcond=None)[0], truth)


def test_one_regional_measurement_adds_one_constraint():
    assert np.linalg.matrix_rank(module.observation_matrix([80, 90], 5, (0,))) == 2


def test_invalid_observations_rejected():
    with pytest.raises(ValueError):
        module.observation_matrix([80], 3, (0, 0))
    with pytest.raises(ValueError):
        module.observation_matrix([np.nan], 3)
