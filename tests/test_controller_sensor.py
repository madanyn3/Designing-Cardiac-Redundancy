import pytest
from cardiac_redundancy.controller_sensor import TrailingObservation


def test_no_fabricated_history_and_exact_constant():
    s = TrailingObservation(1)
    assert s.update(0, 80, 5) is None
    assert s.update(.4, 80, 5) is None
    assert s.update(1, 80, 5) == pytest.approx((80, 5))


def test_irregular_samples_linear_ramp_and_boundary_interpolation():
    s = TrailingObservation(1)
    for t in (0, .4, .9, 1.4): result = s.update(t, 10*t, 2*t)
    assert result == pytest.approx((9, 1.8))
    assert s.update(1.7, 17, 3.4) == pytest.approx((12, 2.4))


def test_identical_prefix_is_causal_and_invalid_input_rejected():
    a, b = TrailingObservation(.5), TrailingObservation(.5)
    for t in (0, .2, .5, .6):
        assert a.update(t, 70, 3) == b.update(t, 70, 3)
    a.update(.7, 0, 0)
    assert b.update(.7, 70, 3) == pytest.approx((70, 3))
    with pytest.raises(ValueError): b.update(.7, 80, 4)
    with pytest.raises(ValueError): b.update(.8, float('nan'), 4)
