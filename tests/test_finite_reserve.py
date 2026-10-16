import numpy as np
import pytest
from cardiac_redundancy.finite_reserve import cover_deficit


def test_zero_reserve_covers_nothing():
    t=np.linspace(0,2,51)
    b,c,e=cover_deficit(t,np.ones((len(t),1)),[0])
    np.testing.assert_array_equal(c,0)
    assert e==0


def test_release_limited_exponential_and_conservation():
    t=np.linspace(0,2,51)
    b,c,e=cover_deficit(t,np.ones((len(t),1))*10,[1],tau=2)
    np.testing.assert_allclose(b[:,0],np.exp(-t/2),rtol=1e-7)
    np.testing.assert_allclose(b+c,1,atol=1e-12)


def test_large_reserve_covers_only_requested_deficit():
    t=np.linspace(0,2,51)
    b,c,e=cover_deficit(t,np.ones((len(t),1))*.2,[10])
    assert np.isclose(c[-1,0],.4)
    assert np.isclose(b[-1,0],9.6)


def test_negative_inventory_is_rejected():
    with pytest.raises(ValueError):
        cover_deficit([0,1],np.ones((2,1)),[-1])
