import numpy as np
import pytest
from cardiac_redundancy.transport import marker_rhs


def test_zero_storage_junction_mixes_inflows_without_creating_marker():
    edges=[('a','A','H'),('b','B','H'),('out','H','C')]
    r=marker_rhs([10,10,10],[10,0,0],edges,[2,1,3],['A','B','C'])
    np.testing.assert_allclose(r,[-2,0,2],atol=1e-12)
    assert sum(r)==pytest.approx(0)


def test_reverse_junction_and_unbalanced_flow():
    edges=[('a','A','H'),('out','H','B')]
    np.testing.assert_allclose(marker_rhs([10,10],[0,10],edges,[-1,-1],['A','B']),[1,-1])
    with pytest.raises(ValueError):marker_rhs([10,10],[0,10],edges,[1,2],['A','B'])
