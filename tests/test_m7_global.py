import numpy as np
import pytest
from cardiac_redundancy.global_effects import trajectories,elementary_effects,summarize

def test_seed_prefix_and_interior_coverage():
    a=trajectories(8,6,7202);b=trajectories(16,6,7202)
    np.testing.assert_array_equal(a,b[:8])
    assert a.min()>=-1e-12 and a.max()<=1+1e-12
    assert np.any((a>0.1)&(a<.9))
    for t in a:
        dx=np.diff(t,axis=0)
        assert np.all((abs(dx)>1e-10).sum(axis=1)==1)
        assert np.all((abs(dx)>1e-10).sum(axis=0)==1)
        np.testing.assert_allclose(abs(dx[abs(dx)>1e-10]),2/3)

def test_linear_effects_have_exact_sign_and_magnitude():
    x=trajectories(8,3,1);weights=np.array([2,-3,0]);y=x@weights+7
    e=elementary_effects(x,y)
    np.testing.assert_allclose(e,np.tile(weights,(8,1)),atol=1e-14)
    s=summarize(e);np.testing.assert_allclose(s['mu_star'],[2,3,0],atol=1e-14)
    assert max(s['sigma'])<1e-14

def test_interaction_effect_depends_on_background():
    x=trajectories(8,2,2);y=x[:,:,0]*x[:,:,1];e=elementary_effects(x,y)
    for k,xx in enumerate(x):
        for j,dx in enumerate(np.diff(xx,axis=0)):
            changed=np.flatnonzero(abs(dx)>1e-10)[0]
            assert e[k,changed]==pytest.approx(xx[j,1-changed])
    assert np.std(e[:,0])>0

def test_incomplete_trajectory_rejected():
    with pytest.raises(ValueError):elementary_effects(np.zeros((2,3,2)),np.zeros((2,3)))

def test_bootstrap_reproducible():
    e=np.arange(24).reshape(8,3)
    assert summarize(e,3,100)==summarize(e,3,100)
