import numpy as np
from scipy.integrate import solve_ivp
from cardiac_redundancy.oxygen_numerics import bounded_jacobian,periodic_rk4


def test_bounded_difference_matches_analytic_with_inert_ledger():
    def rhs(t,y):
        return np.array([-2*y[0]+y[1]**2,-y[1],2*y[0]-y[1]**2])
    y=np.array([3.,.5,1000.])
    expected=np.array([[-2,1,0],[0,-1,0],[2,-1,0]])
    np.testing.assert_allclose(bounded_jacobian(rhs)(0,y),expected,atol=1e-7)


def test_bounded_jacobian_stiff_transfer_matches_analytic_inventory():
    def rhs(t,y):return np.array([-1000*y[0],1000*y[0],0.])
    sol=solve_ivp(rhs,(0,.02),[1.,0.,0.],method='Radau',jac=bounded_jacobian(rhs),rtol=1e-9,atol=1e-11)
    np.testing.assert_allclose(sol.y[:,-1],[np.exp(-20),1-np.exp(-20),0],atol=1e-10)
    np.testing.assert_allclose(sol.y.sum(axis=0),1,atol=1e-12)


def test_periodic_rk4_matches_analytic_forcing_and_refines():
    errors=[]
    for steps in (8,16):
        force=np.cos(2*np.pi*np.linspace(0,1,2*steps+1))
        t,y=periodic_rk4(lambda k,y:np.array([force[k],-force[k]]),
            [1.,2.],1/steps,steps,3,steps//8)
        exact=1+np.sin(2*np.pi*t)/(2*np.pi)
        errors.append(np.max(abs(y[:,0]-exact)))
        np.testing.assert_allclose(y.sum(axis=1),3,atol=1e-12)
        assert t[-1]==3
    assert errors[1]<errors[0]/15 and errors[1]<2e-6
