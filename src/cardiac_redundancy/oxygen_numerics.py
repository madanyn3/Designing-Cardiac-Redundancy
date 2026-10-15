"""Bounded finite-difference Jacobian for oxygen systems with inert ledgers.

SciPy's adaptive numerical differencing can enlarge perturbations indefinitely
for identically zero storage/ledger columns. A fixed, relative perturbation
avoids stepping outside bounded gas/store domains during Jacobian evaluation.
This affects derivative estimation only, never accepted state values.
"""
import numpy as np


def bounded_jacobian(rhs):
    def jac(time,state):
        base=rhs(time,state);out=np.empty((len(state),len(state)))
        for j in range(len(state)):
            step=1e-8*max(1.,abs(state[j]));trial=state.copy();trial[j]+=step
            out[:,j]=(rhs(time,trial)-base)/step
        return out
    return jac


def periodic_rk4(rate_at_index, initial, step_s, steps_per_period, periods, save_every):
    """RK4 on a periodic forcing grid, including precomputed half-step values.

    ``rate_at_index(k, state)`` uses forcing index k in [0, 2*N], where
    N=steps_per_period. Even indices are full steps, odd indices midpoints.
    Only the FORCING repeats; the oxygen state evolves continuously. Aligning
    steps with hydraulic knots avoids adaptive sampling across their corners.
    Saved states are never clipped; the caller must check physical bounds.
    """
    if not np.isfinite(step_s) or step_s<=0 or any(not isinstance(x,int) or x<=0 for x in (steps_per_period,periods,save_every)):
        raise ValueError('Positive step and integer grid counts required')
    total=steps_per_period*periods
    if total%save_every:raise ValueError('Output stride must divide horizon')
    y=np.asarray(initial,dtype=float).copy();states=[y.copy()]
    for j in range(total):
        k=2*(j%steps_per_period)
        a=rate_at_index(k,y)
        b=rate_at_index(k+1,y+step_s*a/2)
        c=rate_at_index(k+1,y+step_s*b/2)
        d=rate_at_index(k+2,y+step_s*c)
        y+=step_s*(a+2*b+2*c+d)/6
        if (j+1)%save_every==0:states.append(y.copy())
    return np.arange(len(states))*step_s*save_every,np.array(states)
