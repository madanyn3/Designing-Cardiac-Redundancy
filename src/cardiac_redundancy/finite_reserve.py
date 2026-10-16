"""Finite, nonrefilling reserve diagnostic on prescribed deficit trajectories.

Units are mL O2-equivalent bookkeeping, not actual extra oxygen. This optimistic
coverage diagnostic cannot certify oxidizability, force or tissue recovery.
"""
import numpy as np
from scipy.integrate import solve_ivp


def cover_deficit(times, deficit, initial, tau=1., rtol=1e-8, max_step=.01):
    t, d, b0 = np.asarray(times, float), np.asarray(deficit, float), np.asarray(initial, float)
    if t.ndim != 1 or len(t)<2 or np.any(np.diff(t)<=0) or d.shape != (len(t), len(b0)) or b0.ndim != 1:
        raise ValueError('Invalid time/organ shapes')
    if any(np.any(~np.isfinite(x)) for x in (t,d,b0)) or np.any(d<0) or np.any(b0<0) or not np.isfinite(tau) or tau<=0:
        raise ValueError('Finite nonnegative inputs and positive release time required')
    m = len(b0)
    def rhs(time, y):
        demand = np.array([np.interp(time, t, d[:,j]) for j in range(m)])
        rate = np.minimum(demand, np.maximum(y[:m], 0)/tau)
        return np.r_[-rate, rate]
    sol = solve_ivp(rhs, (t[0],t[-1]), np.r_[b0,np.zeros(m)], t_eval=t,
                    rtol=rtol, atol=rtol*1e-3, max_step=max_step)
    if not sol.success:
        raise RuntimeError(sol.message)
    remaining, covered = sol.y[:m].T, sol.y[m:].T
    error = float(np.max(np.abs(remaining+covered-b0)))
    if remaining.min() < -1e-9 or covered.min() < -1e-9 or error > 1e-8:
        raise RuntimeError('Finite reserve conservation/domain failure')
    return remaining, covered, error
