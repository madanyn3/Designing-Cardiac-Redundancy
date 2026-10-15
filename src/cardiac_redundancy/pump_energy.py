"""Mechanical energy ledger for the abstract head/resistance/capacity ports.

The capacity ceiling is interpreted as an ideal dissipative flow limiter. This
closes the hydraulic balance; it is not a measured hardware efficiency model.
Signed port work is retained: passive pressure-driven flow can consume energy
from the circulation without creating recoverable input-energy credit.
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from .selective import finite_head_flow

# Conventional mmHg to Pa times mL to m^3: pressure times volume gives joules.
MMHG_ML_TO_J = 133.322387415e-6
POWER_FIELDS = ('head_W', 'port_W', 'resistance_loss_W', 'limiter_loss_W')


def stage_power(pin, pout, head, resistance, capacity, enabled=True):
    """Return one stage's flow, pressure rise and signed power balance in SI."""
    if not isinstance(enabled, bool):
        raise ValueError('enabled must be boolean')
    available_flow = finite_head_flow(pin, pout, head, resistance, capacity)
    q = available_flow if enabled else 0.
    rise = pout-pin
    # For positive flow the remaining pressure drop is zero when uncapped;
    # when capped it belongs to the explicit ideal flow-limiter comparator.
    residual_drop = max(0., head-rise-resistance*q) if q else 0.
    k = MMHG_ML_TO_J
    return dict(flow_mL_s=q, pressure_rise_mmHg=rise,
                head_W=k*q*head, port_W=k*q*rise,
                resistance_loss_W=k*resistance*q*q,
                limiter_loss_W=k*q*residual_drop)


def integrate_power(time_s, records):
    """Integrate one continuous segment; split at discrete control changes.

    Records must include both endpoints. Work is cumulative trapezoidal
    quadrature on the supplied output grid, so grid refinement must be checked
    separately from ODE convergence. Add segment totals across interventions;
    never interpolate across a discontinuous power change.
    """
    t = np.asarray(time_s, dtype=float)
    if t.ndim != 1 or len(t)<2 or not np.all(np.isfinite(t)) or np.any(np.diff(t)<=0):
        raise ValueError('At least two strictly increasing finite times required')
    if len(records)!=len(t):
        raise ValueError('One power record per time required')
    result = {}
    for field in POWER_FIELDS:
        p = np.asarray([r[field] for r in records], dtype=float)
        if not np.all(np.isfinite(p)):
            raise ValueError('Nonfinite power')
        result[field.replace('_W','_J')] = cumulative_trapezoid(p,t,initial=0.)
    return result


def pump_power(model, t, state, active=True):
    """Evaluate both configured stages without changing the hydraulic state."""
    p, _ = model.pressures(t, state, active=active)
    c = model.control_at(t)
    return {
        side: stage_power(p[a],p[b],c[side+'_head'],
                          getattr(model.pump_parameters,side).resistance_mmHg_s_mL,
                          getattr(model.pump_parameters,side).capacity_mL_s,
                          c['pump_on'])
        for side,a,b in (('right','VC','PA'),('left','PU','AO'))
    }
