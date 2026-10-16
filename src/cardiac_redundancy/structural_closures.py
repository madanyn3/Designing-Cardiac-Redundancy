"""Exploratory alternative equation forms, not calibrated physiological laws.

These pure functions do not change any frozen simulation. Coupling them into a
circulation/transport experiment requires new conservation and refinement checks.
"""
import numpy as np


def flow_extraction(flow, reference_flow, reference_extraction):
    """Crone-shaped extraction with E(Qref)=Eref, using signed-flow magnitude.

    PS=-Qref*log(1-Eref). At zero flow E tends to one for positive PS,
    but delivery Q*E is zero. This is a transport-form experiment, not a
    claim that haemoglobin-bound oxygen obeys a single tracer compartment.
    """
    q, qr, er = np.broadcast_arrays(np.asarray(flow, float), np.asarray(reference_flow, float), np.asarray(reference_extraction, float))
    if np.any(~np.isfinite(q)) or np.any(~np.isfinite(qr)) or np.any(~np.isfinite(er)) or np.any(qr <= 0) or np.any((er < 0) | (er >= 1)):
        raise ValueError('Finite flow, positive reference flow and extraction in [0,1) required')
    ps = -qr*np.log1p(-er)
    ratio = np.divide(ps, np.abs(q), out=np.full_like(ps, np.inf), where=q != 0)
    return np.where(ps == 0, 0., -np.expm1(-ratio))


def transition_fraction(time, target, duration, form):
    """Explicit demand schedule; supply does not automatically reduce demand.

    Linear duration reaches the target; exponential duration is its time constant.
    Their integrated demands intentionally differ. Neither time is empirical.
    """
    t = np.asarray(time, float)
    if np.any(~np.isfinite(t)) or np.any(t < 0) or not np.isfinite(target) or not 0 < target <= 1 or not np.isfinite(duration) or duration < 0:
        raise ValueError('Invalid demand-transition inputs')
    if form not in ('step', 'linear', 'exponential'):
        raise ValueError('Unknown transition form')
    if form == 'step' or duration == 0:
        return np.full_like(t, target)
    remaining = np.maximum(0., 1-t/duration) if form == 'linear' else np.exp(-t/duration)
    return target+(1-target)*remaining
