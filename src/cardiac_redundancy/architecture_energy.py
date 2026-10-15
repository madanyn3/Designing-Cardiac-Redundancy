"""Energy-source accounting at identical prescribed hydraulic boundaries.

This is not a realizable mechanical-device model. Positive and returned work
are separate; negative port power never finances positive work without an
explicit storage/regeneration model. Each stage is clipped separately.
"""
import numpy as np

def external_input(time_s,stage_port_W,efficiency,overhead_W=0.):
    t=np.asarray(time_s,dtype=float);p=np.asarray(stage_port_W,dtype=float)
    if t.ndim!=1 or len(t)<2 or np.any(~np.isfinite(t)) or np.any(np.diff(t)<=0):
        raise ValueError('Strictly increasing finite times required')
    if p.ndim!=2 or p.shape[0]!=len(t) or p.shape[1]!=2 or np.any(~np.isfinite(p)):
        raise ValueError('Two finite stage power histories required')
    if not np.isfinite(efficiency) or not 0<efficiency<=1 or not np.isfinite(overhead_W) or overhead_W<0:
        raise ValueError('Efficiency in (0,1], nonnegative overhead required')
    positive=float(np.trapezoid(np.maximum(p,0).sum(axis=1),t))
    returned=float(np.trapezoid(np.maximum(-p,0).sum(axis=1),t))
    signed=float(np.trapezoid(p.sum(axis=1),t))
    return dict(positive_port_J=positive,returned_port_J=returned,signed_port_J=signed,
                external_input_J=positive/efficiency+overhead_W*(t[-1]-t[0]),
                efficiency=efficiency,overhead_W=overhead_W)
