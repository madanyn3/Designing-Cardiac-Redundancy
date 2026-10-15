"""Historical two-sphere arithmetic, not a myocardial constitutive model.

Exact-shell mode changes geometry only. It deliberately retains the historical
thin-wall stress prescription to expose its effect, not to validate mechanics.
"""
from dataclasses import dataclass
import math

MMHG_TO_PA = 133.322387415

@dataclass(frozen=True)
class PrimerInputs:
    external_mL_min: float = 700.
    rate_per_min: float = 60.
    ejection_fraction: float = .45
    wall_pressure_mmHg: float = 90.
    stress_kPa: float = 22.
    density_g_mL: float = 1.05
    self_flow_mL_min_g: float = 250/300

    def __post_init__(self):
        if any(not math.isfinite(v) or v<=0 for v in vars(self).values()):
            raise ValueError('Positive finite historical inputs required')
        if self.ejection_fraction>1:
            raise ValueError('Ejection fraction cannot exceed one')

def coefficients(p, exact_shell=False):
    alpha=p.wall_pressure_mmHg*MMHG_TO_PA/1000/(2*p.stress_kPa)
    ratio=1+alpha+alpha*alpha/3
    k=6*p.density_g_mL*alpha/(p.rate_per_min*p.ejection_fraction)
    if exact_shell:k*=ratio
    return alpha,ratio,k,p.self_flow_mL_min_g*k

def geometry(flow,p,exact_shell=False):
    if not math.isfinite(flow) or flow<=0:raise ValueError('Positive flow required')
    alpha,ratio,k,c=coefficients(p,exact_shell)
    edv=flow/(p.rate_per_min*p.ejection_fraction)
    r=(3*edv/(4*math.pi))**(1/3)
    return dict(flow_mL_min=flow,stroke_volume_mL=flow/p.rate_per_min,
                EDV_each_mL=edv,radius_cm=r,wall_cm=alpha*r,
                h_over_r=alpha,two_shell_mass_g=k*flow,
                self_flow_mL_min=c*flow,feedback_gain=c,
                exact_to_thin_same_dimensions=ratio)

def fixed_point(p,exact_shell=False):
    c=coefficients(p,exact_shell)[3]
    if c>=1:raise ValueError('No positive finite fixed point for positive external load')
    return geometry(p.external_mL_min/(1-c),p,exact_shell)

def iterate(p,seed,atol=1e-6,exact_shell=False,max_updates=10000):
    if not math.isfinite(seed) or seed<0 or not math.isfinite(atol) or atol<=0:
        raise ValueError('Finite nonnegative seed and positive tolerance required')
    c=coefficients(p,exact_shell)[3]
    if c>=1:raise ValueError('Noncontractive iteration')
    values=[float(seed)]
    for _ in range(max_updates):
        residual=abs(p.external_mL_min+c*values[-1]-values[-1])
        # For this linear contraction residual/(1-c) is the exact error.
        if residual/(1-c)<=atol:return dict(values=values,error_bound_mL_min=residual/(1-c),updates=len(values)-1)
        values.append(p.external_mL_min+c*values[-1])
    raise RuntimeError('Iteration budget exhausted')

def steady_power_W(flow_mL_min,head_mmHg):
    if not all(math.isfinite(v) for v in (flow_mL_min,head_mmHg)) or min(flow_mL_min,head_mmHg)<0:
        raise ValueError('Nonnegative steady flow/head required')
    return flow_mL_min*1e-6/60*head_mmHg*MMHG_TO_PA
