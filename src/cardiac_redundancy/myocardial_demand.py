"""Normal-reference-normalized myocardial demand sensitivity, not a PVA model.

Fractions and tissue mass ratios are explicit assumptions. Positive transmural
boundary-work rate is a workload surrogate; it is not total contractile energy.
The closure reproduces the chosen normal demand by construction, not validation.
"""
from dataclasses import dataclass
import math

@dataclass(frozen=True)
class DemandAssumptions:
    maintenance_fraction: float=.15
    activation_fraction: float=.25
    work_fraction: float=.60
    mass_ratio: float=.25
    work_cost_multiplier: float=1.
    contractility_exponent: float=1.
    def __post_init__(self):
        vals=vars(self)
        if any(isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) for x in vals.values()):raise ValueError('Finite numeric assumptions required')
        if min(self.maintenance_fraction,self.activation_fraction,self.work_fraction)<0 or not math.isclose(self.maintenance_fraction+self.activation_fraction+self.work_fraction,1.,abs_tol=1e-12):raise ValueError('Fractions must sum to one')
        if self.mass_ratio<=0 or self.work_cost_multiplier<=0 or self.contractility_exponent<0:raise ValueError('Positive mass/cost and nonnegative exponent required')

def side_demand(normal_demand_mL_min,positive_work_W,normal_positive_work_W,rate_bpm,contractility,assumptions=DemandAssumptions(),normal_rate_bpm=80.):
    values=(normal_demand_mL_min,positive_work_W,normal_positive_work_W,rate_bpm,contractility,normal_rate_bpm)
    if any(not math.isfinite(x) for x in values) or min(normal_demand_mL_min,positive_work_W,rate_bpm)<0 or min(normal_positive_work_W,contractility,normal_rate_bpm)<=0:raise ValueError('Invalid normal reference or operating point')
    a=assumptions
    basal=normal_demand_mL_min*a.maintenance_fraction*a.mass_ratio
    activation=normal_demand_mL_min*a.activation_fraction*a.mass_ratio*(rate_bpm/normal_rate_bpm)*contractility**a.contractility_exponent
    work=normal_demand_mL_min*a.work_fraction*(positive_work_W/normal_positive_work_W)*a.work_cost_multiplier
    return dict(maintenance_mL_min=basal,activation_mL_min=activation,work_mL_min=work,total_mL_min=basal+activation+work)
