"""Outcome policy is independent of tissue demand and hydraulic restriction."""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class PreservationChoice:
    demand_multiplier: float = 1.
    conductance_factor: float = 1.
    outcome_required: bool = True

    def __post_init__(self):
        if not math.isfinite(self.demand_multiplier) or self.demand_multiplier<0:raise ValueError('Invalid demand multiplier')
        if not math.isfinite(self.conductance_factor) or not 0<=self.conductance_factor<=1:raise ValueError('Invalid conductance factor')
        if not isinstance(self.outcome_required,bool):raise ValueError('Explicit outcome boolean required')


def apply_policy(organs,choices):
    """Return explicit demand vector, outcome mask and separate control requests.

    Controls are requests, not automatically applied to existing hydraulic
    data. Existing limb names denote pooled limbs, never lower limbs alone.
    """
    names={o.name for o in organs}
    if set(choices)!=names:raise ValueError('Every modeled organ requires an explicit choice')
    return dict(demand_mL_min=[o.demand_mL_min*choices[o.name].demand_multiplier for o in organs],
        outcome_required=[choices[o.name].outcome_required for o in organs],
        requested_conductance={o.name:choices[o.name].conductance_factor for o in organs})
