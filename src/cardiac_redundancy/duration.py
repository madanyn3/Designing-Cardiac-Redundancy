"""First-order substrate-limited supported duration.

The project's headline endpoint is how long organ perfusion can be sustained,
not whether a single instant satisfies the screens. This module gives the
closed-form bound that sets the time scale, so a long simulation is verified
against an expectation instead of being the only source of the number.

The structure of the limit: oxidative demand must be met by glucose, hepatic
export is capped by a rate rather than by the size of the glycogen store, and
the difference is drawn from the circulating glucose pool. Duration is
therefore pool divided by net drain, and the store never empties -- the
reference liver holds 330 mmol, which at the 0.21 mmol/min cap would last
over a day. The cap, not the store, is what binds.

This is an estimate, not a result. It assumes constant demand, export already
at its cap, and a fully accessible pool with no declared lower glucose limit.
Transport lags (release_time_s, blood_target_time_s), demand transients during
startup and any physiological floor on blood glucose all shorten it. It must
be verified against a coupled long-horizon replay before being reported, and
it is a substrate bound only: it assumes ventilation continues to supply
oxygen, so it says nothing about an apneic or gas-exchange-limited case.

No claim about survival, recovery or injury is made. A duration here is the
interval over which the model's declared delivery screens can be met.
"""
from dataclasses import dataclass
import math
from .fuel import O2_ML_STPD_PER_MMOL_GLUCOSE


@dataclass(frozen=True)
class SubstrateBudget:
    """Declared substrate accounting for one scenario.

    Defaults are the M4 reference parameters with the M3 physical inventory.
    """
    release_cap_mmol_min: float = .21
    plasma_glucose_mmol_L: float = 5.
    plasma_to_whole_blood_ratio: float = 1.11
    blood_volume_L: float = 5.25
    glycogen_mmol_L_liver: float = 229.
    liver_volume_L: float = 1.44
    accessible_fraction: float = 1.

    def __post_init__(self):
        for name, value in vars(self).items():
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f'{name} must be a finite number')
            if value < 0:
                raise ValueError(f'{name} must be nonnegative')
        if min(self.plasma_to_whole_blood_ratio, self.blood_volume_L) <= 0:
            raise ValueError('Positive blood volume and plasma ratio required')
        if not 0 <= self.accessible_fraction <= 1:
            raise ValueError('Accessible fraction must lie in [0, 1]')

    @property
    def circulating_pool_mmol(self):
        """Whole-blood glucose available, from the plasma-equivalent reading."""
        return (self.plasma_glucose_mmol_L/self.plasma_to_whole_blood_ratio
                * self.blood_volume_L*self.accessible_fraction)

    @property
    def glycogen_total_mmol(self):
        return self.glycogen_mmol_L_liver*self.liver_volume_L*self.accessible_fraction


def supported_duration(demand_mL_O2_min, budget=SubstrateBudget()):
    """Minutes until the circulating glucose pool cannot cover the shortfall.

    Returns a record rather than a bare number, because the interesting cases
    are the ones with no finite answer: when hepatic export alone covers
    demand, the model is not substrate-limited at all and the duration is
    bounded by something this calculation does not describe.
    """
    if (isinstance(demand_mL_O2_min, bool) or not isinstance(demand_mL_O2_min, (int, float))
            or not math.isfinite(demand_mL_O2_min) or demand_mL_O2_min < 0):
        raise ValueError('Finite nonnegative oxidative demand required')
    if not isinstance(budget, SubstrateBudget):
        raise ValueError('Validated substrate budget required')
    glucose_demand = demand_mL_O2_min/O2_ML_STPD_PER_MMOL_GLUCOSE
    drain = glucose_demand-budget.release_cap_mmol_min
    pool = budget.circulating_pool_mmol
    record = dict(demand_mL_O2_min=float(demand_mL_O2_min),
                  glucose_demand_mmol_min=glucose_demand,
                  release_cap_mmol_min=budget.release_cap_mmol_min,
                  export_covers_mL_O2_min=budget.release_cap_mmol_min*O2_ML_STPD_PER_MMOL_GLUCOSE,
                  net_drain_mmol_min=drain, circulating_pool_mmol=pool,
                  glycogen_total_mmol=budget.glycogen_total_mmol,
                  store_limited_minutes=(None if budget.release_cap_mmol_min == 0
                                         else budget.glycogen_total_mmol/budget.release_cap_mmol_min),
                  basis='first-order: constant demand, export at cap, no glucose floor',
                  interpretation='interval over which the declared delivery screens can be met; '
                                 'not survival, recovery or injury')
    if drain <= 0:
        record.update(minutes=None, limited_by='not substrate-limited: hepatic export covers demand',
                      margin_mmol_min=-drain)
        return record
    record.update(minutes=pool/drain, limited_by='circulating glucose pool against hepatic export cap')
    return record


def duration_sensitivity(demand_mL_O2_min, budget=SubstrateBudget(), factor=.1):
    """How a fractional change in demand or in export moves the duration.

    Reports the two elasticities side by side because they are the quantities a
    design and a physiology respectively control: a pump can change demand, but
    only the declared scenario changes the export cap.
    """
    base = supported_duration(demand_mL_O2_min, budget)
    if base['minutes'] is None:
        return dict(base=base, demand_elasticity=None, export_elasticity=None)
    from dataclasses import replace
    higher = supported_duration(demand_mL_O2_min*(1+factor), budget)
    exported = supported_duration(demand_mL_O2_min,
                                  replace(budget, release_cap_mmol_min=budget.release_cap_mmol_min*(1+factor)))
    elasticity = lambda other: ((other['minutes']-base['minutes'])/base['minutes']/factor
                                if other['minutes'] is not None else None)
    return dict(base=base, factor=factor,
                demand_elasticity=elasticity(higher), export_elasticity=elasticity(exported),
                interpretation='fractional change in duration per fractional change in each input')
