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

A floor on circulating glucose is a declared axis rather than a fixed choice.
With no floor, duration means time to complete exhaustion of the pool, which
is not a physiological endpoint. The evidence does not supply a single cutoff:
cerebral glucose transport begins to limit utilisation near 3 mmol/L, the
ADA/EASD joint statement sets 3.0 mmol/L (54 mg/dL) as the level to report in
trials, measured symptom thresholds cluster at 3 mmol/L, and cognitive
impairment appears around 3.1-3.4 mmol/L, while transport half-saturation
estimates range from about 1.5 to 6.2 mmol/L with wide uncertainty. So the
duration is reported against the floor as a curve, with 3.0 mmol/L available as
an evidence-anchored reference point, and the default of zero preserves the
model's own behaviour.

This is an estimate, not a result. It assumes constant demand and export
already at its cap.
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
    #: Plasma glucose below which the pool is treated as unavailable. Zero
    #: reproduces the model's own behaviour (duration to exhaustion); 3.0 is
    #: the evidence-anchored reference. Declared, never tuned to a result.
    glucose_floor_mmol_L: float = 0.
    #: Hepatic export is capped at release_cap TIMES this fraction, as in
    #: GlucoseTransport.evaluate. Every M5, M7 and M8 run has used 0.5, so the
    #: effective cap is 0.105 mmol/min, not 0.21.
    liver_function_fraction: float = .5
    #: Share of oxidative demand charged to glucose. The model charges all of
    #: it, which is why this defaults to 1.0, but resting metabolism runs near
    #: a respiratory quotient of 0.8, so roughly a third of energy comes from
    #: carbohydrate and the rest from fat. The model has no fatty-acid or
    #: ketone pathway, so this is reported as a declared axis rather than
    #: implemented: below about 0.13 hepatic export alone covers demand and
    #: the circulating pool never drains at all.
    glucose_fraction: float = 1.

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
        if self.glucose_floor_mmol_L >= self.plasma_glucose_mmol_L:
            raise ValueError('Glucose floor must lie below the initial plasma concentration')
        if not 0 <= self.liver_function_fraction <= 1:
            raise ValueError('Liver function fraction must lie in [0, 1]')
        if not 0 < self.glucose_fraction <= 1:
            raise ValueError('Glucose fraction must lie in (0, 1]')

    @property
    def circulating_pool_mmol(self):
        """Whole-blood glucose above the declared floor, from plasma-equivalent."""
        usable = self.plasma_glucose_mmol_L-self.glucose_floor_mmol_L
        return usable/self.plasma_to_whole_blood_ratio*self.blood_volume_L*self.accessible_fraction

    @property
    def effective_export_mmol_min(self):
        """What the release law actually permits: cap times liver function."""
        return self.release_cap_mmol_min*self.liver_function_fraction

    @property
    def total_pool_mmol(self):
        """The whole pool, ignoring the floor; duration to exhaustion uses this."""
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
    glucose_demand = budget.glucose_fraction*demand_mL_O2_min/O2_ML_STPD_PER_MMOL_GLUCOSE
    export = budget.effective_export_mmol_min
    drain = glucose_demand-export
    pool = budget.circulating_pool_mmol
    record = dict(demand_mL_O2_min=float(demand_mL_O2_min),
                  glucose_demand_mmol_min=glucose_demand,
                  release_cap_mmol_min=budget.release_cap_mmol_min,
                  liver_function_fraction=budget.liver_function_fraction,
                  effective_export_mmol_min=export, glucose_fraction=budget.glucose_fraction,
                  export_covers_mL_O2_min=export*O2_ML_STPD_PER_MMOL_GLUCOSE,
                  net_drain_mmol_min=drain, circulating_pool_mmol=pool,
                  glycogen_total_mmol=budget.glycogen_total_mmol,
                  store_limited_minutes=(None if export == 0
                                         else budget.glycogen_total_mmol/export),
                  glucose_floor_mmol_L=budget.glucose_floor_mmol_L,
                  total_pool_mmol=budget.total_pool_mmol,
                  basis='first-order: constant demand, export at cap, declared glucose floor',
                  interpretation='interval over which the declared delivery screens can be met; '
                                 'not survival, recovery or injury')
    if drain <= 0:
        record.update(minutes=None, limited_by='not substrate-limited: hepatic export covers demand',
                      margin_mmol_min=-drain)
        return record
    record.update(minutes=pool/drain, minutes_to_exhaustion=budget.total_pool_mmol/drain,
                  limited_by=('circulating glucose pool against hepatic export cap'
                              + ('' if budget.glucose_floor_mmol_L == 0 else
                                 ', down to the declared %g mmol/L floor' % budget.glucose_floor_mmol_L)))
    return record


#: Floors worth reporting. 3.0 mmol/L is the ADA/EASD trial-reporting level, the
#: measured symptom threshold in several studies, and roughly where cerebral
#: glucose transport begins to limit utilisation.
REFERENCE_FLOORS_MMOL_L = (0., 2.0, 2.5, 3.0, 3.5)


def duration_curve(demand_mL_O2_min, budget=SubstrateBudget(), floors=REFERENCE_FLOORS_MMOL_L):
    """Supported duration against the declared glucose floor.

    Reported as a curve because the evidence supplies a range rather than a
    cutoff, and because the choice moves the answer by a factor of three. A
    single number here would hide the most consequential assumption in it.
    """
    from dataclasses import replace
    rows = []
    for floor in floors:
        if floor >= budget.plasma_glucose_mmol_L:
            continue
        r = supported_duration(demand_mL_O2_min, replace(budget, glucose_floor_mmol_L=floor))
        rows.append(dict(glucose_floor_mmol_L=floor, minutes=r['minutes'],
                         usable_pool_mmol=r['circulating_pool_mmol']))
    return dict(demand_mL_O2_min=float(demand_mL_O2_min), rows=rows,
                anchor='3.0 mmol/L: ADA/EASD trial-reporting level and measured symptom '
                       'threshold; cerebral transport begins to limit near this value',
                interpretation='the floor is a declared axis, not a fitted parameter')


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
