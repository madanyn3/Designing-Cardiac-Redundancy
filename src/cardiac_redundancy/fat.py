"""Arm B: perfusion-limited fatty-acid supply and the obligate-glucose share.

Arm A charges every joule of oxidative demand to glucose, which is why its
supported duration is short. Real resting metabolism draws most of its oxidative
ATP from fatty acid, and the heart in particular prefers it. Arm B asks what
changes when tissue may oxidise fat, and the answer turns out to depend on
circulation rather than on fuel stores: adipose tissue must be perfused well
enough to carry the fatty acid away.

Two limits apply, and which one binds is the finding rather than an input:

  * **lipolytic capacity** -- how fast adipose can release fatty acid at all;
  * **perfusion transport** -- adipose blood flow multiplied by the
    arterio-venous concentration difference the blood can sustain.

Both are declared assumptions. The transport term is the one that makes adipose
conductance a design variable, because the auxiliary pump chooses how much flow
the depot receives.

One tissue is excluded on physiological grounds rather than by assumption: the
brain cannot oxidise long-chain fatty acid, which does not cross the blood-brain
barrier in useful quantity. It can run substantially on ketone bodies, but
ketogenesis takes days of fasting to raise plasma ketones enough to matter, and
this study's horizon is minutes to hours. So the brain's demand stays obligately
glucose here, and a ketone contribution is a separate long-horizon axis, not a
term in this calculation.

Nothing here models fatty-acid transport as a conserved species. It is a supply
controller in the same style as the hepatic glucose export law, which the M4
specification likewise describes as bookkeeping rather than hormone kinetics.

`FuelOxygenTransport` already accepts a per-organ `glucose_share` and records
`other_substrates_unassessed` whenever any share falls below one. Arm A leaves
every share at one. Arm B sets them from a declared and measured fat supply,
which turns that flag from an unassessed gap into an assessed one: the
non-glucose fraction is now backed by a perfusion-limited supply rather than
assumed available.
"""
from dataclasses import dataclass
import math

#: Complete oxidation of palmitate consumes 23 mmol O2 per mmol, so one mmol of
#: fatty acid carries 23 x 22.414 mL O2 STPD of oxidative capacity.
O2_ML_STPD_PER_MMOL_FFA = 23*22.414


@dataclass(frozen=True)
class LipolysisParameters:
    """Declared fatty-acid supply assumptions. None is a measurement here."""
    #: Whole-body lipolytic rate. Reported palmitate turnover in fasting humans
    #: spans roughly 1.4 to 2.7 umol/kg/min; the upper end is taken as capacity.
    whole_body_capacity_umol_kg_min: float = 2.7
    body_mass_kg: float = 70.
    #: Arterio-venous fatty-acid difference the depot can sustain. Plasma
    #: concentrations run about 0.5 to 1.0 mmol/L fasting and reach 1.5 to 2
    #: in prolonged fasting, so a 1.0 mmol/L difference is a mid assumption.
    #: This parameter wants a source before any result rests on it.
    arteriovenous_difference_mmol_L: float = 1.
    #: Fraction of the depot's lipolytic capacity that is accessible.
    accessible_fraction: float = 1.

    def __post_init__(self):
        for name, value in vars(self).items():
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f'{name} must be a finite number')
            if value < 0:
                raise ValueError(f'{name} must be nonnegative')
        if self.body_mass_kg <= 0:
            raise ValueError('Positive body mass required')
        if not 0 <= self.accessible_fraction <= 1:
            raise ValueError('Accessible fraction must lie in [0, 1]')

    @property
    def lipolytic_capacity_mmol_min(self):
        return (self.whole_body_capacity_umol_kg_min*self.body_mass_kg/1000.
                * self.accessible_fraction)


def fat_supply(adipose_flow_mL_min, parameters=LipolysisParameters()):
    """Fatty acid the depot can actually deliver, and which limit binds.

    `adipose_flow_mL_min` is the perfusion the circulation delivers to the
    depot, which in the supported state is the branch's reference flow scaled
    by whatever conductance factor its role receives.
    """
    if (isinstance(adipose_flow_mL_min, bool)
            or not isinstance(adipose_flow_mL_min, (int, float))
            or not math.isfinite(adipose_flow_mL_min) or adipose_flow_mL_min < 0):
        raise ValueError('Finite nonnegative adipose flow required')
    if not isinstance(parameters, LipolysisParameters):
        raise ValueError('Validated lipolysis parameters required')
    lipolytic = parameters.lipolytic_capacity_mmol_min
    transport = adipose_flow_mL_min*parameters.arteriovenous_difference_mmol_L/1000.
    supply = min(lipolytic, transport)
    return dict(adipose_flow_mL_min=float(adipose_flow_mL_min),
                lipolytic_capacity_mmol_min=lipolytic,
                perfusion_transport_mmol_min=transport,
                supply_mmol_min=supply,
                supply_mL_O2_eq_min=supply*O2_ML_STPD_PER_MMOL_FFA,
                binding=('adipose perfusion' if transport < lipolytic
                         else 'lipolytic capacity'),
                interpretation='declared supply controller, not a conserved fatty-acid '
                               'species and not hormone kinetics')


def fuel_split(total_demand_mL_O2_min, brain_demand_mL_O2_min, supply):
    """Split oxidative demand into the part fat can cover and the rest.

    The brain's share stays on glucose because long-chain fatty acid does not
    reach it; every other tissue may use fat up to what the depot delivers.
    Returns the glucose fraction the substrate budget should then be charged.
    """
    for name, value in (('total demand', total_demand_mL_O2_min),
                        ('brain demand', brain_demand_mL_O2_min)):
        if (isinstance(value, bool) or not isinstance(value, (int, float))
                or not math.isfinite(value) or value < 0):
            raise ValueError(f'Finite nonnegative {name} required')
    if brain_demand_mL_O2_min > total_demand_mL_O2_min:
        raise ValueError('Brain demand cannot exceed total demand')
    available = supply['supply_mL_O2_eq_min'] if isinstance(supply, dict) else float(supply)
    non_brain = total_demand_mL_O2_min-brain_demand_mL_O2_min
    from_fat = min(available, non_brain)
    glucose_demand = brain_demand_mL_O2_min+(non_brain-from_fat)
    return dict(total_demand_mL_O2_min=float(total_demand_mL_O2_min),
                brain_demand_mL_O2_min=float(brain_demand_mL_O2_min),
                non_brain_demand_mL_O2_min=non_brain,
                fat_available_mL_O2_eq_min=available,
                met_by_fat_mL_O2_eq_min=from_fat,
                non_brain_covered_fraction=(from_fat/non_brain) if non_brain > 0 else 1.,
                glucose_demand_mL_O2_min=glucose_demand,
                glucose_fraction=(glucose_demand/total_demand_mL_O2_min
                                  if total_demand_mL_O2_min > 0 else 1.),
                floor_fraction=(brain_demand_mL_O2_min/total_demand_mL_O2_min
                                if total_demand_mL_O2_min > 0 else 1.),
                interpretation='the brain share is obligately glucose on this horizon; '
                               'ketone adaptation takes days and is a separate axis')


#: Tissues that cannot oxidise long-chain fatty acid in this model. Red cells
#: and renal medulla are also obligate glucose consumers physiologically, but
#: neither is represented as a separate organ here, so the brain is the only
#: one the model can hold to the constraint.
OBLIGATE_GLUCOSE_ORGANS = ('brain',)


def glucose_shares(organ_names, demands_mL_O2_min, supply,
                   obligate=OBLIGATE_GLUCOSE_ORGANS):
    """Per-organ fraction of oxidative demand charged to glucose.

    Fat is allocated across the non-obligate organs in proportion to their
    demand, so every one of them receives the same fractional coverage. A
    priority ordering would be a further assumption; proportional allocation
    makes none beyond the supply itself.

    Returns the array `FuelOxygenTransport` expects, plus the accounting that
    produced it, so a run records why each share has the value it does.
    """
    names = list(organ_names)
    demands = [float(d) for d in demands_mL_O2_min]
    if len(names) != len(demands):
        raise ValueError('One demand per organ required')
    if any(d < 0 or not math.isfinite(d) for d in demands):
        raise ValueError('Finite nonnegative organ demands required')
    available = supply['supply_mL_O2_eq_min'] if isinstance(supply, dict) else float(supply)
    if not math.isfinite(available) or available < 0:
        raise ValueError('Finite nonnegative fat supply required')
    obligate = set(obligate)
    non_obligate = sum(d for n, d in zip(names, demands) if n not in obligate)
    covered = min(available, non_obligate)
    fraction = (covered/non_obligate) if non_obligate > 0 else 0.
    shares = [1. if n in obligate else 1.-fraction for n in names]
    return dict(shares=shares, organ_names=names,
                obligate_organs=sorted(obligate & set(names)),
                non_obligate_demand_mL_O2_min=non_obligate,
                fat_available_mL_O2_eq_min=available,
                met_by_fat_mL_O2_eq_min=covered,
                covered_fraction=fraction,
                allocation='proportional to demand across non-obligate organs',
                interpretation='shares are assessed against a perfusion-limited supply, '
                               'not assumed; the obligate organs stay wholly on glucose')
