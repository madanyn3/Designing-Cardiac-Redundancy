"""T8.1 registry: optimizable design controls versus uncertain physiology.

Two separations are enforced here rather than left to each runner.

First, a search may choose only the controls in DesignPoint. Quantities in
ScenarioPoint describe uncertain or speculative physiology and are varied as
declared scenarios; letting a search choose them would lower the requirement
instead of meeting it.

Second, size is not free. The auxiliary cavity scale sets the pump's capacity,
so it must also set its myocardial mass and therefore its oxygen cost, and the
beat rate must enter the activation cost it is charged for. Both couplings are
explicit assumptions, with their own declared uncertainty, not measurements.

Bounds are search and scenario assumptions inherited from the M7 registry and
the verified parameter validators. They are not measured biological limits,
and admissibility here is never adequacy.
"""
from dataclasses import dataclass, fields
import math
from .biological_pump import BiologicalParameters
from .myocardial_demand import DemandAssumptions

#: Declared search/scenario intervals with provenance. Keys match the field
#: names below so a protocol file and the validators cannot drift apart.
BOUNDS = {
    'left_scale': ((.1, 1.), 'BiologicalParameters validator; capacity proxy, not anatomical mass'),
    'right_scale': ((.1, 1.), 'BiologicalParameters validator; capacity proxy, not anatomical mass'),
    'septal_scale': ((.1, 1.), 'BiologicalParameters validator'),
    'pericardial_scale': ((.1, 1.), 'BiologicalParameters validator'),
    'period_s': ((.4, 1.5), 'BiologicalParameters validator; 40-150 beats/min'),
    'left_contractility': ((.5, 2.), 'BiologicalParameters validator'),
    'right_contractility': ((.5, 2.), 'BiologicalParameters validator'),
    'septal_contractility': ((.5, 2.), 'BiologicalParameters validator'),
    'coronary_conductance': ((.5, 2.), 'M7 global factor bounds; anatomical/design hypothesis'),
    'preservation_conductance': ((.1, .5), 'M7 global factor bounds; vascular allocation control'),
    'extraction_limit': ((.25, .75), 'M7 global factor bounds; leading deficit effect'),
    'mass_coupling': ((.8, 1.6), 'declared uncertainty on specific myocardial mass per cavity scale'),
    'work_cost_multiplier': ((.75, 1.5), 'M7 global factor bounds'),
    'preservation_fraction': ((.1, .5), 'M7 global sampled range; speculative suppressed state'),
    'glycogen_mmol_L_liver': ((150., 400.), 'M1 inventory range around the 229 reference'),
    'release_cap_mmol_min': ((.1, .4), 'M1 release knob around the 0.21 reference'),
    'liver_function_fraction': ((.25, 1.), 'M5/M7 runs used 0.5; export capability is uncertain'),
}


def _check(obj):
    for f in fields(obj):
        x = getattr(obj, f.name)
        if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x):
            raise ValueError(f'{f.name} must be a finite number')
        lo, hi = BOUNDS[f.name][0]
        if not lo <= x <= hi:
            raise ValueError(f'{f.name} must be within [{lo}, {hi}]')


@dataclass(frozen=True)
class DesignPoint:
    """Controls a search may choose. Defaults are the M5/M7 reference fixture."""
    left_scale: float = .25
    right_scale: float = .25
    septal_scale: float = .25
    pericardial_scale: float = .25
    period_s: float = .75
    left_contractility: float = 1.
    right_contractility: float = 1.
    septal_contractility: float = 1.
    coronary_conductance: float = 1.
    preservation_conductance: float = .25

    def __post_init__(self):
        _check(self)
        BiologicalParameters(**self.biological_fields())

    def biological_fields(self):
        return {f.name: getattr(self, f.name) for f in fields(self)
                if f.name not in ('coronary_conductance', 'preservation_conductance')}

    def biological(self):
        return BiologicalParameters(**self.biological_fields())

    @property
    def rate_bpm(self):
        """Beats per minute implied by the chosen period; charged in demand."""
        return 60./self.period_s

    @property
    def cavity_scale(self):
        """Mean free-wall scale used as the capacity and mass proxy."""
        return (self.left_scale+self.right_scale)/2


@dataclass(frozen=True)
class ScenarioPoint:
    """Uncertain or speculative physiology. Declared, never optimized."""
    extraction_limit: float = .5
    mass_coupling: float = 1.
    work_cost_multiplier: float = 1.
    preservation_fraction: float = .25
    glycogen_mmol_L_liver: float = 229.
    release_cap_mmol_min: float = .21
    liver_function_fraction: float = .5

    def __post_init__(self):
        _check(self)


def mass_ratio(design, scenario):
    """Auxiliary myocardial mass relative to the native reference.

    Cavity volume scales with the chosen scale, so mass is charged in
    proportion, with mass_coupling carrying the uncertainty in specific mass.
    The reference fixture (scale 0.25, coupling 1.0) returns 0.25, reproducing
    the M5/M7 demand parameter exactly.
    """
    if not isinstance(design, DesignPoint) or not isinstance(scenario, ScenarioPoint):
        raise ValueError('Validated design and scenario required')
    value = design.cavity_scale*scenario.mass_coupling
    lo, hi = BOUNDS['mass_coupling'][0]
    if not .1*lo <= value <= 1.*hi:
        raise ValueError('Coupled mass ratio outside the declared domain')
    return value


def demand_assumptions(design, scenario, reference):
    """DemandAssumptions with the coupled mass ratio and declared work cost.

    `reference` supplies the frozen maintenance/activation/work split from
    configs/m5-myocardial-demand.json; its mass and cost entries are replaced
    by the coupled values so a larger or costlier pump cannot be free.
    """
    entries = {k: v for k, v in dict(reference).items()
               if k not in ('mass_ratio', 'work_cost_multiplier')}
    return DemandAssumptions(mass_ratio=mass_ratio(design, scenario),
                             work_cost_multiplier=scenario.work_cost_multiplier,
                             **entries)


def contractility(design, side):
    if side not in ('left', 'right'):
        raise ValueError('Side must be left or right')
    return getattr(design, side+'_contractility')


def registry():
    """Machine-readable two-column registry for the frozen protocol file."""
    design = {f.name for f in fields(DesignPoint)}
    scenario = {f.name for f in fields(ScenarioPoint)}
    if design & scenario or design | scenario != set(BOUNDS):
        raise ValueError('Registry columns must partition the declared bounds')
    return dict(
        design_controls={k: dict(bounds=list(BOUNDS[k][0]), provenance=BOUNDS[k][1],
                                 reference=getattr(DesignPoint(), k)) for k in sorted(design)},
        scenario_factors={k: dict(bounds=list(BOUNDS[k][0]), provenance=BOUNDS[k][1],
                                  reference=getattr(ScenarioPoint(), k)) for k in sorted(scenario)},
        couplings=dict(
            mass_ratio='cavity_scale * mass_coupling; charged in auxiliary oxygen demand',
            activation_cost='scales with rate_bpm = 60/period_s in side_demand',
            coronary_capacity='auxiliary coronary resistance already divides by cavity scale '
                              'in BiologicalPumpModel.flows; coronary_conductance multiplies it'),
        interpretation='search assumptions and declared scenarios; not measured limits')
