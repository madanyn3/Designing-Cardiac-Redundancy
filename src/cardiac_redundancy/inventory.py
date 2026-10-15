"""Physical blood accounting for the project circuit, separate from source coordinates.

V is physical blood volume; U is the unstressed offset. V-U determines pressure
in a linear vascular reservoir, but all of V participates in blood conservation
and later oxygen mixing. The septal geometric coordinate is not a blood pool.
This module does not add regional flow laws, pumps or event localization.
"""
from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np

from .source import VOLUMES, SourceModel
from .mechanics import ModelDomainError


@dataclass(frozen=True)
class Reservoir:
    name: str
    initial_mL: float
    unstressed_mL: float
    compliance_mL_mmHg: float | None
    external_mmHg: float


class PhysicalSourceReference:
    """Pre-split 5,250 mL reference: exact constant translation of source states.

    Recover the before-split VC offset by returning all limb-pool U to VC.
    No regional beds exist in this fixture. V-U and all four valve flows equal
    their source values, so its ODE solution must be identical after translation.
    This is a verification fixture, not the final eleven-reservoir circuit.
    """

    def __init__(self, source_path, project_path):
        inventory = PhysicalInventory(source_path, project_path)
        self.source = SourceModel(source_path)
        self.settings = self.source.settings
        offsets = inventory.unstressed_volumes[:6].copy()
        offsets[VOLUMES.index("VC")] += inventory.unstressed_volumes[6:].sum()
        self.translation = np.r_[offsets,np.zeros(4)]
        self.initial_state = self.source.initial_state+self.translation
        if not np.isclose(self.initial_state[:6].sum(),inventory.total_mL,rtol=0,atol=1e-9):
            raise ValueError("Pre-split physical reference does not match M2 total")

    def pressures(self,t,state,active=True):
        return self.source.pressures(t,np.asarray(state)-self.translation,active=active)

    def rhs(self,t,state,open_valves,active=True):
        # U is constant: d(V-U)/dt = dV/dt. No source RHS modification needed.
        return self.source.rhs(t,np.asarray(state)-self.translation,open_valves,active=active)


class PhysicalInventory:
    """Validated eleven-reservoir layout from the unchanged M2 configuration.

Central order matches SourceModel, followed by the five named limb pools.
LV/RV have physical cavity volumes but no linear vascular compliance: their
nonlinear wall/pericardium laws must be evaluated by the cardiac model.
"""

    def __init__(self, source_path, project_path):
        source = json.loads(Path(source_path).read_text(encoding="utf-8-sig"))
        project = json.loads(Path(project_path).read_text(encoding="utf-8-sig"))
        self.total_mL = float(project["total_blood_volume_mL"])
        if not np.isfinite(self.total_mL) or self.total_mL <= 0:
            raise ValueError("Positive finite physical inventory required")
        if not np.isclose(self.total_mL, project["target_blood_volume_mL"], rtol=0, atol=1e-9):
            raise ValueError("Physical inventory targets disagree")
        initial = project["initial_physical_volumes_mL"]
        offsets = project["unstressed_offsets_mL"]
        pools = project["limb_pools"]
        self.names = VOLUMES + tuple(pools)
        if len(set(self.names)) != len(self.names) or set(initial) != set(self.names) or set(offsets) != set(VOLUMES):
            raise ValueError("Physical reservoir names/offsets disagree")
        components = {"AO":"aorta", "VC":"vena_cava", "PA":"pulmonary_artery", "PU":"pulmonary_vein"}
        reservoirs = []
        for name in VOLUMES:
            compliance = None
            if name in components:
                elastance = source["source_parameters"][components[name]]["E_es_"+name.lower()]
                if not np.isfinite(elastance) or elastance <= 0:
                    raise ValueError("Positive finite vascular elastance required")
                compliance = 1/elastance
            elif offsets[name] != 0:
                raise ValueError("M2 native cavity offsets must remain zero")
            external = source["pressure_external_mmHg"].get(name, 0)
            reservoirs.append(Reservoir(name, initial[name], offsets[name], compliance, external))
        for name, pool in pools.items():
            if not np.isclose(initial[name], pool["volume_initial_mL"], rtol=0, atol=1e-9):
                raise ValueError("Duplicate limb initial volumes disagree")
            reservoirs.append(Reservoir(name, initial[name], pool["unstressed_mL"],
                                        pool["compliance_mL_mmHg"], pool["external_pressure_mmHg"]))
        for r in reservoirs:
            if not np.all(np.isfinite([r.initial_mL,r.unstressed_mL,r.external_mmHg])) or r.unstressed_mL < 0:
                raise ValueError("Invalid reservoir volume/offset/pressure")
            if r.compliance_mL_mmHg is not None and (not np.isfinite(r.compliance_mL_mmHg) or r.compliance_mL_mmHg <= 0):
                raise ValueError("Positive finite compliance required")
        self.reservoirs = tuple(reservoirs)
        self.validate(self.initial_volumes)
        if not np.isclose(self.initial_volumes.sum(), self.total_mL, rtol=0, atol=1e-9):
            raise ValueError("Physical blood inventory does not sum to target")

    @property
    def initial_volumes(self):
        """Fresh array, so integration cannot mutate configuration state."""
        return np.array([r.initial_mL for r in self.reservoirs],dtype=float)

    @property
    def unstressed_volumes(self):
        return np.array([r.unstressed_mL for r in self.reservoirs],dtype=float)

    def coordinate_volumes(self, physical):
        """Translate volume coordinates; never remove U from mass accounting."""
        physical = np.asarray(physical,dtype=float)
        if physical.shape != (len(self.names),) or not np.all(np.isfinite(physical)):
            raise ModelDomainError("Wrong physical volume shape or nonfinite value")
        return physical-self.unstressed_volumes

    def domain_margins(self, physical):
        """Boundary distances in mL: cavity V; vascular V-U.

        Algebraic margins also accept out-of-domain trial values for future
        event localization. Negative margins are not silently clipped.
        """
        coordinates = self.coordinate_volumes(physical)
        return dict(zip(self.names, map(float,coordinates)))

    def validate(self, physical):
        margins = self.domain_margins(physical)
        if np.any(np.asarray(physical) <= 0) or min(margins.values()) < 0:
            raise ModelDomainError("Nonpositive physical volume or negative vascular transmural pressure")

    def vascular_pressures(self, physical):
        """Absolute vascular pressures; LV/RV require the cardiac wall laws.

        Evaluation is algebraic, allowing trial states outside the validity
        domain. Call validate at initialization; integration must localize exits.
        """
        v = self.coordinate_volumes(physical)
        return {r.name:r.external_mmHg+v[i]/r.compliance_mL_mmHg
                for i,r in enumerate(self.reservoirs) if r.compliance_mL_mmHg is not None}

    def incidence(self, edges):
        """Each directed edge removes and adds exactly the same signed flow.

        Closing a route means Q=0, not deleting its reservoir. Negative Q
        reverses transfer. Zero-storage junctions must first be eliminated by
        their flow-balance equations; they are not physical reservoirs here.
        """
        matrix = np.zeros((len(self.names),len(edges)))
        indices = {name:i for i,name in enumerate(self.names)}
        for j,(donor,recipient) in enumerate(edges):
            if donor not in indices or recipient not in indices or donor == recipient:
                raise ValueError("Distinct known physical endpoints required")
            matrix[indices[donor],j] = -1
            matrix[indices[recipient],j] = 1
        return matrix
