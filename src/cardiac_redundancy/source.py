"""Source-only six-volume/four-flow model, excluding project extensions.

Valve modes are inputs to rhs. Hybrid switching and periodic integration are
intentionally separate: rhs alone is not a complete healthy-valve simulator.
"""
import json
import math
from pathlib import Path

import numpy as np

from .mechanics import (ModelDomainError, Wall, driver, inertial_derivative,
                        resistive_flow, septal_root, vascular_pressure)

# These tuples define array order everywhere: changing them changes the state
# contract, incidence rows and saved output interpretation. Configuration dict
# insertion order is intentionally not used to construct the state vector.
VOLUMES = ("LV", "RV", "AO", "VC", "PA", "PU")
VALVES = ("MT", "AV", "TC", "PV")
# Columns 0..3 are mitral/aortic/tricuspid/pulmonary valves; the last two are
# systemic and pulmonary vascular resistances. All endpoints are reservoirs.
EDGES = (("PU", "LV"), ("LV", "AO"), ("VC", "RV"), ("RV", "PA"),
         ("AO", "VC"), ("PA", "PU"))
INCIDENCE = np.zeros((6, 6))
# One edge removes Q from its donor and adds Q to its recipient. Every column
# sums to zero, so sum(dV/dt)=0 for arbitrary signed flows, not just at rest.
for column, (donor, recipient) in enumerate(EDGES):
    INCIDENCE[VOLUMES.index(donor), column] = -1
    INCIDENCE[VOLUMES.index(recipient), column] = 1
INCIDENCE.flags.writeable = False


class SourceModel:
    """Assemble the source-only model from the frozen JSON-subset YAML file.

    State order:
        [V_LV, V_RV, V_AO, V_VC, V_PA, V_PU, Q_MT, Q_AV, Q_TC, Q_PV].
    The first six entries are mL, the last four mL/s. The source's 1,500 mL
    coordinate sum is not the project's 5,250 mL whole-body blood inventory.
    Four open/closed valve flags live outside this continuous state vector.
    """

    def __init__(self, path):
        """Load source parameters, wall laws, valve coefficients and initial state.

        Loading checks the wall/valve parameters and volume sum used here;
        it is not a general schema validator for every possible configuration.
        `source_parameters.heart_parameters.V_tot` is deliberately unused, as
        resolved in M2. No unit conversion is applied to the corrected values.
        """
        self.config = json.loads(Path(path).read_text(encoding="utf-8-sig"))
        self.parameters = p = self.config["source_parameters"]
        self.settings = self.config["solver_targets"]
        self.walls = []  # Order must match septal_root(left, right, septum).
        for component, suffix in (("lvf_calculator", "lvf"), ("rvf_calculator", "rvf"), ("septum", "spt")):
            w = p[component]
            self.walls.append(Wall(w["E_es_" + suffix], w["V_d_" + suffix],
                                   w["P_0_" + suffix], w["lambda_" + suffix], w["V_0_" + suffix]))
        h = p["heart_parameters"]
        self.resistances = np.array([h["R_" + v.lower()] for v in VALVES])
        self.inertances = np.array([h["L_" + v.lower()] for v in VALVES])
        if not np.all(np.isfinite(np.r_[self.resistances, self.inertances])) or min(*self.resistances, *self.inertances) <= 0:
            raise ModelDomainError("Invalid valve parameters")
        self.initial_state = np.array([self.config["initial_volumes_mL"][v] for v in VOLUMES]
                                      + [self.config["initial_valve_flows_mL_s"][v] for v in VALVES], dtype=float)
        expected = self.config["circulating_coordinate_sum_mL"]
        if not math.isclose(self.initial_state[:6].sum(), expected, abs_tol=1e-9, rel_tol=0):
            raise ModelDomainError("Initial volume inventory mismatch")

    def pressures(self, t, state, active=True, septal_hint=0.0):
        """Compute six reservoir pressures and diagnostics without advancing time.

        Returns (pressure_by_node, diagnostics). Pressures/residual are mmHg,
        septal coordinate is mL, and activation is dimensionless. `active=False`
        removes the common active driver from both free walls AND the septum;
        passive wall pressure, pericardium and vascular pressure remain.

        Inputs are not mutated. `septal_hint` is optional search guidance,
        not a persistent septal state. Nonpositive reservoir/cavity volumes
        cause an exception here; this is not event localization of a boundary.
        """
        state = np.asarray(state)
        if state.shape != (10,) or not np.all(np.isfinite(state)) or np.any(state[:6] <= 0):
            raise ModelDomainError("Nonfinite state, wrong state shape or nonpositive cavity volume")
        d = self.config["driver"]
        activation = driver(t, d["period_s"], d["A"], d["B_per_s2"], d["C_s"]) if active else 0.0
        lv, rv = state[:2]
        septal, residual = septal_root(lv, rv, activation, *self.walls,
                                      tolerance=self.settings["septal_residual_tolerance_mmHg"], hint=septal_hint)
        # Both chambers share an external pericardial pressure depending on
        # their combined cavity volume; this is direct mechanical interaction.
        pcd = self.parameters["pericardium"]
        ext = self.config["pressure_external_mmHg"]
        try:
            peri = ext["pericardium"] + pcd["P_0_pcd"] * math.expm1(pcd["lambda_pcd"] * (lv + rv - pcd["V_0_pcd"]))
        except OverflowError as exc:
            raise ModelDomainError("Pericardial exponential overflow") from exc
        # Wall laws give pressure relative to the pericardium. Add the common
        # surrounding pressure once to obtain chamber pressures for flow laws.
        pressures = {"LV": self.walls[0].pressure(lv - septal, activation) + peri,
                     "RV": self.walls[1].pressure(rv + septal, activation) + peri}
        # The source has zero unstressed offsets. PA/PU include the thoracic
        # pressure offset; AO/VC do not. Project offsets are not loaded here.
        for node, component in (("AO", "aorta"), ("VC", "vena_cava"), ("PA", "pulmonary_artery"), ("PU", "pulmonary_vein")):
            pars = self.parameters[component]
            pressures[node] = vascular_pressure(state[VOLUMES.index(node)], pars["E_es_" + node.lower()], ext[node], pars["V_d_" + node.lower()])
        if not all(math.isfinite(x) for x in pressures.values()):
            raise ModelDomainError("Nonfinite pressure")
        return pressures, {"activation": activation, "septal_coordinate_mL": septal,
                           "septal_residual_mmHg": residual, "pericardial_pressure_mmHg": peri}

    def rhs(self, t, state, open_valves, active=True):
        """Return the ten derivatives dy/dt for a fixed valve mode.

        The solver calls this many times at trial states. Entries 0..5 are
        dV/dt [mL/s]; 6..9 are dQ/dt [mL/s^2]. `open_valves` must contain four
        flags in VALVES order. This method neither detects a valve event nor
        changes the supplied state/modes; integration.py owns those actions.
        """
        p, _ = self.pressures(t, state, active=active)
        modes = np.asarray(open_valves, dtype=bool)
        if modes.shape != (4,):
            raise ModelDomainError("Four explicit valve modes required")
        # Closed flow coordinates are constrained inactive DOFs. Implicit solvers
        # perturb them for numerical Jacobians; they must not transfer blood.
        # This evaluates the mode's vector field, not a reset of stored states.
        flows = np.where(modes, np.asarray(state[6:]), 0.0)
        h = self.parameters["heart_parameters"]
        # np.r_ concatenates arrays: four dynamic valve flows plus two flows
        # computed instantly from vascular pressure differences.
        network_flows = np.r_[flows, resistive_flow(p["AO"], p["VC"], h["R_sys"]),
                              resistive_flow(p["PA"], p["PU"], h["R_pul"])]
        dq = np.zeros(4)
        for i, (a, b) in enumerate(EDGES[:4]):
            if modes[i]:
                dq[i] = inertial_derivative(p[a] - p[b], flows[i], self.resistances[i], self.inertances[i])
        # @ is matrix multiplication. This assembles all six inflow-minus-
        # outflow equations with shared edge values, avoiding lost/doubled blood.
        return np.r_[INCIDENCE @ network_flows, dq]
