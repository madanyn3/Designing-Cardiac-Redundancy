"""Reproducible source-model knob experiments and periodic-cycle measurements.

Every knob operates on a deep copy. The frozen source remains untouched. A case
must warm up separately; none of these choices are physiological safe ranges.
"""
from copy import deepcopy

import numpy as np

from .integration import integrate_source
from .mechanics import ModelDomainError
from .source import VOLUMES


def apply_knobs(original, knobs):
    """Translate user-level knobs into authoritative source configuration fields.

    HR scales the Gaussian duration with T (B*T^2 and C/T stay fixed). Changing
    compliance divides elastance. A thoracic change updates all three actual
    pressure offsets. Initial redistribution conserves the coordinate inventory.
    """
    config = deepcopy(original)
    p = config["source_parameters"]
    factors = {
        "systemic_resistance_factor": ("heart_parameters", "R_sys"),
        "pulmonary_resistance_factor": ("heart_parameters", "R_pul"),
        "lv_elastance_factor": ("lvf_calculator", "E_es_lvf"),
        "rv_elastance_factor": ("rvf_calculator", "E_es_rvf"),
        "septal_elastance_factor": ("septum", "E_es_spt"),
        "lv_passive_stiffness_factor": ("lvf_calculator", "lambda_lvf"),
        "rv_passive_stiffness_factor": ("rvf_calculator", "lambda_rvf"),
        "pericardial_stiffness_factor": ("pericardium", "lambda_pcd"),
    }
    compliances = {"aortic": ("aorta", "ao"), "venous": ("vena_cava", "vc"),
                   "pulmonary_arterial": ("pulmonary_artery", "pa"), "pulmonary_venous": ("pulmonary_vein", "pu")}
    for name, value in knobs.items():
        if not isinstance(value, (float,int)) or not np.isfinite(value):
            raise ValueError(f"Nonfinite/non-numeric knob: {name}")
        if name in factors:
            if value <= 0: raise ValueError(f"Positive factor required: {name}")
            component, field = factors[name]
            p[component][field] *= value
        elif name.endswith("_compliance_factor") and name.removesuffix("_compliance_factor") in compliances:
            if value <= 0: raise ValueError(f"Positive compliance factor: {name}")
            component, suffix = compliances[name.removesuffix("_compliance_factor")]
            p[component]["E_es_"+suffix] /= value
        elif name in ("valve_resistance_factor", "valve_inertance_factor"):
            if value <= 0: raise ValueError(f"Positive valve factor: {name}")
            prefix = "R_" if name == "valve_resistance_factor" else "L_"
            for suffix in ("mt", "av", "tc", "pv"):
                p["heart_parameters"][prefix+suffix] *= value
        elif name == "heart_rate_bpm":
            if value <= 0: raise ValueError("Positive heart rate required")
            d = config["driver"]
            old, new = d["period_s"], 60/value
            d["period_s"], d["C_s"], d["B_per_s2"] = new, d["C_s"]*new/old, d["B_per_s2"]*(old/new)**2
        elif name == "activation_amplitude":
            if not 0 < value <= 1: raise ValueError("Periodic study amplitude must be in (0,1]")
            config["driver"]["A"] = value
        elif name == "thoracic_pressure_mmHg":
            for node in ("PA", "PU", "pericardium"):
                config["pressure_external_mmHg"][node] = value
        elif name == "initial_VC_to_PU_mL":
            config["initial_volumes_mL"]["VC"] -= value
            config["initial_volumes_mL"]["PU"] += value
        else:
            raise ValueError(f"Unknown knob: {name}")
    if min(config["initial_volumes_mL"].values()) <= 0:
        raise ValueError("Initial transfer depleted a source reservoir")
    # Preserve source_parameters as the imported provenance record where fields
    # are duplicated. driver/pressure_external/initial_volumes are authoritative.
    return config


def periodic_error(current, previous):
    """Worst same-phase normalized state change; numerical units follow M2."""
    if current.shape != previous.shape: raise ValueError("Phase grids must match")
    return float(np.max(np.abs(current-previous)/np.maximum(1.,np.abs(previous))))


def sampled_checks(model, trajectory, inventory):
    """Sampled invariants; report separately from physiological adequacy."""
    states = trajectory.states
    result = {"relative_volume_drift": float(np.max(np.abs(states[:6].sum(axis=0)-inventory))/inventory),
              "minimum_volume_mL": float(states[:6].min()),
              "minimum_valve_flow_mL_s": float(states[6:].min())}
    if result["relative_volume_drift"] >= model.settings["mass_drift_relative_limit"]:
        raise RuntimeError(f"Conservation failed: {result}")
    if result["minimum_volume_mL"] <= 0:
        raise ModelDomainError("Nonpositive sampled cavity/reservoir volume")
    if result["minimum_valve_flow_mL_s"] < -model.settings["absolute_flow_tolerance_mL_s"]:
        raise RuntimeError(f"Reverse healthy-valve flow: {result}")
    return result


def warm_up(model, protocol, solver, progress=None):
    """Require three consecutive passing phase-grid comparisons, not endpoints only.

    Cycles are restarted at an integer phase boundary using the preceding end
    state. All states and the full original inventory are retained. No clipping.
    """
    period = model.config["driver"]["period_s"]
    y = model.initial_state.copy()
    inventory = y[:6].sum()
    previous, streak, history = None, 0, []
    for cycle in range(1,model.settings["warmup_max_cycles"]+1):
        tr = integrate_source(model, period, initial_state=y, start=(cycle-1)*period,
                              sample_step=period/protocol["warmup_phase_intervals"], **solver)
        check = sampled_checks(model,tr,inventory)
        err = None if previous is None else periodic_error(tr.states,previous)
        streak = streak+1 if err is not None and err < model.settings["periodic_state_relative_tolerance"] else 0
        history.append({"cycle":cycle,"phase_error":err,"passing_streak":streak,**check})
        if progress: progress(history)
        y, previous = tr.states[:,-1].copy(), tr.states.copy()
        if streak >= protocol["consecutive_periodic_comparisons"]:
            return y, cycle*period, history
    raise RuntimeError("Periodic warm-up did not converge within configured cycle limit")


def waveform_summary(times, values):
    """Trapezoidal mean and resolved sampled extrema, with explicit time weights."""
    times, values = np.asarray(times), np.asarray(values)
    if times.ndim != 1 or values.shape != times.shape or len(times)<2 or np.any(np.diff(times)<=0):
        raise ValueError("Matching samples and strictly increasing time required")
    if not np.all(np.isfinite(times)) or not np.all(np.isfinite(values)):
        raise ValueError("Finite samples required")
    return {"mean": float(np.trapezoid(values,times)/(times[-1]-times[0])),
            "minimum":float(values.min()),"maximum":float(values.max())}


def measure(model, tr, cycles):
    """Measure each complete cycle; averages are over the requested ten cycles.

    This computes output metrics from dense-solution samples. Resolution is
    checked with independent integration/output refinements in the study suite.
    """
    pressure, residual = [], []
    for t,y in zip(tr.time,tr.states.T):
        p,m = model.pressures(t,y)
        pressure.append([p[node] for node in VOLUMES])
        residual.append(abs(m["septal_residual_mmHg"]))
    pressure = np.asarray(pressure)
    period = model.config["driver"]["period_s"]
    rows = []
    for cycle in range(cycles):
        left, right = tr.time[0]+cycle*period, tr.time[0]+(cycle+1)*period
        # Explicitly add exact cycle boundaries even when the sampling grid and
        # period are incommensurate; interpolate only these boundary samples.
        inner = tr.time[(tr.time>left)&(tr.time<right)]
        times = np.r_[left,inner,right]
        y = np.array([np.interp(times,tr.time,row) for row in tr.states])
        p = np.array([np.interp(times,tr.time,row) for row in pressure.T])
        ao,pa = waveform_summary(times,p[2]),waveform_summary(times,p[4])
        row = {"MAP_mmHg":ao["mean"],"SBP_mmHg":ao["maximum"],"DBP_mmHg":ao["minimum"],
               "PA_systolic_mmHg":pa["maximum"],"PA_diastolic_mmHg":pa["minimum"],
               "CO_L_min":waveform_summary(times,y[7])["mean"]*.06,
               "RV_output_L_min":waveform_summary(times,y[9])["mean"]*.06,
               "LV_SV_mL":float(np.ptp(y[0])),"RV_SV_mL":float(np.ptp(y[1])),
               "LV_EDV_mL":float(y[0].max()),"LV_ESV_mL":float(y[0].min()),
               "RV_EDV_mL":float(y[1].max()),"RV_ESV_mL":float(y[1].min()),
               "VC_mean_mmHg":waveform_summary(times,p[3])["mean"],
               "PU_mean_mmHg":waveform_summary(times,p[5])["mean"]}
        rows.append(row)
    return {"mean":{key:float(np.mean([r[key] for r in rows])) for key in rows[0]},
            "cycles":rows,"maximum_sampled_septal_residual_mmHg":max(residual)}, pressure
