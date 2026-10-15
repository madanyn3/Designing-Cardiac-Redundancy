"""Run a source-only transient smoke case; not a periodic reproduction claim."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cardiac_redundancy.source import SourceModel, VOLUMES
from cardiac_redundancy.integration import integrate_source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cycles", type=int, default=2)
    parser.add_argument("--method", choices=["Radau","BDF","DOP853"], default="Radau")
    args = parser.parse_args()
    if args.cycles < 1: parser.error("cycles must be positive")
    config = ROOT / "configs/baseline.yaml"
    model = SourceModel(config)
    trajectory = integrate_source(model, args.cycles*model.config["driver"]["period_s"], method=args.method)
    out = ROOT / "reports/m3-source-smoke"
    out.mkdir(parents=True, exist_ok=True)
    pressures, residuals = [], []
    for t,y in zip(trajectory.time, trajectory.states.T):
        p,m = model.pressures(t,y)
        pressures.append([p[node] for node in VOLUMES])
        residuals.append(abs(m["septal_residual_mmHg"]))
    pressures = np.asarray(pressures)
    drift = float(np.max(np.abs(trajectory.states[:6].sum(axis=0)-1500))/1500)
    min_flow = float(trajectory.states[6:].min())
    if drift >= model.settings["mass_drift_relative_limit"] or min_flow < -model.settings["absolute_flow_tolerance_mL_s"]:
        raise RuntimeError(f"Smoke conservation/valve check failed: drift={drift}, minimum flow={min_flow}")
    manifest = {
        "classification": "source transient smoke run; no periodicity, reproduction or convergence claim",
        "python": platform.python_version(),
        "packages": {name:importlib.metadata.version(name) for name in ["numpy","scipy","matplotlib","pytest"]},
        "method": args.method, "cycles": args.cycles,
        "solver_settings": model.settings, "sample_step_s": .001,
        "config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
        "code_sha256": {str(path.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(path.read_bytes()).hexdigest()
                        for path in sorted((ROOT/"src/cardiac_redundancy").glob('*.py'))},
        "relative_volume_drift": drift, "minimum_sampled_valve_flow_mL_s": min_flow,
        "maximum_sampled_septal_residual_mmHg": max(residuals),
        "nfev": trajectory.nfev, "segments": trajectory.segments,
        "events": trajectory.events,
    }
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    np.savez_compressed(out/"trajectory.npz",time_s=trajectory.time,states=trajectory.states,pressures_mmHg=pressures)
    fig,axes = plt.subplots(1,3,figsize=(13,3.8))
    for node in ["AO","PA","LV","RV"]:
        axes[0].plot(trajectory.time,pressures[:,VOLUMES.index(node)],label=node)
    axes[0].set(xlabel="Time (s)",ylabel="Pressure (mmHg)")
    axes[0].legend()
    for node in ["LV","RV"]:
        i=VOLUMES.index(node)
        axes[1].plot(trajectory.time,trajectory.states[i],label=node)
        axes[2].plot(trajectory.states[i],pressures[:,i],label=node)
    axes[1].set(xlabel="Time (s)",ylabel="Cavity volume (mL)")
    axes[2].set(xlabel="Cavity volume (mL)",ylabel="Pressure (mmHg)")
    axes[1].legend(); axes[2].legend()
    fig.suptitle("Source initialization transient — not a reproduced periodic baseline")
    fig.tight_layout()
    fig.savefig(out/"source-transient.png",dpi=150)
    plt.close(fig)
    print(json.dumps({key:manifest[key] for key in ["classification","relative_volume_drift","minimum_sampled_valve_flow_mL_s","maximum_sampled_septal_residual_mmHg","nfev","segments"]},indent=2))


if __name__ == "__main__":
    main()
