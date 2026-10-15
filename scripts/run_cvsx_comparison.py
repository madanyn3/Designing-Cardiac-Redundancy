"""Execute unmodified pinned CVSX equations, with independent JAX/Diffrax numerics.

Run with .venv-cvsx/Scripts/python.exe. No project RHS is imported. CVSX uses
litres, mmHg and seconds. Its clipped fixed-state valves differ numerically from
our event-reset valves; convergence comparisons must quantify that distinction.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
UPSTREAM = ROOT / "references/m3/cvsx-upstream"
sys.path.insert(0, str(UPSTREAM))
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
import numpy as np
import diffrax
from cvsx import models, parameters, components, cardiac_drivers


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pth", type=float, default=-4)
    parser.add_argument("--pericardium-only", action="store_true")
    parser.add_argument("--rtol", type=float, default=1e-8)
    parser.add_argument("--driver-center", type=float, default=.375)
    parser.add_argument("--warmup-cycles", type=int, default=150)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)
    config = json.loads((ROOT / "configs/baseline.yaml").read_text())
    params = parameters.build_parameter_tree("revie", inertial=True)
    # Only cross-implementation parameter mismatch: upstream rounds L_mt to
    # 7.6967e-5; match the pinned CellML 7.6968e-5, explicitly (mL -> L).
    params["mt"] = components.Valve(
        r=params["mt"].r, l=config["source_parameters"]["heart_parameters"]["L_mt"]*1000,
        inertial=True)
    params["p_pl"] = args.pth
    params["cd"] = cardiac_drivers.GaussianCardiacDriver(
        a=jnp.array([1.]), b=jnp.array([80.]), c=jnp.array([args.driver_center]), hr=jnp.array(80.))
    model = models.SmithCVS(**params, p_pl_affects_pu_and_pa=not args.pericardium_only)
    y0 = {"v_"+k.lower(): jnp.array(v/1000) for k,v in config["initial_volumes_mL"].items()}
    y0.update({"q_"+k.lower(): jnp.array(v/1000) for k,v in config["initial_valve_flows_mL_s"].items()})
    if args.warmup_cycles < 2:
        parser.error("At least two warmup cycles required")
    warmup = args.warmup_cycles*.75
    end = warmup+7.5
    manifest = {"python": platform.python_version(), "arguments": vars(args),
        "packages": {n: importlib.metadata.version(n) for n in
                     ("jax", "jaxlib", "equinox", "diffrax", "numpy", "scipy")},
        "source_hashes": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in sorted((UPSTREAM/"cvsx").glob("*.py"))},
        "adapter_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "parameter_override": {"L_mt_mmHg_s2_mL": 7.6968e-5},
        "warmup_s": warmup, "measurement_s": 7.5, "sample_step_s": 0.0005}
    (out/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    (out/"adapter.py").write_bytes(Path(__file__).read_bytes())
    # Save preceding cycles too, to audit same-phase periodicity independently.
    ts = jnp.linspace(warmup-1.5, end, 18001)
    root_solver = diffrax.NewtonNonlinearSolver(rtol=1e-11, atol=1e-12, max_steps=100)

    @jax.jit
    def solve():
        return diffrax.diffeqsolve(diffrax.ODETerm(model), diffrax.Tsit5(),
            t0=0., t1=end, dt0=1e-5, y0=y0, args=(root_solver,),
            stepsize_controller=diffrax.PIDController(rtol=args.rtol, atol=args.rtol*0.001, dtmax=0.001),
            max_steps=2000000, saveat=diffrax.SaveAt(ts=ts))

    start = time.perf_counter()
    solution = solve()
    jax.block_until_ready(solution.ys)
    _, outputs = jax.jit(lambda t,y: model(t,y,(root_solver,),return_outputs=True))(ts,solution.ys)
    times = np.asarray(ts)
    states = {k: np.asarray(v) for k,v in solution.ys.items()}
    outputs = {k: np.asarray(v) for k,v in outputs.items()}
    mask = times >= warmup
    co = float(np.trapz(outputs["q_av"][mask], times[mask])/7.5*60)
    # Normalize in project units (mL and mL/s); avoid hiding near-zero flow error.
    phase_error = max(float(np.max(np.abs(v[-1501:]-v[-3001:-1500])*1000 /
                    np.maximum(1, np.abs(v[-3001:-1500])*1000))) for v in states.values())
    inventory = sum(v for k,v in states.items() if k.startswith("v_"))
    metrics = {"CO_L_min": co, "MAP_mmHg": float(np.trapz(outputs["p_ao"][mask],times[mask])/7.5),
        "SBP_mmHg": float(outputs["p_ao"][mask].max()), "DBP_mmHg": float(outputs["p_ao"][mask].min()),
        "phase_error_last_cycle": phase_error,
        "relative_volume_drift": float(np.max(np.abs(inventory-1.5))/1.5),
        "elapsed_s": time.perf_counter()-start,
        "steps": {k:int(v) for k,v in solution.stats.items()}}
    metrics["cycle_CO_L_min"] = [float(np.trapz(outputs["q_av"][i:i+1501],times[i:i+1501])/.75*60)
                                for i in range(0,18000,1500)]
    np.savez_compressed(out/"trajectory.npz", time_s=times,
                        **{"state_"+k:v for k,v in states.items()}, **outputs)
    (out/"result.json").write_text(json.dumps(metrics,indent=2)+"\n")
    print(json.dumps(metrics), flush=True)


if __name__ == "__main__":
    main()
