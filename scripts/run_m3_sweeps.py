"""Run frozen-source parameter experiments. Each case saves its effective config."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
# Avoid multiplying BLAS threads by the number of independent case workers.
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
from cardiac_redundancy.source import SourceModel
from cardiac_redundancy.integration import integrate_source
from cardiac_redundancy.experiments import apply_knobs, warm_up, measure, sampled_checks


def write_json(path, data):
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n",encoding="utf-8")


def run_case(case, protocol, out):
    start_clock=time.perf_counter()
    directory=Path(out)/case["name"]
    directory.mkdir(parents=True,exist_ok=True)
    result={"case":case,"status":"running"}
    write_json(directory/"result.json",result)
    try:
        base=json.loads((ROOT/protocol.get("base_config", "configs/baseline.yaml")).read_text())
        config=apply_knobs(base,case.get("knobs",{}))
        # Source-form diagnostics specify offsets explicitly; these are not
        # user-facing physiological knobs or changes to the frozen baseline.
        for node, value in case.get("external_pressure_overrides_mmHg", {}).items():
            if node not in config["pressure_external_mmHg"] or not np.isfinite(value):
                raise ValueError(f"Invalid external pressure override: {node}")
            config["pressure_external_mmHg"][node] = value
        write_json(directory/"effective-config.json",config)
        model=SourceModel(directory/"effective-config.json")
        solver={"method":"Radau",**case.get("solver",{})}
        y,t,history=warm_up(model,protocol,solver,
            progress=lambda h:write_json(directory/"warmup.json",h))
        step=case.get("sample_step_s",protocol["sample_step_s"])
        tr=integrate_source(model,protocol["measurement_cycles"]*config["driver"]["period_s"],
            initial_state=y,start=t,sample_step=step,**solver)
        checks=sampled_checks(model,tr,model.initial_state[:6].sum())
        metrics,pressures=measure(model,tr,protocol["measurement_cycles"])
        if metrics["maximum_sampled_septal_residual_mmHg"] >= model.settings["septal_residual_tolerance_mmHg"]:
            raise RuntimeError("Sampled septal residual exceeds threshold")
        np.savez_compressed(directory/"trajectory.npz",time_s=tr.time,states=tr.states,pressures_mmHg=pressures)
        write_json(directory/"events.json",tr.events)
        result.update(status="passed_numerical_checks",warmup_cycles=len(history),
            final_phase_error=history[-1]["phase_error"],metrics=metrics,checks=checks,
            measurement_start_s=t,solver=solver,sample_step_s=step)
    except Exception as exc:
        result.update(status="failed",error_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc())
    result["elapsed_s"]=time.perf_counter()-start_clock
    write_json(directory/"result.json",result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol",default="configs/m3-sweeps.json")
    parser.add_argument("--cases",nargs="*",help="Case names; omit for all")
    parser.add_argument("--workers",type=int,default=1)
    parser.add_argument("--out",default="reports/m3-sweeps")
    parser.add_argument("--knob",action="append",default=[],metavar="NAME=VALUE",
                        help="Run a single custom case; repeat for combined knobs")
    args=parser.parse_args()
    protocol_path=ROOT/args.protocol
    protocol=json.loads(protocol_path.read_text())
    names=[c["name"] for c in protocol["cases"]]
    if len(names)!=len(set(names)): parser.error("Duplicate case names")
    if args.cases and set(args.cases)-set(names): parser.error("Unknown case name")
    cases=[c for c in protocol["cases"] if not args.cases or c["name"] in args.cases]
    if args.knob:
        if args.cases: parser.error("Use either --cases or --knob")
        try:
            knobs=dict((key,float(value)) for key,value in (item.split("=",1) for item in args.knob))
            apply_knobs(json.loads((ROOT/"configs/baseline.yaml").read_text()),knobs)
        except (ValueError,TypeError) as exc:
            parser.error(str(exc))
        cases=[{"name":"custom","knobs":knobs}]
    if args.workers<1: parser.error("workers must be positive")
    out=ROOT/args.out
    out.mkdir(parents=True,exist_ok=True)
    paths=[ROOT/protocol.get("base_config", "configs/baseline.yaml"),protocol_path,Path(__file__),*sorted((ROOT/"src/cardiac_redundancy").glob("*.py"))]
    manifest={"protocol":protocol,"selected_cases":[c["name"] for c in cases],
              "python":platform.python_version(),"workers":args.workers,
              "packages":{n:importlib.metadata.version(n) for n in ("numpy","scipy","matplotlib","pytest")},
              "hashes":{str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    # Save the exact small implementation/configuration snapshot alongside its
    # hashes. Later source edits must not erase which code generated a result.
    for path in paths:
        destination=out/"implementation-snapshot"/path.relative_to(ROOT)
        destination.parent.mkdir(parents=True,exist_ok=True)
        destination.write_bytes(path.read_bytes())
    write_json(out/"run-manifest.json",manifest)
    failures=0
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run_case,c,protocol,str(out)) for c in cases]
        for f in as_completed(futures):
            r=f.result()
            print(r["case"]["name"],r["status"],r.get("warmup_cycles"),r.get("error",""),flush=True)
            failures+=r["status"]=="failed"
    return int(failures>0)


if __name__=="__main__":
    raise SystemExit(main())
