"""Generate a source-sweep report and figures from saved results, including failures."""
import argparse
import json
import os
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR",str(ROOT/".matplotlib"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COMPARE_KEYS=["CO_L_min","RV_output_L_min","LV_SV_mL","RV_SV_mL",
              "MAP_mmHg","SBP_mmHg","DBP_mmHg","PA_systolic_mmHg","PA_diastolic_mmHg"]


def event_comparison(directory, reference, candidate):
    def first_events(result):
        path=Path(result["_directory"])
        events=json.loads((path/"events.json").read_text())
        config=json.loads((path/"effective-config.json").read_text())
        start=result["measurement_start_s"]
        period=config["driver"]["period_s"]
        return [(e["valve"],e["action"],e["time_s"]-start) for e in events
                if start-1e-10<=e["time_s"]<start+period-1e-10]
    ref,other=first_events(reference),first_events(candidate)
    if [(v,a) for v,a,t in ref]!=[(v,a) for v,a,t in other]:
        return {"passed":False,"reason":"Different valve/action sequences"}
    delta=max((abs(r[2]-o[2]) for r,o in zip(ref,other)),default=0)
    return {"passed":delta<.001,"maximum_time_difference_ms":delta*1000}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory",default="reports/m3-sweeps")
    parser.add_argument("--repair-directory",nargs="+",default=[],help="Overlay separately archived reruns; retain original failures")
    args=parser.parse_args()
    directory=ROOT/args.directory
    manifest=json.loads((directory/"run-manifest.json").read_text())
    results=[]
    for name in manifest["selected_cases"]:
        path=directory/name/"result.json"
        r=json.loads(path.read_text()) if path.exists() else {"case":{"name":name},"status":"missing"}
        r["_directory"]=str(path.parent)
        results.append(r)
    repair_notes=[]
    for repair_dir in args.repair_directory:
        repair_root=ROOT/repair_dir
        repair_manifest=json.loads((repair_root/"run-manifest.json").read_text())
        for name in repair_manifest["selected_cases"]:
            path=repair_root/name/"result.json"
            if not path.exists(): continue
            repaired=json.loads(path.read_text())
            repaired["_directory"]=str(path.parent)
            for i,old in enumerate(results):
                if old["case"]["name"]==name:
                    repair_notes.append({"case":name,"original_status":old["status"],"original_error":old.get("error"),"rerun_status":repaired["status"]})
                    results[i]=repaired
    by_name={r["case"]["name"]:r for r in results}
    baseline=by_name.get("baseline")
    lines=["# M3 parameter sweeps and numerical verification", "",
           "2026-10-04. Generated from saved source-only runs. Exploration factors are assumptions, not physiological normal/safe ranges. This is not organ-support or clinical validation.","",
           "Protocol: [predeclared methods](../../docs/m3-sweep-protocol.md). Explanations and commands: [model knobs and numerical methods](../../docs/m3-model-knobs.md). Each case folder contains its effective configuration, warm-up history and result; successful cases also contain arrays and event records.","",
           "## All requested cases", "",
           "Status `pass` below means periodic initialization and sampled numerical checks passed, not that a physiological target was met.","",
           "| Case | Status | Warm-up cycles | MAP | SBP/DBP | CO (L/min) | LV SV (mL) |", "|---|---|---:|---:|---:|---:|---:|"]
    for r in results:
        name=r["case"]["name"]
        result_link=os.path.relpath(Path(r["_directory"])/"result.json",directory).replace('\\','/')
        if r["status"]=="passed_numerical_checks":
            m=r["metrics"]["mean"]
            lines.append(f'| [{name}]({result_link}) | pass | {r["warmup_cycles"]} | {m["MAP_mmHg"]:.2f} | {m["SBP_mmHg"]:.2f}/{m["DBP_mmHg"]:.2f} | {m["CO_L_min"]:.4f} | {m["LV_SV_mL"]:.3f} |')
        else:
            lines.append(f'| [{name}]({result_link}) | {r["status"]} | — | — | — | — | — |')
    for r in results:
        if r["status"]!="passed_numerical_checks":
            lines += ['',f'Failure `{r["case"]["name"]}`: {r.get("error_type","")} — {r.get("error","No complete result")}']
    if repair_notes:
        lines += ['', '## Archived failed attempts and repairs', '',
                  'The original and repair directories each retain their own run manifest and code hashes. Successful original cases are retained. Reruns repair an empty-output-grid sampling error: short valve segments remain integrated but skip an empty interpolation request. No physiological parameter, equation or acceptance threshold changed. Original failures remain in their original folders. The sparse-output regression test covers the repair.', '',
                  '| Case | Original status/error | Rerun status |','|---|---|---|']
        for note in repair_notes:
            lines.append(f'| {note["case"]} | {note["original_status"]}: {note["original_error"]} | {note["rerun_status"]} |')
    evaluation={"source_screens":[],"comparisons":{},"repairs":repair_notes,"failures":[r["case"]["name"] for r in results if r["status"]!="passed_numerical_checks"]}
    if baseline and baseline["status"]=="passed_numerical_checks":
        b=baseline["metrics"]["mean"]
        lines += ["","## Frozen source reproduction screens","",
                  "The CO target is derived from 70 mL × 80/min. These coarse paper comparisons are not independent human validation.","",
                  "| Metric | Measured | Reference | Allowed absolute deviation | Result |","|---|---:|---:|---:|---|"]
        for key,target,tolerance in [("LV_SV_mL",70,7),("CO_L_min",5.6,.56),("SBP_mmHg",120,10),
                                      ("DBP_mmHg",80,10),("PA_systolic_mmHg",25,5),("PA_diastolic_mmHg",8,5)]:
            passed=abs(b[key]-target)<=tolerance
            evaluation["source_screens"].append({"metric":key,"value":b[key],"reference":target,"absolute_tolerance":tolerance,"passed":passed})
            lines.append(f'| {key} | {b[key]:.4f} | {target} | {tolerance} | {"PASS" if passed else "FAIL"} |')
        th=by_name.get("thoracic_zero")
        if th and th["status"]=="passed_numerical_checks":
            decrease=100*(1-th["metrics"]["mean"]["CO_L_min"]/b["CO_L_min"])
            evaluation["thoracic_CO_decrease_pct"]={"value":decrease,"passed":4<=decrease<=14}
            lines += ["",f'Thoracic pressure −4 → 0 mmHg: CO decreases **{decrease:.3f}%**; predeclared screen 9±5 percentage points: **{"PASS" if 4<=decrease<=14 else "FAIL"}**.']
            if not 4<=decrease<=14:
                lines += ['', 'This discrepancy prevents a full source-perturbation reproduction claim. The original Smith manuscript abstract and Results/Figure 8 were rechecked: approximately 9% is the stated result. Runtime offsets match the pinned PMR equations and the default CVSX pressure placement (PA, PU and pericardium); a dedicated unit test checks that the offset enters each affected pressure exactly once. The exact Figure 8 parameter/protocol equivalence with this later reconstruction has not been established. Do not attribute the discrepancy definitively to a particular cause or tune it away.', '', 'Follow-up numerical diagnostics: [thoracic investigation](../m3-thoracic-diagnostic/diagnosis.md). Source: [Smith author manuscript](https://ir.canterbury.ac.nz/server/api/core/bitstreams/bdd17d64-212e-46af-ac38-6eee0f28f06b/content), Results on manuscript p.18 and Figure 8.']
        lines += ["","## Numerical and initialization comparisons","",
                  "Percent change is relative to the baseline; the largest change among the nine declared flow/SV/pressure metrics must be <0.5%. Valve times are compared after aligning the first measured cycle, with a <1 ms threshold.","",
                  "| Case | Largest metric difference (%) | Largest first-cycle valve-time difference (ms) | Joint result |",
                  "|---|---:|---:|---|"]
        for name in ["baseline_refined","baseline_BDF","baseline_repeat","baseline_output_refined","initial_transfer_50","initial_transfer_minus50"]:
            r=by_name.get(name)
            if not r or r["status"]!="passed_numerical_checks": continue
            differences={key:100*abs(r["metrics"]["mean"][key]-b[key])/abs(b[key]) for key in COMPARE_KEYS}
            ev=event_comparison(directory,baseline,r)
            maximum=max(differences.values())
            passed=maximum<.5 and ev["passed"]
            evaluation["comparisons"][name]={"differences_pct":differences,"events":ev,"passed":passed}
            timing=ev.get("maximum_time_difference_ms")
            timing_text=f'{timing:.6f}' if timing is not None else 'sequence mismatch'
            lines.append(f'| {name} | {maximum:.6f} | {timing_text} | {"PASS" if passed else "FAIL"} |')
        good=[r for r in results if r["status"]=="passed_numerical_checks"]
        worst_drift=max(r["checks"]["relative_volume_drift"] for r in good)
        worst_root=max(r["metrics"]["maximum_sampled_septal_residual_mmHg"] for r in good)
        lines += ["",f'Largest measured-window relative volume drift across successful cases: **{worst_drift:.3e}**. Largest sampled septal residual: **{worst_root:.3e} mmHg**. Warm-up histories retain their separate checks.',
                  "",f'Baseline mean left/right output difference: **{abs(b["CO_L_min"]-b["RV_output_L_min"]):.3e} L/min**.',
                  "","## Reading the parameter responses","",
                  "These are model-derived responses at the declared settings, not measured human responses. Systemic resistance changes both pressure and output; it does not divide a fixed cardiac output. Contractility alters ejection while the closed circulation redistributes filling. Rate changes should be read alongside stroke volume rather than assuming unchanged volume per beat."]
        for name in ["systemic_R_half","systemic_R_double","LV_contractility_half","rate_60","rate_100"]:
            r=by_name.get(name)
            if r and r["status"]=="passed_numerical_checks":
                m=r["metrics"]["mean"]
                lines += ["",f'- **{name}:** MAP {m["MAP_mmHg"]:.2f} mmHg, CO {m["CO_L_min"]:.3f} L/min, LV SV {m["LV_SV_mL"]:.2f} mL.']
        for name in ["aortic_compliance_half","aortic_compliance_double"]:
            r=by_name.get(name)
            if r and r["status"]=="passed_numerical_checks":
                m=r["metrics"]["mean"]
                lines += ["",f'- **{name}:** pulse pressure {m["SBP_mmHg"]-m["DBP_mmHg"]:.2f} mmHg, versus baseline {b["SBP_mmHg"]-b["DBP_mmHg"]:.2f}; CO {m["CO_L_min"]:.3f} L/min. Pressure pulsatility changes much more than mean output in these two selected cases.']
        lines += ["","## Parameter response figures","","![One-at-a-time source response](parameter-response.png)",
                  "","![Baseline final measured cycle](baseline-cycle.png)"]
        physical=[r for r in good if r["case"]["name"]=="baseline" or r["case"].get("knobs") and not r["case"]["name"].startswith("initial_")]
        fig,axes=plt.subplots(1,3,figsize=(13,7),sharey=True)
        labels=[r["case"]["name"].replace('_',' ') for r in physical]
        y=np.arange(len(physical))
        for ax,key,label in zip(axes,["CO_L_min","MAP_mmHg","LV_SV_mL"],["Cardiac output (L/min)","MAP (mmHg)","LV stroke volume (mL)"]):
            ax.barh(y,[r["metrics"]["mean"][key] for r in physical],color=['#1d4ed8' if r["case"]["name"]=='baseline' else '#63a7c4' for r in physical])
            ax.axvline(b[key],color='#1d4ed8',linestyle='--',linewidth=1)
            ax.set_xlabel(label); ax.grid(axis='x',alpha=.2)
        axes[0].set_yticks(y,labels); axes[0].invert_yaxis()
        fig.suptitle("Source parameter experiments — independent periodic warm-up; no organ-support claim")
        fig.tight_layout();fig.savefig(directory/"parameter-response.png",dpi=160);plt.close(fig)
        with np.load(directory/"baseline/trajectory.npz") as data:
            t,y,p=data['time_s'],data['states'],data['pressures_mmHg']
            mask=t>=t[-1]-.75-1e-10
            phase=t[mask]-(t[-1]-.75)
            fig,axes=plt.subplots(1,3,figsize=(13,4))
            for node,index in [("AO",2),("PA",4),("LV",0),("RV",1)]: axes[0].plot(phase,p[mask,index],label=node)
            for name,index in [("MT",6),("AV",7),("TC",8),("PV",9)]: axes[1].plot(phase,y[index,mask],label=name)
            for name,index in [("LV",0),("RV",1)]: axes[2].plot(y[index,mask],p[mask,index],label=name)
            axes[0].set(xlabel="Phase time (s)",ylabel="Pressure (mmHg)")
            axes[1].set(xlabel="Phase time (s)",ylabel="Valve flow (mL/s)")
            axes[2].set(xlabel="Cavity volume (mL)",ylabel="Pressure (mmHg)")
            for ax in axes: ax.legend();ax.grid(alpha=.2)
            fig.suptitle("Unchanged source: final cycle after periodic initialization")
            fig.tight_layout();fig.savefig(directory/"baseline-cycle.png",dpi=160);plt.close(fig)
    lines += ["","## What is still unestablished","",
              "The subsequent [independent pinned CVSX comparison](../m3-source-audit/report.md) passes for baseline and the thoracic-pressure perturbation (maximum difference across 15 metrics: 0.003206%). The original journal experiment remains unresolved; see the [source investigation](../m3-source-audit/investigation.md). Project organ branches, auxiliary hydraulics, collapse-domain events and oxygen/respiratory coupling are not implemented by these sweeps. Simultaneous-event and longer failure-event stress tests remain open. Selected one-at-a-time cases do not validate every parameter or combination, and source reproduction is not proof of human physiological accuracy.","",
              "To regenerate this report: `python scripts/report_m3_sweeps.py" + (" --repair-directory " + " ".join(args.repair_directory) if args.repair_directory else "") + "` using the project environment. Run commands and custom knobs are in the linked guide. Failed cases remain part of this report; do not remove them to improve the apparent pass rate.",""]
    (directory/"evaluation.json").write_text(json.dumps(evaluation,indent=2)+"\n",encoding='utf-8')
    (directory/"report.md").write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps({"reported":len(results),"failed_or_incomplete":evaluation['failures']},indent=2))


if __name__=="__main__": main()
