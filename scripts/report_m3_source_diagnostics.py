"""Summarize rejected explanations and transient timing without target fitting."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"reports/m3-source-audit"


def load(path):
    return json.loads((ROOT/path/"result.json").read_text())


def main():
    b=load("reports/m3-cvsx/baseline")
    z=load("reports/m3-cvsx/thoracic_zero")
    d=load("reports/m3-cvsx/driver027_baseline")
    dz=load("reports/m3-cvsx/driver027_thoracic_zero")
    early=load("reports/m3-cvsx/first14_baseline")
    earlyz=load("reports/m3-cvsx/first14_thoracic_zero")
    only=load("reports/m3-source-form-diagnostic/pericardium_only_baseline")["metrics"]["mean"]
    change=load("reports/m3-source-form-diagnostic/pericardium_only_change")["metrics"]["mean"]
    drop=lambda base,pert:100*(1-pert["CO_L_min"]/base["CO_L_min"])
    results={"CVSX_settled_drop_pct":drop(b,z),"driver027_settled_drop_pct":drop(d,dz),
        "pericardium_only_formulation_drop_pct":drop(only,z),
        "pericardium_only_intervention_drop_pct":drop(b,change),
        "beats5_to14_mean_drop_pct":drop(early,earlyz),
        "CVSX_refinement_CO_difference_pct":{
            name:100*abs(load("reports/m3-cvsx/"+name+"_refined")["CO_L_min"]-v["CO_L_min"])/v["CO_L_min"]
            for name,v in (("baseline",b),("thoracic_zero",z))}}
    drops=100*(1-np.array(earlyz["cycle_CO_L_min"])/np.array(early["cycle_CO_L_min"]))
    results["beats3_to14_drop_pct"]=dict(zip(map(str,range(3,15)),map(float,drops)))
    (OUT/"diagnostics.json").write_text(json.dumps(results,indent=2)+"\n")
    fig,ax=plt.subplots(figsize=(8,4),layout="constrained")
    ax.plot(range(3,15),drops,"o-",label="Individual early-cycle CO reduction")
    ax.axhline(results["CVSX_settled_drop_pct"],color="C1",label="Settled reduction")
    ax.axhline(9,color="gray",linestyle="--",label="Published approximate reduction")
    ax.set(xlabel="Cycle from inherited initial conditions",ylabel="CO reduction (%)",
           title="The measurement window matters; early cycles are not periodic")
    ax.legend()
    fig.savefig(OUT/"transient-window.png",dpi=160)
    plt.close(fig)
    print(json.dumps(results,indent=2))


if __name__=="__main__": main()
