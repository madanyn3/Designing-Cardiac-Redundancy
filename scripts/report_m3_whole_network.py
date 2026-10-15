"""Numerical verification and achieved regional flows for whole-network controls."""
import json,os
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'reports/m3-whole-network'
os.environ.setdefault('MPLCONFIGDIR',str(OUT/'mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def read(p):return json.loads(p.read_text())
def load(n,r):return [np.load(OUT/n/f"segment-{s['index']}.npz") for s in r['stages']]

def main():
    protocol=read(OUT/'manifest.json')['protocol'];results={c['name']:read(OUT/c['name']/'result.json') for c in protocol['cases']}
    comparisons={};passed=True;summaries={}
    for n,r in results.items():
        if r['status']=='failed':comparisons[n]={'passed':False,'error':r['error']};passed=False;continue
        if n.endswith('-nominal'):continue
        refname=n.rsplit('-',1)[0]+'-nominal';ref=results[refname]
        if ref['status']=='failed':passed=False;continue
        aa=load(refname,ref);bb=load(n,r);err=0.
        for a,b in zip(aa,bb):
            ta=a['time']-ref['arrest_s'];tb=b['time']-r['arrest_s'];mask=ta<=min(ta[-1],tb[-1])
            for key in ('states','pressures','flows'):
                xa=a[key] if key=='states' else a[key].T;xb=b[key] if key=='states' else b[key].T
                interp=np.array([np.interp(ta[mask],tb,row) for row in xb]);scale=np.maximum(1,np.max(abs(xa[:,mask]),axis=1))
                err=max(err,float(np.max(abs(interp-xa[:,mask])/scale[:,None])))
        ea=[e for s in ref['stages'] for e in s['valves']];eb=[e for s in r['stages'] for e in s['valves']]
        same=[(e['valve'],e['action']) for e in ea]==[(e['valve'],e['action']) for e in eb]
        dt=max([abs((a['time_s']-ref['arrest_s'])-(b['time_s']-r['arrest_s']))*1000 for a,b in zip(ea,eb)],default=0) if same else None
        term=(ref['termination'] or {}).get('node')==(r['termination'] or {}).get('node')
        good=bool(err<.005 and same and dt<1 and term and abs(ref['observed_s']-r['observed_s'])<.001 and len(aa)==len(bb))
        comparisons[n]=dict(max_error_pct=100*err,valve_ms=dt,sequence_matches=same,termination_matches=term,passed=good);passed &= good
    fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
    for n,r in results.items():
        if not n.endswith('-nominal') or r['status']=='failed':continue
        ds=load(n,r);t=np.hstack([d['time']-r['arrest_s'] for d in ds]);q=np.vstack([d['flows'] for d in ds]);p=np.vstack([d['pressures'] for d in ds]);fn=r['flow_names']
        mask=(t>=2)&(t<=3);summary={'status':r['status'],'observed_s':r['observed_s']}
        if np.count_nonzero(mask)>1:
            def avg(x):return float(np.trapezoid(x[mask],t[mask])/(t[mask][-1]-t[mask][0]))
            summary['flows_mL_min']={name:avg(q[:,i])*60 for i,name in enumerate(fn)};summary['AO_mmHg']=avg(p[:,2])
        summaries[n]=summary
        label=n.removesuffix('-nominal');axes[0].plot(t,p[:,2],label=label)
        axes[1].plot(t,(q[:,fn.index('kidney_protected')]+q[:,fn.index('kidney_other')])*60,label=label)
    for a in axes:a.axvline(3,color='grey',ls=':');a.axvline(5,color='grey',ls=':')
    axes[0].set_ylabel('Aortic pressure (mmHg)');axes[0].legend(fontsize=8)
    axes[1].set_ylabel('Kidney-pair flow (mL/min)');axes[1].set_xlabel('Time after arrest / regional control (s)')
    fig.suptitle('Whole-network controls; restore conductance from 3 to 5 seconds');fig.tight_layout();fig.savefig(OUT/'controls.png',dpi=150);plt.close(fig)
    rows=['# Whole-network regional controls — M3','',
      '2026-10-06. Implements independent conductance controls beyond limbs, with a validated partition of tissue-supply edges into active support and preservation roles. Roles do not prescribe adequate flow or confer recoverability. No oxygen/fuel, anticoagulation, tissue injury or native-repair model is added.','',
      '## Objective mapping','',
      'Active supply edges: brain, lungs (nutritive supply), respiratory_muscle (surrogate), liver_arterial. Pulmonary exchange and hepatic outflow are explicit dependent paths. Native coronary, both kidney, all five digestive/portal, all limb and remaining-tissue supplies belong to preservation. Every supply edge must be assigned exactly once. All native tissues and 5,250 mL blood remain accounted for.','',
      '## Implementation and experiment','',
      'WholeNetworkModel extends the verified selective storage/B1 fixture. Independent factors multiply branch conductance; cerebral/coronary pressure laws and finite blood stores remain. The hepatic junction is solved with the effective controlled conductances before flows are calculated. All-open RHS is exactly equal to the previous model, permitting hash-checked reuse of its periodic seed. Independent numerical comparisons share that same verified seed.','',
      'Five scenarios each use nominal Radau, refined Radau and BDF: all open; preservation supply conductances at .25; preservation supply conductances at zero; digestive/portal supply at zero with hepatic artery retained; and only the other kidney closed. Quarter and zero are numerical stress fixtures, not preservation minima. In particular the zero case stops native coronary supply and cannot claim recoverability. Zero-flow comparator inclusion is not a recommendation to run the preservation set without blood.','',
      'Arrest/restriction occurs at phase .375 s. Frozen B1 heads ramp to 20/100 mmHg in one second, Rp=.1 and capacity=120 mL/s per stage. Returns remain open; limb collaterals remain at the previous zero comparator setting. All regional conductances restore linearly during seconds 3–5, and the observation ends at 8 s or a domain stop. Native activation remains off: this is gate restoration, not repaired-heart handover.','',
      '## Achieved flow during seconds 2–3','',
      'Transient averages before reopening. Conductance factors are not imposed flow fractions. Liver arterial and portal contributions are shown separately; their combined serial outflow is not extra cardiac output.','',
      '| Case | AO mmHg | Brain mL/min | Kidney pair mL/min | Native coronary mL/min | Hepatic artery mL/min | Portal contributors mL/min | Left support L/min |','|---|---:|---:|---:|---:|---:|---:|---:|']
    gut=['stomach_oesophagus','small_intestine','large_intestine','pancreas','spleen']
    for n,s in summaries.items():
        if 'flows_mL_min' not in s:rows.append(f"| {n} | no valid 2–3 s window | | | | | | |");continue
        q=s['flows_mL_min'];rows.append(f"| {n} | {s['AO_mmHg']:.3f} | {q['brain']:.2f} | {q['kidney_protected']+q['kidney_other']:.2f} | {q['coronary_left']+q['coronary_right']:.2f} | {q['liver_arterial']:.2f} | {sum(q[k] for k in gut):.2f} | {q['support_left']/1000:.4f} |")
    rows += ['', '![Regional control trajectories](controls.png)','', '## Verification','',
      '| Comparison | Max scaled difference % | Valve timing ms | Pass |','|---|---:|---:|---|']
    for n,c in comparisons.items():rows.append(f"| {n} | {c.get('max_error_pct')} | {c.get('valve_ms')} | {c['passed']} |")
    rows += ['',f'Overall declared comparisons pass: **{passed}**.',
      f"Maximum relative blood drift: {max(r.get('relative_blood_drift',0) for r in results.values()):.3e}.",
      f"Maximum hepatic zero-storage residual: {max(r.get('maximum_hepatic_residual',0) for r in results.values()):.3e} mL/s.",
      f"Maximum state discontinuity at segment restart: {max(r.get('state_jump',0) for r in results.values()):.3e}.",
      '', 'Limits: maximum coordinate-scaled trajectory difference <0.5%, using max(1 native unit, reference coordinate peak); matched valve sequences and times <1 ms; matched domain node/end time. Unit tests also cover exact all-open regression, independent kidney/coronary controls, forward/reverse hepatic balance, ramps, disconnected-H rejection and exhaustive role validation. [75-test record](../m3-whole-network-tests.xml).','',
      '## Reproduce and remaining scope','',
      'Run `.venv/Scripts/python.exe scripts/run_m3_whole_network.py`, then `.venv/Scripts/python.exe scripts/report_m3_whole_network.py`. [Protocol and role map](../../configs/m3-whole-network-protocol.json), [model](../../src/cardiac_redundancy/whole_network.py). Each run retains states, pressures, all branch flows, event records and code/configuration/seed snapshots.','',
      'This removes the limb-only software limitation. Anatomical gates and collateral topology, additional non-limb regional blood storage, preservation-flow/duration thresholds, liver glycogen/glucose and oxygen/respiratory chemistry remain unresolved. Pump resistance/capacity are still fixture constants; energy/size maps and native-heart handover are next hydraulic/interface work. No repair-window or full-recovery result is established.','']
    (OUT/'report.md').write_text('\n'.join(rows),encoding='utf-8')
    (OUT/'verification.json').write_text(json.dumps(dict(passed=passed,comparisons=comparisons,scenarios=summaries),indent=2)+'\n')
    print(json.dumps(dict(passed=passed,comparisons=comparisons),indent=2));return 0 if passed else 1

if __name__=='__main__':raise SystemExit(main())
