"""Routing verification and conservative passive-marker postprocessing."""
import hashlib,json,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'reports/m3-selective'
os.environ.setdefault('MPLCONFIGDIR',str(OUT/'matplotlib-cache'))
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cardiac_redundancy.transport import marker_rhs

def read(p):return json.loads(p.read_text())
def load(name,r):
    return [np.load(OUT/name/f"segment-{s['index']}.npz") for s in r['segments']]

def marker(name,r,segments):
    effective=read(OUT/name/'effective.json');names=r['volume_names'];edges=effective['edges']
    amount=np.zeros(16);amount[6:]=segments[0]['states'][6:16,0]
    initial=sum(amount);history=[];drift=0.;minimum=1.;maximum=0.
    for d in segments:
        ts=d['time_s'];v=d['states'][:16];q=d['flows'].T
        def rhs(t,a):
            vi=np.array([np.interp(t,ts,row) for row in v]);qi=np.array([np.interp(t,ts,row) for row in q])
            return marker_rhs(vi,a,edges,qi,names)
        if ts[-1]==ts[0]:continue
        sol=solve_ivp(rhs,(ts[0],ts[-1]),amount,rtol=1e-8,atol=1e-9,max_step=.02,t_eval=ts)
        if not sol.success:raise RuntimeError('Marker solve failed')
        amount=sol.y[:,-1];concentration=sol.y/v
        minimum=min(minimum,float(concentration.min()));maximum=max(maximum,float(concentration.max()))
        drift=max(drift,float(np.max(abs(sol.y.sum(axis=0)-initial)))/initial)
        history.append(np.vstack((ts-r['arrest_s'],sol.y[6:].sum(axis=0))))
    result=dict(initial_marker_mL=initial,final_marker_in_limb_mL=float(sum(amount[6:])),
                relative_marker_drift=drift,minimum_concentration=minimum,maximum_concentration=maximum,
                passed=bool(drift<1e-8 and minimum>=-1e-6 and maximum<=1.001))
    (OUT/name/'marker.json').write_text(json.dumps(result,indent=2)+'\n')
    return result,np.hstack(history)

def main():
    manifest=read(OUT/'manifest.json');cases=manifest['protocol']['cases']
    results={c['name']:read(OUT/c['name']/'result.json') for c in cases};comparisons={};ok=True
    for name,r in results.items():
        if r['status']=='failed':comparisons[name]={'passed':False,'error':r['error']};ok=False;continue
        if not name.endswith(('-refined','-BDF')):continue
        base_name=name.rsplit('-',1)[0]+'-nominal';ref=results[base_name]
        if ref['status']=='failed':ok=False;continue
        aa=load(base_name,ref);bb=load(name,r);error=0.
        for a,b in zip(aa,bb):
            ta=a['time_s']-ref['arrest_s'];tb=b['time_s']-r['arrest_s'];mask=ta<=min(ta[-1],tb[-1])
            for key in ('states','pressures','flows'):
                xa=a[key] if key=='states' else a[key].T;xb=b[key] if key=='states' else b[key].T
                interp=np.array([np.interp(ta[mask],tb,row) for row in xb]);scale=np.maximum(1,np.max(abs(xa[:,mask]),axis=1))
                error=max(error,float(np.max(abs(xa[:,mask]-interp)/scale[:,None])))
        ea=[e for s in ref['segments'] for e in s['valves']];eb=[e for s in r['segments'] for e in s['valves']]
        same=[(e['valve'],e['action']) for e in ea]==[(e['valve'],e['action']) for e in eb]
        dt=max([abs((a['time_s']-ref['arrest_s'])-(b['time_s']-r['arrest_s']))*1000 for a,b in zip(ea,eb)],default=0) if same else None
        enddiff=abs((aa[-1]['time_s'][-1]-ref['arrest_s'])-(bb[-1]['time_s'][-1]-r['arrest_s']))*1000
        term=(ref['termination'] or {}).get('node')==(r['termination'] or {}).get('node')
        passed=bool(error<.005 and same and dt<1 and enddiff<1 and term and len(aa)==len(bb));ok &= passed
        comparisons[name]=dict(max_scaled_error_pct=100*error,valve_sequence=same,valve_ms=dt,endpoint_ms=enddiff,termination_matches=term,passed=passed)
    nominal={n:r for n,r in results.items() if not n.endswith(('-refined','-BDF')) and r['status']!='failed'}
    markers={};summaries={};fig,axes=plt.subplots(3,1,figsize=(10,9),sharex=True)
    for name,r in nominal.items():
        ds=load(name,r);a=np.hstack([d['time_s']-r['arrest_s'] for d in ds]);v=np.hstack([d['states'][:16] for d in ds])
        pressures=np.vstack([d['pressures'] for d in ds]);flows=np.vstack([d['flows'] for d in ds])
        total=v[6:].sum(axis=0);fn=r['flow_names']
        limbq=sum(flows[:,i] for i,n in enumerate(fn) if n.startswith('limb_') and n.endswith('_exchange'))*60
        marker_result,history=marker(name,r,ds);markers[name]=marker_result;ok &= marker_result['passed']
        # Event logs carry identical before/after states; trajectory endpoints
        # additionally verify continuity between actual solver segments.
        jumps=[float(np.max(abs(x['states'][:,-1]-y['states'][:,0]))) for x,y in zip(ds,ds[1:])]
        events=read(OUT/name/'events.json');jumps += [float(np.max(abs(np.array(e['state_before'])-e['state_after']))) for e in events]
        jump=max(jumps,default=0);ok &= jump<1e-9
        before=a<=4+1e-9
        summary=dict(initial_limb_mL=float(total[0]),limb_mL_at_reopen=float(total[np.where(before)[0][-1]]),
                     final_limb_mL=float(total[-1]),minimum_limb_mL=float(total.min()),maximum_limb_mL=float(total.max()),
                     final_AO_mmHg=float(pressures[-1,2]),final_limb_exchange_mL_min=float(limbq[-1]),
                     minimum_AO_after_reopening_mmHg=float(pressures[a>=4,2].min()) if np.any(a>=4) else None,
                     segment_state_jump=jump,observed_s=float(a[-1]),status=r['status'])
        summaries[name]=summary
        if name.endswith('-nominal'):
            label=name.removesuffix('-nominal');axes[0].plot(a,total,label=label);axes[1].plot(a,pressures[:,2],label=label);axes[2].plot(history[0],history[1]/marker_result['initial_marker_mL'],label=label)
    for ax in axes:ax.axvline(4,color='grey',ls=':',lw=1);ax.axvline(6,color='grey',ls=':',lw=1)
    axes[0].set_ylabel('Regional blood (mL)');axes[0].legend(ncol=3,fontsize=8)
    axes[1].set_ylabel('Aortic pressure (mmHg)');axes[2].set_ylabel('Original regional marker\nremaining in region (fraction)');axes[2].set_xlabel('Seconds after arrest and restriction')
    fig.suptitle('Exploratory storage and B1 pump fixture; reopen at 4 s (ramp ends at 6 s)');fig.tight_layout();fig.savefig(OUT/'routing.png',dpi=150);plt.close(fig)
    lines=['# Selective routing, retained blood and reopening — M3 verification','',
      '2026-10-06. Implements the approved exploratory routing proposal. These are numerical fixtures, not measured limb anatomy, preservation-flow prescriptions or evidence of full recovery. No anticoagulant, oxygen uptake, cell injury, coagulation, edema or no-reflow kinetics are simulated.','',
      '## Configuration and boundaries','',
      'Each original limb pool is split into proximal arterial/microvascular and distal venous storage, conserving existing physical and unstressed blood. Numerical proximal fractions 0.1/0.2/0.3 are independently warmed to periodic states. A 0.1/0.9 gate/exchange resistance split preserves the original series resistance at the reference construction pressure drop; compliance is constructed to match reference proximal/distal pressures. These are explicit engineering sensitivity choices, not literature-derived regional anatomy. The five tissue pools still aggregate all limbs.','',
      'All gates act together in the main cases. Default restriction means gate conductance 0.25, not 25% achieved tissue flow. Because the gate is only 10% of original inlet resistance, this setting makes the combined gate/exchange resistance 1.3 times its original value at fixed pressures. The residual collateral conductance 0.02 is a test choice relative to open gate conductance, not a measured leakage rate. No collateral topology or independent tissue-selective hardware is validated.','',
      'Native activation is removed at phase 0.375 s. B1 right/left heads ramp to 20/100 mmHg over one second, with Rp=.1 mmHg s/mL and capacity 120 mL/s per stage. These are the frozen M2 verification seeds, not optimized pump settings. Inlet/return controls switch at arrest; reopening starts at 4 s and is either immediate or a two-second conductance ramp. The observation ends at 10 s or a localized domain boundary. All intervals are engineering test durations, not ischemic tolerance estimates.','',
      '## Outcomes','',
      '| Case | Initial regional mL | At reopening mL | Final regional mL | Final AO mmHg | Observed s | Status |',
      '|---|---:|---:|---:|---:|---:|---|']
    for n,s in summaries.items():lines.append(f"| {n} | {s['initial_limb_mL']:.3f} | {s['limb_mL_at_reopen']:.3f} | {s['final_limb_mL']:.3f} | {s['final_AO_mmHg']:.3f} | {s['observed_s']:.4f} | {s['status']} |")
    lines += ['', 'If a case terminates before 4 s, its “at reopening” column is the last valid observation; no reopening occurred. Domain stops are not successful preservation or proof of biological infeasibility.','',
      'Arterial exclusion with open return drains only part of the regional inventory. Closing both ends preserves combined regional volume during isolation; internal arterial-to-venous redistribution still occurs. Venous obstruction with continued inlet flow accumulates blood. Reopening the trapped case produces transient blood-volume overshoot; this is intravascular filling, not simulated edema. The original-blood marker keeps changing after volume begins to recover, demonstrating that refilling and washout are different quantities. These observations apply to the selected fixtures.','',
      '![Routing and retained blood](routing.png)','',
      'The passive marker labels blood initially in regional stores. Its global amount is conserved while it moves and mixes, including reversed flows and zero-storage hepatic mixing. Remaining regional marker is not identical to regional blood volume: replacement blood can refill a region. The marker is not viable-cell count, oxygen debt or anticoagulant concentration. It is reconstructed from saved hydraulic samples using an independent transport solve; local mixing is idealized.','',
      '## Verification','',
      '| Comparison | Scaled max error % | Valve timing ms | Pass |','|---|---:|---:|---|']
    for n,c in comparisons.items():lines.append(f"| {n} | {c.get('max_scaled_error_pct')} | {c.get('valve_ms')} | {c['passed']} |")
    lines += ['',f'Overall declared checks passed: **{bool(ok)}**.',
      'Trajectory criterion: per-coordinate maximum absolute difference divided by max(1 in its native unit, nominal peak magnitude), below 0.5%; matching valve sequences and domain endpoints within 1 ms. Storage sensitivity and venous-obstruction cases use nominal settings only; they are not independently refined/BDF-verified.',
      f"Maximum relative hydraulic blood drift: {max((r.get('relative_volume_drift',0) for r in results.values()),default=0):.3e}.",
      f"Maximum passive-marker amount drift: {max((r['relative_marker_drift'] for r in markers.values()),default=0):.3e}.",
      '', 'Unit tests cover storage/unstressed-inventory partition, reference resistance equivalence, trapping, leakage, pump shutoff/saturation, conductance ramps, analytic marker mixing, signed donor reversal and hepatic mixing. [Test results](../m3-selective-tests.xml).','',
      '## Reproduce and next work','',
      'Run `.venv/Scripts/python.exe scripts/run_m3_selective.py`, then `.venv/Scripts/python.exe scripts/report_m3_selective.py`. [Protocol](../../configs/m3-selective-protocol.json), [model](../../src/cardiac_redundancy/selective.py), [transport](../../src/cardiac_redundancy/transport.py), [research basis and limitations](../../docs/m3-selective-perfusion-reperfusion.md). The run snapshot preserves the hydraulic implementation and inputs; report metadata hashes the later passive-marker code.','',
      'Next obtain anatomical storage/collateral inputs and replace or sweep the fixture choices; extend pump delay/mismatch/domain tests and organ controls. M4 adds oxygen and chemistry, with injury/coagulation limitations explicit. No full-reperfusion or functional-recovery claim follows from these hydraulic/marker checks. M3 remains open.','']
    (OUT/'report.md').write_text('\n'.join(lines),encoding='utf-8')
    output=dict(passed=bool(ok),comparisons=comparisons,outcomes=summaries,markers=markers,
                report_code_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),ROOT/'src/cardiac_redundancy/transport.py')})
    for p in (Path(__file__),ROOT/'src/cardiac_redundancy/transport.py'):
        target=OUT/'report-implementation-snapshot'/p.relative_to(ROOT)
        target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(p.read_bytes())
    (OUT/'verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2));return 0 if ok else 1

if __name__=='__main__':raise SystemExit(main())
