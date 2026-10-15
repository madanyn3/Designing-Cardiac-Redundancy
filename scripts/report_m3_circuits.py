"""Compare declared circuit trajectories, events and localized domain exits."""
import json,os,sys,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'reports/m3-circuits'
sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.pump_energy import stage_power, integrate_power, MMHG_ML_TO_J
os.environ.setdefault('MPLCONFIGDIR',str(OUT/'mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read(p):return json.loads(p.read_text())
def load(name,r):
    records=[]
    for s in r['segments']:
        with np.load(OUT/name/f"segment-{s['index']}.npz") as z:records.append({k:z[k].copy() for k in z.files})
    return records

def energy_records(r,ds):
    """Carry the mechanical ledger forward using each saved segment's controls."""
    case=r['case'];architecture=case['architecture'];totals={};flow_error=0.
    if architecture=='BIO':
        # Boundary pressure-volume work includes active and passive mechanics;
        # it is not auxiliary metabolic energy or the ideal head-port ledger.
        for chamber in ('aLV','aRV'):
            i=r['volume_names'].index(chamber)
            totals[chamber+'_boundary_work_J']=sum(float(-MMHG_ML_TO_J*np.trapezoid(d['pressures'][:,i],d['states'][i])) for d in ds)
        return totals,flow_error
    edges=([('support_right','VC','PA',20.,case.get('right_capacity',120)),('support_left','PU','AO',100.,case.get('left_capacity',120))]
           if architecture=='B1' else [('oxygenator','VC','AO',100.,120.)]+([('support_left','PU','AO',100.,case.get('left_capacity',120))] if architecture=='E2' else []))
    for name,a,b,head,capacity in edges:
        totals[name]={}
        for d in ds:
            relative=d['time']-r['arrest_s'];records=[]
            # At the end of the preceding segment the port is still released;
            # its isolation belongs to the next segment, preserving both sides.
            enabled=bool(relative[0]>=case.get('delay',0)-1e-10 and not(case.get('restart') and relative[0]>=7-1e-10))
            for j,t in enumerate(relative):
                factor=np.clip(t-case.get('delay',0),0.,1.)
                if case.get('restart'):factor*=1-np.clip(t-6,0.,1.)
                p=d['pressures'][j]
                record=stage_power(p[r['volume_names'].index(a)],p[r['volume_names'].index(b)],head*factor,.1,capacity,enabled)
                flow_error=max(flow_error,abs(record['flow_mL_s']-d['flows'][j,r['flow_names'].index(name)]))
                records.append(record)
            if len(relative)>1:
                for k,v in integrate_power(relative,records).items():totals[name][k]=totals[name].get(k,0.)+float(v[-1])
    return totals,flow_error

def main():
    protocol=read(OUT/'manifest.json')['protocol']
    names=[c['name']+'-'+s['label'] for c in protocol['cases'] for s in protocol['solvers']]
    results={n:read(OUT/n/'result.json') for n in names}
    comparisons={};passed=True
    for name,r in results.items():
        if r['status'] not in ('completed','model_domain_limit'):
            comparisons[name]=dict(passed=False,error=r.get('error',r['status']));passed=False;continue
        if name.endswith('-nominal'):continue
        refname=name.rsplit('-',1)[0]+'-nominal';ref=results[refname]
        if ref['status'] not in ('completed','model_domain_limit'):passed=False;continue
        aa=load(refname,ref);bb=load(name,r);error=0.
        for a,b in zip(aa,bb):
            ta=a['time']-ref['arrest_s'];tb=b['time']-r['arrest_s']
            mask=(ta>=tb[0]-1e-10)&(ta<=min(ta[-1],tb[-1]))
            for key in ('states','pressures','flows'):
                xa=a[key] if key=='states' else a[key].T;xb=b[key] if key=='states' else b[key].T
                interpolated=np.array([np.interp(ta[mask],tb,row) for row in xb])
                scale=np.maximum(1,np.max(abs(xa[:,mask]),axis=1))
                error=max(error,float(np.max(abs(interpolated-xa[:,mask])/scale[:,None])))
        ea=[e for s in ref['segments'] for e in s['events']];eb=[e for s in r['segments'] for e in s['events']]
        sequence=[(e['valve'],e['action']) for e in ea]==[(e['valve'],e['action']) for e in eb]
        dt=max([abs((a['time_s']-ref['arrest_s'])-(b['time_s']-r['arrest_s']))*1000 for a,b in zip(ea,eb)],default=0.) if sequence else None
        term=(ref['termination'] or {}).get('node')==(r['termination'] or {}).get('node')
        end_delta=abs((aa[-1]['time'][-1]-ref['arrest_s'])-(bb[-1]['time'][-1]-r['arrest_s']))*1000
        good=error<.005 and sequence and dt<1 and term and end_delta<1 and len(aa)==len(bb)
        comparisons[name]=dict(passed=bool(good),max_scaled_difference_pct=100*error,event_ms=dt,sequence_matches=sequence,
                               domain_matches=term,end_time_difference_ms=end_delta)
        passed &= good
    energies={};max_energy_flow_error=0.;max_work_grid_pct=0.
    for name,r in results.items():
        if r['status'] in ('completed','model_domain_limit'):
            ds=load(name,r)
            energies[name],err=energy_records(r,ds);max_energy_flow_error=max(max_energy_flow_error,err)
            coarse=[]
            for d in ds:
                ix=np.unique(np.r_[np.arange(0,len(d['time']),2),len(d['time'])-1])
                coarse.append({k:v[:,ix] if k=='states' else v[ix] for k,v in d.items()})
            reduced,_=energy_records(r,coarse)
            def flattened(e):
                return np.array([v for x in e.values() for v in (x.values() if isinstance(x,dict) else [x])])
            a=flattened(energies[name]);b=flattened(reduced)
            max_work_grid_pct=max(max_work_grid_pct,float(np.max(abs(a-b)/np.maximum(1,abs(a))))*100)
    passed &= max_energy_flow_error<1e-8
    passed &= max_work_grid_pct<.5
    summaries={};fig,axes=plt.subplots(2,1,figsize=(11,8),sharex=True)
    for name,r in results.items():
        if not name.endswith('-nominal') or r['status'] not in ('completed','model_domain_limit'):continue
        ds=load(name,r);t=np.hstack([d['time']-r['arrest_s'] for d in ds]);p=np.vstack([d['pressures'] for d in ds]);q=np.vstack([d['flows'] for d in ds])
        summaries[name]=dict(status=r['status'],observed_s=float(t[-1]),domain_node=(r['termination'] or {}).get('node'),
                             final_AO_mmHg=float(p[-1,2]),minimum_AO_mmHg=float(p[:,2].min()),
                             maximum_PU_mmHg=float(p[:,5].max()),maximum_VC_mmHg=float(p[:,3].max()))
        axes[0].plot(t,p[:,2],label=name.removesuffix('-nominal'))
        axes[1].plot(t,p[:,5],label=name.removesuffix('-nominal'))
    axes[0].set_ylabel('Aortic pressure (mmHg)');axes[0].legend(fontsize=8,ncol=3)
    axes[1].set_ylabel('Pulmonary venous pressure (mmHg)');axes[1].set_xlabel('Time after arrest (s)')
    for a in axes:
        for t in (4,6,7):a.axvline(t,color='grey',ls=':',alpha=.5)
    fig.tight_layout();fig.savefig(OUT/'circuits.png',dpi=150);plt.close(fig)
    metrics=dict(passed=bool(passed),comparisons=comparisons,scenarios=summaries,
                 maximum_relative_blood_drift=max(r.get('relative_blood_drift',0) for r in results.values()),
                 maximum_septal_residual_mmHg=max(r.get('root_residual_mmHg',0) for r in results.values()),
                 maximum_state_jump=max(r.get('state_jump',0) for r in results.values()),
                 maximum_work_grid_difference_pct=max_work_grid_pct,
                 maximum_energy_reconstruction_flow_error_mL_s=max_energy_flow_error,energy=energies)
    write=lambda p,x:p.write_text(json.dumps(x,indent=2)+'\n')
    write(OUT/'verification.json',metrics)
    provenance=[Path(__file__),ROOT/'src/cardiac_redundancy/pump_energy.py',OUT/'manifest.json',OUT/'warm-BIO/seed.npz',OUT/'warm-BIO/history.json']
    for n in names:
        provenance.append(OUT/n/'result.json')
        provenance.extend(sorted((OUT/n).glob('segment-*.npz')))
    write(OUT/'analysis-manifest.json',dict(hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in provenance}))
    lines=['# Combined M3 circuit verification','',
           'Handover, delayed/unequal support, depletion and E1/E2/BIO hydraulic fixtures. Prescribed restart is not modeled tissue recovery. Eight seconds is a numerical observation horizon.',
           '',f'Overall declared comparisons pass: **{passed}**.','',
           '| Scenario | Outcome | Observed s | Domain node | Final AO mmHg | Peak PU mmHg |',
           '|---|---|---:|---|---:|---:|']
    for n,s in summaries.items():lines.append(f"| {n} | {s['status']} | {s['observed_s']:.6f} | {s['domain_node'] or '-'} | {s['final_AO_mmHg']:.3f} | {s['maximum_PU_mmHg']:.3f} |")
    lines+=['','![Circuit pressure trajectories](circuits.png)','',
            '| Comparison | Scaled trajectory difference % | Valve timing ms | Pass |','|---|---:|---:|---|']
    for n,s in comparisons.items():lines.append(f"| {n} | {s.get('max_scaled_difference_pct','-')} | {s.get('event_ms','-')} | {s['passed']} |")
    lines+=['',f"Maximum relative blood drift {metrics['maximum_relative_blood_drift']:.3e}; septal residual {metrics['maximum_septal_residual_mmHg']:.3e} mmHg; state jump {metrics['maximum_state_jump']:.3e}.",
            '', 'Acceptance: trajectory differences <0.5% using max(1 native unit, reference peak magnitude) per coordinate; same valve-event sequence, event timing <1 ms; same domain node and terminal time within 1 ms. A matched domain stop passes localization verification, not support feasibility.',
            '', 'The BIO network is independently warmed with cavities isolated and auxiliary coronary loads open. Warm-up history and seed are archived. Native-seed reuse checks original source/configuration hashes. Manifests retain implementation/configuration snapshots and runtime versions.',
            '', '[Protocol and code explanation](../../docs/m3-combined-circuits.md). [Automated tests](../m3-circuit-tests.xml). Rerun `scripts/run_m3_circuits.py`, then `scripts/report_m3_circuits.py` with the project Python environment.',
            '', f'Mechanical-energy records are included per scenario in verification.json. Maximum work change from approximately 5 to 10 ms sampling is {max_work_grid_pct:.6f}% (max(1 J, reference magnitude) scaling; <0.5% screen). Abstract ports retain the head/port/resistance/limiter ledger; reconstruction is checked against saved flows. BIO records signed chamber boundary work -integral(P dV), which includes passive mechanics and is not metabolic input energy. Values end at the actual domain stop when one occurs; do not compare unequal-duration totals as equal-duration budgets.',
            '', 'Limitations: no gas exchange in the ideal oxygenator; no injury/recovery prediction, restoration controller, thrombosis, anatomical isolation validation or finite fuel model. Auxiliary alpha=.25 and stress capacities are fixtures, not derived designs. Biological withdrawal suppresses activation while its passive valves remain connected; it does not clamp nonzero inertial flow.']
    (OUT/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(metrics,indent=2));return 0 if passed else 1

if __name__=='__main__':raise SystemExit(main())
