"""Read saved pressure/oxygen trajectories; classify screens without viability claims."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
from scipy.interpolate import interp1d
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.constraints import pressure_windows,oxygen_screen,combined_status
from cardiac_redundancy.oxygen import OxygenParameters
from cardiac_redundancy.organ_oxygen import OrganOxygenTransport,map_demands
from cardiac_redundancy.whole_network import WholeNetworkModel


def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))


def main():
    policy=read('configs/m4-constraint-policy.json'); contract=read(policy['pressure_contract'])['pressure']
    source=ROOT/'reports/m4-demand-transitions'; prior=read('reports/m4-demand-transitions/result.json')
    assert prior['status']=='passed'
    h=WholeNetworkModel(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml',arterial_fraction=.2,inlet_resistance_fraction=.1)
    organs=map_demands(read('data/m1-demand-budget.json')['baseline_mL_min'],read('configs/m4-organ-mapping.json')['mapping'],h.flow_names,.75)
    m=OrganOxygenTransport(h.volume_names,h.edges,organs,OxygenParameters(**read('configs/m4-oxygen-reference.json')['parameters']))
    assert prior['organ_names']==[o.name for o in organs] and prior['volume_names']==list(h.volume_names)
    n=len(h.volume_names); count=len(organs); used=[]
    def hydraulics(label):
        directory=ROOT/f'reports/m3-whole-network/preservation_quarter-{label}'
        chunks=[]
        for i in range(4):
            path=directory/f'segment-{i}.npz';used.append(path);chunks.append(dict(np.load(path)))
        time=np.concatenate([s['time'] for s in chunks]); unique,rev=np.unique(time[::-1],return_index=True); idx=len(time)-1-rev
        # Retain the right-hand event sample; interpolation is a sampled-data diagnostic.
        arrays={k:np.concatenate([s[k].T if k=='states' else s[k] for s in chunks])[idx] for k in ('states','pressures','flows')}
        return unique-unique[0],unique[0],arrays
    ht,origin,hyd=hydraulics('nominal');rt,_,ref=hydraulics('refined')
    ps=interp1d(ht,hyd['pressures'],axis=0);vf=interp1d(ht,hyd['states'][:,:n],axis=0);qf=interp1d(ht,hyd['flows'],axis=0)
    lower=contract['proposed_lower_exploratory_screens_mmHg'];thresholds=dict(mean=lower['MAP'],maximum=lower['SBP'],minimum=lower['DBP'])
    pressure=pressure_windows(ht,hyd['pressures'][:,h.volume_names.index('AO')],thresholds)
    refined_pressure=pressure_windows(rt,ref['pressures'][:,h.volume_names.index('AO')],thresholds)
    discrepancy=max(abs(a['metrics_mmHg'][k]-b['metrics_mmHg'][k]) for a,b in zip(pressure['windows'],refined_pressure['windows']) for k in thresholds)
    pressure_agrees=all(a['status']==b['status'] for a,b in zip(pressure['windows'],refined_pressure['windows']))
    if not pressure_agrees:pressure['status']='unassessed'
    sensitivity=[]
    for mean in contract['mean_AO_cases_mmHg']:
        for maximum in contract['stricter_SBP_screen_cases_mmHg']:
            for minimum in contract['DBP_screen_cases_mmHg']:
                screen=pressure_windows(ht,hyd['pressures'][:,h.volume_names.index('AO')],dict(mean=mean,maximum=maximum,minimum=minimum))
                sensitivity.append(dict(mean=mean,maximum=maximum,minimum=minimum,status=screen['status']))
    out=ROOT/'reports/m4-constraints';out.mkdir(parents=True,exist_ok=True); results=[]
    for case in [c for c in prior['cases'] if c['solver']=='refined']:
        name=case['case']['name']; runs={}
        for solver in ('nominal','refined'):
            path=source/f'{name}-{solver}.npz';used.append(path);data=dict(np.load(path))
            t,idx=np.unique(data['time_s'],return_index=True);state=data['oxygen_states'][idx];demand=data['demand_mL_min'][idx]
            diagnostics=[m.evaluate(vf(ti),yi,qf(ti),4690/60,di)[1] for ti,yi,di in zip(t,state,demand)]
            flux={k:np.array([d[k] for d in diagnostics]) for k in ('demand','delivered','extractable','consumed','unmet')}
            runs[solver]=(t,state,flux)
        t,state,flux=runs['refined'];tn,_,fn=runs['nominal']
        errors=np.max(np.abs(flux['unmet']-interp1d(tn,fn['unmet'],axis=0)(t)),axis=0)
        screens=[];pvalues=ps(t);flow=qf(t);regional=[]
        for j,o in enumerate(organs):
            tol=max(policy['oxygen_tolerance']['absolute_floor_mL_s'],10*errors[j])
            result=oxygen_screen(t,flux['demand'][:,j],flux['consumed'][:,j],flux['delivered'][:,j],state[:,n+2+count+j],o.demand_mL_min,tol)
            result.update(organ=o.name,role='active' if o.parent_pool in policy['active_parent_pools'] else 'preservation',paired_rate_error_mL_s=float(errors[j]))
            result['tolerance_sensitivity']={str(mult):oxygen_screen(t,flux['demand'][:,j],flux['consumed'][:,j],flux['delivered'][:,j],state[:,n+2+count+j],o.demand_mL_min,max(1e-8,mult*errors[j]))['status'] for mult in (1,10,100)}
            screens.append(result)
            _,a,b=next(e for e in h.edges if e[0]==o.edge);nodes=(*h.volume_names,'H')
            pin=pvalues[:,nodes.index(a)];pout=pvalues[:,nodes.index(b)]
            effective=pin-pout
            if o.edge=='brain':effective=pin-np.maximum(pout,h.base.closure_settings['ICP_mmHg']+h.base.closure_settings['tone_offset_mmHg'])
            if o.edge.startswith('coronary_'):
                pc=[]
                for ti,yi in zip(t,interp1d(ht,hyd['states'],axis=0)(t)):
                    p,_=h.pressures(origin+ti,yi,active=False);pc.append(h.base.closing_pressures(p)[o.edge])
                effective=pin-np.maximum(pout,pc)
            regional.append(dict(organ=o.name,edge=o.edge,canonical_inlet=a,canonical_outlet=b,
                inlet_range_mmHg=[float(pin.min()),float(pin.max())],outlet_range_mmHg=[float(pout.min()),float(pout.max())],
                forward_effective_gradient_range_mmHg=[float(effective.min()),float(effective.max())],
                reverse_flow_observed=bool(np.any(flow[:,h.flow_names.index(o.edge)]<0)),pressure_adequacy='unassessed',
                note='edge endpoints, not tissue capillary pressures; liver edge begins at mixed hepatic junction'))
        numerical='valid' if np.all(np.isfinite(state)) and np.min(state[:,:n+1])>=0 and np.max(np.abs(m.balance(state)-m.balance(state[0])))<1e-7 else 'numerical_failure'
        combined=combined_status(numerical,pressure['status'],[s['status'] for s in screens],policy['unresolved_requirements'])
        results.append(dict(case=name,assessment=combined,oxygen=screens,regional_pressures=regional))
        np.savez_compressed(out/f'{name}-diagnostics.npz',time_s=t,**flux,central_and_reservoir_pressures_mmHg=pvalues)
    inputs=[Path(__file__),*used,source/'result.json',ROOT/'configs/m4-constraint-policy.json',ROOT/policy['pressure_contract'],
        *sorted((ROOT/'src/cardiac_redundancy').glob('*.py')),
        *[ROOT/p for p in ('data/m1-demand-budget.json','configs/m4-organ-mapping.json','configs/m4-oxygen-reference.json','configs/baseline.yaml','configs/m3-evidence-regional.yaml','configs/m2-frozen-overlay.yaml')]]
    payload=dict(status='assessment_completed',pressure=pressure,pressure_refinement_max_difference_mmHg=discrepancy,
        pressure_classification_agrees=pressure_agrees,pressure_sensitivity=sensitivity,organ_names=[o.name for o in organs],
        pressure_node_order=[*h.volume_names,'H'],flux_units='mL O2 STPD/s',cases=results,
        hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs})
    (out/'result.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    lines=['# M4 joint constraint assessment — 2026-10-07','','**Research screens only. Biological recovery remains unassessed for every run.**','',
        f"Central AO lower-screen result: **{pressure['status']}**. Eight complete one-second windows; mean/maximum/minimum thresholds 65/90/50 mmHg. Maximum nominal/refined window-metric difference {discrepancy:.6g} mmHg; classifications agree: {pressure_agrees}.",'',
        '| Case | Active branches with oxygen violations | Preservation branches with oxygen violations | Overall |','|---|---:|---:|---|']
    for r in results:
        a=sum(s['role']=='active' and s['status']=='violated' for s in r['oxygen']);p=sum(s['role']=='preservation' and s['status']=='violated' for s in r['oxygen'])
        lines.append(f"| {r['case']} | {a} | {p} | {r['assessment']['status']} |")
    lines+=['','A branch count is not a severity or survival score. Liver is assessed after depleted portal mixing. All preservation branches retain requirements. Deficits are never reset; their baseline-demand-equivalent minutes are not tolerated ischemic minutes.', '',
        'The saved JSON lists each organ, actual paired-rate tolerance, peak and accumulated unmet demand, first sampled crossing bracket, extraction undefined at zero supply, regional pressure ranges and all unresolved requirements. Pressure sensitivity covers the 40 M2 combinations. Incomplete windows are unassessed; failures of numerics/domain are classified separately by the assessment API.', '',
        'No eight-second result predicts a repair window. Artificial initial oxygenation, prescribed ventilation and hypothetical demand ramps remain as in the transition report. Regional pressure minima, functional recovery, fuel delivery, CO2, microvascular distribution and blood/reperfusion compatibility are unresolved.', '',
        'Reproduce: `.venv\\Scripts\\python.exe scripts/assess_m4_constraints.py`. [Detailed results](result.json); [policy and sources](../../docs/m4-constraints.md).']
    (out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(dict(cases=len(results),pressure_status=pressure['status'],pressure_refinement_mmHg=discrepancy)))


if __name__=='__main__':main()
