"""Joint BIO conductance/workload/demand experiment, with archived inputs.

Every conductance change reruns hydraulics. Metabolic demand is then recalculated
from that trajectory, rather than copied from the old reference pump. The replay
remains offline: shortages do not automatically reduce prescribed contraction.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.biological_pump import BiologicalPumpModel,BiologicalParameters
from cardiac_redundancy.biological_energy import biological_power
from cardiac_redundancy.circuit_integration import integrate_circuit
from cardiac_redundancy.myocardial_demand import DemandAssumptions,side_demand
from cardiac_redundancy.pump_energy import MMHG_ML_TO_J
import run_m5_organ_preservation as replay
OUT=ROOT/'reports/m5-joint-preservation'
ARGS=(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml')
def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))
def write(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def load_joint(case):
    path=f"reports/m5-joint-preservation/{case['name']}/hydraulic-refined.npz"
    a=dict(np.load(ROOT/path));a['time']=a.pop('time_s');a['flows']=a.pop('flows_mL_s')
    return [a],[path,f"reports/m5-joint-preservation/{case['name']}/hydraulics.json"]
def run(case):
    cfg=read('configs/m5-joint-preservation.json');duration=cfg['duration_s']
    source=read(f"reports/m5-biological/{case['pump']}/result.json")
    params=source['parameters'];seedpath=f"reports/m5-biological/{case['pump']}/seed.npz";seed=np.load(ROOT/seedpath)
    roles=read('configs/m3-whole-network-protocol.json')
    normal=read('reports/m5-normal-comparison/verification.json')['runs'][-1]['metrics']['cycle']
    cost=read('configs/m5-myocardial-demand.json');native=read('data/m1-demand-budget.json')['baseline_mL_min']['heart']
    weights=read('configs/m4-organ-mapping.json')['mapping']['heart']
    dest=OUT/case['name'];dest.mkdir(parents=True,exist_ok=True)
    write(dest/'status.json',dict(status='hydraulics_running',case=case))
    runs=[];arrays=[]
    for label,tol,step,sample in [('nominal',1e-7,.01,.001),('refined',1e-9,.005,.0005)]:
        m=BiologicalPumpModel(*ARGS,biological=BiologicalParameters(**params),arterial_fraction=.2,inlet_resistance_fraction=.1)
        pre=integrate_circuit(m,.375,start=float(seed['time_s']),initial_state=seed['state'],rtol=tol,max_step=step,sample_step=sample)
        if pre.termination:raise RuntimeError('Pre-arrest domain exit')
        arrest=float(pre.time[-1]);m.release_auxiliary(arrest)
        m.native.set_regional_controls(time=arrest,factors=dict.fromkeys(roles['preservation_supply_edges'],case['preservation_conductance']))
        pieces=[];state=pre.states[:,-1];start=arrest;events=[]
        for end in (1.,duration):
            tr=integrate_circuit(m,arrest+end-start,start=start,initial_state=state,active=False,rtol=tol,max_step=step,sample_step=sample)
            if tr.termination:
                write(dest/'status.json',dict(status='hydraulic_domain_exit',case=case,termination=tr.termination));raise RuntimeError('Support domain exit')
            pieces.append(tr);events+=tr.events;state=tr.states[:,-1];start=float(tr.time[-1])
        absolute,idx=np.unique(np.concatenate([x.time for x in pieces]),return_index=True)
        y=np.concatenate([x.states for x in pieces],axis=1)[:,idx];t=absolute-arrest
        p=[];q=[];records=[];root=0.
        for ti,yi in zip(absolute,y.T):
            pressure,meta=m.pressures(ti,yi,False);flow=m.flows(ti,yi,[True]*8,False,pressure)
            p.append([pressure[n] for n in m.volume_names]);q.append([flow[n] for n in m.flow_names]);records.append(biological_power(m,ti,yi))
            root=max(root,abs(meta['septal_residual_mmHg']),abs(meta['aux_septal_residual_mmHg']))
        p=np.array(p);q=np.array(q)
        drift=float(np.max(abs(y[:18].sum(axis=0)-5250)))
        assert drift<1e-5 and root<1e-8 and y[18:].min()>-1e-7
        sides={};demand={}
        for side,j in [('left',16),('right',17)]:
            powers={key:np.array([r[side][key] for r in records]) for key in records[0][side]}
            work={key.replace('_W','_J'):float(np.trapezoid(v,t)) for key,v in powers.items() if key.endswith('_W')}
            pv=float(-MMHG_ML_TO_J*np.sum(.5*(p[1:,j]+p[:-1,j])*np.diff(y[j])))
            inertial=abs(work['valve_storage_rate_J']-(powers['valve_stored_J'][-1]-powers['valve_stored_J'][0]))
            assert abs(pv-work['chamber_J'])<.005 and inertial<.005 and abs(powers['balance_residual_W']).max()<1e-10
            positive=float(np.trapezoid(np.maximum(powers['transmural_W'],0),t)/duration)
            sides[side]=dict(work=work,positive_transmural_mean_W=positive,PV_error_J=abs(pv-work['chamber_J']),storage_error_J=inertial)
            demand[side]=side_demand(native*weights['coronary_'+side],positive,normal['sides'][side]['positive_transmural_mean_W'],60/params['period_s'],max(params[side+'_contractility'],params['septal_contractility']),DemandAssumptions(**cost['reference']),cost['normal_rate_bpm'])
        mask=t>=duration-2;ao=p[mask,m.volume_names.index('AO')];window=t[mask]
        mean=float(np.trapezoid(ao,window)/(window[-1]-window[0]))
        runs.append(dict(label=label,blood_error_mL=drift,root_residual_mmHg=root,sides=sides,demand=demand,events=events,mean_AO_mmHg=mean,min_AO_mmHg=float(ao.min()),max_AO_mmHg=float(ao.max()),pressure_screen=bool(mean>=65 and ao.min()>=50 and ao.max()<=90),auxiliary_left_L_min=float(np.trapezoid(q[mask,m.flow_names.index('aAV')],window)/(window[-1]-window[0])*.06)))
        # Retain fine output for mechanical quadrature; replay the same 5 ms
        # grid as prior T5.2 experiments to keep transport comparisons matched.
        np.savez_compressed(dest/f'fine-{label}.npz',time_s=t,states=y,pressures_mmHg=p,flows_mL_s=q)
        keep=np.unique(np.r_[np.arange(0,len(t),round(.005/sample)),len(t)-1])
        np.savez_compressed(dest/f'hydraulic-{label}.npz',time_s=t[keep],states=y[:,keep],pressures_mmHg=p[keep],flows_mL_s=q[keep])
        arrays.append((t,y.T,p,q))
        print(case['name'],label,'hydraulics/energy verified',flush=True)
    errors=[]
    for j in (1,2,3):
        a=arrays[0][j];b=np.array([np.interp(arrays[0][0],arrays[1][0],v) for v in arrays[1][j].T]).T
        errors.append(float(np.max(abs(a-b)/np.maximum(1,np.max(abs(a),axis=0)))))
    assert max(errors)<.005
    e0,e1=[r['events'] for r in runs]
    assert [(e['valve'],e['action']) for e in e0]==[(e['valve'],e['action']) for e in e1]
    assert max([abs(a['time_s']-b['time_s']) for a,b in zip(e0,e1)],default=0)<.001
    assert max(abs(runs[0]['sides'][s]['work'][k]-runs[1]['sides'][s]['work'][k]) for s in ('left','right') for k in runs[0]['sides'][s]['work'])<.005
    profile=next(c for c in read('configs/m5-organ-preservation.json')['cases'] if c['name']==cfg['demand_case'])
    metabolic=dict(case,architecture='BIO',sides=runs[-1]['demand'],profile=profile['profile'],liver_function_fraction=profile['liver_function_fraction'],auxiliary_mL_min=sum(d['total_mL_min'] for d in runs[-1]['demand'].values()))
    write(dest/'hydraulics.json',dict(status='verified',case=case,parameters=params,volume_names=m.volume_names,flow_names=m.flow_names,runs=runs,scaled_errors=errors,inputs=[seedpath,f"reports/m5-biological/{case['pump']}/result.json"]))
    replay.OUT=OUT;replay.hydraulics=load_joint
    result=replay.run(metabolic)
    write(dest/'status.json',dict(status='verified',case=case))
    return dict(case=case,hydraulics=runs[-1],metabolic=result['records'][-1])
if __name__=='__main__':
    cases=read('configs/m5-joint-preservation.json')['cases']
    if len(sys.argv)>1:cases=[c for c in cases if c['name'] in sys.argv[1:]]
    with ProcessPoolExecutor(max_workers=3) as pool:results=list(pool.map(run,cases))
    write(OUT/('verification.json' if len(sys.argv)==1 else 'selected-verification.json'),dict(status='verified',results=results))
