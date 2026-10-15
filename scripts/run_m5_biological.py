"""T5.1 BIO control characterization. All limits are exploratory assumptions."""
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.biological_pump import BiologicalPumpModel,BiologicalParameters
from cardiac_redundancy.circuit_integration import integrate_circuit
from run_m5_dynamic import write,read
OUT=ROOT/'reports/m5-biological'
ARGS=(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml')
CONFIG=read('configs/m5-biological-pump.json')
CASES=CONFIG['cases']

def run(name):
    pars=BiologicalParameters(**(CONFIG['reference'] | CASES[name]));m=BiologicalPumpModel(*ARGS,biological=pars,arterial_fraction=.2,inlet_resistance_fraction=.1)
    # Independently warm every passive fixture, including altered chamber sizes.
    y=m.initial_state.copy();last=None;streak=0;history=[]
    for cycle in range(1,201):
        tr=integrate_circuit(m,.75,start=(cycle-1)*.75,initial_state=y,sample_step=.0075)
        if tr.termination:raise RuntimeError(str(tr.termination))
        err=None if last is None else float(np.max(abs(tr.states-last)/np.maximum(1,abs(last))))
        streak=streak+1 if err is not None and err<1e-5 else 0
        history.append(err);last=tr.states;y=tr.states[:,-1]
        if streak>=3:break
    else:raise RuntimeError('Warm-up did not converge')
    seed=y.copy();origin=float(tr.time[-1]);dest=OUT/name;dest.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(dest/'seed.npz',state=seed,time_s=origin)
    records=[];arrays=[]
    for label,tol,step in [('nominal',1e-7,.01),('refined',1e-9,.005)]:
        m=BiologicalPumpModel(*ARGS,biological=pars,arterial_fraction=.2,inlet_resistance_fraction=.1)
        pre=integrate_circuit(m,.375,start=origin,initial_state=seed,rtol=tol,max_step=step,sample_step=.005)
        if pre.termination:raise RuntimeError('Pre-arrest failure')
        t=float(pre.time[-1]);arrest=t;y=pre.states[:,-1];m.release_auxiliary(t)
        m.native.set_regional_controls(time=t,factors=dict.fromkeys(read('configs/m3-whole-network-protocol.json')['preservation_supply_edges'],.25))
        rows=[];times=[];ps=[];qs=[];events=[];root=0.;termination=None
        for end in (1.,8.):
            tr=integrate_circuit(m,arrest+end-t,start=t,initial_state=y,active=False,rtol=tol,max_step=step,sample_step=.005)
            for ti,yi in zip(tr.time,tr.states.T):
                p,meta=m.pressures(ti,yi,False);q=m.flows(ti,yi,[True]*8,False,p)
                rows.append(yi);times.append(ti-arrest);ps.append([p[n] for n in m.volume_names+('H',)]);qs.append([q[n] for n in m.flow_names])
                root=max(root,abs(meta['aux_septal_residual_mmHg']),abs(meta['septal_residual_mmHg']))
                assert abs((m.incidence@np.array(list(q[n] for n in m.flow_names)))[-1])<1e-8
            events+=tr.events;y=tr.states[:,-1];t=float(tr.time[-1]);termination=tr.termination
            if termination:break
        time,idx=np.unique(times,return_index=True);z=np.array(rows)[idx];p=np.array(ps)[idx];q=np.array(qs)[idx]
        drift=float(np.max(abs(z[:,:18].sum(axis=1)-5250)))
        if drift>1e-5 or root>1e-8 or z[:,18:].min() < -1e-7:raise RuntimeError('Invariant failure')
        np.savez_compressed(dest/f'{label}.npz',time_s=time,states=z.T,pressures_mmHg=p,flows_mL_s=q)
        mask=time>=time[-1]-2
        records.append(dict(label=label,rtol=tol,max_step_s=step,termination=termination,blood_error_mL=drift,root_residual_mmHg=root,events=events,arrest_s=arrest,
            final_two_second_MAP_mmHg=float(np.trapezoid(p[mask,2],time[mask])/(time[mask][-1]-time[mask][0])),
            auxiliary_left_L_min=float(np.trapezoid(q[mask,m.flow_names.index('aAV')],time[mask])/(time[mask][-1]-time[mask][0])*.06)))
        arrays.append((time,z,p,q))
    errors=[]
    for j in (1,2,3):
        x=arrays[0][j];b=np.array([np.interp(arrays[0][0],arrays[1][0],row) for row in arrays[1][j].T]).T
        errors.append(float(np.max(abs(x-b)/np.maximum(1,np.max(abs(x),axis=0)))))
    e0,e1=[r['events'] for r in records]
    sequence=[(e['valve'],e['action']) for e in e0]==[(e['valve'],e['action']) for e in e1]
    event_error=max([abs(a['time_s']-b['time_s']) for a,b in zip(e0,e1)],default=0)
    if max(errors)>.005 or not sequence or event_error>.001 or any(r['termination'] for r in records):raise RuntimeError('Refinement/domain failure')
    result=dict(status='verified',parameters=asdict(pars),warm_cycles=cycle,warm_error_history=history,records=records,scaled_errors=errors,event_difference_s=event_error,volume_names=m.volume_names,flow_names=m.flow_names)
    write(dest/'result.json',result);print(name,'verified',flush=True);return result

if __name__=='__main__':
    with ProcessPoolExecutor(max_workers=3) as pool:results=list(pool.map(run,CASES))
    write(OUT/'verification.json',dict(status='verified',cases=CASES,results=results))
