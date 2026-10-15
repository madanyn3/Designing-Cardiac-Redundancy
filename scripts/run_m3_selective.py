"""Versioned storage/routing/B1 hydraulic experiments; no recovery claims."""
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib,json,platform,sys,traceback
from pathlib import Path
import numpy as np
import scipy
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.selective import SelectiveModel
from cardiac_redundancy.integration import integrate_source
OUT=ROOT/'reports/m3-selective'

def write(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8')

def model(f):
    return SelectiveModel(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml',arterial_fraction=f,inlet_resistance_fraction=.1)

def warm(f):
    dest=OUT/f'warm-{f}';m=model(f);y=m.initial_state.copy();previous=None;streak=0;history=[]
    for cycle in range(1,201):
        tr=integrate_source(m,.75,start=(cycle-1)*.75,initial_state=y,sample_step=.0075)
        if tr.termination:raise RuntimeError(f'Warm-up domain exit {tr.termination}')
        error=None if previous is None else float(np.max(abs(tr.states-previous)/np.maximum(1,abs(previous))))
        streak=streak+1 if error is not None and error<1e-5 else 0
        history.append({'cycle':cycle,'error':error});write(dest/'warmup.json',history)
        previous=tr.states;y=tr.states[:,-1].copy()
        if streak>=3:break
    else:raise RuntimeError('No periodic seed')
    np.savez_compressed(dest/'seed.npz',state=y,time_s=tr.time[-1])
    write(dest/'result.json',dict(cycles=cycle,error=error,fraction=f))
    return f,cycle

def run(case):
    dest=OUT/case['name'];write(dest/'result.json',{'status':'running'})
    try:
        m=model(case['fraction']);seed=np.load(OUT/f"warm-{case['fraction']}"/'seed.npz')
        t=float(seed['time_s']);y=seed['state'].copy();solver=case['solver'];traces=[];events=[]
        pre=integrate_source(m,.375,start=t,initial_state=y,sample_step=.005,**solver)
        if pre.termination:raise RuntimeError('Pre-arrest domain exit')
        y=pre.states[:,-1].copy();t=pre.time[-1];arrest=t
        before=m.pressures(t,y,active=True)[0]
        m.set_controls(time=t,ramp_s=1,pump_on=case.get('support',True),right_head=20,left_head=100)
        # Vascular restriction is immediate; support head ramps separately by
        # retaining its original time ramp while setting vascular endpoints.
        c=case['controls'];start,duration,old,new=m.ramp
        for k,v in c.items():old[k]=v;new[k]=v
        events.append(dict(kind='arrest_and_restriction',time_s=t,state_before=y.tolist(),state_after=y.tolist(),pressure_before=before,pressure_after=m.pressures(t,y,active=False)[0],controls=c))
        termination=None;root=0.;drift=0.;minmargin=float('inf');minflow=float('inf')
        # End exactly at ramp knots and the reopening event. No solver step
        # straddles a discontinuous control change.
        for index,(end_offset,reopen) in enumerate([(1,False),(4,False),(6,True),(10,False)]):
            if reopen:
                before_controls=m.control_at(t)
                m.set_controls(time=t,ramp_s=case['reopen_ramp_s'],inlet=1,outlet=1,collateral=0)
                events.append(dict(kind='reopen',time_s=t,state_before=y.tolist(),state_after=y.tolist(),before=before_controls,after=m.controls))
            tr=integrate_source(m,arrest+end_offset-t,start=t,initial_state=y,active=False,sample_step=.005,**solver)
            ps=[];qs=[]
            for ti,yi in zip(tr.time,tr.states.T):
                p,meta=m.pressures(ti,yi,active=False);q=m.flows(ti,yi,[True]*4,active=False,pressures=p)
                ps.append([p[n] for n in m.volume_names+('H',)]);qs.append([q[n] for n in m.flow_names])
                root=max(root,abs(meta['septal_residual_mmHg']))
                minmargin=min(minmargin,min(m.domain_margins(yi).values()))
            drift=max(drift,float(np.max(abs(tr.states[:16].sum(axis=0)-5250)))/5250)
            minflow=min(minflow,float(tr.states[16:].min()))
            np.savez_compressed(dest/f'segment-{index}.npz',time_s=tr.time,states=tr.states,pressures=np.array(ps),flows=np.array(qs))
            traces.append(dict(index=index,start_s=float(tr.time[0]),end_s=float(tr.time[-1]),valves=tr.events))
            t=tr.time[-1];y=tr.states[:,-1].copy()
            if tr.termination:termination=tr.termination;break
        if drift>=1e-8 or root>=1e-8 or minflow < -m.settings['absolute_flow_tolerance_mL_s']:raise RuntimeError('Numerical invariant failure')
        result=dict(status='model_domain_limit' if termination else 'completed',case=case,arrest_s=arrest,
            termination=termination,segments=traces,relative_volume_drift=drift,root_residual=root,
            minimum_domain_margin=minmargin,minimum_valve_flow=minflow,final_AO_mmHg=ps[-1][2],
            final_limb_volume_mL=float(sum(y[6:16])),volume_names=m.volume_names,flow_names=m.flow_names)
        write(dest/'events.json',events);write(dest/'effective.json',dict(storage_fraction=m.fraction,inlet_resistance_fraction=m.inlet_fraction,U=m.U.tolist(),C=m.C,edges=m.edges))
    except Exception as e:result=dict(status='failed',case=case,error=str(e),traceback=traceback.format_exc())
    write(dest/'result.json',result);return case['name'],result['status'],result.get('error','')

def main():
    protocol=json.loads((ROOT/'configs/m3-selective-protocol.json').read_text())
    paths=[Path(__file__),ROOT/'configs/m3-selective-protocol.json',*sorted((ROOT/'src/cardiac_redundancy').glob('*.py')),
           *[ROOT/n for n in ('configs/baseline.yaml','configs/m3-evidence-regional.yaml','configs/m2-frozen-overlay.yaml')]]
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    for p in paths:
        d=OUT/'implementation-snapshot'/p.relative_to(ROOT);d.parent.mkdir(parents=True,exist_ok=True);d.write_bytes(p.read_bytes())
    write(OUT/'manifest.json',dict(protocol=protocol,hashes=hashes,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    with ProcessPoolExecutor(max_workers=3) as pool:
        for future in as_completed([pool.submit(warm,f) for f in sorted(set(c['fraction'] for c in protocol['cases']))]):
            print('periodic seed',future.result(),flush=True)
    with ProcessPoolExecutor(max_workers=4) as pool:
        results=[]
        for future in as_completed([pool.submit(run,c) for c in protocol['cases']]):
            r=future.result();results.append(r);print(*r,flush=True)
    return int(any(r[1]=='failed' for r in results))

if __name__=='__main__':raise SystemExit(main())
