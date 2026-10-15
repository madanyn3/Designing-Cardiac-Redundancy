"""Final M3 inventory/compliance sensitivity; numerical fixtures, not anatomy."""
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib,json,sys,traceback
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.whole_network import WholeNetworkModel
from cardiac_redundancy.integration import integrate_source
OUT=ROOT/'reports/m3-storage-sensitivity'

def write(p,r):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(r,indent=2)+'\n')
def model(f):
    return WholeNetworkModel(ROOT/'configs/baseline.yaml',OUT/f'factor-{f}/project.json',ROOT/'configs/m2-frozen-overlay.yaml',arterial_fraction=.2,inlet_resistance_fraction=.1)

def prepare(f):
    c=json.loads((ROOT/'configs/m3-evidence-regional.yaml').read_text());dv=du=0.
    for name,p in c['limb_pools'].items():
        dv+=(f-1)*p['volume_initial_mL'];du+=(f-1)*p['unstressed_mL']
        for k in ('volume_initial_mL','unstressed_mL','compliance_mL_mmHg'):p[k]*=f
        c['initial_physical_volumes_mL'][name]=p['volume_initial_mL']
    c['initial_physical_volumes_mL']['VC']-=dv;c['unstressed_offsets_mL']['VC']-=du
    c['status']='M3 sensitivity fixture: limb V/U/C scaled together; compensate VC V/U; no resistance recalibration'
    write(OUT/f'factor-{f}/project.json',c)
    m=model(f);y=m.initial_state.copy();previous=None;streak=0;history=[]
    for cycle in range(1,201):
        tr=integrate_source(m,.75,start=(cycle-1)*.75,initial_state=y,sample_step=.0075)
        if tr.termination:raise RuntimeError('Storage warm-up domain exit')
        error=None if previous is None else float(np.max(abs(tr.states-previous)/np.maximum(1,abs(previous))))
        streak=streak+1 if error is not None and error<1e-5 else 0
        history.append(dict(cycle=cycle,error=error));write(OUT/f'factor-{f}/warmup.json',history)
        previous=tr.states;y=tr.states[:,-1].copy()
        if cycle%20==0:print('storage',f,'warm',cycle,flush=True)
        if streak>=3:
            np.savez_compressed(OUT/f'factor-{f}/seed.npz',state=y,time_s=tr.time[-1]);return f,cycle
    raise RuntimeError('No periodic seed')

def run(f,label,solver):
    dest=OUT/f'factor-{f}'/label
    try:
        m=model(f);seed=np.load(OUT/f'factor-{f}/seed.npz');t=float(seed['time_s']);y=seed['state'].copy()
        pre=integrate_source(m,.375,start=t,initial_state=y,sample_step=.005,**solver)
        if pre.termination:raise RuntimeError('Pre-arrest domain exit')
        y=pre.states[:,-1].copy();t=float(pre.time[-1]);start=t
        roles=json.loads((ROOT/'configs/m3-whole-network-protocol.json').read_text())
        m.set_regional_controls(time=t,factors=dict.fromkeys(roles['preservation_supply_edges'],.25))
        m.set_controls(time=t,ramp_s=1,pump_on=True,right_head=20,left_head=100)
        drift=0.;residual=0.;jump=0.;minq=0.;stages=[];termination=None
        dest.mkdir(parents=True,exist_ok=True)
        for i,end in enumerate((1,4,6,8)):
            if i==2:m.set_regional_controls(time=t,ramp_s=2,factors=dict.fromkeys(m.regional_factors,1.))
            tr=integrate_source(m,start+end-t,start=t,initial_state=y,active=False,sample_step=.005,**solver)
            jump=max(jump,float(np.max(abs(tr.states[:,0]-y))));ps=[];qs=[]
            for ti,yi in zip(tr.time,tr.states.T):
                p,meta=m.pressures(ti,yi,False);q=m.flows(ti,yi,[True]*4,False,p)
                ps.append([p[n] for n in m.volume_names+('H',)]);qs.append([q[n] for n in m.flow_names])
                residual=max(residual,abs(meta['septal_residual_mmHg']))
            drift=max(drift,float(np.max(abs(tr.states[:16].sum(axis=0)-5250)))/5250);minq=min(minq,float(tr.states[16:].min()))
            np.savez_compressed(dest/f'segment-{i}.npz',time=tr.time,states=tr.states,pressures=ps,flows=qs)
            stages.append(dict(index=i,events=tr.events));t=float(tr.time[-1]);y=tr.states[:,-1].copy()
            if tr.termination:termination=tr.termination;break
        if drift>=1e-8 or residual>=1e-8 or jump>1e-9 or minq < -1e-6:raise RuntimeError('Invariant failure')
        r=dict(status='model_domain_limit' if termination else 'completed',factor=f,solver=solver,start_s=start,termination=termination,stages=stages,
               relative_blood_drift=drift,septal_residual=residual,state_jump=jump,volume_names=m.volume_names,flow_names=m.flow_names)
    except Exception as e:r=dict(status='failed',error=str(e),traceback=traceback.format_exc())
    write(dest/'result.json',r);return f,label,r['status']

def main():
    factors=(.5,1.,1.5)
    settings={'nominal':dict(method='Radau',rtol=1e-7,max_step=.001),'refined':dict(method='Radau',rtol=1e-8,max_step=.0005),'BDF':dict(method='BDF',rtol=1e-7,max_step=.001)}
    paths=[Path(__file__),*sorted((ROOT/'src/cardiac_redundancy').glob('*.py')),ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml']
    for p in paths:
        d=OUT/'implementation-snapshot'/p.relative_to(ROOT);d.parent.mkdir(parents=True,exist_ok=True);d.write_bytes(p.read_bytes())
    write(OUT/'manifest.json',dict(factors=factors,solvers=settings,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}))
    with ProcessPoolExecutor(max_workers=3) as pool:
        for fut in as_completed([pool.submit(prepare,f) for f in factors]):print('seed',fut.result(),flush=True)
        results=[]
        for fut in as_completed([pool.submit(run,f,n,s) for f in factors for n,s in settings.items()]):
            r=fut.result();results.append(r);print(*r,flush=True)
    return int(any(r[-1]=='failed' for r in results))

if __name__=='__main__':raise SystemExit(main())
