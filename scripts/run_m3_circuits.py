"""Combined M3 handover, stress and circuit protocol with archived provenance."""
from concurrent.futures import ProcessPoolExecutor, as_completed
import json, hashlib, sys, traceback, platform
from pathlib import Path
import numpy as np
import scipy
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.circuit import CircuitModel
from cardiac_redundancy.backup_pump import PumpParameters, PumpStage
from cardiac_redundancy.circuit_integration import integrate_circuit
OUT=ROOT/'reports/m3-circuits'
ARGS=(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml')

def write(p,r):
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(r,indent=2)+'\n')

def model(case):
    return CircuitModel(*ARGS,architecture=case['architecture'],arterial_fraction=.2,inlet_resistance_fraction=.1,
        pump_parameters=PumpParameters(PumpStage(.1,case.get('right_capacity',120)),PumpStage(.1,case.get('left_capacity',120))))

def warm_bio():
    m=model(dict(architecture='BIO'));y=m.initial_state.copy();last=None;streak=0;history=[]
    for cycle in range(1,201):
        tr=integrate_circuit(m,.75,start=(cycle-1)*.75,initial_state=y,sample_step=.0075)
        if tr.termination:raise RuntimeError('Auxiliary warm-up domain exit: '+str(tr.termination))
        error=None if last is None else float(np.max(abs(tr.states-last)/np.maximum(1,abs(last))))
        streak=streak+1 if error is not None and error<1e-5 else 0
        history.append(dict(cycle=cycle,error=error));write(OUT/'warm-BIO/history.json',history)
        y=tr.states[:,-1].copy();last=tr.states
        if cycle%10==0:print('BIO warm cycle',cycle,'error',error,flush=True)
        if streak>=3:
            np.savez_compressed(OUT/'warm-BIO/seed.npz',state=y,time_s=tr.time[-1]);return
    raise RuntimeError('BIO failed periodic criterion within 200 cycles')

def run(case,setting,protocol):
    name=case['name']+'-'+setting['label'];dest=OUT/name
    write(dest/'result.json',dict(status='running'))
    try:
        m=model(case);solver={k:v for k,v in setting.items() if k!='label'}
        seed=OUT/'warm-BIO/seed.npz' if case['architecture']=='BIO' else ROOT/'reports/m3-selective/warm-0.2/seed.npz'
        with np.load(seed) as z:y=z['state'].copy();t=float(z['time_s'])
        pre=integrate_circuit(m,.375,start=t,initial_state=y,sample_step=.005,**solver)
        if pre.termination:raise RuntimeError('Pre-arrest domain exit')
        y=pre.states[:,-1].copy();t=float(pre.time[-1]);arrest=t
        # Synthetic finite-inventory stress, not a patient initialization.
        if case.get('deplete'):
            i=m.volume_names.index(case['deplete']);j=m.volume_names.index('VC')
            excess=y[i]-(m.native.U[i]+.1);y[i]-=excess;y[j]+=excess
        roles=json.loads((ROOT/'configs/m3-whole-network-protocol.json').read_text())
        m.native.set_regional_controls(time=t,factors=dict.fromkeys(roles['preservation_supply_edges'],.25))
        delay=case.get('delay',0);active=False;events=[];segments=[];drift=0.;root=0.;jump=0.;minflow=0.;termination=None
        knots=sorted(set([0,delay,delay+1,4,6,7,protocol['duration_s']]))
        for idx,offset in enumerate(knots[:-1]):
            before=y.copy();p_before=m.pressures(t,y,active)[0]
            actions=[]
            if offset==delay:
                if m.architecture=='BIO':m.release_auxiliary(t)
                else:m.native.set_controls(time=t,ramp_s=1,pump_on=True,right_head=20,left_head=100)
                actions.append('support_release')
            if offset==4:
                m.native.set_regional_controls(time=t,ramp_s=2,factors=dict.fromkeys(m.native.regional_factors,1.))
                actions.append('regional_reopening')
                if case.get('restart')==4:active=True;actions.append('prescribed_native_restart')
            if offset==6 and case.get('restart'):
                if m.architecture=='BIO':m.withdrawal=t
                else:m.native.set_controls(time=t,ramp_s=1,right_head=0,left_head=0)
                actions.append('backup_withdrawal')
            if offset==7 and case.get('restart') and m.architecture!='BIO':
                m.native.set_controls(time=t,pump_on=False);actions.append('support_isolation')
            if actions:events.append(dict(time_s=t,actions=actions,state_jump=float(np.max(abs(y-before))),pressure_before=p_before,pressure_after=m.pressures(t,y,active)[0]))
            tr=integrate_circuit(m,knots[idx+1]-offset,start=t,initial_state=y,active=active,sample_step=.005,**solver)
            jump=max(jump,float(np.max(abs(tr.states[:,0]-y))))
            ps=[];qs=[]
            for ti,yi in zip(tr.time,tr.states.T):
                p,meta=m.pressures(ti,yi,active);q=m.flows(ti,yi,[True]*len(m.valve_names),active,p)
                ps.append([p[n] for n in m.volume_names+('H',)]);qs.append([q[n] for n in m.flow_names])
                root=max(root,abs(meta['septal_residual_mmHg']),abs(meta.get('aux_septal_residual_mmHg',0)))
            drift=max(drift,float(np.max(abs(tr.states[:m.volume_count].sum(axis=0)-5250)))/5250)
            minflow=min(minflow,float(tr.states[m.volume_count:].min()))
            np.savez_compressed(dest/f'segment-{idx}.npz',time=tr.time,states=tr.states,pressures=ps,flows=qs)
            segments.append(dict(index=idx,active=active,events=tr.events,start_s=float(tr.time[0]),end_s=float(tr.time[-1])))
            t=float(tr.time[-1]);y=tr.states[:,-1].copy()
            if tr.termination:termination=tr.termination;break
        if drift>=1e-8 or root>=1e-8 or jump>1e-9 or minflow < -1e-6:raise RuntimeError('Invariant check failed')
        r=dict(status='model_domain_limit' if termination else 'completed',case=case,solver=setting,arrest_s=arrest,termination=termination,
               segments=segments,events=events,volume_names=m.volume_names,flow_names=m.flow_names,relative_blood_drift=drift,
               root_residual_mmHg=root,state_jump=jump,minimum_valve_flow=minflow,final_AO_mmHg=ps[-1][2])
    except Exception as e:r=dict(status='failed',error=str(e),traceback=traceback.format_exc(),case=case,solver=setting)
    write(dest/'result.json',r);return name,r['status'],r.get('error','')

def main():
    protocol=json.loads((ROOT/'configs/m3-circuit-protocol.json').read_text())
    # Native seed is valid because historical source/configurations are intact.
    old=json.loads((ROOT/'reports/m3-selective/manifest.json').read_text())
    for n,h in old['hashes'].items():
        if Path(n).parts[0] in ('src','configs') and hashlib.sha256((ROOT/n).read_bytes()).hexdigest()!=h:raise RuntimeError('Stale native seed: '+n)
    files=[Path(__file__),ROOT/'configs/m3-circuit-protocol.json',*sorted((ROOT/'src/cardiac_redundancy').glob('*.py')),
           *ARGS,ROOT/'reports/m3-selective/warm-0.2/seed.npz']
    for p in files:
        target=OUT/'implementation-snapshot'/p.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(p.read_bytes())
    write(OUT/'manifest.json',dict(protocol=protocol,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                                 python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    # Run independent abstract circuits while the biological seed warms.
    results=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        warm=pool.submit(warm_bio)
        futures=[pool.submit(run,c,s,protocol) for c in protocol['cases'] if c['architecture']!='BIO' for s in protocol['solvers']]
        for f in as_completed(futures):
            r=f.result();results.append(r);print(*r,flush=True)
        warm.result()
        futures=[pool.submit(run,c,s,protocol) for c in protocol['cases'] if c['architecture']=='BIO' for s in protocol['solvers']]
        for f in as_completed(futures):
            r=f.result();results.append(r);print(*r,flush=True)
    return int(any(r[1]=='failed' for r in results))

if __name__=='__main__':raise SystemExit(main())
