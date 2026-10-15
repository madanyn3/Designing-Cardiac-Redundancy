"""Verify whole-network controls with explicit active/preservation roles."""
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib,json,sys,traceback
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.whole_network import WholeNetworkModel,validate_roles
from cardiac_redundancy.integration import integrate_source
OUT=ROOT/'reports/m3-whole-network';SEED=ROOT/'reports/m3-selective/warm-0.2/seed.npz'

def write(p,r):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(r,indent=2)+'\n')
def model():return WholeNetworkModel(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml',arterial_fraction=.2,inlet_resistance_fraction=.1)

def run(case,protocol):
    dest=OUT/case['name'];write(dest/'result.json',{'status':'running'})
    try:
        m=model();validate_roles(m,protocol);seed=np.load(SEED);t=float(seed['time_s']);y=seed['state'].copy()
        solver=case['solver'];pre=integrate_source(m,.375,start=t,initial_state=y,sample_step=.005,**solver)
        if pre.termination:raise RuntimeError('Pre-arrest domain exit')
        y=pre.states[:,-1].copy();t=float(pre.time[-1]);arrest=t
        m.set_controls(time=t,ramp_s=1,pump_on=True,right_head=20,left_head=100)
        m.set_regional_controls(time=t,factors=case['factors'])
        events=[dict(kind='arrest_and_regional_controls',time_s=t,factors=m.regional_at(t),state=y.tolist())]
        stages=[];drift=0.;root=0.;hepatic=0.;minflow=0.;termination=None;continuity=0.
        for i,finish in enumerate((1,3,5,8)):
            if i==2:
                m.set_regional_controls(time=t,ramp_s=2,factors=dict.fromkeys(m.regional_factors,1.))
                events.append(dict(kind='reopen',time_s=t,initial_factors=m.regional_at(t),target=m.regional_factors,state=y.tolist()))
            tr=integrate_source(m,arrest+finish-t,start=t,initial_state=y,active=False,sample_step=.005,**solver)
            continuity=max(continuity,float(np.max(abs(tr.states[:,0]-y))))
            ps=[];qs=[]
            for ti,yi in zip(tr.time,tr.states.T):
                p,meta=m.pressures(ti,yi,active=False);q=m.flows(ti,yi,[True]*4,active=False,pressures=p)
                ps.append([p[n] for n in m.volume_names+('H',)]);qs.append([q[n] for n in m.flow_names])
                root=max(root,abs(meta['septal_residual_mmHg']))
                hepatic=max(hepatic,abs(sum((1 if b.recipient=='H' else -1)*q[b.name] for b in m.base.hepatic_branches)))
            drift=max(drift,float(np.max(abs(tr.states[:16].sum(axis=0)-5250)))/5250)
            minflow=min(minflow,float(tr.states[16:].min()))
            np.savez_compressed(dest/f'segment-{i}.npz',time=tr.time,states=tr.states,pressures=ps,flows=qs)
            stages.append(dict(index=i,valves=tr.events));y=tr.states[:,-1].copy();t=float(tr.time[-1])
            if tr.termination:termination=tr.termination;break
        if drift>=1e-8 or root>=1e-8 or hepatic>1e-9 or minflow < -1e-6 or continuity>1e-9:raise RuntimeError('Invariant failure')
        r=dict(status='model_domain_limit' if termination else 'completed',case=case,arrest_s=arrest,observed_s=t-arrest,
               termination=termination,stages=stages,volume_names=m.volume_names,flow_names=m.flow_names,
               relative_blood_drift=drift,maximum_hepatic_residual=hepatic,maximum_septal_residual=root,state_jump=continuity)
        write(dest/'events.json',events)
    except Exception as e:r=dict(status='failed',case=case,error=str(e),traceback=traceback.format_exc())
    write(dest/'result.json',r);return case['name'],r['status'],r.get('error','')

def main():
    protocol=json.loads((ROOT/'configs/m3-whole-network-protocol.json').read_text());validate_roles(model(),protocol)
    # All-open RHS regression is exact. Reuse only the matching periodic seed.
    old=json.loads((ROOT/'reports/m3-selective/manifest.json').read_text())
    for name,digest in old['hashes'].items():
        p=Path(name)
        if p.parts[0] in ('src','configs') and hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=digest:raise RuntimeError('Stale seed: '+name)
    paths=[Path(__file__),SEED,ROOT/'configs/m3-whole-network-protocol.json',*sorted((ROOT/'src/cardiac_redundancy').glob('*.py')),
           *[ROOT/n for n in ('configs/baseline.yaml','configs/m3-evidence-regional.yaml','configs/m2-frozen-overlay.yaml')]]
    for p in paths:
        target=OUT/'implementation-snapshot'/p.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(p.read_bytes())
    write(OUT/'manifest.json',dict(protocol=protocol,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}))
    with ProcessPoolExecutor(max_workers=4) as pool:
        results=[]
        for f in as_completed([pool.submit(run,c,protocol) for c in protocol['cases']]):
            r=f.result();results.append(r);print(*r,flush=True)
    return int(any(r[1]=='failed' for r in results))

if __name__=='__main__':raise SystemExit(main())
