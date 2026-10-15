"""Verify F0 from the recorded periodic pressure-aware baseline; no support."""
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import platform
import sys
import traceback
import numpy as np
import scipy

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.pressure_regional import PressureRegionalModel
from cardiac_redundancy.interventions import integrate_arrest

OUT=ROOT/'reports/m3-arrest'
SEED=ROOT/'reports/m3-pressure-regional/baseline'


def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')


def run(case,protocol):
    dest=OUT/case['name']; dest.mkdir(parents=True,exist_ok=True)
    result={'case':case,'status':'running'}; write(dest/'result.json',result)
    try:
        m=PressureRegionalModel(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml')
        seed=np.load(SEED/'trajectory.npz'); start=float(seed['time_s'][-1]); y=seed['states'][:,-1]
        arrest=start+case['phase_s']
        run=integrate_arrest(m,initial_state=y,start=start,arrest_time=arrest,
                             end=arrest+protocol['post_arrest_duration_s'],sample_step=.01,**case['solver'])
        write(dest/'interventions.json',run.interventions)
        valve_events=[e for _,tr in run.segments for e in tr.events]
        write(dest/'valves.json',valve_events)
        drift=max(float(np.max(abs(tr.states[:11].sum(axis=0)-sum(y[:11]))))/sum(y[:11]) for _,tr in run.segments)
        jump=float(np.max(abs(np.array(run.interventions[0]['before']['state'])-run.interventions[0]['after']['state']))) if run.interventions else None
        root=0.; minimum=float('inf'); minq=float('inf')
        for i,(active,tr) in enumerate(run.segments):
            ps=[]; qs=[]
            for t,state in zip(tr.time,tr.states.T):
                p,meta=m.pressures(t,state,active=active)
                q=m.flows(t,state,[True]*4,active=active,pressures=p)
                root=max(root,abs(meta['septal_residual_mmHg']))
                minimum=min(minimum,min(m.domain_margins(state).values()))
                ps.append([p[n] for n in m.volume_names+('H',)])
                qs.append([q[n] for n in m.flow_names])
            minq=min(minq,float(tr.states[11:].min()))
            np.savez_compressed(dest/f'segment-{i}.npz',time_s=tr.time,states=tr.states,
                                pressures_mmHg=np.array(ps),flows_mL_s=np.array(qs),active=active)
        if drift>=1e-8 or jump!=0 or root>=1e-8 or minq < -m.settings['absolute_flow_tolerance_mL_s']:
            raise RuntimeError('Conservation, continuity, valve or algebraic check failed')
        post=run.segments[-1][1]
        p,_=m.pressures(post.time[-1],post.states[:,-1],active=False)
        q=m.flows(post.time[-1],post.states[:,-1],[True]*4,active=False,pressures=p)
        result.update(status='model_domain_limit' if run.termination else 'completed',
            termination=run.termination,arrest_time_s=arrest,post_segment=len(run.segments)-1,
            final_elapsed_s=float(post.time[-1]-arrest),relative_volume_drift=drift,
            state_jump=jump,maximum_septal_residual_mmHg=root,minimum_domain_margin_mL=minimum,
            minimum_valve_flow_mL_s=minq,final_AO_mmHg=p['AO'],final_VC_mmHg=p['VC'],
            final_brain_flow_mL_min=q['brain']*60,final_AV_flow_mL_min=q['AV']*60)
    except Exception as exc:
        result.update(status='failed',error=str(exc),traceback=traceback.format_exc())
    write(dest/'result.json',result)
    return result


def main():
    protocol=json.loads((ROOT/'configs/m3-arrest-protocol.json').read_text())
    previous=json.loads((ROOT/'reports/m3-pressure-regional/manifest.json').read_text())
    # A cached periodic state is reusable only when its constitutive model and
    # integration implementation still match. New orchestration is separate.
    for name,digest in previous['hashes'].items():
        if name.startswith(('src/','configs/')) and name!='configs/m3-pressure-regional-protocol.json':
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                raise RuntimeError(f'Periodic seed no longer matches: {name}')
    seed_result=json.loads((SEED/'result.json').read_text())
    if seed_result['status']!='passed_numerical_checks': raise RuntimeError('Unverified seed')
    files=[Path(__file__),ROOT/'configs/m3-arrest-protocol.json',SEED/'trajectory.npz',SEED/'result.json',
           *sorted((ROOT/'src/cardiac_redundancy').glob('*.py')),
           *[ROOT/p for p in ('configs/baseline.yaml','configs/m3-evidence-regional.yaml','configs/m2-frozen-overlay.yaml')]]
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    for p in files:
        target=OUT/'implementation-snapshot'/p.relative_to(ROOT)
        target.parent.mkdir(parents=True,exist_ok=True); target.write_bytes(p.read_bytes())
    write(OUT/'manifest.json',dict(protocol=protocol,hashes=hashes,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures=[pool.submit(run,c,protocol) for c in protocol['cases']]
        results=[]
        for f in as_completed(futures):
            r=f.result(); results.append(r)
            print(r['case']['name'],r['status'],r.get('error',''),flush=True)
    return int(any(r['status']=='failed' for r in results))


if __name__=='__main__':raise SystemExit(main())
