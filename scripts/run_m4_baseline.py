"""Periodic unsuppressed oxygen reference on verified native M3 hydraulics.

Compares no-credit and finite-store assumptions; does not fit observations.
"""
import json,sys,hashlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import root
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.whole_network import WholeNetworkModel
from cardiac_redundancy.integration import integrate_source
from cardiac_redundancy.organ_oxygen import map_demands
from cardiac_redundancy.oxygen import OxygenParameters,oxygen_pressure
from cardiac_redundancy.respiratory import RespiratoryOxygenTransport
from cardiac_redundancy.oxygen_numerics import bounded_jacobian


def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))


def solve_store(task):
    store,resume=task
    out=ROOT/'reports/m4-baseline'
    h=WholeNetworkModel(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml',arterial_fraction=.2,inlet_resistance_fraction=.1)
    z=np.load(out/'native-cycle.npz');t=z['time_s'];n=len(h.volume_names)
    vf=interp1d(t,z['states'][:n],axis=1);qf=interp1d(t,z['flows'],axis=0)
    organs=map_demands(read('data/m1-demand-budget.json')['baseline_mL_min'],read('configs/m4-organ-mapping.json')['mapping'],h.flow_names,.75)
    records=[]
    m=RespiratoryOxygenTransport(h.volume_names,h.edges,organs,OxygenParameters(**read('configs/m4-oxygen-reference.json')['parameters']),store_seconds=store)
    y=m.initial_state(vf(0));initial=m.balance(y);b=m.base_size;physical=np.r_[np.arange(n+1),b,np.arange(b+3,len(y))]
    # Solve a mean-flow balance only to initialize the periodic calculation.
    # The final accepted reference retains the full pulsatile hydraulic cycle.
    av=np.trapezoid(z['states'][:n],t,axis=1)/.75
    aq=np.trapezoid(z['flows'],t,axis=0)/.75
    average_seed=m.initial_state(av)
    def average_residual(x):
        trial=average_seed.copy();trial[:n+1]=x
        return m.rhs(av,trial,aq)[:n+1]
    steady=root(average_residual,average_seed[:n+1],tol=1e-9)
    if np.max(abs(steady.fun))>1e-8:raise RuntimeError('Mean-flow initializer failed')
    y[:n]=steady.x[:n]/av*vf(0);y[n]=steady.x[n]
    seed_path=out/f'store-{store}-initial-seed.npz'
    if resume and (out/f'store-{store}-refined.npz').exists():
        y=np.load(out/f'store-{store}-refined.npz')['oxygen_states'][0].copy()
        np.savez_compressed(seed_path,state=y)
    warm_history=[]
    for step in ((.00125,.000625) if resume else (.0025,.00125,.000625)):
        grid=np.linspace(0,.75,round(.75/step)+1);mid=(grid[1:]+grid[:-1])/2
        vv=vf(grid).T;qq=qf(grid);vm=vf(mid).T;qm=qf(mid);streak=0
        for cycle in range(1,1001):
            y[n+1:b]=0.;y[b+2]=0.;previous=y.copy()
            for k in range(len(grid)-1):
                k1=m.rhs(vv[k],y,qq[k]);k2=m.rhs(vm[k],y+step*k1/2,qm[k])
                k3=m.rhs(vm[k],y+step*k2/2,qm[k]);k4=m.rhs(vv[k+1],y+step*k3,qq[k+1])
                y+=step*(k1+2*k2+2*k3+k4)/6
            error=float(np.max(abs(y[physical]-previous[physical])/np.maximum(1,abs(y[physical]))))
            warm_history.append(dict(step_s=step,cycle=cycle,maximum_scaled_change=error))
            streak=streak+1 if error<1e-7 else 0
            if cycle%50==0:print('baseline store',store,'RK4 step',step,'cycle',cycle,'periodic error',error,flush=True)
            if streak>=2:break
        if streak<2:raise RuntimeError('Periodic fixed-step iteration failed')
    chunk=cycle
    # Ledger reset only defines the reporting interval; physical inventories persist.
    physical_seed=y.copy();y[n+1:b]=0.;y[b+2]=0.
    pairs={};metrics=[]
    for label,rtol,step in [('nominal',1e-8,.00125),('refined',1e-10,.000625)]:
        rhs=lambda ti,x:m.rhs(vf(min(ti,.75)),x,qf(min(ti,.75)))
        sol=solve_ivp(rhs, (0,.75),y,method='Radau',jac=bounded_jacobian(rhs),rtol=rtol,atol=rtol*.01,max_step=step,t_eval=t)
        if not sol.success:raise RuntimeError(sol.message)
        states=sol.y.T;pairs[label]=states
        diagnostics=[m.evaluate(vf(ti),yi,qf(ti))[1] for ti,yi in zip(t,states)]
        uptake=np.array([d['consumed'] for d in diagnostics]);exchange=np.array([d['blood_extracted'] for d in diagnostics])
        unmet=np.array([d['unmet'] for d in diagnostics]);demand=np.array([d['demand'] for d in diagnostics])
        error=float(np.max(abs(m.balance(states)-m.balance(y))))
        store_values=states[:,b+3:]
        if error>1e-7 or np.min(states[:,:n+1])<0 or store_values.min() < -1e-7 or np.max(store_values-m.store_capacity)>1e-7:raise RuntimeError('Baseline accounting/domain failure')
        total=np.trapezoid(uptake,t,axis=0)/.75*60
        ao=states[:,h.volume_names.index('AO')]/vf(t)[h.volume_names.index('AO')]
        vc=states[:,h.volume_names.index('VC')]/vf(t)[h.volume_names.index('VC')]
        r=dict(store_seconds=store,solver=label,rtol=rtol,max_step_s=step,initialization='mean-flow root or hashed resume seed; refined fixed-step periodic iteration',
            consumed_mL_min=float(total.sum()),demand_mL_min=float(sum(o.demand_mL_min for o in organs)),
            unmet_mL_min=float(np.trapezoid(unmet.sum(axis=1),t)/.75*60),
            per_organ_consumed_mL_min=total.tolist(),AO_content_mL_dL=float(np.trapezoid(ao,t)/.75*100),
            VC_content_mL_dL=float(np.trapezoid(vc,t)/.75*100),AO_PO2_mmHg=float(np.trapezoid([oxygen_pressure(x) for x in ao],t)/.75),
            maximum_balance_error_mL=error,final_physical_period_change=float(np.max(abs(states[-1,physical]-states[0,physical])/np.maximum(1,abs(states[0,physical])))))
        metrics.append(r);records.append(r)
        np.savez_compressed(out/f'store-{store}-{label}.npz',time_s=t,oxygen_states=states,consumed=uptake,unmet=unmet,blood_extracted=exchange,demand=demand)
    difference=float(np.max(abs(pairs['nominal']-pairs['refined'])))
    if difference>.001 or metrics[-1]['final_physical_period_change']>1e-6:raise RuntimeError(f'Baseline oxygen refinement/periodicity failure: difference={difference}, periodic={metrics[-1]["final_physical_period_change"]}')
    (out/f'store-{store}-warmup.json').write_text(json.dumps(dict(history=warm_history,refinement_difference_mL=difference),indent=2)+'\n')
    print('baseline complete',store,metrics[-1],flush=True)
    return records


def main():
    out=ROOT/'reports/m4-baseline';out.mkdir(parents=True,exist_ok=True)
    h=WholeNetworkModel(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml',arterial_fraction=.2,inlet_resistance_fraction=.1)
    seed=np.load(ROOT/'reports/m3-selective/warm-0.2/seed.npz');start=float(seed['time_s'])
    tr=integrate_source(h,.75,start=start,initial_state=seed['state'],rtol=1e-8,max_step=.0005,sample_step=.0025)
    if tr.termination:raise RuntimeError('Normal hydraulic domain exit')
    t=tr.time-start;n=len(h.volume_names);qs=[];ps=[]
    for ti,yi in zip(tr.time,tr.states.T):
        p,_=h.pressures(ti,yi);q=h.flows(ti,yi,[True]*4,pressures=p)
        qs.append([q[k] for k in h.flow_names]);ps.append([p[k] for k in (*h.volume_names,'H')])
    qs=np.array(qs);ps=np.array(ps)
    np.savez_compressed(out/'native-cycle.npz',time_s=t,states=tr.states,flows=qs,pressures=ps)
    vf=interp1d(t,tr.states[:n],axis=1);qf=interp1d(t,qs,axis=0)
    with ProcessPoolExecutor(max_workers=3) as pool:
        records=[r for group in pool.map(solve_store,[(k,'--resume' in sys.argv) for k in (0,5,15)]) for r in group]
    organs=map_demands(read('data/m1-demand-budget.json')['baseline_mL_min'],read('configs/m4-organ-mapping.json')['mapping'],h.flow_names,.75)
    paths=[Path(__file__),ROOT/'reports/m3-selective/warm-0.2/seed.npz',*sorted((ROOT/'src/cardiac_redundancy').glob('*.py')),
        *[ROOT/p for p in ('configs/baseline.yaml','configs/m3-evidence-regional.yaml','configs/m2-frozen-overlay.yaml','configs/m4-oxygen-reference.json','configs/m4-organ-mapping.json','data/m1-demand-budget.json')]]
    if '--resume' in sys.argv:paths.extend(sorted(out.glob('*-initial-seed.npz')))
    result=dict(status='verified_periodic_reference',volume_names=h.volume_names,flow_names=h.flow_names,organ_names=[o.name for o in organs],
        cases=records,hydraulic_periodic_volume_difference_mL=float(np.max(abs(tr.states[:n,-1]-tr.states[:n,0]))),
        CO_L_min=float(np.trapezoid(qs[:,h.flow_names.index('AV')],t)/.75*.06),
        validation='independent matched-subject validation unavailable; comparison is contextual, not a pass/fail calibration',
        hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
