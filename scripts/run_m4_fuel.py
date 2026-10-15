"""Paired coupled arrest replays and maintained-flow finite-fuel controls.

Long controls prescribe the balanced mean native flow indefinitely. They test
fuel accounting ONLY, not auxiliary-pump endurance or recoverable repair time.
"""
from dataclasses import replace
from pathlib import Path
import hashlib,json,sys
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.circuit import CircuitModel
from cardiac_redundancy.organ_oxygen import map_demands
from cardiac_redundancy.oxygen import OxygenParameters
from cardiac_redundancy.respiratory import RespiratoryOxygenTransport
from cardiac_redundancy.oxygen_numerics import bounded_jacobian
from cardiac_redundancy.fuel import (FuelParameters,GlucoseTransport,FuelOxygenTransport,
    O2_ML_STPD_PER_MMOL_GLUCOSE,GLUCOSE_MG_PER_MMOL,GLUCOSYL_MG_PER_MMOL)
from run_m4_respiratory import load_segments
OUT=ROOT/'reports/m4-fuel'


def read(path):return json.loads((ROOT/path).read_text(encoding='utf-8'))
def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def network(architecture='B1',scenario='unsuppressed'):
    h=CircuitModel(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',
        ROOT/'configs/m2-frozen-overlay.yaml',architecture=architecture,arterial_fraction=.2,inlet_resistance_fraction=.1)
    organs=map_demands(read('data/m1-demand-budget.json')['scenarios'][scenario]['organ_mL_min'],
        read('configs/m4-organ-mapping.json')['mapping'],h.flow_names,.75,10 if architecture=='BIO' else None)
    return h,organs


def short(case,p):
    write(OUT/case['name']/'result.json',dict(status='running',case=case,parameters=vars(p)))
    h,organs=network(case['architecture']);n=len(h.volume_names);m=len(organs)
    o=RespiratoryOxygenTransport(h.volume_names,h.edges,organs,
        OxygenParameters(**read('configs/m4-oxygen-reference.json')['parameters']),store_seconds=5)
    g=GlucoseTransport(h.volume_names,h.edges,organs,p,extraction_limit=.5)
    model=FuelOxygenTransport(o,g);k=model.oxygen_size
    segments,paths=load_segments(case);origin=float(segments[0]['time'][0]);vol=segments[0]['states'][:n,0]
    y0=model.initial_state(vol)
    # Use the verified respiratory run's initial state, not its later trajectory:
    # substrate-limited O2 uptake and respiratory feedback are recomputed here.
    seed_path=f"reports/m4-respiratory/{case['seed']}/refined.npz"
    seed=np.load(ROOT/seed_path);y0[:k]=seed['oxygen_states'][0]
    records=[];states=[];last_flux=None
    settings=[('nominal',1e-11,.005),('refined',1e-12,.0025)] if case.get('tight_solver') else [('nominal',1e-8,.01),('refined',1e-10,.005)]
    for label,tol,step in settings:
        y=y0.copy();times=[];rows=[];volumes=[];flows=[]
        for segment in segments:
            tt=segment['time'];vf=interp1d(tt,segment['states'][:n],axis=1);qf=interp1d(tt,segment['flows'],axis=0)
            controls=dict(mode='feedback',sweep_mL_s_STPD=100 if case['architecture']=='E1' else 0,
                          liver_function_fraction=case.get('liver_function_fraction',1.))
            sol=solve_ivp(lambda t,z:model.rhs(vf(t),z,qf(t),**controls),(tt[0],tt[-1]),y,
                method='RK45',rtol=tol,atol=tol*.001,max_step=step,t_eval=tt)
            if not sol.success:raise RuntimeError(sol.message)
            y=sol.y[:,-1];times.extend(sol.t-origin);rows.extend(sol.y.T)
            volumes.extend(vf(tt).T);flows.extend(qf(tt))
        time,idx=np.unique(times,return_index=True);rows=np.array(rows)[idx];volumes=np.array(volumes)[idx];flows=np.array(flows)[idx]
        diags=[model.evaluate(v,z,q,**controls)[1] for v,z,q in zip(volumes,rows,flows)]
        ob,gb=model.balances(rows);initial_ob,initial_gb=model.balances(y0)
        error_o=float(np.max(abs(ob-initial_ob)));error_g=float(np.max(abs(gb-initial_gb)))
        go=rows[:,k:];request=np.array([v.demand_mL_min for v in organs])/60
        identity_o=float(np.max(abs(rows[:,n+2:n+2+m]+rows[:,n+2+m:n+2+2*m]-time[:,None]*request)))
        identity_g=float(np.max(abs(go[:,n+1:n+1+m]+go[:,n+1+m:n+1+2*m]-time[:,None]*request/O2_ML_STPD_PER_MMOL_GLUCOSE)))
        stoich=float(np.max(abs(rows[:,n+2:n+2+m]-go[:,n+1:n+1+m]*O2_ML_STPD_PER_MMOL_GLUCOSE)))
        physical_o=np.c_[rows[:,:n+1],rows[:,o.base_size+1],rows[:,o.base_size+3:k]]
        if max(error_o,identity_o,stoich)>1e-6 or max(error_g,identity_g)>1e-8 or go[:,:n+1].min() < -1e-8 or physical_o.min() < -1e-7:
            raise RuntimeError('Coupled accounting/domain check failed')
        if np.max(rows[:,o.base_size+3:k]-o.store_capacity)>1e-7 or np.max(rows[:,o.base_size+1])>o.respiratory.oxygenator_gas_volume_mL_STPD+1e-7:
            raise RuntimeError('Oxygen capacity exceeded')
        deficit=np.array([d['unmet'] for d in diags]);consumed=np.array([d['consumed'] for d in diags])
        ceiling=np.array([d['glucose_ceiling_mL_s'] for d in diags])
        dest=OUT/case['name'];dest.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(dest/f'{label}.npz',time_s=time,states=rows,volumes_mL=volumes,
                            oxygen_consumed_mL_s=consumed,oxygen_unmet_mL_s=deficit,glucose_ceiling_mL_s=ceiling)
        records.append(dict(label=label,solver='RK45',rtol=tol,atol=tol*.001,max_step_s=step,
            oxygen_balance_error_mL=error_o,glucose_balance_error_mmol=error_g,
            oxygen_demand_identity_error_mL=identity_o,glucose_demand_identity_error_mmol=identity_g,
            oxidation_stoichiometry_error_mL=stoich,unmet_oxygen_mL=float(rows[-1,n+2+m:n+2+2*m].sum()),
            glycogen_remaining_mmol=float(go[-1,n]),glucose_remaining_mmol=float(go[-1,:n].sum()),
            released_mmol=float(go[-1,-3]),minimum_plasma_equivalent_mmol_L=float(np.min(go[:,:n]/volumes*1000*p.plasma_to_whole_blood_ratio))))
        states.append(rows);last_flux=(deficit,consumed)
    do=float(np.max(abs(states[0][:,:k]-states[1][:,:k])))
    dg=float(np.max(abs(states[0][:,k:]-states[1][:,k:])))
    if do>.002 or dg>1e-5:raise RuntimeError(f'Coupled refinement failed: {do}, {dg}')
    result=dict(status='verified',scope='coupled saved-hydraulic replay; no recoverable-duration inference',duration_s=float(time[-1]),case=case,parameters=vars(p),
        volume_names=h.volume_names,organ_names=[v.name for v in organs],records=records,
        paired_max_oxygen_state_difference=do,paired_max_glucose_state_difference_mmol=dg,
        inputs=paths+[seed_path])
    write(OUT/case['name']/'result.json',result);return result


def long_control(name,p,scenario='unsuppressed',liver_fraction=1.,extraction=.5):
    write(OUT/name/'result.json',dict(status='running',name=name,parameters=vars(p),scenario=scenario))
    h,organs=network(scenario=scenario);n=len(h.volume_names);m=len(organs)
    native=np.load(ROOT/'reports/m4-baseline/native-cycle.npz');t=native['time_s']
    v=np.trapezoid(native['states'][:n],t,axis=1)/(t[-1]-t[0])
    q=np.trapezoid(native['flows'],t,axis=0)/(t[-1]-t[0])
    # Frozen mean volumes need EXACTLY balanced flows. Project the tiny sampled
    # periodic residual onto the incidence nullspace, including hepatic H.
    raw=q.copy();a=h.incidence;q=q-a.T@np.linalg.pinv(a@a.T)@(a@q)
    if np.max(abs(a@q))>1e-9:raise RuntimeError('Mean-flow projection failed')
    g=GlucoseTransport(h.volume_names,h.edges,organs,p,extraction_limit=extraction)
    y0=g.initial_state(v);d=np.array([o.demand_mL_min for o in organs])/60/O2_ML_STPD_PER_MMOL_GLUCOSE
    sample=np.arange(0,21600.1,10);states=[];records=[]
    for label,tol,step in [('nominal',1e-7,30.),('refined',1e-9,15.)]:
        rhs=lambda t,y:g.evaluate(v,y,q,d,liver_function_fraction=liver_fraction,numerical_trial=True)[0]
        sol=solve_ivp(rhs,(0,21600),y0,method='Radau',jac=bounded_jacobian(rhs),rtol=tol,atol=tol*.001,max_step=step,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        if sol.y[:n+1].min() < -1e-8:raise RuntimeError('Accepted solver step has negative fuel inventory')
        y=sol.sol(sample).T;balance=float(np.max(abs(g.balance(y)-g.balance(y0))))
        identity=float(np.max(abs(y[:,n+1:n+1+m]+y[:,n+1+m:n+1+2*m]-sample[:,None]*d)))
        if balance>1e-7 or identity>1e-7 or y[:,:n+1].min() < -1e-8:raise RuntimeError('Long control accounting failed')
        deficits=np.array([g.evaluate(v,row,q,d,liver_function_fraction=liver_fraction)[1]['unmet_mmol_s'] for row in y])
        bad=np.flatnonzero(np.any(deficits>1e-8,axis=1))
        first=None if not len(bad) else [float(sample[max(0,bad[0]-1)]),float(sample[bad[0]])]
        dest=OUT/name;dest.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(dest/f'{label}.npz',time_s=sample,states=y,volumes_mL=v,flows_mL_s=q,unmet_mmol_s=deficits)
        records.append(dict(label=label,solver='Radau',rtol=tol,atol=tol*.001,max_step_s=step,
            balance_error_mmol=balance,demand_identity_error_mmol=identity,first_sampled_glucose_shortfall_bracket_s=first,
            final_glycogen_mmol=float(y[-1,n]),final_blood_glucose_mmol=float(y[-1,:n].sum()),
            total_glucose_uptake_mmol=float(y[-1,n+1:n+1+m].sum()),total_unmet_glucose_mmol=float(y[-1,n+1+m:n+1+2*m].sum()),
            total_nonoxidative_uptake_mmol=float(y[-1,-2]),released_mmol=float(y[-1,-3])))
        states.append(y)
    difference=float(np.max(abs(states[0]-states[1])))
    if difference>1e-4:raise RuntimeError(f'Long refinement failed: {difference}')
    result=dict(status='verified',scope='6-hour maintained mean native flow, glucose-only necessary-supply control; oxygen and pump energy not simulated',
        name=name,parameters=vars(p),scenario=scenario,liver_function_fraction=liver_fraction,extraction_limit=extraction,
        volume_names=h.volume_names,organ_names=[o.name for o in organs],demand_mmol_min=(d*60).tolist(),
        maximum_flow_projection_mL_s=float(np.max(abs(q-raw))),records=records,paired_max_state_difference_mmol=difference,
        inputs=['reports/m4-baseline/native-cycle.npz'])
    write(OUT/name/'result.json',result);return result


def main():
    cfg=read('configs/m4-fuel-reference.json');p=FuelParameters(**cfg['parameters']);results=[]
    cases=[
        (dict(name='native_cycle',architecture='B1',hydraulic='periodic',seed='native120_maintained'),p),
        (dict(name='B1_reference',architecture='B1',hydraulic='whole',seed='B1_feedback_store5'),p),
        (dict(name='B1_low_glucose',architecture='B1',hydraulic='whole',seed='B1_feedback_store5'),replace(p,initial_plasma_glucose_mmol_L=.001)),
        (dict(name='B1_liver_loss',architecture='B1',hydraulic='whole',seed='B1_feedback_store5',liver_function_fraction=0,tight_solver=True),p),
        (dict(name='handover',architecture='B1',hydraulic='handover',seed='handover_feedback',tight_solver=True),p),
        (dict(name='E1_fed',architecture='E1',hydraulic='E1',seed='E1_fed',tight_solver=True),p),
        (dict(name='BIO',architecture='BIO',hydraulic='BIO',seed='BIO_feedback',tight_solver=True),p)]
    controls=[('mean_reference',p,{}),('mean_no_release',replace(p,release_cap_mmol_min=0),{}),
        ('mean_total_output_scale',replace(p,release_cap_mmol_min=1.9*75/GLUCOSE_MG_PER_MMOL),{}),
        ('mean_high_release',replace(p,release_cap_mmol_min=2.1),{}),
        ('mean_export_capacity',replace(p,release_cap_mmol_min=2.1,export_increment_mmol_L=3.),{}),
        ('mean_low_store',replace(p,glycogen_mmol_L_liver=160,release_cap_mmol_min=2.1),{}),
        ('mean_high_store',replace(p,glycogen_mmol_L_liver=274,release_cap_mmol_min=2.1),{}),
        ('mean_reduced_demand',replace(p,release_cap_mmol_min=2.1),dict(scenario='stronger_reduction')),
        ('mean_nonoxidative',replace(p,release_cap_mmol_min=2.1,nonoxidative_mmol_min=.1),{})]
    for case,params in cases:
        r=read(f"reports/m4-fuel/{case['name']}/result.json") if '--manifest-only' in sys.argv else short(case,params)
        if r['status']!='verified' or r['case']!=case or r['parameters']!=vars(params):raise RuntimeError('Saved coupled case does not match protocol')
        results.append(r);print(case['name'],r['status'],flush=True)
    for name,params,kwargs in controls:
        r=read(f'reports/m4-fuel/{name}/result.json') if '--manifest-only' in sys.argv else long_control(name,params,**kwargs)
        if r['status']!='verified' or r['parameters']!=vars(params) or r['scenario']!=kwargs.get('scenario','unsuppressed'):raise RuntimeError('Saved mean-flow case does not match protocol')
        results.append(r);print(name,r['status'],flush=True)
    paths={Path(__file__),ROOT/'configs/m4-fuel-reference.json',ROOT/'configs/m4-organ-mapping.json',ROOT/'configs/m4-oxygen-reference.json',
        ROOT/'data/m1-demand-budget.json',ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml',
        ROOT/'scripts/run_m4_respiratory.py',*list((ROOT/'src/cardiac_redundancy').glob('*.py'))}
    for r in results:paths.update(ROOT/path for path in r['inputs'])
    source_hashes={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(paths)}
    outputs={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(OUT.glob('*/*')) if path.is_file()}
    write(OUT/'manifest.json',dict(status='verified',cases=len(results),paired_runs=2*len(results),
        software=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__),hashes=source_hashes,output_hashes=outputs,
        initial=dict(glycogen_mmol=p.glycogen_mmol,polymer_g=p.glycogen_mmol*GLUCOSYL_MG_PER_MMOL/1000,
                     glucose_equivalent_g=p.glycogen_mmol*GLUCOSE_MG_PER_MMOL/1000,
                     blood_glucose_mmol=5250*p.blood_concentration_mmol_mL),results=results))


if __name__=='__main__':main()
