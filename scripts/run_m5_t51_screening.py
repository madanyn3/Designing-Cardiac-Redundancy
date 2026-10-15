"""Paired M5 joint metabolic replays of verified B1 and BIO trajectories.

Hydraulic trajectories are prescribed inputs: neither glucose deficits nor pump
self-demand reduce the hydraulic driver here. This tests conditional delivery,
not closed-loop energetic survival or recoverable support duration.
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
OUT=ROOT/'reports/m5-t51-screening'


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
    segments,paths=hydraulics(case);origin=float(segments[0]['time'][0]);vol=segments[0]['states'][:n,0]
    y0=model.initial_state(vol)
    # Transfer baseline CONTENT, not inventory, to each candidate blood volume.
    # Additional auxiliary blood uses the donor's baseline content; its stores
    # start full as an explicit initialization assumption. Ledgers start at zero.
    base=np.load(ROOT/'reports/m4-baseline/native-cycle.npz')
    seed=np.load(ROOT/'reports/m4-baseline/store-5-refined.npz')
    baseline=read('reports/m4-baseline/result.json')
    bv=np.array([np.interp(.375,base['time_s'],row) for row in base['states'][:16]])
    bz=np.array([np.interp(.375,seed['time_s'],row) for row in seed['oxygen_states'].T])
    for j,name in enumerate(h.volume_names):
        donor={'aLV':'PU','aRV':'VC'}.get(name,name);i=baseline['volume_names'].index(donor)
        y0[j]=vol[j]*bz[i]/bv[i]
    y0[n]=bz[16]
    for j,organ in enumerate(organs):
        if organ.name in baseline['organ_names']:
            i=baseline['organ_names'].index(organ.name)
            y0[o.base_size+3+j]=bz[16+5+2*len(baseline['organ_names'])+i]
    seed_path='reports/m4-baseline/store-5-refined.npz'
    records=[];states=[];last_flux=None
    settings=[('nominal',1e-9,.005),('refined',1e-11,.0025)]
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
        records.append(dict(per_organ_unmet_oxygen_mL=dict(zip([v.name for v in organs],rows[-1,n+2+m:n+2+2*m].tolist())),label=label,solver='RK45',rtol=tol,atol=tol*.001,max_step_s=step,
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
        inputs=paths+[seed_path,'reports/m4-baseline/native-cycle.npz','reports/m4-baseline/result.json'])
    write(OUT/case['name']/'result.json',result);return result



def hydraulics(case):
    if case['architecture']=='BIO':
        path=f"reports/m5-biological/{case['trajectory']}/refined.npz"
        result=read(f"reports/m5-biological/{case['trajectory']}/result.json")
        a=dict(np.load(ROOT/path));a['time']=a.pop('time_s');a['flows']=a.pop('flows_mL_s')
        segments=[a];paths=[path]
    else:
        directory=f"reports/m5-dynamic/{case['trajectory']}-refined"
        result=read(directory+'/result.json');segments=[];paths=[]
        for seg in result['segments']:
            path=directory+f"/segment-{seg['index']}.npz";a=dict(np.load(ROOT/path));keep=a['time_s']<=8.+1e-9
            if not keep.any():continue
            a['time']=a['time_s'][keep];a['states']=a['states'][:,keep];a['flows']=a['flows_mL_s'][keep]
            segments.append(a);paths.append(path)
    h,_=network(case['architecture'])
    if list(h.volume_names)!=result['volume_names'] or list(h.flow_names)!=result['flow_names']:raise ValueError('Trajectory ordering mismatch')
    return segments,paths

def run(case):
    result=short(case,FuelParameters(**read('configs/m4-fuel-reference.json')['parameters']))
    print(case['name'],result['status'],flush=True);return result

if __name__=='__main__':
    from concurrent.futures import ProcessPoolExecutor
    cases=[dict(name='B1_'+name,architecture='B1',trajectory=name) for name in ('transmural_gate','left_weak')]
    cases += [dict(name='BIO_'+name,architecture='BIO',trajectory=name) for name in ('reference','faster','stronger_left','asymmetric')]
    if len(sys.argv)>1:cases=[c for c in cases if c['name'] in sys.argv[1:]]
    with ProcessPoolExecutor(max_workers=3) as pool:results=list(pool.map(run,cases))
    write(OUT/('verification.json' if len(sys.argv)==1 else 'selected-verification.json'),dict(status='verified',results=results))
