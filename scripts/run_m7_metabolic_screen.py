"""M7 conditional metabolic sensitivity replay, derived from the frozen M6 runner.

Only declared metabolic knobs vary; hydraulic geometry and absolute oxygen
stores remain fixed. The mass ratio is a demand-model parameter in this study,
not a physical resizing experiment. Every case uses paired numerical settings.

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
from cardiac_redundancy.bilateral_support import BilateralSupport, BilateralRespiration, LungAllocation
from cardiac_redundancy.activation_timing import ActivationTiming
from cardiac_redundancy.organ_preservation import OrganPreservation
from cardiac_redundancy.organ_oxygen import map_demands
from cardiac_redundancy.oxygen import OxygenParameters
from cardiac_redundancy.respiratory import RespiratoryOxygenTransport
from cardiac_redundancy.oxygen_numerics import bounded_jacobian
from cardiac_redundancy.fuel import (FuelParameters,GlucoseTransport,FuelOxygenTransport,
    O2_ML_STPD_PER_MMOL_GLUCOSE,GLUCOSE_MG_PER_MMOL,GLUCOSYL_MG_PER_MMOL)
from run_m4_respiratory import load_segments
OUT=ROOT/'reports/m7-metabolic-screen'


def read(path):return json.loads((ROOT/path).read_text(encoding='utf-8'))
def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def network(case):
    cfg=read('configs/m5-final-comparison.json')
    h=BilateralSupport(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',
        ROOT/'configs/m2-frozen-overlay.yaml',device=case['device'],allocation=LungAllocation(**cfg['lung_allocation']),arterial_fraction=.2,inlet_resistance_fraction=.1)
    mapping=read('configs/m4-organ-mapping.json')['mapping']
    share=cfg['lung_allocation']['supported_share']
    mapping['lungs']={'lung_supported':share,'lung_preserved':1-share}
    organs=map_demands(read('data/m1-demand-budget.json')['baseline_mL_min'],mapping,h.flow_names,.75,10 if case['device']=='BIO' else None)
    if case['device']=='BIO':
        organs=tuple(replace(o,demand_mL_min=case['sides'][o.name.removeprefix('aux_coronary_')]['total_mL_min']) if o.name.startswith('aux_coronary_') else o for o in organs)
    return h,organs

def ramp_integral(t,timing):
    x=max(0.,t-timing.latency_s);u=min(x,timing.ramp_s)
    return u*u/(2*timing.ramp_s)+max(0.,x-timing.ramp_s)


def short(case,p):
    write(OUT/case['name']/'result.json',dict(status='running',case=case,parameters=vars(p)))
    h,organs=network(case);n=len(h.volume_names);m=len(organs)
    o=BilateralRespiration(h.volume_names,h.edges,organs,allocation=h.allocation,
        parameters=OxygenParameters(**read('configs/m4-oxygen-reference.json')['parameters']),store_seconds=5)
    if case['device']=='BIO':o.store_capacity[-2:]=read('configs/m5-myocardial-demand.json')['auxiliary_store_capacity_mL']
    g=GlucoseTransport(h.volume_names,h.edges,organs,p,extraction_limit=case['knobs']['extraction_limit'])
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
    y0[o.base_size+1]=bz[16]*(1-h.allocation.supported_share)
    for j,organ in enumerate(organs):
        if organ.name in ('lung_supported','lung_preserved'):
            i=baseline['organ_names'].index('lungs');fraction=h.allocation.supported_share if organ.name=='lung_supported' else 1-h.allocation.supported_share
            y0[o.base_size+3+j]=fraction*bz[16+5+2*len(baseline['organ_names'])+i]
        elif organ.name in baseline['organ_names']:
            i=baseline['organ_names'].index(organ.name)
            y0[o.base_size+3+j]=bz[16+5+2*len(baseline['organ_names'])+i]
    seed_path='reports/m4-baseline/store-5-refined.npz'
    entries={k:dict(v) for k,v in read('configs/m5-organ-preservation.json')['cases'][1]['profile'].items()}
    entries.pop('lungs');entries['lung_supported']=dict(role='active',fraction=1.,ramp_s=0.)
    entries['lung_preserved']=dict(role='preservation',fraction=.25,ramp_s=0.)
    for name in ('coronary_left','coronary_right','kidney_protected'):entries[name]['role']='active'
    if case['device']=='mechanical':
        for name in ('aux_coronary_left','aux_coronary_right'):entries.pop(name)
    for entry in entries.values():
        if entry['role']=='preservation':entry['fraction']=case['knobs']['preservation_fraction']
    profile=OrganPreservation(organs,entries)
    trajectory=np.load(ROOT/'reports/m6-coupled/observed_arrest/hydraulic-refined.npz')
    knots=trajectory['time_s'];amplitude=trajectory['aux_amplitude']
    duration=float(knots[-1]);mean_amp=float(np.trapezoid(amplitude,knots)/duration)
    table=np.array([profile.demand(float(t)) for t in knots])
    if case.get('demand_restore') is not None:
        fraction=np.clip((knots-case['demand_restore'])/2,0,1)
        for j,organ in enumerate(organs):
            if organ.parent_pool!='auxiliary_heart':
                table[:,j]+=(organ.demand_mL_min-table[:,j])*fraction
    aux_idx=[i for i,o in enumerate(organs) if o.parent_pool=='auxiliary_heart']
    for i in aux_idx:
        v=case['sides'][organs[i].name.removeprefix('aux_coronary_')]
        table[:,i]=v['maintenance_mL_min']+amplitude*(v['activation_mL_min']+(v['work_mL_min']/mean_amp if mean_amp>0 else 0))
    cumulative=np.vstack([np.zeros(m),np.cumsum((table[:-1]+table[1:])*np.diff(knots)[:,None]/120,axis=0)])
    def demand_at(t,integral=False):
        t=float(np.clip(t,knots[0],knots[-1]));i=min(len(knots)-2,max(0,np.searchsorted(knots,t,side='right')-1))
        dt=t-knots[i];slope=(table[i+1]-table[i])/(knots[i+1]-knots[i])
        return cumulative[i]+(table[i]*dt+.5*slope*dt*dt)/60 if integral else table[i]+slope*dt
    def controls_at(t):
        return dict(mode='apnea' if case.get('apnea') is not None and t>=case['apnea'] else 'feedback',
                    sweep_mL_s_STPD=0.,liver_function_fraction=case['liver_function_fraction'],demand_mL_min=demand_at(t))
    records=[];states=[];last_flux=None
    settings=[('nominal',1e-9,.005),('refined',1e-11,.0025)]
    for label,tol,step in settings:
        y=y0.copy();times=[];rows=[];volumes=[];flows=[]
        for segment in segments:
            tt=segment['time'];vf=interp1d(tt,segment['states'][:n],axis=1);qf=interp1d(tt,segment['flows'],axis=0)
            controls=dict(mode='feedback',sweep_mL_s_STPD=100 if case['architecture']=='E1' else 0,
                          liver_function_fraction=case.get('liver_function_fraction',1.))
            sol=solve_ivp(lambda t,z:model.rhs(vf(t),z,qf(t),**controls_at(t)),(tt[0],tt[-1]),y,
                method='RK45',rtol=tol,atol=tol*.001,max_step=step,t_eval=tt)
            if not sol.success:raise RuntimeError(sol.message)
            y=sol.y[:,-1];times.extend(sol.t-origin);rows.extend(sol.y.T)
            volumes.extend(vf(tt).T);flows.extend(qf(tt))
        time,idx=np.unique(times,return_index=True);rows=np.array(rows)[idx];volumes=np.array(volumes)[idx];flows=np.array(flows)[idx]
        diags=[model.evaluate(v,z,q,**controls_at(t+origin))[1] for t,v,z,q in zip(time,volumes,rows,flows)]
        ob,gb=model.balances(rows);initial_ob,initial_gb=model.balances(y0)
        error_o=float(np.max(abs(ob-initial_ob)));error_g=float(np.max(abs(gb-initial_gb)))
        go=rows[:,k:];request=np.array([v.demand_mL_min for v in organs])/60
        expected=np.array([demand_at(t,True) for t in time])
        identity_o=float(np.max(abs(rows[:,n+2:n+2+m]+rows[:,n+2+m:n+2+2*m]-expected)))
        identity_g=float(np.max(abs(go[:,n+1:n+1+m]+go[:,n+1+m:n+1+2*m]-expected/O2_ML_STPD_PER_MMOL_GLUCOSE)))
        stoich=float(np.max(abs(rows[:,n+2:n+2+m]-go[:,n+1:n+1+m]*O2_ML_STPD_PER_MMOL_GLUCOSE)))
        physical_o=np.c_[rows[:,:n+1],rows[:,o.base_size+1],rows[:,o.base_size+3:k]]
        if max(error_o,identity_o,stoich)>1e-6 or max(error_g,identity_g)>1e-8 or go[:,:n+1].min() < -1e-8 or physical_o.min() < -1e-7:
            raise RuntimeError('Coupled accounting/domain check failed')
        if np.max(rows[:,o.base_size+3:k]-o.store_capacity)>1e-7 or np.max(rows[:,o.base_size+1])>o.parameters.alveolar_volume_mL_STPD*(1-h.allocation.supported_share)+1e-7:
            raise RuntimeError('Oxygen capacity exceeded')
        supported_gas=rows[:,n]-rows[:,o.base_size+1]
        if supported_gas.min() < -1e-7 or supported_gas.max()>o.parameters.alveolar_volume_mL_STPD*h.allocation.supported_share+1e-7:raise RuntimeError('Supported lung gas outside domain')
        deficit=np.array([d['unmet'] for d in diags]);consumed=np.array([d['consumed'] for d in diags])
        ceiling=np.array([d['glucose_ceiling_mL_s'] for d in diags])
        dest=OUT/case['name'];dest.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(dest/f'metabolic-{label}.npz',time_s=time,states=rows,volumes_mL=volumes,
                            oxygen_consumed_mL_s=consumed,oxygen_unmet_mL_s=deficit,glucose_ceiling_mL_s=ceiling,requested_O2_mL_s=np.array([d['demand'] for d in diags]),lung_gas_O2_mL=np.array([d['lung_gas_O2_mL'] for d in diags]),lung_fractions=np.array([d['lung_fractions'] for d in diags]),lung_transfer_mL_s=np.array([d['lung_transfer_mL_s'] for d in diags]))
        aux_bad=np.flatnonzero(np.any(deficit[:,aux_idx]>1e-8,axis=1))
        first_aux=None if not len(aux_bad) else [float(time[max(0,aux_bad[0]-1)]),float(time[aux_bad[0]])]
        group_metrics={}
        for role,indices in profile.groups().items():
            bad=np.flatnonzero(np.any(deficit[:,indices]>1e-8,axis=1))
            group_metrics[role]=dict(
                target_final_mL_min=float(demand_at(time[-1])[indices].sum()),
                requested_O2_mL=float(expected[-1,indices].sum()),
                unmet_O2_mL=float(rows[-1,n+2+m+np.array(indices,dtype=int)].sum()),
                first_shortfall_bracket_s=None if not len(bad) else [float(time[max(0,bad[0]-1)]),float(time[bad[0]])])
        records.append(dict(groups=group_metrics,first_sampled_auxiliary_shortfall_bracket_s=first_aux,prescribed_driver_status='unsupported_under_declared_demand' if len(aux_bad) else 'no_auxiliary_shortfall_observed',per_organ_unmet_oxygen_mL=dict(zip([v.name for v in organs],rows[-1,n+2+m:n+2+2*m].tolist())),label=label,solver='RK45',rtol=tol,atol=tol*.001,max_step_s=step,
            oxygen_balance_error_mL=error_o,glucose_balance_error_mmol=error_g,
            oxygen_demand_identity_error_mL=identity_o,glucose_demand_identity_error_mmol=identity_g,
            oxidation_stoichiometry_error_mL=stoich,unmet_oxygen_mL=float(rows[-1,n+2+m:n+2+2*m].sum()),
            glycogen_remaining_mmol=float(go[-1,n]),glucose_remaining_mmol=float(go[-1,:n].sum()),
            released_mmol=float(go[-1,-3]),minimum_plasma_equivalent_mmol_L=float(np.min(go[:,:n]/volumes*1000*p.plasma_to_whole_blood_ratio))))
        states.append(rows);last_flux=(deficit,consumed)
    do=float(np.max(abs(states[0][:,:k]-states[1][:,:k])))
    dg=float(np.max(abs(states[0][:,k:]-states[1][:,k:])))
    if do>.002 or dg>1e-5:raise RuntimeError(f'Coupled refinement failed: {do}, {dg}')
    result=dict(status='verified',scope='coupled saved-hydraulic replay; no recoverable-duration inference',duration_s=float(time[-1]),case=case,parameters=vars(p),profile=entries,
        volume_names=h.volume_names,organ_names=[v.name for v in organs],records=records,
        paired_max_oxygen_state_difference=do,paired_max_glucose_state_difference_mmol=dg,
        inputs=paths+[seed_path,'reports/m4-baseline/native-cycle.npz','reports/m4-baseline/result.json'])
    write(OUT/case['name']/'result.json',result);return result



def hydraulics(case):
    path="reports/m6-coupled/observed_arrest/hydraulic-refined.npz"
    a=dict(np.load(ROOT/path));a['time']=a.pop('time_s');a['flows']=a.pop('flows_mL_s')
    if case.get('apnea') is None:return [a],[path]
    cut=case['apnea'];segments=[]
    for mask in (a['time']<=cut,a['time']>=cut):
        segments.append(dict(time=a['time'][mask],states=a['states'][:,mask],flows=a['flows'][mask]))
    return segments,[path]


def run(case):
    import time
    from cardiac_redundancy.myocardial_demand import side_demand,DemandAssumptions
    start=time.perf_counter()
    try:
        h=read('reports/m6-coupled/observed_arrest/hydraulics.json')
        normal=read('reports/m5-normal-comparison/verification.json')['runs'][-1]['metrics']['cycle']
        reference=read('configs/m5-myocardial-demand.json')['reference']
        for key in ('mass_ratio','work_cost_multiplier'):reference[key]=case['knobs'][key]
        native=read('data/m1-demand-budget.json')['baseline_mL_min']['heart']
        weights=read('configs/m4-organ-mapping.json')['mapping']['heart']
        sides={s:side_demand(native*weights['coronary_'+s],h['runs'][-1]['positive_transmural_mean_W'][j],normal['sides'][s]['positive_transmural_mean_W'],80,1,DemandAssumptions(**reference)) for j,s in enumerate(('left','right'))}
        case=dict(case,device='BIO',architecture='BIO',liver_function_fraction=.5,sides=sides)
        pars=read('configs/m4-fuel-reference.json')['parameters']
        for key in ('glycogen_mmol_L_liver','release_cap_mmol_min'):pars[key]=case['knobs'][key]
        result=short(case,FuelParameters(**pars))
        result['runtime_s']=time.perf_counter()-start
        write(OUT/case['name']/'result.json',result)
        print(case['name'],result['records'][-1]['unmet_oxygen_mL'],flush=True)
        return result
    except Exception as exc:
        write(OUT/case['name']/'result.json',dict(status='error',case=case,error=repr(exc),runtime_s=time.perf_counter()-start))
        raise


def cases():
    config=read('configs/m7-metabolic-screen.json')
    nominal={f['id']:f['nominal'] for f in config['factors']}
    result=[dict(name='nominal',knobs=nominal)]
    for f in config['factors']:
        if not f['low'] < f['nominal'] < f['high']:raise ValueError('Invalid bounds')
        for end in ('low','high'):
            result.append(dict(name=f['id']+'_'+end,knobs=dict(nominal,**{f['id']:f[end]})))
    return result


if __name__=='__main__':
    from concurrent.futures import ProcessPoolExecutor
    selected=cases()
    if len(sys.argv)>1:selected=[c for c in selected if c['name'] in sys.argv[1:]]
    with ProcessPoolExecutor(max_workers=3) as pool:list(pool.map(run,selected))
