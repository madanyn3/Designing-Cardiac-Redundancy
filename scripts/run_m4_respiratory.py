"""Combined T4.4/T4.5 replay experiments with baseline-derived oxygen inventories."""
from concurrent.futures import ProcessPoolExecutor,as_completed
from dataclasses import replace
import hashlib,json,sys,traceback
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.circuit import CircuitModel
from cardiac_redundancy.organ_oxygen import map_demands
from cardiac_redundancy.oxygen import OxygenParameters
from cardiac_redundancy.respiratory import RespiratoryOxygenTransport,RespiratoryParameters
from cardiac_redundancy.oxygen_numerics import bounded_jacobian,periodic_rk4
from cardiac_redundancy.constraints import oxygen_screen,combined_status,pressure_windows
from cardiac_redundancy.preservation import PreservationChoice,apply_policy
OUT=ROOT/'reports/m4-respiratory'


def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))
def write(p,data):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def load_segments(case):
    if case['hydraulic']=='periodic':
        z=dict(np.load(ROOT/'reports/m4-baseline/native-cycle.npz'))
        return [dict(time=z['time_s'],states=z['states'],flows=z['flows'],pressures=z['pressures'])],['reports/m4-baseline/native-cycle.npz']
    if case['hydraulic']=='whole':directory='reports/m3-whole-network/preservation_quarter-nominal';indices=range(4)
    else:
        directory=f"reports/m3-circuits/{case['hydraulic']}-nominal"
        result=read(directory+'/result.json')
        if result['status']!='completed':raise ValueError('Hydraulic case incomplete')
        indices=[s['index'] for s in result['segments']]
    paths=[f'{directory}/segment-{i}.npz' for i in indices]
    return [dict(np.load(ROOT/p)) for p in paths],paths


def run(case):
    name=case['name'];dest=OUT/name
    try:
        protocol=read('configs/m4-respiratory-protocol.json')
        h=CircuitModel(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml',architecture=case['architecture'],arterial_fraction=.2,inlet_resistance_fraction=.1)
        organs=map_demands(read('data/m1-demand-budget.json')['baseline_mL_min'],read('configs/m4-organ-mapping.json')['mapping'],h.flow_names,.75,case.get('auxiliary_total_mL_min'))
        parameters=OxygenParameters(**read('configs/m4-oxygen-reference.json')['parameters'])
        parameters=replace(parameters,exchange_fraction=case.get('exchange_fraction',1.))
        rp=RespiratoryParameters(reference_ventilation_mL_s_BTPS=protocol['respiratory']['reference_ventilation_mL_min_BTPS']/60,
            response_time_s=protocol['respiratory']['response_time_s'],
            oxygenator_gas_volume_mL_STPD=protocol['oxygenator']['gas_volume_mL_STPD'],
            oxygenator_transfer_limit_mL_s=protocol['oxygenator']['transfer_limit_mL_s'])
        m=RespiratoryOxygenTransport(h.volume_names,h.edges,organs,parameters,respiratory=rp,store_seconds=case['store_seconds'])
        n=len(h.volume_names);b=m.base_size;count=len(organs)
        segments,hydraulic_paths=load_segments(case)
        periodic=case['hydraulic']=='periodic'
        first=segments[0];origin=float(first['time'][0]);initial_volumes=first['states'][:n,0]
        baseline=read('reports/m4-baseline/result.json');base_names=baseline['volume_names'];base_organs=baseline['organ_names']
        if case['hydraulic'] in ('periodic','whole'):
            expected=baseline
        else:expected=read(f"reports/m3-circuits/{case['hydraulic']}-nominal/result.json")
        if list(h.flow_names)!=expected['flow_names'] or list(h.volume_names)!=expected['volume_names']:
            raise ValueError('Saved hydraulic ordering does not match oxygen graph')
        z=np.load(ROOT/f"reports/m4-baseline/store-{case['store_seconds']}-refined.npz")
        native=np.load(ROOT/'reports/m4-baseline/native-cycle.npz');phase=0. if periodic else .375
        bz=interp1d(z['time_s'],z['oxygen_states'],axis=0)(phase)
        bv=interp1d(native['time_s'],native['states'][:len(base_names)],axis=1)(phase)
        y0=m.initial_state(initial_volumes)
        for j,node in enumerate(h.volume_names):
            reference='PU' if node=='aLV' else 'VC' if node=='aRV' else node
            k=base_names.index(reference);y0[j]=initial_volumes[j]*bz[k]/bv[k]
        y0[n]=bz[len(base_names)]
        base_size=len(base_names)+2+2*len(base_organs)
        for j,o in enumerate(organs):
            if o.name in base_organs:y0[b+3+j]=bz[base_size+3+base_organs.index(o.name)]
        controls={key:case[key] for key in ('mode','ventilation_fraction','drive_fraction','sweep_mL_s_STPD') if key in case}
        required_choices={o.name:PreservationChoice() for o in organs}
        optional_choices={o.name:PreservationChoice(outcome_required=not o.parent_pool.startswith('limb_')) for o in organs}
        required=apply_policy(organs,required_choices);optional=apply_policy(organs,optional_choices)
        saved={};records=[]
        for setting in protocol['solvers']:
            setting=dict(setting)
            setting['rtol']*=case.get('solver_rtol_scale',1.)
            setting['max_step_s']=min(setting['max_step_s'],case.get('solver_max_step_s',setting['max_step_s']))
            y=y0.copy();ts=[];ys=[];vls=[];qls=[];pls=[]
            duration=120. if periodic else float(segments[-1]['time'][-1]-origin)
            for segment in segments:
                tt=segment['time'];vf=interp1d(tt,segment['states'][:n],axis=1);qf=interp1d(tt,segment['flows'],axis=0);pf=interp1d(tt,segment['pressures'],axis=0)
                if periodic:
                    start,end=0.,duration
                    v=lambda t:vf(t%.75);q=lambda t:qf(t%.75);p=lambda t:pf(t%.75)
                    sample=np.linspace(0,duration,1201)
                else:
                    start,end=float(tt[0]),float(tt[-1]);v=vf;q=qf;p=pf;sample=tt
                rhs=lambda t,x:m.rhs(v(t),x,q(t),**controls)
                if periodic:
                    step=protocol['periodic_solver'][setting['label']+'_step_s']
                    per=round(.75/step);fg=np.linspace(0,.75,2*per+1)
                    vgrid=vf(fg).T;qgrid=qf(fg)
                    time,rows=periodic_rk4(lambda k,x:m.rhs(vgrid[k],x,qgrid[k],**controls),
                        y,step,per,round(duration/.75),round(.1/step))
                    y=rows[-1].copy();ts.extend(time);ys.extend(rows)
                    sample=time
                else:
                    sol=solve_ivp(rhs,(start,end),y,method='Radau',jac=bounded_jacobian(rhs),rtol=setting['rtol'],atol=setting['rtol']*.01,max_step=setting['max_step_s'],t_eval=sample)
                    if not sol.success:raise RuntimeError(sol.message)
                    y=sol.y[:,-1].copy();ts.extend(sol.t-origin);ys.extend(sol.y.T)
                vls.extend([v(ti) for ti in sample]);qls.extend([q(ti) for ti in sample]);pls.extend([p(ti) for ti in sample])
            times=np.array(ts);states=np.array(ys);volumes=np.array(vls);flows=np.array(qls);pressures=np.array(pls)
            unique,idx=np.unique(times,return_index=True);times=unique;states=states[idx];volumes=volumes[idx];flows=flows[idx];pressures=pressures[idx]
            diags=[m.evaluate(vi,yi,qi,**controls)[1] for vi,yi,qi in zip(volumes,states,flows)]
            flux={k:np.array([d[k] for d in diags]) for k in ('delivered','extractable','consumed','demand','unmet','blood_extracted')}
            ventilation=np.array([d['ventilation_mL_s_BTPS'] for d in diags]);external=np.array([d['oxygenator_external_mL_s'] for d in diags])
            balance=float(np.max(abs(m.balance(states)-m.balance(y0))))
            demand_error=float(np.max(abs(states[:,n+2:n+2+count]+states[:,n+2+count:b]-times[:,None]*m.reference_demand/60)))
            if balance>1e-6 or demand_error>1e-6 or np.min(states[:,:n+1])<0 or states[:,b+1].min() < -1e-7 or states[:,b+1].max()>rp.oxygenator_gas_volume_mL_STPD+1e-7 or states[:,b+3:].min() < -1e-7 or np.max(states[:,b+3:]-m.store_capacity)>1e-7:
                raise RuntimeError(f'Accounting/domain failure balance={balance}, demand={demand_error}, blood/gas_min={states[:,:n+1].min()}, oxygenator_range={states[:,b+1].min(),states[:,b+1].max()}, tissue_min={states[:,b+3:].min()}, tissue_excess={np.max(states[:,b+3:]-m.store_capacity)}')
            dest.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(dest/f"{setting['label']}.npz",time_s=times,oxygen_states=states,volumes_mL=volumes,
                pressures_mmHg=pressures,ventilation_mL_s_BTPS=ventilation,**flux)
            saved[setting['label']]=(states,flux,times)
            screens=[oxygen_screen(times,flux['demand'][:,j],flux['consumed'][:,j],flux['delivered'][:,j],states[:,n+2+count+j],o.demand_mL_min,1e-6,blood_extraction=flux['blood_extracted'][:,j]) for j,o in enumerate(organs)]
            th=dict(mean=65,maximum=90,minimum=50);pressure=pressure_windows(times,pressures[:,h.volume_names.index('AO')],th)
            unknown=read('configs/m4-constraint-policy.json')['unresolved_requirements']
            assessments={key:combined_status('valid',pressure['status'],[s['status'] for s,flag in zip(screens,mask) if flag],unknown)
                for key,mask in [('all_tissues',required['outcome_required']),('pooled_limbs_optional',optional['outcome_required'])]}
            actual_solver=dict(method='fixed RK4',step_s=step,label=setting['label']) if periodic else dict(setting,method='Radau')
            r=dict(solver=actual_solver,oxygen_conservation_error_mL=balance,demand_accounting_error_mL=demand_error,
                total_unmet_mL=float(states[-1,n+2+count:b].sum()),total_consumed_mL=float(states[-1,n+2:n+2+count].sum()),
                final_native_gas_oxygen_mL=float(states[-1,n]),native_net_external_oxygen_mL=float(states[-1,n+1]),
                initial_total_physical_oxygen_mL=float(m.balance(y0)),final_oxygenator_oxygen_mL=float(states[-1,b+1]),
                oxygenator_net_external_oxygen_mL=float(states[-1,b+2]),final_tissue_store_mL=float(states[-1,b+3:].sum()),
                final_AO_content_mL_dL=float(states[-1,h.volume_names.index('AO')]/volumes[-1,h.volume_names.index('AO')]*100),
                final_ventilation_L_min_BTPS=float(ventilation[-1]*.06),organ_screens=screens,pressure=pressure,assessments=assessments)
            records.append(r)
        delta=abs(saved['nominal'][0]-saved['refined'][0])
        difference=float(np.max(np.delete(delta,b,axis=1)))
        ventilation_difference=float(np.max(delta[:,b]))
        rate_error=np.max(abs(saved['nominal'][1]['unmet']-saved['refined'][1]['unmet']),axis=0)
        if difference>.002 or ventilation_difference>.002:raise RuntimeError(f'Refinement exceeds tolerance: oxygen={difference} mL, ventilation={ventilation_difference} mL BTPS/s')
        # Reassess supply against paired numerical uncertainty, never physiological grace.
        for r in records:
            state,flux,time=saved[r['solver']['label']]
            r['organ_screens']=[oxygen_screen(time,flux['demand'][:,j],flux['consumed'][:,j],flux['delivered'][:,j],state[:,n+2+count+j],o.demand_mL_min,max(1e-8,10*rate_error[j]),blood_extraction=flux['blood_extracted'][:,j]) for j,o in enumerate(organs)]
            r['assessments']={key:combined_status('valid',r['pressure']['status'],[s['status'] for s,flag in zip(r['organ_screens'],mask) if flag],unknown)
                for key,mask in [('all_tissues',required['outcome_required']),('pooled_limbs_optional',optional['outcome_required'])]}
        result=dict(status='verified',case=case,hydraulic_inputs=hydraulic_paths,volume_names=h.volume_names,organ_names=[o.name for o in organs],
            store_capacity_mL=m.store_capacity.tolist(),records=records,refinement_difference_mL=difference,
            refinement_ventilation_difference_mL_s_BTPS=ventilation_difference,
            physiological_validation='unavailable; response, stores and membrane parameters are assumptions',
            initialization='periodic B1 oxygen contents; BIO cavity contents transferred from PU/VC and auxiliary stores assumed full; no BIO metabolic warm-up')
    except Exception as e:result=dict(status='failed',case=case,error=str(e),traceback=traceback.format_exc())
    write(dest/'result.json',result);return name,result['status'],result.get('error','')


def main():
    protocol=read('configs/m4-respiratory-protocol.json')
    if '--manifest-only' not in sys.argv:
        with ProcessPoolExecutor(max_workers=3) as pool:
            for future in as_completed([pool.submit(run,c) for c in protocol['cases']]):
                print(*future.result(),flush=True)
    results=[]
    for c in protocol['cases']:
        r=read(f"reports/m4-respiratory/{c['name']}/result.json")
        if r['case']!=c:raise RuntimeError('Saved case differs from current protocol')
        results.append((c['name'],r['status'],r.get('error','')))
    paths=[Path(__file__),*sorted((ROOT/'src/cardiac_redundancy').glob('*.py')),
        *[ROOT/p for p in ('data/m1-demand-budget.json','configs/baseline.yaml','configs/m3-evidence-regional.yaml','configs/m2-frozen-overlay.yaml')],
        *sorted((ROOT/'configs').glob('m4-*.json')),
        ROOT/'reports/m4-baseline/native-cycle.npz',
        *[ROOT/f'reports/m4-baseline/store-{k}-refined.npz' for k in (0,5,15)],ROOT/'reports/m4-baseline/result.json']
    for c in protocol['cases']:
        for p in load_segments(c)[1]:paths.append(ROOT/p)
    write(OUT/'manifest.json',dict(status='verified' if all(r[1]=='verified' for r in results) else 'incomplete',results=results,
        software=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__),
        reconstruction='Manifest can consolidate individually rerun cases; each saved case must match the protocol exactly.',
        output_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for c in protocol['cases'] for p in sorted((OUT/c['name']).glob('*')) if p.is_file()},
        hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}))
    return int(any(r[1]!='verified' for r in results))


if __name__=='__main__':raise SystemExit(main())
