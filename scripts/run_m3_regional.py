"""Run the step-2 regional comparator with independent warm-up per case."""
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import traceback

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
import scipy
from cardiac_redundancy.regional import RegionalModel
from cardiac_redundancy.pressure_regional import PressureRegionalModel
from cardiac_redundancy.integration import integrate_source
from cardiac_redundancy.experiments import periodic_error


def write(path,data):
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def checks(model,tr):
    if tr.termination is not None:
        raise RuntimeError('Model domain exit: '+str(tr.termination))
    drift=float(np.max(abs(tr.states[:11].sum(axis=0)-5250))/5250)
    minimum=float((tr.states[:11]-model.inventory.unstressed_volumes[:,None]).min())
    minq=float(tr.states[11:].min())
    if drift>=1e-8 or minimum<0 or minq < -model.settings['absolute_flow_tolerance_mL_s']:
        raise RuntimeError(f'Regional checks failed: {drift}, {minimum}, {minq}')
    return {'relative_volume_drift':drift,'minimum_domain_margin_mL':minimum,'minimum_valve_flow_mL_s':minq}


def run_case(case,protocol,out):
    out=Path(out)/case['name']; out.mkdir(parents=True,exist_ok=True)
    result={'case':case,'status':'running'}; write(out/'result.json',result)
    try:
        paths=(ROOT/'configs/baseline.yaml',ROOT/protocol.get('project_config','configs/project-extension.yaml'),ROOT/'configs/m2-frozen-overlay.yaml')
        model=(PressureRegionalModel(*paths,closure=case.get('closure',{}))
               if protocol.get('pressure_aware',False) else RegionalModel(*paths))
        transfer=case.get('initial_VC_to_PU_mL',0)
        model.initial_state[3]-=transfer; model.initial_state[5]+=transfer
        write(out/'effective-model.json',{'variant':protocol['variant'],
            'volume_names':model.volume_names,'initial_state':model.initial_state.tolist(),
            'branches':[vars(b) for b in model.branches],
            'closure':getattr(model,'closure_settings',None)})
        y=model.initial_state.copy(); previous=None; streak=0; history=[]
        period=model.config['driver']['period_s']
        for cycle in range(1,model.settings['warmup_max_cycles']+1):
            tr=integrate_source(model,period,initial_state=y,start=(cycle-1)*period,
                sample_step=period/protocol['warmup_phase_intervals'],**case['solver'])
            check=checks(model,tr)
            error=None if previous is None else periodic_error(tr.states,previous)
            streak=streak+1 if error is not None and error<model.settings['periodic_state_relative_tolerance'] else 0
            history.append({'cycle':cycle,'phase_error':error,'passing_streak':streak,**check})
            write(out/'warmup.json',history)
            y=tr.states[:,-1].copy(); previous=tr.states
            if streak>=protocol['consecutive_periodic_comparisons']: break
        else: raise RuntimeError('Regional periodic warm-up did not converge')
        tr=integrate_source(model,period*protocol['measurement_cycles'],initial_state=y,start=cycle*period,
            sample_step=case.get('sample_step_s',protocol['sample_step_s']),**case['solver'])
        check=checks(model,tr)
        ps=[]; qs=[]; root=[]
        for t,y in zip(tr.time,tr.states.T):
            p,m=model.pressures(t,y)
            # Saved valve states already obey event semantics; use them directly.
            q=model.flows(t,y,[True]*4,pressures=p)
            ps.append([p[n] for n in model.volume_names+('H',)])
            qs.append(list(q.values())); root.append(abs(m['septal_residual_mmHg']))
        ps=np.array(ps); qs=np.array(qs)
        average=lambda a:float(np.trapezoid(a,tr.time)/(tr.time[-1]-tr.time[0]))
        meanflows={n:average(qs[:,i]) for i,n in enumerate(model.flow_names)}
        ao=ps[:,2]
        metrics={'MAP_mmHg':average(ao),'SBP_mmHg':float(ao.max()),'DBP_mmHg':float(ao.min()),
            'CO_L_min':meanflows['AV']*.06,'RV_output_L_min':meanflows['PV']*.06,
            'PA_systolic_mmHg':float(ps[:,4].max()),'PA_diastolic_mmHg':float(ps[:,4].min()),
            'LV_SV_mL':float(np.ptp(tr.states[0,-int(round(period/(tr.time[1]-tr.time[0])))-1:])),
            'LV_EDV_mL':float(tr.states[0].max()),'LV_ESV_mL':float(tr.states[0].min())}
        balance=meanflows['AV']-meanflows['PV']-meanflows['lungs']
        hepatic=float(np.max(abs(sum((1 if b.recipient=='H' else -1)*qs[:,model.flow_names.index(b.name)] for b in model.hepatic_branches))))
        if max(root)>=model.settings['septal_residual_tolerance_mmHg'] or hepatic>1e-9:
            raise RuntimeError('Regional algebraic residual failure')
        result.update(status='passed_numerical_checks',warmup_cycles=cycle,metrics=metrics,
            mean_branch_flows_mL_s=meanflows,checks=check,maximum_septal_residual_mmHg=max(root),
            maximum_hepatic_balance_residual_mL_s=hepatic,
            LV_minus_RV_minus_lung_supply_mL_s=balance,
            measurement_start_s=float(tr.time[0]))
        np.savez_compressed(out/'trajectory.npz',time_s=tr.time,states=tr.states,pressures_mmHg=ps,flows_mL_s=qs)
        write(out/'events.json',tr.events)
    except Exception as exc:
        result.update(status='failed',error=str(exc),traceback=traceback.format_exc())
    write(out/'result.json',result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--out',default='reports/m3-regional')
    parser.add_argument('--protocol',default='configs/m3-regional-protocol.json')
    args=parser.parse_args()
    protocol=json.loads((ROOT/args.protocol).read_text())
    out=ROOT/args.out; out.mkdir(parents=True,exist_ok=True)
    paths=[Path(__file__),ROOT/args.protocol,ROOT/'configs/baseline.yaml',
           ROOT/protocol.get('project_config','configs/project-extension.yaml'),ROOT/'configs/m2-frozen-overlay.yaml',
           *sorted((ROOT/'src/cardiac_redundancy').glob('*.py'))]
    manifest={'protocol':protocol,'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
        'hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    for p in paths:
        dest=out/'implementation-snapshot'/p.relative_to(ROOT); dest.parent.mkdir(parents=True,exist_ok=True); dest.write_bytes(p.read_bytes())
    write(out/'manifest.json',manifest)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run_case,c,protocol,str(out)) for c in protocol['cases']]
        results=[]
        for future in as_completed(futures):
            r=future.result(); results.append(r)
            print(r['case']['name'],r['status'],r.get('warmup_cycles'),r.get('error',''),flush=True)
    return int(any(r['status']=='failed' for r in results))


if __name__=='__main__': raise SystemExit(main())
