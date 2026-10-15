"""Paired dynamic B1 stage experiments; keep hydraulic and metabolic claims apart."""
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import hashlib,json,sys,traceback
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.backup_pump import BackupPumpModel,PumpParameters,PumpStage
from cardiac_redundancy.bounded_pump import PreloadPumpModel,PreloadGate
from cardiac_redundancy.integration import integrate_source
from cardiac_redundancy.constraints import pressure_windows
OUT=ROOT/'reports/m5-dynamic'
SEED=ROOT/'reports/m3-whole-network/preservation_quarter-nominal/segment-0.npz'


def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))
def write(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def run(case,setting):
    label=setting['label'];dest=OUT/(case['name']+'-'+label)
    write(dest/'result.json',dict(status='running',case=case,solver=setting))
    try:
        protocol=read('configs/m5-dynamic-protocol.json');args=(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml')
        kwargs=dict(arterial_fraction=.2,inlet_resistance_fraction=.1,
            pump_parameters=PumpParameters(PumpStage(.1,case['right_capacity']),PumpStage(.1,case['left_capacity'])))
        if case['gated']:
            model=PreloadPumpModel(*args,**kwargs,right_gate=PreloadGate(),left_gate=PreloadGate(external_pressure_mmHg=case['left_external_mmHg']))
        else:model=BackupPumpModel(*args,**kwargs)
        seed=np.load(SEED);y=seed['states'][:,0].copy();origin=float(seed['time'][0]);n=len(model.volume_names)
        expected=read('reports/m3-whole-network/preservation_quarter-nominal/result.json')
        if list(model.volume_names)!=expected['volume_names'] or list(model.flow_names)!=expected['flow_names']:raise ValueError('Saved seed ordering mismatch')
        transfer=0.
        if 'PU_margin_mL' in case:
            j=model.volume_names.index('PU');v=model.volume_names.index('VC')
            transfer=float(y[j]-model.U[j]-case['PU_margin_mL']);y[j]-=transfer;y[v]+=transfer
        if min(model.domain_margins(y).values())<=0:raise ValueError('Invalid initial domain')
        model.set_controls(time=origin,ramp_s=protocol['activation_ramp_s'],pump_on=True,
            right_head=protocol['heads_mmHg']['right'],left_head=protocol['heads_mmHg']['left'])
        roles=read('configs/m3-whole-network-protocol.json')
        model.set_regional_controls(time=origin,factors=dict.fromkeys(roles['preservation_supply_edges'],protocol['preservation_conductance_factor']))
        settings={k:v for k,v in setting.items() if k!='label'};segments=[];t=origin;blood_error=root_error=junction_error=pool_error=jump=0.;termination=None
        for i,end in enumerate((protocol['activation_ramp_s'],case.get('duration_s',protocol['duration_s']))):
            tr=integrate_source(model,origin+end-t,start=t,initial_state=y,active=False,sample_step=protocol['sample_step_s'],**settings)
            jump=max(jump,float(np.max(abs(tr.states[:,0]-y))))
            ps=[];qs=[];gates=[];net=[]
            for ti,yi in zip(tr.time,tr.states.T):
                p,meta=model.pressures(ti,yi,active=False)
                q=model.flows(ti,yi,[True]*4,active=False,pressures=p)
                ordered=np.array([q[name] for name in model.flow_names]);dv=model.incidence@ordered
                # Pulmonary vascular inventory includes PA and PU. Pulmonary
                # exchange cancels internally; bronchial/nutritive return and
                # residual native valve flows prohibit blindly setting Qr=Ql.
                pulmonary=q['support_right']+q['PV']+q['lungs']-q['support_left']-q['MT']
                pool_error=max(pool_error,abs(dv[model.volume_names.index('PA')]+dv[model.volume_names.index('PU')]-pulmonary))
                junction_error=max(junction_error,abs(dv[-1]));root_error=max(root_error,abs(meta['septal_residual_mmHg']))
                ps.append([p[name] for name in model.volume_names+('H',)]);qs.append(ordered);net.append(pulmonary)
                gates.append([model.right_gate.fraction(p['VC']),model.left_gate.fraction(p['PU'])] if case['gated'] else [1.,1.])
            blood_error=max(blood_error,float(np.max(abs(tr.states[:n].sum(axis=0)-5250))))
            np.savez_compressed(dest/f'segment-{i}.npz',time_s=tr.time-origin,states=tr.states,pressures_mmHg=ps,flows_mL_s=qs,preload_fraction=gates,pulmonary_net_mL_s=net)
            segments.append(dict(index=i,valves=tr.events));y=tr.states[:,-1];t=float(tr.time[-1]);termination=tr.termination
            if termination:break
        if blood_error>1e-5 or max(root_error,junction_error,pool_error,jump)>1e-8:raise RuntimeError('Dynamic invariant failure')
        result=dict(status='model_domain_limit' if termination else 'completed',case=case,solver=setting,
            arrest_s=origin,observed_s=t-origin,termination=termination,segments=segments,
            volume_names=model.volume_names,flow_names=model.flow_names,initial_PU_to_VC_transfer_mL=transfer,
            blood_conservation_error_mL=blood_error,septal_residual_mmHg=root_error,hepatic_residual_mL_s=junction_error,
            pulmonary_balance_identity_error_mL_s=pool_error,state_jump_mL=jump,
            metabolic_adequacy='not assessed in hydraulic characterization; M4 coupling required',recovery='unassessed')
    except Exception as e:result=dict(status='failed',case=case,solver=setting,error=str(e),traceback=traceback.format_exc())
    write(dest/'result.json',result);return case['name'],label,result['status'],result.get('error','')


def assemble():
    protocol=read('configs/m5-dynamic-protocol.json');comparisons={};summaries={}
    for case in protocol['cases']:
        runs=[read(f"reports/m5-dynamic/{case['name']}-{s['label']}/result.json") for s in protocol['solvers']]
        if any(r['status']=='failed' for r in runs):raise RuntimeError('Failed run: '+case['name'])
        for r,setting in zip(runs,protocol['solvers']):
            if r['case']!=case or r['solver']!=setting or r['status'] not in ('completed','model_domain_limit'):raise RuntimeError('Protocol mismatch')
        arrays=[[dict(np.load(OUT/f"{case['name']}-{r['solver']['label']}"/f"segment-{s['index']}.npz")) for s in r['segments']] for r in runs]
        errors={key:0. for key in ('states','pressures_mmHg','flows_mL_s')}
        if len(arrays[0])!=len(arrays[1]):raise RuntimeError('Segment mismatch')
        for a,b in zip(*arrays):
            mask=a['time_s']<=min(a['time_s'][-1],b['time_s'][-1])
            for key in errors:
                x=a[key] if key=='states' else a[key].T;z=b[key] if key=='states' else b[key].T
                z=np.array([np.interp(a['time_s'][mask],b['time_s'],row) for row in z])
                scale=np.maximum(1,np.max(abs(x[:,mask]),axis=1))
                errors[key]=max(errors[key],float(np.max(abs(x[:,mask]-z)/scale[:,None])))
        events=[[e for s in r['segments'] for e in s['valves']] for r in runs]
        sequence=[(e['valve'],e['action']) for e in events[0]]==[(e['valve'],e['action']) for e in events[1]]
        event_difference=max([abs(a['time_s']-b['time_s']) for a,b in zip(*events)],default=0.)
        domain=(runs[0]['termination'] or {}).get('node')==(runs[1]['termination'] or {}).get('node')
        duration_difference=abs(runs[0]['observed_s']-runs[1]['observed_s'])
        passed=max(errors.values())<.005 and sequence and event_difference<.001 and domain and duration_difference<.001
        comparisons[case['name']]=dict(passed=passed,scaled_errors=errors,matching_valve_sequence=sequence,
            valve_time_difference_s=event_difference,matching_domain=domain,stop_time_difference_s=duration_difference)
        r=runs[1];a=arrays[1];times=np.concatenate([z['time_s'] for z in a]);time,idx=np.unique(times,return_index=True)
        states=np.concatenate([z['states'].T for z in a])[idx];p=np.concatenate([z['pressures_mmHg'] for z in a])[idx];q=np.concatenate([z['flows_mL_s'] for z in a])[idx]
        fn=r['flow_names'];vn=r['volume_names'];last=time>=max(1.,time[-1]-2)
        mean=lambda x:None if sum(last)<2 else float(np.trapezoid(x[last],time[last])/(time[last][-1]-time[last][0]))
        pulmonary=states[:,vn.index('PA')]+states[:,vn.index('PU')]
        summary=dict(status=r['status'],observed_s=r['observed_s'],domain_node=(r['termination'] or {}).get('node'),
            final_AO_mmHg=float(p[-1,vn.index('AO')]),final_PU_mmHg=float(p[-1,vn.index('PU')]),final_VC_mmHg=float(p[-1,vn.index('VC')]),
            pulmonary_volume_change_mL=float(pulmonary[-1]-pulmonary[0]),
            tail_window_s=[float(time[last][0]),float(time[-1])] if sum(last)>1 else None,
            tail_mean_AO_mmHg=mean(p[:,vn.index('AO')]),tail_right_mL_s=mean(q[:,fn.index('support_right')]),tail_left_mL_s=mean(q[:,fn.index('support_left')]),
            regional_tail_flow_mL_min={name:None if mean(q[:,j]) is None else 60*mean(q[:,j]) for j,name in enumerate(fn)},
            pressure_screen=pressure_windows(time,p[:,vn.index('AO')],protocol['pressure_screen_mmHg']))
        summaries[case['name']]=summary
    paths=[Path(__file__),ROOT/'configs/m5-dynamic-protocol.json',ROOT/'configs/m3-whole-network-protocol.json',SEED,
        ROOT/'reports/m3-whole-network/preservation_quarter-nominal/result.json',
        *[ROOT/p for p in ('configs/baseline.yaml','configs/m3-evidence-regional.yaml','configs/m2-frozen-overlay.yaml')],
        *sorted((ROOT/'src/cardiac_redundancy').glob('*.py'))]
    result=dict(status='verified' if all(c['passed'] for c in comparisons.values()) else 'failed_refinement',protocol=protocol,
        comparisons=comparisons,scenarios=summaries,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        output_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.glob('*/*') if p.is_file()})
    write(OUT/'verification.json',result);return result['status']


def main():
    for p,h in read('data/m2-freeze-manifest.json')['sha256'].items():
        if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h:raise RuntimeError('Frozen source changed')
    protocol=read('configs/m5-dynamic-protocol.json')
    if '--report-only' not in sys.argv:
        with ProcessPoolExecutor(max_workers=3) as pool:
            for f in as_completed([pool.submit(run,c,s) for c in protocol['cases'] for s in protocol['solvers']]):print(*f.result(),flush=True)
    status=assemble();print(status);return int(status!='verified')


if __name__=='__main__':raise SystemExit(main())
