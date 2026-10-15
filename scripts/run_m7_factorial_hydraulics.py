"""M7 paired hydraulic factorial runner derived from frozen M6 mechanics.

Only explicit conductance and contractility controls vary; initial volume
inventory, sensor history, observer and horizon are kept identical.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.sensitivity_circulation import SensitivityCirculation
from cardiac_redundancy.biological_pump import BiologicalParameters
from cardiac_redundancy.support_controller import SupportController,ControllerParameters
from cardiac_redundancy.controller_sensor import TrailingObservation
from cardiac_redundancy.circuit_integration import integrate_circuit
from cardiac_redundancy.biological_energy import biological_power
from cardiac_redundancy.pump_energy import MMHG_ML_TO_J
OUT=ROOT/'reports/m7-factorial'
def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))
def write(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def run(case):
    cfg=read('configs/m6-coupled-protocol.json');cc=read('configs/m6-controller.json')
    cfg['sample_period_s']=case.get('sample_period_s',cfg['sample_period_s'])
    dest=OUT/case['name'];write(dest/'status.json',dict(status='running'))
    records=[];arrays=[]
    for label,tol,step in [('nominal',1e-7,.01),('refined',1e-9,.005)]:
        m=SensitivityCirculation(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml',arterial_fraction=.2,inlet_resistance_fraction=.1,coronary_conductance=case['knobs']['coronary_conductance'],biological=BiologicalParameters(**{k:case['knobs']['contractility'] for k in ('left_contractility','right_contractility','septal_contractility')}))
        y=np.load(ROOT/'reports/m5-joint-preservation/reference_g25/fine-refined.npz')['states'][:,0].copy()
        m.apply_support_controls(0,case['knobs']['preservation_conductance'])
        if case.get('deplete_PU') or case.get('domain_boundary'):
            delta=m.domain_margins(y)['PU']-(0 if case.get('domain_boundary') else 1.)
            y[5]-=delta;y[3]+=delta
        if case.get('inlet_multiplier'):m.valve_R[4]*=case['inlet_multiplier']
        controller=SupportController(ControllerParameters(**cc['parameters']))
        sensor=TrailingObservation(cfg['filter_window_s']);rng=np.random.default_rng(cfg['seed'])
        # Declared synthetic observer priming; no physical inventory is added.
        # Every policy receives the same initial pressure/history convention.
        initial_pressure=m.pressures(0,y)[0]['AO']
        for t in np.round(np.arange(-1,0,cfg['sample_period_s']),8):sensor.update(float(t),initial_pressure,5.213513)
        ts=[];ys=[];ps=[];qs=[];amps=[];powers=[];events=[];termination=None;root=0.;jumps=0.
        release=.5 if case['policy']=='fixed' else None
        if release is not None:m.aux_envelope.set(release,1,1)
        ticks=np.round(np.arange(0,cfg['duration_s']+.0001,cfg['sample_period_s']),8)
        regional_names=list(m.native.regional_factors)
        for tick,end in zip(ticks[:-1],ticks[1:]):
            tick=float(tick);end=float(end)
            for key,target,ramp in [('restart',1,1),('rearrest',0,0)]:
                if case.get(key)==tick:m.native_envelope.set(tick,target,ramp)
            if case.get('withdraw')==tick:m.aux_envelope.set(tick,0,case['withdrawal_s'])
            if case.get('reopen')==tick:
                # Reopen systemic regional paths; retain declared one-lung
                # allocation so pulmonary capacity is not silently restored.
                names=[n for n in regional_names if n not in ('lungs','pulmonary_exchange')]
                m.native.set_regional_controls(time=tick,ramp_s=2,factors=dict.fromkeys(names,1.))
            p,_=m.pressures(tick,y);obs_p=p['AO'];obs_q=y[19]*.06
            if case.get('noise'):obs_p+=rng.normal(0,3);obs_q+=rng.normal(0,.1)
            if case.get('stuck_high'):obs_p=85.
            obs=sensor.update(tick,obs_p,obs_q)
            if case['policy']=='observed':
                before=len(controller.events)
                controller.update(tick,*obs,sensor_valid=not ('sensor_loss' in case and tick>=case['sensor_loss']),
                    model_valid=min(m.domain_margins(y).values())>0,
                    adequate_for_recovery=False)
                # Recovery is not automatically certified from pressure/flow.
                # Scripted handovers are separate scenarios with replay checks.
                for event in controller.events[before:]:
                    if event['after']=='ACTIVATE':
                        release=controller.release_s
                        m.aux_envelope.set(controller.ramp_start,1,controller.p.ramp_s)
                    if event['after']=='FAULT':
                        # Cancel any pending future activation before holding.
                        held=m.aux_envelope.value(tick)
                        m.aux_envelope.events=[e for e in m.aux_envelope.events if e[0]<=tick]
                        m.aux_envelope.set(tick,held,0);release=None
            bounds=sorted(set([tick,end]+([release] if release is not None and tick<=release<=end else [])))
            for start,stop in zip(bounds[:-1],bounds[1:]):
                if release is not None and m.trigger is None and abs(start-release)<1e-9:m.release_auxiliary(start)
                tr=integrate_circuit(m,stop-start,start=start,initial_state=y,rtol=tol,max_step=step,sample_step=.005)
                jumps=max(jumps,float(np.max(abs(tr.states[:,0]-y))))
                for ti,yi in zip(tr.time,tr.states.T):
                    if ts and ti<=ts[-1]+1e-12:continue
                    p,meta=m.pressures(ti,yi);q=m.flows(ti,yi,[True]*8,pressures=p)
                    power=biological_power(m,ti,yi)
                    ts.append(ti);ys.append(yi);ps.append([p[n] for n in m.volume_names]);qs.append([q[n] for n in m.flow_names]);amps.append(m.aux_envelope.value(ti))
                    powers.append([[power[s][k] for k in ('chamber_W','transmural_W','balance_residual_W')] for s in ('left','right')])
                    root=max(root,abs(meta['septal_residual_mmHg']),abs(meta['aux_septal_residual_mmHg']))
                events+=tr.events;y=tr.states[:,-1].copy()
                if tr.termination:termination=tr.termination;break
            if termination:break
        t=np.array(ts);v=np.array(ys).T;p=np.array(ps);q=np.array(qs);power=np.array(powers)
        drift=float(np.max(abs(v[:18].sum(axis=0)-5250)))
        assert drift<1e-5 and jumps<1e-10 and root<1e-8 and v[18:].min()>-1e-7
        assert abs(power[:,:,2]).max()<1e-10
        work=np.trapezoid(power[:,:,0],t,axis=0).tolist()
        pv=[float(-MMHG_ML_TO_J*np.sum(.5*(p[1:,j]+p[:-1,j])*np.diff(v[j]))) for j in (16,17)]
        # 5 ms saved quadrature is independently checked, not solver error.
        pv_error=max(abs(a-b) for a,b in zip(work,pv));assert pv_error<.05
        mask=t>=max(0,t[-1]-2);win=t[mask]
        mean=lambda z:float(np.trapezoid(z[mask],win)/(win[-1]-win[0])) if len(win)>1 else None
        record=dict(label=label,termination=termination,release_s=m.trigger,controller_events=controller.events,
            blood_error_mL=drift,state_jump=jumps,root_residual_mmHg=root,PV_quadrature_error_J=pv_error,
            final_MAP_mmHg=mean(p[:,2]),aux_left_L_min=mean(q[:,m.flow_names.index('aAV')]*.06),
            chamber_work_J=work,positive_transmural_mean_W=np.trapezoid(np.maximum(power[:,:,1],0),t,axis=0).tolist(),
            final_controller_state=controller.state if case['policy']=='observed' else 'prescribed',events=events)
        if len(t)>1:record['positive_transmural_mean_W']=[w/(t[-1]-t[0]) for w in record['positive_transmural_mean_W']]
        np.savez_compressed(dest/f'hydraulic-{label}.npz',time_s=t,states=v,pressures_mmHg=p,flows_mL_s=q,aux_amplitude=np.array(amps))
        records.append(record);arrays.append((t,v,p,q))
        print(case['name'],label,'hydraulics done',flush=True)
    t,a,p,q=arrays[0];tb,b,pb,qb=arrays[1]
    errors=[]
    for x,z in ((a.T,b.T),(p,pb),(q,qb)):
        interpolated=np.array([np.interp(t,tb,col) for col in z.T]).T
        errors.append(float(np.max(abs(x-interpolated)/np.maximum(1,np.max(abs(x),axis=0)))))
    assert max(errors)<.005
    assert records[0]['release_s']==records[1]['release_s']
    assert [(e['after'],e['time_s']) for e in records[0]['controller_events']]==[(e['after'],e['time_s']) for e in records[1]['controller_events']]
    result=dict(status='verified',case=case,runs=records,scaled_errors=errors,volume_names=m.volume_names,flow_names=m.flow_names)
    write(dest/'hydraulics.json',result);write(dest/'status.json',dict(status='hydraulics_verified'))
    return result

