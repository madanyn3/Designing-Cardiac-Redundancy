"""T5.2 first increment: work from verified T5.1 BIO trajectories."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.biological_pump import BiologicalPumpModel,BiologicalParameters
from cardiac_redundancy.biological_energy import biological_power
from cardiac_redundancy.circuit_integration import integrate_circuit
from cardiac_redundancy.pump_energy import MMHG_ML_TO_J
OUT=ROOT/'reports/m5-biological-energy'
def main():
    OUT.mkdir(parents=True,exist_ok=True);results={};inputs=[]
    for name in ('reference','faster','stronger_left','asymmetric'):
        base=ROOT/'reports/m5-biological'/name;r=json.loads((base/'result.json').read_text());inputs.append(base/'result.json');runs=[]
        for label in ('nominal','refined'):
            m=BiologicalPumpModel(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml',biological=BiologicalParameters(**r['parameters']),arterial_fraction=.2,inlet_resistance_fraction=.1)
            # Refine output sampling as well as ODE tolerances: PV quadrature
            # across fast ejection needs a finer grid than T5.1's 5 ms plots.
            seedpath=base/'seed.npz';inputs.append(seedpath);seed=np.load(seedpath)
            tol,step,sample=(1e-7,.01,.001) if label=='nominal' else (1e-9,.005,.0005)
            m.trigger=None;m.valve_enabled[4:]=[False]*4
            pre=integrate_circuit(m,.375,start=float(seed['time_s']),initial_state=seed['state'],rtol=tol,max_step=step,sample_step=sample)
            if pre.termination:raise RuntimeError('Pre-arrest domain')
            arrest=float(pre.time[-1]);m.release_auxiliary(arrest)
            roles=json.loads((ROOT/'configs/m3-whole-network-protocol.json').read_text())
            m.native.set_regional_controls(time=arrest,factors=dict.fromkeys(roles['preservation_supply_edges'],.25))
            pieces=[];state=pre.states[:,-1];start=arrest
            for end in (1.,8.):
                tr=integrate_circuit(m,arrest+end-start,start=start,initial_state=state,active=False,rtol=tol,max_step=step,sample_step=sample)
                if tr.termination:raise RuntimeError('Support domain')
                pieces.append(tr);state=tr.states[:,-1];start=float(tr.time[-1])
            absolute,idx=np.unique(np.concatenate([tr.time for tr in pieces]),return_index=True)
            y=np.concatenate([tr.states for tr in pieces],axis=1)[:,idx];t=absolute-arrest
            assert np.max(abs(y[:18].sum(axis=0)-5250))<1e-5
            pressures=[]
            for ti,yi in zip(absolute,y.T):
                pressure=m.pressures(ti,yi,False)[0];pressures.append([pressure[node] for node in m.volume_names])
            pressures=np.array(pressures)
            records=[biological_power(m,ti,yi) for ti,yi in zip(absolute,y.T)];sides={}
            np.savez_compressed(OUT/f'{name}-{label}.npz',time_s=t,states=y,pressures_mmHg=pressures)
            for side,j in (('left',16),('right',17)):
                powers={key:np.array([x[side][key] for x in records]) for key in records[0][side]}
                work={key.replace('_W','_J'):float(np.trapezoid(v,t)) for key,v in powers.items() if key.endswith('_W')}
                delta=float(powers['valve_stored_J'][-1]-powers['valve_stored_J'][0])
                storage_error=abs(work['valve_storage_rate_J']-delta)
                coarse=np.arange(0,len(t),2);coarse=np.unique(np.r_[coarse,len(t)-1])
                quadrature=max(abs(work[k.replace('_W','_J')]-float(np.trapezoid(v[coarse],t[coarse]))) for k,v in powers.items() if k.endswith('_W'))
                pc=pressures[:,j];pv=float(-MMHG_ML_TO_J*np.sum(.5*(pc[1:]+pc[:-1])*np.diff(y[j])))
                pv_error=abs(pv-work['chamber_J'])
                if max(storage_error,quadrature,pv_error)>.005 or abs(powers['balance_residual_W']).max()>1e-10:raise RuntimeError((name,side,'energy check',storage_error,quadrature,pv_error))
                sides[side]=dict(work=work,valve_storage_change_J=delta,storage_integral_error_J=storage_error,coarse_grid_difference_J=quadrature,PV_work_J=pv,PV_quadrature_difference_J=pv_error,maximum_balance_residual_W=float(abs(powers['balance_residual_W']).max()))
            runs.append(dict(label=label,rtol=tol,max_step_s=step,sample_step_s=sample,sides=sides))
            print(name,label,'verified',flush=True)
        difference=max(abs(runs[0]['sides'][s]['work'][k]-runs[1]['sides'][s]['work'][k]) for s in ('left','right') for k in runs[0]['sides'][s]['work'])
        assert difference<.005
        results[name]=dict(status='verified',duration_s=float(t[-1]),runs=runs,paired_work_difference_J=difference)
    files=inputs+[Path(__file__),ROOT/'src/cardiac_redundancy/biological_energy.py',ROOT/'src/cardiac_redundancy/biological_pump.py']
    result=dict(status='verified',scope='mechanical ledger only; metabolic calibration remains open',absolute_energy_tolerance_J=.005,results=results,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    (OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    rows=['# T5.2 BIO mechanical-energy first increment','','Eight-second transient work; signed chamber work, not steady-cycle PVA or total metabolic energy.','','| BIO case | Left chamber J | Right chamber J | Total chamber J | Total port J | Valve loss J |','|---|---:|---:|---:|---:|---:|']
    for name,r in results.items():
        s=r['runs'][-1]['sides'];l=s['left']['work'];v=s['right']['work'];rows.append(f"| {name} | {l['chamber_J']:.5f} | {v['chamber_J']:.5f} | {l['chamber_J']+v['chamber_J']:.5f} | {l['port_J']+v['port_J']:.5f} | {l['valve_loss_J']+v['valve_loss_J']:.5f} |")
    rows+=['','Checks: pointwise energy identity, endpoint valve inertial storage versus integrated storage rate, direct pressure-volume quadrature, paired hydraulic inputs and 0.5/1 ms and 1/2 ms output-grid comparisons. Absolute energy tolerance: 0.005 J. See verification.json for achieved errors.','','No loop is artificially closed across changing cavity inventory. Initial elastic/inertial energy and external-pressure exchange prohibit treating transient net work as total active muscle expenditure. Auxiliary coronary supply is assessed separately in T5.1; this report does not clear those deficits.','','[T5.2 evidence and execution plan](../../docs/m5-t52-energy.md).']
    (OUT/'report.md').write_text('\n'.join(rows)+'\n')
if __name__=='__main__':main()
