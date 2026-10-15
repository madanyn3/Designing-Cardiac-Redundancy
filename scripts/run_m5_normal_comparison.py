"""Normal-heart mechanical comparator using the same BIO energy equations."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.whole_network import WholeNetworkModel
from cardiac_redundancy.integration import integrate_source
from cardiac_redundancy.biological_energy import chamber_power
from cardiac_redundancy.pump_energy import MMHG_ML_TO_J
OUT=ROOT/'reports/m5-normal-comparison'
def main():
    OUT.mkdir(parents=True,exist_ok=True);seed=np.load(ROOT/'reports/m4-baseline/native-cycle.npz');runs=[]
    for label,tol,step,sample in [('nominal',1e-7,.01,.001),('refined',1e-9,.005,.0005)]:
        m=WholeNetworkModel(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml',arterial_fraction=.2,inlet_resistance_fraction=.1)
        tr=integrate_source(m,8.375,initial_state=seed['states'][:,0],active=True,rtol=tol,max_step=step,sample_step=sample)
        assert tr.termination is None
        n=len(m.volume_names);R=m.base.source.resistances;L=m.base.source.inertances;rows=[];pressures=[]
        for t,y in zip(tr.time,tr.states.T):
            p,_=m.pressures(t,y,True);side={}
            for name,node,inn,out,j in [('left','LV','PU','AO',0),('right','RV','VC','PA',2)]:
                side[name]=chamber_power(p[inn],p[node],p[out],p['PERI'],y[n+j],y[n+j+1],R[j],R[j+1],L[j],L[j+1])
            rows.append(side);pressures.append([p[node] for node in m.volume_names])
        pressure=np.array(pressures);metrics={}
        for window,lo,hi in [('cycle',0,.75),('matched_8s',.375,8.375)]:
            mask=(tr.time>=lo-1e-9)&(tr.time<=hi+1e-9);t=tr.time[mask];idx=np.flatnonzero(mask);sides={}
            for name,j in [('left',0),('right',1)]:
                w={k.replace('_W','_J'):float(np.trapezoid([rows[i][name][k] for i in idx],t)) for k in rows[0][name] if k.endswith('_W')}
                raw=np.array([rows[i][name]['transmural_W'] for i in idx]);positive=float(np.trapezoid(np.maximum(raw,0),t)/(hi-lo))
                pv=float(-MMHG_ML_TO_J*np.sum(.5*(pressure[idx[1:],j]+pressure[idx[:-1],j])*np.diff(tr.states[j,idx])))
                coarse=idx[::2];coarse=np.unique(np.r_[coarse,idx[-1]])
                quad=abs(w['chamber_J']-float(np.trapezoid([rows[i][name]['chamber_W'] for i in coarse],tr.time[coarse])))
                delta=rows[idx[-1]][name]['valve_stored_J']-rows[idx[0]][name]['valve_stored_J']
                assert max(abs(pv-w['chamber_J']),quad,abs(delta-w['valve_storage_rate_J']))<.005
                sides[name]=dict(work=w,positive_transmural_mean_W=positive,PV_error_J=abs(pv-w['chamber_J']),coarse_error_J=quad)
            metrics[window]=dict(duration_s=hi-lo,sides=sides,CO_L_min=float(np.trapezoid(tr.states[n+1,mask],t)/(hi-lo)*.06),mean_AO_mmHg=float(np.trapezoid(pressure[mask,2],t)/(hi-lo)),minimum_AO_mmHg=float(pressure[mask,2].min()),maximum_AO_mmHg=float(pressure[mask,2].max()))
        drift=float(np.max(abs(tr.states[:n].sum(axis=0)-5250)));assert drift<1e-5
        np.savez_compressed(OUT/f'{label}.npz',time_s=tr.time,states=tr.states,pressures_mmHg=pressure)
        runs.append(dict(label=label,rtol=tol,max_step_s=step,sample_step_s=sample,blood_error_mL=drift,metrics=metrics));print(label,'verified',flush=True)
    err=max(abs(runs[0]['metrics'][w]['sides'][side]['work'][k]-runs[1]['metrics'][w]['sides'][side]['work'][k]) for w in metrics for side in ('left','right') for k in metrics[w]['sides'][side]['work']);assert err<.005
    files=[Path(__file__),ROOT/'reports/m4-baseline/native-cycle.npz',*list((ROOT/'src/cardiac_redundancy').glob('*.py'))]
    result=dict(status='verified',runs=runs,paired_work_error_J=err,reference='project normal resting model; not measured individual human',hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    (OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
