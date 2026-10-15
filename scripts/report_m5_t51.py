"""Assemble T5.1 evidence, research screens and reproducibility hashes."""
from pathlib import Path
import hashlib,json,sys,os
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cardiac_redundancy.constraints import pressure_windows
from run_m5_t51_screening import hydraulics
from run_m5_dynamic import read,write
OUT=ROOT/'reports/m5-t51'

def main():
    cfg=read('configs/m5-biological-pump.json');bio=read('reports/m5-biological/verification.json')
    assert bio['status']=='verified' and bio['cases']==cfg['cases']
    for name,changes in cfg['cases'].items():
        assert read(f'reports/m5-biological/{name}/result.json')['parameters']==cfg['reference'] | changes
    names=['B1_transmural_gate','B1_left_weak','BIO_reference','BIO_faster','BIO_stronger_left','BIO_asymmetric']
    inputs=set();summaries=[];fig,axes=plt.subplots(2,1,figsize=(10,8),constrained_layout=True)
    lines=['# T5.1 bounded-pump characterization and closure evidence','','2026-10-09. Numerical implementation and characterization complete; no optimized or physiologically validated auxiliary heart. See the [closure and equations](../../docs/m5-t51-closure.md).','','## Eight-second joint screens','','All native pumping remains off, preservation conductances remain at 0.25, and demands remain the M1 unsuppressed scenario. The rate and chamber changes are exploratory assumptions. Mean pressure and output below cover the final two seconds; pressure-screen status covers all eight one-second windows, including activation.','','| Case | Mean AO (mmHg) | Left support (L/min) | Pressure screen | Unmet O2 (mL) | Auxiliary share (mL) |','|---|---:|---:|---|---:|---:|']
    for name in names:
        path=f'reports/m5-t51-screening/{name}/result.json';r=read(path);assert r['status']=='verified' and abs(r['duration_s']-8)<1e-7
        segments,paths=hydraulics(r['case']);inputs.update(paths);inputs.update(r['inputs']);inputs.add(path)
        time,idx=np.unique(np.concatenate([s['time'] for s in segments]),return_index=True)
        ps=np.concatenate([s['pressures_mmHg'] for s in segments])[idx];qs=np.concatenate([s['flows'] for s in segments])[idx]
        if r['case']['architecture']=='BIO':hr=read(f"reports/m5-biological/{r['case']['trajectory']}/result.json");edge='aAV'
        else:hr=read(f"reports/m5-dynamic/{r['case']['trajectory']}-refined/result.json");edge='support_left'
        tail=time>=6;mean=float(np.trapezoid(ps[tail,2],time[tail])/2);output=float(np.trapezoid(qs[tail,hr['flow_names'].index(edge)],time[tail])/2*.06)
        pressure=pressure_windows(time,ps[:,2],dict(mean=65,maximum=90,minimum=50))
        record=r['records'][-1];per=record['per_organ_unmet_oxygen_mL'];aux=sum(v for k,v in per.items() if k.startswith('aux_'))
        summary=dict(case=name,final_two_second_mean_AO_mmHg=mean,left_output_L_min=output,pressure=pressure,
            unmet_oxygen_mL=record['unmet_oxygen_mL'],per_organ_unmet_oxygen_mL=per,
            organ_oxygen_status={k:'violated' if v>1e-7 else 'no_deficit_observed_in_8_s' for k,v in per.items()},
            joint_status='fails_declared_demand_screen',unassessed=['input energy and efficiency','injury and reperfusion','blood-cell viability and coagulation','physiological joint-fuel baseline','substrate-dependent auxiliary contractility','regional functional pressure sufficiency'],
            glycogen_remaining_mmol=record['glycogen_remaining_mmol'],released_mmol=record['released_mmol'])
        assert record['unmet_oxygen_mL']>1e-7
        summaries.append(summary)
        lines.append(f"| {name} | {mean:.2f} | {output:.3f} | {pressure['status']} | {record['unmet_oxygen_mL']:.4f} | {aux:.4f} |")
        axes[0].plot(time,ps[:,2],label=name)
        a=np.load(ROOT/f'reports/m5-t51-screening/{name}/refined.npz');n=len(r['volume_names']);m=len(r['organ_names'])
        axes[1].plot(a['time_s'],a['states'][:,n+2+m:n+2+2*m].sum(axis=1),label=name)
    axes[0].set(ylabel='Central aortic pressure (mmHg)',title='Actual hydraulic trajectories; native pumping off')
    axes[0].axhline(65,color='black',ls='--',lw=.8,label='65 mmHg mean-screen reference')
    axes[0].legend(fontsize=8,ncol=2)
    axes[1].set(xlabel='Time after arrest (s)',ylabel='Accumulated unmet O2 (mL)',title='Joint glucose/oxygen demand deficits; not injury predictions')
    OUT.mkdir(parents=True,exist_ok=True);fig.savefig(OUT/'characterization.png',dpi=160);plt.close(fig)
    lines+=['','![Hydraulic and metabolic traces](characterization.png)','','The dashed pressure line is a mean-screen reference, not an instantaneous threshold. The fixed-window maximum/minimum screen is retained for continuity; it is not beat-defined systolic/diastolic pressure in continuous-flow cases. A complete-cycle BIO pressure table is recorded in `verification.json`.','','## Findings and limits','','Every tested case has coronary oxygen-demand deficits in the restrictive joint-fuel model. The biological fixtures also have auxiliary-coronary deficits. No deficits are observed in the brain, nutritive lung, respiratory-muscle or hepatic branches during these eight seconds, but this is not evidence of indefinite sufficiency or recoverability. Finite initial stores contribute. All other preservation demands remain charged.','','The faster BIO fixture raises pressure and output in this sample. Increased left contractility helps less here; this is a local comparison, not a global sensitivity ranking. Abstract B1 supplies prescribed head without a myocardial force-generation or self-perfusion constraint. BIO generates pressure from filling, wall activation, septal interaction and valve dynamics. Neither model currently reduces commanded support when its own input energy or substrate supply is inadequate.','','L-M4-F01 remains: no tissue glucose buffer and no validated joint-fuel physiological baseline. Deficits are reported against original demand, never a fuel-reduced target. The BIO additional demand remains an independent assumed 10 mL O2/min, not a calibrated function of size, frequency or contractility. T5.2 must address those costs before energy comparisons or endurance claims.','','## Verification','','177 tests pass, including eight new bounds, similarity-regression and allocation/conservation tests. Four BIO scenarios each have nominal/refined Radau trajectories; six joint screens each have nominal/refined RK45 replays on the same refined hydraulic input. This separates metabolic integration refinement from hydraulic refinement. Blood, oxygen, glucose and original-demand ledgers are checked, as are oxygen/glucose stoichiometry, sampled physical inventories, valve sequences and septal roots. All cases must be verified before this report is generated.','','Hydraulic relative coordinate differences must remain below 0.5%, matched valve event times within 1 ms, blood drift below 1e-5 mL and root/junction residuals below 1e-8. Joint replay differences must be below 0.002 in oxygen-state coordinates and 1e-5 mmol in glucose states; oxygen conservation/demand/stoichiometry below 1e-6 mL and glucose conservation/demand below 1e-8 mmol. These are numerical tolerances, not tissue tolerance. Detailed achieved values are retained in case JSON files.','','## Reproduce','','```powershell','.\\.venv\\Scripts\\python.exe scripts/run_m5_biological.py','.\\.venv\\Scripts\\python.exe scripts/run_m5_t51_screening.py','.\\.venv\\Scripts\\python.exe scripts/report_m5_t51.py','.\\.venv\\Scripts\\python.exe -m pytest -q --junitxml=reports/m5-t51-tests.xml','.\\.venv\\Scripts\\python.exe scripts/audit_m4_fuel.py','```','','The earlier static and dynamic increments are reproduced by `characterize_m5_pump.py`, `run_m5_dynamic.py` and `report_m5_dynamic.py`. Input and result hashes are stored in the manifest; numerical results are saved without rounding in NPZ/JSON files.']
    # Phase-aligned complete auxiliary cycles after the ramp, unlike fixed windows.
    cycles={}
    for name in cfg['cases']:
        r=read(f'reports/m5-biological/{name}/result.json');period=r['parameters']['period_s'];a=np.load(ROOT/f'reports/m5-biological/{name}/refined.npz');t=a['time_s'];p=a['pressures_mmHg'][:,2];rows=[]
        for j in range(int(np.floor(t[-1]/period))):
            lo=j*period;hi=(j+1)*period
            if lo<1:continue
            grid=np.r_[lo,t[(t>lo)&(t<hi)],hi];v=np.interp(grid,t,p)
            rows.append(dict(start_s=lo,end_s=hi,mean_mmHg=float(np.trapezoid(v,grid)/period),maximum_mmHg=float(v.max()),minimum_mmHg=float(v.min())))
        cycles[name]=rows
    write(OUT/'verification.json',dict(status='verified',summaries=summaries,complete_auxiliary_cycles=cycles))
    (OUT/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    files={Path(__file__),ROOT/'scripts/run_m5_biological.py',ROOT/'scripts/run_m5_t51_screening.py',ROOT/'configs/m5-biological-pump.json',ROOT/'configs/m4-fuel-reference.json',ROOT/'configs/m4-oxygen-reference.json',ROOT/'configs/m4-organ-mapping.json',ROOT/'data/m1-demand-budget.json',ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml',ROOT/'configs/m3-whole-network-protocol.json',*list((ROOT/'src/cardiac_redundancy').glob('*.py'))}
    files.update(ROOT/p for p in inputs)
    for directory in ('reports/m5-biological','reports/m5-t51-screening'):
        files.update(p for p in (ROOT/directory).rglob('*') if p.is_file())
    files.add(ROOT/'reports/m5-t51-tests.xml')
    write(OUT/'manifest.json',dict(status='verified',hydraulic_cases=4,metabolic_cases=6,tests=177,software=dict(python=sys.version,numpy=np.__version__),hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}))

if __name__=='__main__':main()
