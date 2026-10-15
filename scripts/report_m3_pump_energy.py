"""Reanalyse verified saved hydraulics; no new ODE trajectories are claimed."""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.pump_energy import stage_power, integrate_power, POWER_FIELDS

SOURCE=ROOT/'reports/m3-whole-network'
OUT=ROOT/'reports/m3-pump-energy'


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    results={};hashes={};max_balance=0.;max_flow=0.;max_grid=0.
    for path in sorted(SOURCE.glob('*/result.json')):
        r=json.loads(path.read_text());name=path.parent.name
        if r['status']!='completed':raise RuntimeError('Incomplete source case: '+name)
        hashes[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
        sums={side:{f.replace('_W','_J'):0. for f in POWER_FIELDS} for side in ('right','left')}
        coarse={side:dict.fromkeys(sums[side],0.) for side in sums}
        peaks={side:dict.fromkeys(POWER_FIELDS,0.) for side in sums}
        for seg in sorted(path.parent.glob('segment-*.npz')):
            hashes[str(seg.relative_to(ROOT))]=hashlib.sha256(seg.read_bytes()).hexdigest()
            with np.load(seg) as z:
                t=z['time'];p=z['pressures'];q=z['flows'];saved={'time_s':t}
                for side,a,b,head in (('right','VC','PA',20.),('left','PU','AO',100.)):
                    records=[]
                    for i,ti in enumerate(t):
                        h=head*np.clip(ti-r['arrest_s'],0.,1.)
                        power=stage_power(p[i,r['volume_names'].index(a)],p[i,r['volume_names'].index(b)],h,.1,120.)
                        max_flow=max(max_flow,abs(power['flow_mL_s']-q[i,r['flow_names'].index('support_'+side)]))
                        max_balance=max(max_balance,abs(power['head_W']-sum(power[k] for k in POWER_FIELDS[1:])))
                        records.append(power)
                    work=integrate_power(t,records)
                    indices=np.unique(np.r_[np.arange(0,len(t),2),len(t)-1])
                    reduced=integrate_power(t[indices],[records[i] for i in indices])
                    for field in work:
                        sums[side][field]+=float(work[field][-1]);coarse[side][field]+=float(reduced[field][-1])
                        saved[side+'_'+field]=work[field]
                    for field in POWER_FIELDS:
                        saved[side+'_'+field]=np.array([x[field] for x in records])
                        peaks[side][field]=max(peaks[side][field],float(saved[side+'_'+field].max()))
            dest=OUT/name;dest.mkdir(exist_ok=True)
            np.savez_compressed(dest/seg.name,**saved)
        max_grid=max(max_grid,max(abs(sums[s][k]-coarse[s][k])/max(1.,abs(sums[s][k]))*100 for s in sums for k in sums[s]))
        results[name]=dict(work_J=sums,peak_W=peaks)
    comparisons={}
    for name,r in results.items():
        base=name.rsplit('-',1)[0];ref=results[base+'-nominal']
        comparisons[name]=max(abs(r['work_J'][s][k]-ref['work_J'][s][k])/max(1.,abs(ref['work_J'][s][k]))*100 for s in ('right','left') for k in ref['work_J'][s])
    passed=max_balance<1e-12 and max_flow<1e-9 and max_grid<.5 and max(comparisons.values())<.5
    metrics=dict(passed=passed,max_balance_residual_W=max_balance,max_saved_flow_difference_mL_s=max_flow,
                 max_grid_work_difference_pct=max_grid,max_solver_work_difference_pct=max(comparisons.values()),cases=results)
    (OUT/'metrics.json').write_text(json.dumps(metrics,indent=2))
    for f in (Path(__file__),ROOT/'src/cardiac_redundancy/pump_energy.py'):
        hashes[str(f.relative_to(ROOT))]=hashlib.sha256(f.read_bytes()).hexdigest()
    (OUT/'manifest.json').write_text(json.dumps(dict(input_hashes=hashes,scope='Postprocessing saved eight-second hydraulic fixtures; original 20/100 mmHg heads, R=.1, capacity=120.'),indent=2))
    lines=['# Pump energy accounting verification','',
           'Mechanical ledger on the saved whole-network study; no new ODE runs. All eight seconds include restoration of conductance from seconds 3 to 5. These are not fixed-routing survival or energy-budget estimates.','',
           '| Case | Right head work J | Left head work J | Total port work J | Total resistance loss J | Total limiter loss J |',
           '|---|---:|---:|---:|---:|---:|']
    for name,r in results.items():
        if not name.endswith('-nominal'):continue
        w=r['work_J'];total=lambda k:sum(w[s][k] for s in w)
        lines.append(f"| {name} | {w['right']['head_J']:.6f} | {w['left']['head_J']:.6f} | {total('port_J'):.6f} | {total('resistance_loss_J'):.6f} | {total('limiter_loss_J'):.6f} |")
    lines += ['',f'Checks passed: **{passed}**. Maximum instantaneous balance residual {max_balance:.3e} W; reconstructed flow difference {max_flow:.3e} mL/s.',
              f'Maximum work change from approximately 5 to 10 ms sampling: {max_grid:.6f}%; maximum nominal/refined/BDF work difference: {max(comparisons.values()):.6f}%. Differences use max(1 J, absolute reference work); acceptance <0.5%.',
              '', 'Each saved integration segment is integrated separately by trapezoidal quadrature, preserving one-sided values at intervention boundaries. Segment energy arrays start at zero; metrics.json sums segment totals. Peak values are sampled peaks, not guaranteed continuous-time extrema.',
              '', 'Head work = signed port work + resistance loss + idealized limiter loss. The limiter interpretation is a modeling assumption. No efficiency, idle power, biological metabolism, finite fuel depletion or recovery claim is included. Negative port work is not credited as energy regeneration.',
              '', '[Model and equations](../../docs/m3-pump-energy.md). [Automated test record](../m3-pump-energy-tests.xml). Reproduce with `.venv/Scripts/python.exe scripts/report_m3_pump_energy.py`.']
    (OUT/'report.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:v for k,v in metrics.items() if k!='cases'},indent=2))
    return 0 if passed else 1


if __name__=='__main__':raise SystemExit(main())
