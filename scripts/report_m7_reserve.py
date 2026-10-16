"""Paired finite-inventory coverage bounds, never an energetic circulation run."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.finite_reserve import cover_deficit
from run_m7_global_metabolism import network


def main():
    out=ROOT/'reports/m7-finite-reserve';out.mkdir(parents=True,exist_ok=True)
    inputs=[Path(__file__),ROOT/'src/cardiac_redundancy/finite_reserve.py',ROOT/'docs/m7-final-structural-protocol.md']
    rows=[]
    for name in ['reference','glucose_flow','combined']:
        folder=ROOT/'reports/m7-structural-replay'/name
        inputs += [folder/'result.json',folder/'metabolic-refined.npz']
        result=json.loads((folder/'result.json').read_text())
        if result['status']!='verified':raise RuntimeError('Unverified source replay')
        h,organs=network(result['case'])
        data=np.load(folder/'metabolic-refined.npz');t=data['time_s']
        d=np.maximum(data['oxygen_unmet_mL_s'],0)
        sampled=np.trapezoid(d,t,axis=0)
        for seconds in [0,1,5]:
            initial=seconds*np.array([o.demand_mL_min for o in organs])/60
            # Initial 1e-7/.01 versus 1e-9/.005 pair missed the 1e-6
            # coverage-difference gate; retained in rejected-initial-settings.
            coarse=cover_deficit(t,d,initial,rtol=1e-9,max_step=.0025)
            fine=cover_deficit(t,d,initial,rtol=1e-11,max_step=.00125)
            error=float(np.max(np.abs(coarse[1]-fine[1])))
            if error>1e-6 or np.any(fine[1][-1]>sampled+1e-7):raise RuntimeError('Coverage/refinement failure')
            row=dict(source=name,initial_seconds=seconds,initial_inventory_mL_equivalent=float(initial.sum()),
                auxiliary_initial_inventory_mL_equivalent=float(sum(initial[j] for j,o in enumerate(organs) if o.parent_pool=='auxiliary_heart')),
                ledger_deficit_mL=result['records'][-1]['unmet_oxygen_mL'],sampled_deficit_mL=float(sampled.sum()),
                covered_mL=float(fine[1][-1].sum()),residual_mL=float((sampled-fine[1][-1]).sum()),
                balance_error_mL=fine[2],paired_difference_mL=error)
            rows.append(row)
            np.savez_compressed(out/f'{name}-{seconds}s.npz',time_s=t,remaining=fine[0],covered=fine[1],initial=initial)
    summary=dict(scope='Optimistic finite substrate-equivalent deficit coverage; additional oxygen and force unverified',rows=rows,
                 numerical_pairs=[dict(rtol=1e-9,max_step=.0025),dict(rtol=1e-11,max_step=.00125)],
                 hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs})
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    lines=['# Finite reserve coverage diagnostic','',summary['scope'],'',
           '| Source | Reserve seconds | Initial inventory (mL equivalent) | Sampled deficit | Covered | Residual |',
           '|---|---:|---:|---:|---:|---:|']
    for r in rows:
        residual='~0' if abs(r['residual_mL'])<1e-7 else f"{r['residual_mL']:.7g}"
        lines.append(f"| {r['source']} | {r['initial_seconds']} | {r['initial_inventory_mL_equivalent']:.7g} | {r['sampled_deficit_mL']:.7g} | {r['covered_mL']:.7g} | {residual} |")
    lines+=['','Initial inventory is explicitly added and never refilled; it is not credited to the frozen blood/liver model. Covered plus remaining equals initial inventory for every organ. Demand is interpolated from saved samples: sampled integrals and original cumulative ledgers are both retained. Saved-waveform integration differs from the original glucose-flow deficit ledger by about 0.7%; this reconstruction error is distinct from reserve-solver refinement. Residuals below 1e-7 mL are displayed as ~0; signed numerical residuals are retained in JSON, not clipped from states. Zero/one/five seconds are assumptions, not measured available glycogen. A zero residual here would not establish sufficient oxygen, energy-supported pumping, recoverability or survival.']
    (out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('\n'.join(lines))


if __name__=='__main__':main()
