"""Report newly coupled flow-regulation experiments with paired checks."""
import json
import hashlib
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reports/m7-regulation'


def main():
    protocol=json.loads((OUT/'input-manifest.json').read_text())
    for path,digest in protocol['hashes'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:raise RuntimeError(f'Changed input {path}')
    rows=[]
    lines=['# Bounded cerebral regulation: coupled circulation and metabolism','',
           '| Form | Final MAP mmHg | Brain flow mL/s | Auxiliary CO L/min | Chamber work J | Deficit mL O2-equivalent |',
           '|---|---:|---:|---:|---:|---:|']
    for spec in protocol['specs']:
        folder=OUT/spec['name']
        h=json.loads((folder/'hydraulics.json').read_text())
        m=json.loads((folder/'result.json').read_text())
        if m['status']!='verified' or any(r['termination'] is not None for r in h['runs']):raise RuntimeError('Unverified or terminated regulation case')
        a=np.load(folder/'hydraulic-refined.npz');t=a['time_s'];mask=t>=t[-1]-2
        idx=h['flow_names'].index('brain')
        brain=float(np.trapezoid(a['flows_mL_s'][mask,idx],t[mask])/(t[-1]-t[mask][0]))
        r=h['runs'][-1];mr=m['records'][-1]
        row=dict(name=spec['name'],MAP_mmHg=r['final_MAP_mmHg'],brain_flow_mL_s=brain,
                 aux_CO_L_min=r['aux_left_L_min'],chamber_work_J=sum(r['chamber_work_J']),
                 deficit_mL=mr['unmet_oxygen_mL'],release_s=r['release_s'],groups=mr['groups'],
                 paired_deficit_difference_mL=abs(mr['unmet_oxygen_mL']-m['records'][0]['unmet_oxygen_mL']))
        rows.append(row)
        lines.append(f"| {row['name']} | {row['MAP_mmHg']:.6f} | {brain:.6f} | {row['aux_CO_L_min']:.6f} | {row['chamber_work_J']:.6f} | {row['deficit_mL']:.8g} |")
    frozen=json.loads((ROOT/'reports/m7-factorial/hydraulic_nominal/result.json').read_text())
    frozen_h=json.loads((ROOT/'reports/m7-factorial/hydraulic_nominal/hydraulics.json').read_text())
    if abs(rows[0]['deficit_mL']-frozen['records'][-1]['unmet_oxygen_mL'])>1e-7 or abs(rows[0]['MAP_mmHg']-frozen_h['runs'][-1]['final_MAP_mmHg'])>1e-6:
        raise RuntimeError('Passive reference reproduction failed')
    lines+=['','The flow target is the existing 11.2 mL/s construction anchor. Gains are hypothetical, instantaneous and bounded; this is not a physiological autoregulation calibration. Closing pressure remains active. Each alternative has new conservative hydraulics and a new workload-dependent metabolic replay. The observed controller may activate at different times because regulation changes pressure; these are whole-system comparisons under the same controller, not a fixed-trigger isolation of vascular effects.',
            '', 'All final MAPs must be interpreted against the declared 65 mmHg screen. Higher cerebral flow alone does not imply more surviving brain tissue; no injury/recovery model exists. Supply-limited pump contraction remains unimplemented.']
    (OUT/'summary.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
    (OUT/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('\n'.join(lines))


if __name__=='__main__':main()
