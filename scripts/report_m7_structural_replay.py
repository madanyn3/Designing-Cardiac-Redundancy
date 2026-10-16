"""Summarize the paired structural replay without assigning viability labels."""
import hashlib
import json
import os
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'reports/m7-structural-replay'
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'reports/m7-plot-cache'))


def main():
    manifest = json.loads((OUT/'input-manifest.json').read_text())
    for path, expected in manifest['hashes'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != expected:
            raise RuntimeError(f'Changed input: {path}')
    frozen = json.loads((ROOT/'reports/m7-factorial/hydraulic_nominal/result.json').read_text())
    reference = json.loads((OUT/'reference/result.json').read_text())
    if reference['status'] != 'verified' or abs(reference['records'][-1]['unmet_oxygen_mL']-frozen['records'][-1]['unmet_oxygen_mL']) > 1e-7:
        raise RuntimeError('Frozen reference reproduction failed')
    rows = []
    lines = ['# T7.4 extraction/transition structural replay', '',
             'Ten-second fixed hydraulic trajectory. All pressures, pump work and flows are prescribed identically; no survival or energetic-feedback claim.', '',
             '| Form | Requested O2 (mL) | Unmet O2 equivalent (mL) | Active / preservation / pump deficit (mL) | Paired O2-state difference (mL) |',
             '|---|---:|---:|---|---:|']
    for spec in manifest['specs']:
        folder = OUT/spec['name']
        result = json.loads((folder/'result.json').read_text())
        if result['status'] != 'verified':
            raise RuntimeError(f'Unverified case: {spec}')
        record = result['records'][-1]
        data = np.load(folder/'metabolic-refined.npz')
        extra = np.maximum(0, data['glucose_ceiling_mL_s']-data['oxygen_consumed_mL_s'])
        row = dict(name=spec['name'], requested_mL=sum(g['requested_O2_mL'] for g in record['groups'].values()),
                   deficit_mL=record['unmet_oxygen_mL'], groups=record['groups'],
                   additional_oxygen_limitation_integral_mL=float(np.trapezoid(extra.sum(axis=1), data['time_s'])),
                   paired_oxygen_difference_mL=result['paired_max_oxygen_state_difference'],
                   paired_deficit_difference_mL=abs(result['records'][0]['unmet_oxygen_mL']-record['unmet_oxygen_mL']),
                   paired_glucose_difference_mmol=result['paired_max_glucose_state_difference_mmol'])
        rows.append(row)
        groups=' / '.join(f"{record['groups'][role]['unmet_O2_mL']:.6g}" for role in ['active','preservation','pump'])
        lines.append(f"| {spec['name']} | {row['requested_mL']:.7g} | {row['deficit_mL']:.7g} | {groups} | {row['paired_oxygen_difference_mL']:.3g} |")
    lines += ['', 'Finite transitions change the integrated requested demand as well as its time course. Flow-dependent ceilings are calibrated at the nominal final-window mean absolute flow, not measured permeability values. Serial hepatic depletion, inventories and demand ledgers are retained.', '',
              'Scope: this completes one structural experiment block only. Pressure regulation, finite glucose buffering and supply-limited pump contraction remain outside this replay. T7.4/T7.5 stay open until their combined decision report is prepared.']
    (OUT/'summary.json').write_text(json.dumps(rows, indent=2)+'\n', encoding='utf-8')
    (OUT/'report.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    labels = [r['name'] for r in rows]
    for ax, field, title in [(axes[0], 'requested_mL', 'Integrated requested oxygen'),
                             (axes[1], 'deficit_mL', 'Unmet aerobic demand')]:
        ax.barh(labels, [r[field] for r in rows])
        ax.set(title=title, xlabel='mL O2 equivalent over 10 s')
        ax.invert_yaxis()
    fig.suptitle('Equation-form sensitivity on identical prescribed hydraulics')
    fig.tight_layout()
    fig.savefig(OUT/'comparison.png', dpi=160)
    plt.close(fig)
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
