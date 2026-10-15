"""Render the verified finite-fuel results without changing historical M4 reports."""
from pathlib import Path
import json,os
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    out=ROOT/'reports/m4-fuel';manifest=json.loads((out/'manifest.json').read_text())
    if manifest['status']!='verified':raise RuntimeError('Only render verified runs')
    cases=manifest['results'];short=[r for r in cases if 'case' in r];long=[r for r in cases if 'case' not in r]
    lines=['# M4 finite-fuel verification','',
        'Completed 2026-10-09. **Verified at explicit exploratory-model scope. No recoverable support time is established.**',
        '',f"{len(short)} coupled short scenarios and {len(long)} maintained-flow controls; {manifest['paired_runs']} paired integrations. [Equations, sources and limitations](../../docs/m4-fuel-accounting.md). [Machine-readable results and provenance](manifest.json).",'',
        '## Coupled short replays','',
        'Original organ demand is retained. Actual glucose withdrawal is tied to actual oxygen consumption; finite oxygen stores and respiratory feedback remain active. Nearly empty glucose is a numerical stress condition, not a proposed subject state. BIO retains 10 mL O2/min additional assumed working-myocardium demand.','',
        '| Case | Duration (s) | Unmet O2 (mL) | Glucose exported (mmol) | Final glycogen (mmol) | Paired glucose-state difference (mmol) |',
        '|---|---:|---:|---:|---:|---:|']
    for r in short:
        v=r['records'][-1]
        lines.append(f"| {r['case']['name']} | {r.get('duration_s',8):g} | {v['unmet_oxygen_mL']:.6g} | {v['released_mmol']:.6g} | {v['glycogen_remaining_mmol']:.6f} | {r['paired_max_glucose_state_difference_mmol']:.3g} |")
    lines+=['','These are fuel/oxygen results, not new hydraulic success claims: the historical BIO pressure failure remains. Initial normal blood glucose buffers a short interruption in liver release, but does not remove instantaneous glucose-delivery limitations. The native-cycle control explicitly tests this: no tissue glucose buffer is credited, so pulsatile glucose delivery can limit uptake despite the existing oxygen buffer. The oxygen-only baseline remains verified; a fully physiological joint-fuel resting baseline is **not established**. This limitation must accompany M5 searches.','',
        '## Long controls: distinguish release from stock','',
        'These six-hour glucose-only calculations impose fixed, balanced mean native flow and fixed mean blood volumes. They presume oxygen availability for requested glucose use. They do not simulate long-duration backup haemodynamics, respiratory failure, pump energy, injury or reperfusion. Times below bracket the first **sampled numerical glucose-demand shortfall**, using 10-second samples and a 1e-8 mmol/s numerical rate tolerance; they are not clinical thresholds or safe durations.','',
        '| Control | Release cap (mmol/min) | Export increment cap (mmol/L) | Initial glycogen (mmol) | First sampled shortfall bracket (min) | Glycogen left at 6 h (mmol) |',
        '|---|---:|---:|---:|---|---:|']
    for r in long:
        v=r['records'][-1];p=r['parameters'];br=v['first_sampled_glucose_shortfall_bracket_s']
        time='none sampled' if br is None else f'{br[0]/60:.2f}–{br[1]/60:.2f}'
        initial=p['liver_volume_L']*p['glycogen_mmol_L_liver']*p['accessible_fraction']
        lines.append(f"| {r['name']} | {p['release_cap_mmol_min']:.6g} | {p['export_increment_mmol_L']:.6g} | {initial:.3f} | {time} | {v['final_glycogen_mmol']:.6g} |")
    lines+=['','The 0.21 mmol/min scale is transferred from observed healthy fasting glycogenolysis; its use as a maximum is an assumption. The approximately 0.791 scale is transferred from **total** hepatic glucose output and is only a hypothetical glycogen-only comparator. The 2.1 scale is an A stress capacity, ten times the reference. Export increment, extraction and hepatic function are additional assumed limits. Larger stores cannot cure a delivery-rate bottleneck.','',
        '![Finite-fuel control trajectories](fuel-controls.png)','',
        '## Verification and reproducibility','']
    allrecords=[v for r in cases for v in r['records']]
    max_glucose=max(v.get('glucose_balance_error_mmol',v.get('balance_error_mmol',0)) for v in allrecords)
    max_o=max(v.get('oxygen_balance_error_mL',0) for v in allrecords)
    stoich=max(v.get('oxidation_stoichiometry_error_mL',0) for v in allrecords)
    max_short=max(r['paired_max_glucose_state_difference_mmol'] for r in short)
    max_long=max(r['paired_max_state_difference_mmol'] for r in long)
    lines += [f'- Maximum glucose-equivalent conservation error: **{max_glucose:.3g} mmol**.',
        f'- Maximum coupled oxygen conservation error: **{max_o:.3g} mL**; oxidation/glucose stoichiometry error: **{stoich:.3g} mL O2**.',
        f'- Maximum paired glucose-state difference: **{max_short:.3g} mmol** in short replays (limit 1e-5), **{max_long:.3g} mmol** in long controls (limit 1e-4).',
        f"- Maximum mean-flow balancing adjustment: **{max(r['maximum_flow_projection_mL_s'] for r in long):.6g} mL/s**; all corrected node balances below 1e-9 mL/s.",
        '- Original-demand identities and physical state bounds are checked separately; no failed run is counted as a physiological success.',
        '- Earlier RK45 pair settings were rejected: B1 reference glucose difference 1.12876e-5 mmol, and liver-loss difference 1.27589e-5 mmol. Tighter documented settings retained the 1e-5 criterion.',
        '', '```powershell',r'.\.venv\Scripts\python.exe -m pytest -q',r'.\.venv\Scripts\python.exe scripts/run_m4_fuel.py',r'.\.venv\Scripts\python.exe scripts/report_m4_fuel.py','```','',
        'The manifest records parameters, software, input/output SHA-256 hashes, solver settings and both runs. `--manifest-only` consolidates existing matching protocol results; it does not rerun integration. Historical M4 oxygen-only artifacts remain historical and are not relabeled fuel-constrained.',
        '', '## Gate disposition','',
        'The finite glycogen/glucose accounting obligation is complete at the declared restrictive exploratory scope. M4 can hand this tested constraint to M5. Mixed fuels, glucose regulatory physiology, nonoxidative cellular demand, liver function, preservation and reperfusion remain explicit limitations; no optimized pump or recoverable duration follows from this report.']
    (out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    fig,axes=plt.subplots(2,1,figsize=(10,8),sharex=True)
    for r in long:
        if r['name'] not in ('mean_reference','mean_no_release','mean_high_release','mean_export_capacity','mean_reduced_demand'):continue
        z=np.load(out/r['name']/'refined.npz');n=len(r['volume_names']);t=z['time_s']/60
        label=r['name'].removeprefix('mean_').replace('_',' ')
        axes[0].plot(t,z['states'][:,n],label=label)
        j=r['volume_names'].index('AO')
        axes[1].plot(t,z['states'][:,j]/z['volumes_mL'][j]*1000*r['parameters']['plasma_to_whole_blood_ratio'],label=label)
    axes[0].set_ylabel('Accessible liver glycogen (mmol)');axes[0].legend(fontsize=9)
    axes[1].set_ylabel('Arterial plasma-equivalent glucose (mmol/L)');axes[1].set_xlabel('Maintained-flow control time (min)')
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('Finite-fuel controls — prescribed native mean flow, not backup support duration')
    fig.tight_layout();fig.savefig(out/'fuel-controls.png',dpi=160);plt.close(fig)
    print(out/'report.md')


if __name__=='__main__':main()
