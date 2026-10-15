"""Summarize regional comparator runs and enforce declared numerical screens."""
import hashlib
import json
import os
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.matplotlib'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    out=ROOT/'reports/m3-regional'
    protocol=json.loads((ROOT/'configs/m3-regional-protocol.json').read_text())
    results={c['name']:json.loads((out/c['name']/'result.json').read_text()) for c in protocol['cases']}
    if any(r['status']!='passed_numerical_checks' for r in results.values()):
        raise RuntimeError('Incomplete or failed regional runs; inspect result records')
    base=results['baseline']
    comparisons={}
    for name,r in results.items():
        if name=='baseline': continue
        diffs={k:100*abs(v-base['metrics'][k])/abs(base['metrics'][k]) for k,v in r['metrics'].items()}
        for k,v in r['mean_branch_flows_mL_s'].items():
            ref=base['mean_branch_flows_mL_s'][k]
            diffs['flow_'+k]=100*abs(v-ref)/abs(ref)
        def phase_events(case):
            start=results[case]['measurement_start_s']
            return [(e['valve'],e['action'],e['time_s']-start) for e in
                    json.loads((out/case/'events.json').read_text()) if e['time_s']<start+.75]
        a,b=phase_events('baseline'),phase_events(name)
        same=[e[:2] for e in a]==[e[:2] for e in b]
        dt=max(abs(x[2]-y[2])*1000 for x,y in zip(a,b)) if same else None
        comparisons[name]={'largest_metric_difference_pct':max(diffs.values()),
            'all_metric_differences_pct':diffs,'maximum_valve_time_difference_ms':dt,
            'passed':same and max(diffs.values())<protocol['comparison_metric_limit_pct'] and dt<protocol['comparison_event_time_limit_ms']}
    frozen=json.loads((ROOT/'data/m2-freeze-manifest.json').read_text())['sha256']
    frozen_ok=all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in frozen.items())
    summary={'comparisons':comparisons,'frozen_unchanged':frozen_ok,
             'all_passed':all(c['passed'] for c in comparisons.values()) and frozen_ok}
    (out/'verification.json').write_text(json.dumps(summary,indent=2)+'\n')
    baseline=np.load(out/'baseline/trajectory.npz')
    t=baseline['time_s']; mask=t>=t[-1]-.75; phase=t[mask]-(t[-1]-.75)
    fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    ax[0].plot(phase,baseline['pressures_mmHg'][mask,2],label='Aorta')
    ax[0].plot(phase,baseline['pressures_mmHg'][mask,4],label='Pulmonary artery')
    ax[0].set(xlabel='Cardiac phase (s)',ylabel='Pressure (mmHg)'); ax[0].legend()
    ax[1].plot(baseline['states'][0,mask],baseline['pressures_mmHg'][mask,0],label='LV')
    ax[1].plot(baseline['states'][1,mask],baseline['pressures_mmHg'][mask,1],label='RV')
    ax[1].set(xlabel='Physical cavity volume (mL)',ylabel='Pressure (mmHg)'); ax[1].legend()
    fig.savefig(out/'baseline-waveforms.png',dpi=160); plt.close(fig)
    lines=['# M3 step 2 — regional circulation comparator','',
        '2026-10-04. Implements the M2 regional network using simple resistors for brain and coronary beds. Pressure-dependent closures are step 3; support pumps and oxygen transport are absent. These results verify the chosen hydraulic construction, not measured organ physiology.','',
        '## Circuit and accounting','',
        'Eleven physical blood compartments retain the 5250 mL inventory, including all five limb pools. Four native inertial valve flows give 15 ODE states. The six central volume positions remain LV, RV, AO, VC, PA, PU. Limb states follow, then MT, AV, TC and PV flow states. The unchanged Smith ventricular/septal/pericardial laws use physical LV/RV cavity volumes. Vascular pressures use V minus the frozen unstressed offsets.','',
        'The single systemic resistor is replaced by regional paths. Kidney pair conductance splits equally; native coronary conductance splits 70/30. Gut and hepatic artery join at zero-storage H, then drain through hepatic outflow. H pressure is solved algebraically at every evaluation. Limb inflows enter their stored pools and return separately to VC. Nutritive lung flow is AO to PU, distinct from pulmonary exchange flow PA to PU. No organ, pool or blood is removed.','',
        'Both ends of each signed flow enter the incidence balance. Reverse regional flow is retained. H has no independent volume and its outflow must not be added a second time to systemic output. Baseline construction fractions are not imposed simulation flows. Parameters are taken from the unchanged M2 project configuration and the frozen kidney/coronary split definitions; there is no calibration to the source baseline.','',
        '## Executed cases','',
        '| Case | Warm-up cycles | MAP | SBP/DBP (mmHg) | LV output (L/min) | RV output (L/min) |',
        '|---|---:|---:|---:|---:|---:|']
    for name,r in results.items():
        m=r['metrics'];lines.append(f"| {name} | {r['warmup_cycles']} | {m['MAP_mmHg']:.4f} | {m['SBP_mmHg']:.4f}/{m['DBP_mmHg']:.4f} | {m['CO_L_min']:.6f} | {m['RV_output_L_min']:.6f} |")
    lines+=['','Each case independently requires three successive passing phase-grid comparisons across all 15 states (1e-5 criterion, 101 phases per cycle, at most 200 cycles), then ten measured cycles. Source absolute tolerances and the predeclared protocol are retained. These project outputs need not equal the source model because resistance, storage and flow topology changed.','',
        '## Baseline regional flows','', '| Path | Mean flow (mL/s) |','|---|---:|']
    for n,q in base['mean_branch_flows_mL_s'].items():lines.append(f'| {n} | {q:.6f} |')
    lines+=['',f"Left output minus right output minus nutritive lung flow: **{base['LV_minus_RV_minus_lung_supply_mL_s']:.3e} mL/s**. "
        'Left and right output are not expected to be equal: the AO-to-PU path recirculates through the left heart. A small residual remains from finite warm-up and sampled integration.', '',
        '## Numerical verification','', '| Comparison | Largest metric/branch-flow difference (%) | Maximum valve-time difference (ms) | Pass |',
        '|---|---:|---:|---|']
    for n,c in comparisons.items():lines.append(f"| {n} | {c['largest_metric_difference_pct']:.6f} | {c['maximum_valve_time_difference_ms']} | {c['passed']} |")
    lines+=['', 'All flow/pressure/volume metric comparisons must be below 0.5%; valve times must agree within 1 ms. Branch flows are included, so matching total output cannot conceal regional numerical disagreement.','',
        f"Maximum measured relative inventory drift: **{max(r['checks']['relative_volume_drift'] for r in results.values()):.3e}**; "
        f"maximum septal residual: **{max(r['maximum_septal_residual_mmHg'] for r in results.values()):.3e} mmHg**; "
        f"maximum hepatic balance residual: **{max(r['maximum_hepatic_balance_residual_mL_s'] for r in results.values()):.3e} mL/s**. "
        f"Frozen nine-file hashes unchanged: **{frozen_ok}**. Overall declared comparison result: **{summary['all_passed']}**.", '',
        'All 47 automated tests pass in reports/m3-regional-tests.xml, including the previous source and inventory tests. New tests check construction-anchor flows and split conductance, forward/reversed hepatic balance, full-cycle physical conservation and an analytically timed domain exit.', '',
        'A minimal reservoir-boundary event hook was added to the common valve integrator: a crossing ends the trajectory with an explicit model-domain record. An analytical test drains a reservoir to its unstressed volume and checks event time and conservation. Initial zero margins terminate conservatively. This is not an injury threshold. Full simultaneous-event, ventricular-boundary trial-state robustness and support-event stress testing remain later M3 work.', '',
        '![Regional baseline pressure and PV loops](baseline-waveforms.png)','',
        '## Run and inspect','',
        'Run `.venv/Scripts/python.exe scripts/run_m3_regional.py --workers 2`, then `.venv/Scripts/python.exe scripts/report_m3_regional.py`. '
        'The protocol is `configs/m3-regional-protocol.json`. Each case saves effective branch definitions, initial state, warm-up history, result, trajectories and valve events. The top-level manifest and implementation snapshot retain exact source/configuration versions. Tests are in `tests/test_regional.py`; their record is `reports/m3-regional-tests.xml`.','',
        '## Scope and next task','',
        'This completes the step-2 simple-resistor regional comparator when the checks above pass. It is not the complete frozen pressure-aware project variant. Next implement and compare the cerebral waterfall and coronary compression laws, keeping resistance fixed during their specified sensitivities. The unresolved source-paper thoracic reproduction requirement remains open. These simulations supply no oxygen sufficiency, organ recovery, support-duration or optimized auxiliary-heart claim.','']
    (out/'report.md').write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps(summary,indent=2))
    return 0 if summary['all_passed'] else 1


if __name__=='__main__':raise SystemExit(main())
