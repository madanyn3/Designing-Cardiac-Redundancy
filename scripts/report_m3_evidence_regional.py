"""Report the amended regional experiment from its saved run manifest."""
import hashlib
import json
import argparse
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reports/m3-evidence-regional'


def main():
    global OUT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',default='reports/m3-evidence-regional')
    args=parser.parse_args()
    OUT=ROOT/args.out
    manifest=json.loads((OUT/'manifest.json').read_text())
    protocol=manifest['protocol']
    results={c['name']:json.loads((OUT/c['name']/'result.json').read_text()) for c in protocol['cases']}
    if any(r['status']!='passed_numerical_checks' for r in results.values()):
        raise RuntimeError('Run incomplete or failed; inspect individual results')
    base=results['baseline']; comparisons={}
    for name,r in results.items():
        if name=='baseline': continue
        differences=[100*abs(v-base['metrics'][k])/abs(base['metrics'][k]) for k,v in r['metrics'].items()]
        differences += [100*abs(v-base['mean_branch_flows_mL_s'][k])/abs(base['mean_branch_flows_mL_s'][k]) for k,v in r['mean_branch_flows_mL_s'].items()]
        def events(n):
            start=results[n]['measurement_start_s']
            return [(e['valve'],e['action'],e['time_s']-start) for e in json.loads((OUT/n/'events.json').read_text()) if e['time_s']<start+.75]
        a,b=events('baseline'),events(name)
        same=[e[:2] for e in a]==[e[:2] for e in b]
        dt=max(abs(x[2]-y[2])*1000 for x,y in zip(a,b)) if same and a else None
        comparisons[name]={'maximum_difference_pct':max(differences),'event_difference_ms':dt,
            'passed':dt is not None and max(differences)<protocol['comparison_metric_limit_pct'] and dt<protocol['comparison_event_time_limit_ms']}
    frozen=json.loads((ROOT/'data/m2-freeze-manifest.json').read_text())['sha256']
    unchanged=all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in frozen.items())
    verification={'comparisons':comparisons,'frozen_unchanged':unchanged,
                  'passed':unchanged and all(r['passed'] for r in comparisons.values())}
    (OUT/'verification.json').write_text(json.dumps(verification,indent=2)+'\n')
    project=json.loads((OUT/'implementation-snapshot'/protocol['project_config']).read_text())
    lines=['# M3 evidence-amended regional baseline','',
        'Amendment started 2026-10-05; verification completed 2026-10-06. Referenced baseline estimates with explicitly unresolved anatomical partitions. Fixed brain/coronary resistors, unchanged storage and pressure anchors, no oxygen transport or auxiliary pump. Numerical consistency is not physiological validation.','',
        '## Implemented amendment','',
        'Five portal branches replace gut: stomach_oesophagus, small_intestine, large_intestine, pancreas and spleen. All join hepatic arterial flow at zero-storage H. The algebraic pressure solve includes every incident branch, including reversed flows. Source mechanics, 5250 mL inventory, kidney 50/50 and coronary 70/30 conductance partitions remain unchanged.','',
        'Parent fractions follow Leggett–Williams; tissue coefficients use the qualified transfers documented in [regional candidates](../../docs/m3-regional-flow-candidates.md) and [digestive split](../../docs/m3-splanchnic-split.md). Limb bone and limb-other retain unsupported historical allocations explicitly. Remaining flow is the sum of named non-limb tissue complements. M1 sub-organ demand reconciliation and accessory-organ mapping remain open.','',
        '## Baseline results','', '| Metric | Value |','|---|---:|']
    for k,v in base['metrics'].items(): lines.append(f'| {k} | {v:.6f} |')
    lines += ['', '## Construction targets versus solved flow','', '| Branch | Target mL/min | Solved mL/min | Basis |','|---|---:|---:|---|']
    q=base['mean_branch_flows_mL_s']
    for row in project['branches']:
        n=row['name']
        actual=q['coronary_left']+q['coronary_right'] if n=='heart' else q['kidney_protected']+q['kidney_other'] if n=='kidneys' else q[n]
        lines.append(f"| {n} | {row['target_flow_mL_s']*60:.3f} | {actual*60:.3f} | {row['class']} |")
    lines += ['', 'Targets construct resistance; they are not imposed during simulation. Return flows are serial and not additional cardiac output.','',
              '## Verification','', '| Case | Warm-up cycles | Max difference % | Valve difference ms | Pass |','|---|---:|---:|---:|---|']
    for n,c in comparisons.items(): lines.append(f"| {n} | {results[n]['warmup_cycles']} | {c['maximum_difference_pct']:.6g} | {c['event_difference_ms']} | {c['passed']} |")
    lines += ['',f"Baseline warm-up: {base['warmup_cycles']} cycles. Each run measures ten cycles. Limits: 0.5% metrics/branch flows and 1 ms valve times.",
        f"Maximum relative volume drift: {max(r['checks']['relative_volume_drift'] for r in results.values()):.3e}.",
        f"Maximum hepatic balance residual: {max(r['maximum_hepatic_balance_residual_mL_s'] for r in results.values()):.3e} mL/s.",
        f"Historical M2 hashes unchanged: {unchanged}. Overall numerical comparison passed: {verification['passed']}.",'',
        'Run: `.venv/Scripts/python.exe scripts/run_m3_regional.py --protocol configs/m3-evidence-regional-protocol.json --out reports/m3-evidence-regional --workers 2`, then `scripts/report_m3_evidence_regional.py` with the same Python. Build inputs with `scripts/audit_m3_flow_candidates.py` and `scripts/build_m3_evidence_regional.py`. The manifest and implementation-snapshot contain exact simulation inputs/code. Automated tests are recorded in `reports/m3-evidence-regional-tests.xml`.','',
        'The old regional report remains a historical run. The amendment does not close oxygen adequacy, organ recovery, selective routing feasibility, hepatic buffering, or the outstanding source-paper reproduction discrepancy.','']
    if protocol.get('pressure_aware'):
        lines[0]='# M3 pressure-aware regional baseline'
        lines[2]='2026-10-06. Nominal cerebral waterfall and coronary compression laws on the evidence-amended allocation. No oxygen transport or auxiliary pump. Numerical consistency is not physiological validation.'
        lines += pressure_report(base,project,manifest)
        lines=[line.replace('scripts/report_m3_evidence_regional.py`','scripts/report_m3_evidence_regional.py --out reports/m3-pressure-regional`')
                    .replace('--protocol configs/m3-evidence-regional-protocol.json --out reports/m3-evidence-regional',
                             '--protocol configs/m3-pressure-regional-protocol.json --out reports/m3-pressure-regional')
                    .replace('reports/m3-evidence-regional-tests.xml','reports/m3-pressure-regional-tests.xml') for line in lines]
    (OUT/'report.md').write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps(verification,indent=2))
    return 0 if verification['passed'] else 1


def pressure_report(base,project,manifest):
    """Compare nominal closure outputs with the recorded fixed-R comparator.

    This is a model-form comparison, distinct from numerical convergence.
    Neither difference is a clinical accuracy score.
    """
    import os
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.matplotlib'))
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ref_path=ROOT/'reports/m3-evidence-regional/baseline/result.json'
    ref=json.loads(ref_path.read_text())
    effective=json.loads((OUT/'baseline/effective-model.json').read_text())['closure']
    contrast={}
    lines=['','## Pressure-dependent laws and model-form comparison','',
        'Implements the M2 signed quasi-static waterfall and ventricular/pericardial compression hypotheses. These are exploratory closures; there is no active cerebral autoregulation, intramyocardial capacitance or validated minimum-flow criterion. Both flow endpoints remain AO and VC, so donor reversal and conservation are unchanged.','',
        f"Cerebral nominal/effective resistance: {effective['brain_nominal_resistance_mmHg_s_mL']:.9f}/{effective['brain_effective_resistance_mmHg_s_mL']:.9f} mmHg s/mL. Nominal Pc=20 mmHg; calibrated from the amended 11.2 mL/s construction flow. Coronary kappa=0.5. Coronary resistances are not recalibrated to recover a desired output.",'',
        '| Quantity | Fixed resistor | Pressure-aware | Change % |','|---|---:|---:|---:|']
    for k in ('MAP_mmHg','SBP_mmHg','DBP_mmHg','CO_L_min'):
        a,b=ref['metrics'][k],base['metrics'][k]; contrast[k]=100*(b-a)/a
        lines.append(f'| {k} | {a:.6f} | {b:.6f} | {contrast[k]:.4f} |')
    for k in ('brain','coronary_left','coronary_right'):
        a,b=ref['mean_branch_flows_mL_s'][k]*60,base['mean_branch_flows_mL_s'][k]*60
        contrast[k]=100*(b-a)/a
        lines.append(f'| {k} mL/min | {a:.6f} | {b:.6f} | {contrast[k]:.4f} |')
    (OUT/'model-form-comparison.json').write_text(json.dumps({'signed_change_pct':contrast,
        'reference_result_sha256':hashlib.sha256(ref_path.read_bytes()).hexdigest(),
        'scope':'nominal laws only; not a sweep of all coupled physiological sensitivities'},indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for directory,label,style in ((ROOT/'reports/m3-evidence-regional','Fixed R','--'),(OUT,'Pressure-aware','-')):
        tr=np.load(directory/'baseline/trajectory.npz')
        result=json.loads((directory/'baseline/result.json').read_text())
        names=list(result['mean_branch_flows_mL_s'])
        t=tr['time_s']; mask=t>=t[-1]-.75; phase=t[mask]-(t[-1]-.75)
        axes[0].plot(phase,tr['flows_mL_s'][mask,names.index('brain')]*60,style,label=label)
        for organ,color in [('coronary_left','C0'),('coronary_right','C1')]:
            axes[1].plot(phase,tr['flows_mL_s'][mask,names.index(organ)]*60,style,color=color,label=label+' '+organ)
    for ax in axes: ax.set(xlabel='Cycle phase (s)',ylabel='Instantaneous flow (mL/min)'); ax.legend(fontsize=8)
    fig.savefig(OUT/'pressure-closure-flows.png',dpi=160); plt.close(fig)
    lines += ['', '![Pressure-closure flow comparison](pressure-closure-flows.png)','',
        'Unit tests cover threshold crossing against an analytic integrated volume, reverse-flow symmetry, hydraulic passivity, disabled-closure equivalence, revised calibration, fixed resistance across the ICP/tone grid, and active/passive conservation. The four coupled runs establish nominal solver/initialization reliability. Coupled one-factor ICP/tone/resistance/kappa experiments are recorded separately in [the sensitivity report](../m3-pressure-sensitivity/report.md); the nominal protocol itself does not execute them.','',
        'Next: the no-support arrest and auxiliary-support interfaces, plus their event/domain tests. Whole-brain regional adequacy and coronary viability remain unpredicted.']
    return lines


if __name__=='__main__': raise SystemExit(main())
