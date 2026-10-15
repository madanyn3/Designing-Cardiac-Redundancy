"""Render verified T5.1 hydraulic comparisons without declaring metabolic adequacy."""
from pathlib import Path
import json,os
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'reports/m5-dynamic'
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    r=json.loads((OUT/'verification.json').read_text(encoding='utf-8'))
    if r['status']!='verified':raise RuntimeError('Do not render failed comparisons as verified')
    rows=['# T5.1 dynamic bounded-pump characterization','',
        '2026-10-09. **Ten hydraulic scenarios, twenty paired integrations verified. T5.1 and M5 remain open.** [Protocol](../../configs/m5-dynamic-protocol.json), [raw results/provenance](verification.json), [M5 plan](../../docs/m5-plan.md).','',
        '## Experiment and scope','',
        'Native ventricular activation stops at the saved reference phase. Right/left heads ramp to 20/100 mmHg over one second. Active supply branches remain open; preservation conductances remain at 0.25 throughout. No regional reopening or native restart occurs. All blood and native tissue routes remain present. These controls are numerical fixtures, not adequate preservation prescriptions.','',
        'The ungated stage is the M3 comparator. The gated stage caps available flow capacity using an assumed 0–4 mmHg inlet-availability ramp. The gauge comparator uses zero external pressure for both inlets; the transmural comparator uses VC external pressure 0 and PU external pressure −4 mmHg, consistent with the frozen source pressure convention. Negative PU gauge pressure need not imply an empty pulmonary venous reservoir. Gate thresholds are not measured human limits.','',
        'The eight initial cases last 12 seconds or until a localized domain boundary. Two cases extend to 60 seconds. Depleted cases initially transfer blood from PU to VC, leaving PU only 0.1 mL above its model boundary; right support is disabled. This synthetic inventory redistribution preserves 5,250 mL but is not a normal initialization.','',
        '## Dynamic outcomes','',
        'Tail averages use the final two seconds, excluding the first activation second. An early-domain-exit row therefore has a shorter horizon than completed rows. A domain boundary is a model-validity limit, not a tissue-injury time.','',
        '| Case | Observed s | Tail mean AO mmHg | Right support L/min | Left support L/min | PA+PU inventory change mL | Outcome |',
        '|---|---:|---:|---:|---:|---:|---|']
    for name,s in r['scenarios'].items():
        f=lambda value:'—' if value is None else f'{value:.3f}'
        rows.append(f"| {name} | {s['observed_s']:.5f} | {f(s['tail_mean_AO_mmHg'])} | {f(None if s['tail_right_mL_s'] is None else s['tail_right_mL_s']*.06)} | {f(None if s['tail_left_mL_s'] is None else s['tail_left_mL_s']*.06)} | {s['pulmonary_volume_change_mL']:.3f} | {s['domain_node']+' domain limit' if s['domain_node'] else 'completed'} |")
    rows+=['','The ungated run reaches the VC boundary at approximately 10.018 s despite a tail aortic mean near 98 mmHg. Its support law can continue drawing from a poorly filled inlet. The transmural-gated reference completes 60 s, with a final-window mean near 96.50 mmHg and left support near 2.841 L/min. This is a tested hydraulic operating point, not a derived minimum pump.','',
        'Weak left support and the two-small-stage case end near 46 mmHg mean aortic pressure. Weak right support initially preserves aortic pressure by changing blood distribution: PA+PU inventory falls about 407 mL by 12 s and 638 mL by 60 s. The latter still has mean aortic pressure near 89.25 mmHg; pressure alone does not reveal the changed circulation.','',
        '![Pressure, inventory and flow](dynamic.png)','',
        '## Why right and left support flows need not be equal','',
        '`d(V_PA + V_PU)/dt = Q_support_right + Q_PV + Q_lungs − Q_support_left − Q_MT`.','',
        'Pulmonary gas-exchange flow is internal to this pool and cancels. Q_lungs is the separate nutritive lung return into PU; Q_PV and Q_MT are native pulmonary-valve and mitral flows. The native heart has no active contraction but its valves and passive chambers remain connected. It can therefore transmit pressure-driven flow.','',
        'In the 60 s transmural reference, right/left support is about 2.702/2.841 L/min, with approximately 0.139 L/min nutritive lung return making up the difference. In the weak-right 60 s case, right support is 1.200 L/min, native passive pulmonary-valve flow is about 1.223 L/min, and nutritive return about 0.134 L/min; together these approximately supply the 2.557 L/min left output. This is not evidence that an unsupported native right ventricle remains physiologically acceptable. Chamber loading, injury and recovery are unassessed.','',
        'The depleted right-disabled cases also retain passive native pulmonary-valve flow. Their completion cannot be interpreted as proof that the pulmonary circulation needs no driving pressure or that right support is unnecessary. No valve isolation was silently added.','',
        '## Regional supply and pressure interpretation','',
        '| Case | Brain mL/min | Respiratory muscle mL/min | Hepatic arterial mL/min | Portal contributors mL/min |',
        '|---|---:|---:|---:|---:|']
    for name in ('transmural_60','right_weak_60','left_weak','both_small'):
        q=r['scenarios'][name]['regional_tail_flow_mL_min']
        portal=sum(q[k] for k in ('stomach_oesophagus','small_intestine','large_intestine','pancreas','spleen'))
        rows.append(f"| {name} | {q['brain']:.3f} | {q['respiratory_muscle']:.3f} | {q['liver_arterial']:.3f} | {portal:.3f} |")
    rows+=['','These achieved flows do not certify organ metabolism. The raw results retain the historical one-second 65/90/50 mean/maximum/minimum research screen for continuity. For these continuous-flow comparators, the window maximum is not a beat-defined systolic pressure; a maximum-only failure must not be interpreted as an independently established physiological failure. No minimum venous filling/organ recovery constraint is newly inferred.','',
        '## Numerical verification','',
        '| Case | Maximum scaled difference % | Valve timing difference ms | Paired check |','|---|---:|---:|---|']
    for name,c in r['comparisons'].items():rows.append(f"| {name} | {100*max(c['scaled_errors'].values()):.6g} | {1000*c['valve_time_difference_s']:.6g} | {c['passed']} |")
    results=[json.loads(p.read_text()) for p in OUT.glob('*/result.json')]
    rows+=['',f"Maximum blood-inventory error: {max(x['blood_conservation_error_mL'] for x in results):.3g} mL. Maximum pointwise pulmonary balance error: {max(x['pulmonary_balance_identity_error_mL_s'] for x in results):.3g} mL/s. State continuity, hepatic mixing and septal roots pass. Paired criteria retain <0.5% coordinate-scaled differences, matched valve sequences/timing <1 ms and matched domain node/time <1 ms. All nine frozen M2 files are checked before execution.",'',
        '## Remaining T5.1 work','',
        '1. Expose and verify independent Smith-shaped auxiliary chamber/rate/contractility settings. A cavity similarity scale is not myocardial mass.',
        '2. Replay candidate hydraulics through the M4 oxygen/glucose model and report violated/unassessed constraints. Retain L-M4-F01 (no tissue glucose buffers), all preservation tissues, finite hepatic release and explicit additional biological-pump demand.',
        '3. Define the conditional capability envelope across pressure, regional supply, inventory and design variables. Hydraulic completion alone is not a feasible biological design.',
        '', 'T5.2 input-energy/efficiency/self-demand, T5.3 mass arithmetic, T5.4 matched architectures and T5.5 unmodeled biological costs remain open. The inherited source discrepancy R-M3-01 is unchanged.','',
        'Reproduce with `.venv\\Scripts\\python.exe scripts/run_m5_dynamic.py` and `.venv\\Scripts\\python.exe scripts/report_m5_dynamic.py`. `--report-only` on the runner checks and consolidates saved matching runs; it does not integrate them.']
    (OUT/'report.md').write_text('\n'.join(rows)+'\n',encoding='utf-8')
    fig,axes=plt.subplots(3,1,figsize=(10,10),sharex=True,layout='constrained')
    for name in ('ungated','gauge_gate','transmural_60','right_weak_60','left_weak','both_small'):
        result=json.loads((OUT/(name+'-refined')/'result.json').read_text());vn=result['volume_names'];fn=result['flow_names']
        ds=[np.load(OUT/(name+'-refined')/f"segment-{s['index']}.npz") for s in result['segments']]
        t=np.concatenate([d['time_s'] for d in ds]);p=np.concatenate([d['pressures_mmHg'] for d in ds]);v=np.concatenate([d['states'].T for d in ds]);q=np.concatenate([d['flows_mL_s'] for d in ds])
        pool=v[:,vn.index('PA')]+v[:,vn.index('PU')]
        label=name.replace('_',' ')
        axes[0].plot(t,p[:,vn.index('AO')],label=label)
        axes[1].plot(t,pool-pool[0]);axes[2].plot(t,q[:,fn.index('support_left')]*.06)
    axes[0].set(ylabel='Aortic pressure (mmHg)',title='T5.1 hydraulic controls — no metabolic/recovery adequacy inferred');axes[0].legend(fontsize=8,ncol=2)
    axes[1].set_ylabel('PA + PU inventory change (mL)');axes[2].set(ylabel='Left auxiliary flow (L/min)',xlabel='Time after loss of native activation (s)')
    for a in axes:a.grid(alpha=.2)
    fig.savefig(OUT/'dynamic.png',dpi=160);plt.close(fig)
    print(OUT/'report.md')


if __name__=='__main__':main()
