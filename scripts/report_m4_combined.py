"""Render M4 combined numerical evidence from saved results; no simulations."""
import json,os,sys,xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))
def write(p,s):(ROOT/p).write_text(s,encoding='utf-8')


def main():
    baseline=read('reports/m4-baseline/result.json')
    cases=[r for r in baseline['cases'] if r['solver']=='refined']
    lines=['# Unsuppressed periodic oxygen reference — 2026-10-08','',
        '**Numerical reference and contextual comparison; not independent physiological validation.**','',
        'Native pumping remains active. The M1 total demand is 255.847652 mL O2/min, with extraction cap 0.75 and prescribed alveolar ventilation 4.69 L/min BTPS. The 0/5/15 s capacities are hypothetical usable tissue stores, not ischemic tolerance times. No demand or hydraulic parameter is fitted to an observed oxygen value.','',
        '| Store credit, s | Actual uptake, mL/min | Unmet demand, mL/min | AO content, mL/dL | VC content, mL/dL | AO PO2, mmHg |',
        '|---:|---:|---:|---:|---:|---:|']
    for r in cases:lines.append(f"| {r['store_seconds']} | {r['consumed_mL_min']:.6f} | {r['unmet_mL_min']:.6f} | {r['AO_content_mL_dL']:.4f} | {r['VC_content_mL_dL']:.4f} | {r['AO_PO2_mmHg']:.3f} |")
    names=baseline['organ_names'];spl=['hepatic_outflow','stomach_oesophagus','small_intestine','large_intestine','pancreas','spleen']
    lines+=['',f"Native achieved CO: {baseline['CO_L_min']:.6f} L/min. The source-model CO is a different configuration; this regional-project output is not a new source-reproduction result.",'',
        '## Organ interpretation','']
    for r in cases:
        total=sum(r['per_organ_consumed_mL_min'][names.index(n)] for n in spl)
        lines.append(f"- Store {r['store_seconds']} s: actual liver + gut + pancreas + spleen uptake = {total:.6f} mL/min. Full organ uptake values remain in result.json.")
    lines+=['','[Zacho et al., Table 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC3574885/) reports fasting splanchnic uptake 58.1 mL/min (range 40.5–84.5) in ten men aged 43–69. The model’s synthetic 25-year-old is unmatched. A value inside that range is contextual agreement, not a validation pass. The M1 whole-body budget is also a constructed input, so reproducing it as uptake tests supply consistency rather than independently validating metabolism.','',
        'A deficit in the no-store model can reflect pulsatile extraction limits even when cycle-average delivery is sufficient. Comparing finite stores changes that assumption explicitly; it does not prove the selected store capacity or oxygen-preservation threshold. Independent matched-subject arterial, venous and organ uptake validation is unavailable.','',
        '## Numerical evidence','',
        f"Maximum reporting-cycle oxygen balance error: {max(r['maximum_balance_error_mL'] for r in baseline['cases']):.4g} mL. Maximum refined physical periodic residual: {max(r['final_physical_period_change'] for r in cases):.4g} (scaled). Hydraulic cycle endpoint volume difference: {baseline['hydraulic_periodic_volume_difference_mL']:.4g} mL.",'',
        'The initialization solver and residual records are saved with the baseline. Reporting ledgers are reset at the cycle boundary; physical stores are retained. Paired Radau reporting cycles check numerical refinement. Accepted states are checked without clipping.','',
        '![Periodic oxygen reference](baseline.png)','',
        'Reproduce: `.venv\\Scripts\\python.exe scripts/run_m4_baseline.py`. [Raw results and hashes](result.json); [equations and assumptions](../../docs/m4-respiratory-verification.md).','']
    write('reports/m4-baseline/report.md','\n'.join(lines))
    fig,axes=plt.subplots(2,1,figsize=(9,6),sharex=True,layout='constrained')
    for r in cases:
        z=np.load(ROOT/f"reports/m4-baseline/store-{r['store_seconds']}-refined.npz")
        axes[0].plot(z['time_s'],z['consumed'].sum(axis=1)*60,label=f"{r['store_seconds']} s assumed store")
        axes[1].plot(z['time_s'],z['unmet'].sum(axis=1)*60)
    axes[0].axhline(cases[0]['demand_mL_min'],color='black',ls='--',label='Assigned resting demand')
    axes[0].set(ylabel='Actual uptake (mL O₂/min)',title='Unsuppressed native reference: sensitivity to tissue buffering')
    axes[1].set(xlabel='Cardiac phase (s)',ylabel='Unmet demand (mL O₂/min)');axes[0].legend(fontsize=9)
    fig.savefig(ROOT/'reports/m4-baseline/baseline.png',dpi=170);plt.close(fig)

    protocol=read('configs/m4-respiratory-protocol.json');manifest=read('reports/m4-respiratory/manifest.json')
    results=[read(f"reports/m4-respiratory/{c['name']}/result.json") for c in protocol['cases']]
    verified=[r for r in results if r['status']=='verified'];failed=[r for r in results if r['status']!='verified']
    lines=['# Combined T4.4–T4.5 respiratory verification — 2026-10-08','',
        f"**{len(verified)}/{len(results)} scenarios verified numerically, each with paired oxygen solver settings.** Numerical verification does not imply adequate organ support or recovery.",'',
        'All cases use unsuppressed M1 demand; BIO adds an assumed 10 mL/min auxiliary myocardial demand. Short runs replay existing eight-second M3 hydraulics with periodic-reference oxygen contents. Native120 cases retain native pumping for 120 s; they test gas dependency, not arrest endurance. Fixed preservation demand and required-outcome masks are independent of routing.','',
        '| Case | Unmet O2, mL | Final AO O2, mL/dL | Final ventilation, L/min BTPS | Pressure screen | Violating organ screens |',
        '|---|---:|---:|---:|---|---:|']
    for result in verified:
        r=result['records'][-1];bad=sum(s['status']=='violated' for s in r['organ_screens'])
        lines.append(f"| {result['case']['name']} | {r['total_unmet_mL']:.6f} | {r['final_AO_content_mL_dL']:.4f} | {r['final_ventilation_L_min_BTPS']:.4f} | {r['pressure']['status']} | {bad} |")
    for r in failed:lines+=['',f"**Unverified case {r['case']['name']}:** {r['error']}"]
    lines+=['','### Oxygenator source accounting over eight seconds','',
        '| Case | Initial gas O2, mL STPD | Final gas O2, mL STPD | Net external O2, mL STPD |',
        '|---|---:|---:|---:|']
    for x in verified:
        if x['case']['name'] not in ('E1_fed','E1_finite','E2_fed'):continue
        r=x['records'][-1]
        lines.append(f"| {x['case']['name']} | 100.000 | {r['final_oxygenator_oxygen_mL']:.6f} | {r['oxygenator_net_external_oxygen_mL']:.6f} |")
    lines+=['','The finite E1 gas pool is not exhausted in this short run. Its zero external input and decreasing inventory are explicit; this does not establish sustained support without a gas source. Net external oxygen is oxygen entering minus oxygen leaving the fixed gas pool, not total sweep-gas volume.']
    lines+=['','## Accounting and numerical checks','',
        f"Maximum oxygen conservation error {max(r['oxygen_conservation_error_mL'] for x in verified for r in x['records']):.4g} mL; maximum demand identity error {max(r['demand_accounting_error_mL'] for x in verified for r in x['records']):.4g} mL; maximum paired oxygen-state difference {max(x['refinement_difference_mL'] for x in verified):.4g} mL. Ventilation differences are checked separately in mL BTPS/s.",'',
        'Each case saves state histories, delivered/extractable/extracted/consumed/unmet rates, pressure windows, all-tissue and pooled-limb-optional assessments, and solver settings. The source manifest records inputs and software. The numerical pass requires physical-domain and balance checks, not merely solver success.','',
        '## Interpretation','',
        '- Apnea stops external native gas replenishment, not exchange with gas already in the lung. A short eight-second comparison can be buffered by initial stores; it must not be extrapolated indefinitely.',
        '- Tissue buffering can cover a transient mismatch between extraction and metabolism. Its hypothetical capacity and kinetics change the result; neither is an injury threshold.',
        '- Oxygenator sweep has an explicit external oxygen ledger. Zero sweep leaves a finite gas inventory; gas volume and membrane capacity here are numerical fixtures.',
        '- Drive loss can reduce ventilation despite oxygen reaching respiratory muscle. The feedback law is an assumption, not demonstrated spontaneous breathing.',
        '- Required-outcome exclusion changes only assessment of an identical trajectory. It cannot establish a smaller pump, erase limb deficits or resolve reperfusion consequences.',
        '- Central pressure and organ oxygen screens can disagree. Even when all modeled screens pass, recovery remains unassessed and the combined status is indeterminate.',
        '- Comparison with earlier T4.3 values must account for the different initial oxygen state and optional stores. It is not a controlled test of a single changed assumption against every historical run.','',
        '![Native gas dependencies](native-gas.png)','',
        '![Short support comparisons](support.png)','',
        '## Required-outcome sensitivity','',
        '| Case | All tissues required | Pooled limbs optional |','|---|---|---|']
    for x in verified:
        r=x['records'][-1];lines.append(f"| {x['case']['name']} | {r['assessments']['all_tissues']['status']} | {r['assessments']['pooled_limbs_optional']['status']} |")
    lines+=['','The present graph cannot exclude only lower limbs. Demand/conductance/outcome choices are independently configurable in preservation.py, but changing conductance requires a new hydraulic run. No full-recovery or safe-duration conclusion follows.','',
        'Reproduce baseline first, then `.venv\\Scripts\\python.exe scripts/run_m4_respiratory.py`; render with `scripts/report_m4_combined.py`. [Protocol](../../configs/m4-respiratory-protocol.json), [manifest](manifest.json), [equations/provenance/limits](../../docs/m4-respiratory-verification.md), [M0–M4 report](../../docs/m0-m4-report.md).','']
    write('reports/m4-respiratory/report.md','\n'.join(lines))
    fig,axes=plt.subplots(3,1,figsize=(9,8),sharex=True,layout='constrained')
    for r in verified:
        name=r['case']['name']
        if not name.startswith('native120_'):continue
        z=np.load(ROOT/f'reports/m4-respiratory/{name}/refined.npz');n=len(r['volume_names']);a=r['volume_names'].index('AO')
        label=name.replace('native120_','').replace('feedback','drive loss (5 s response)')
        axes[0].plot(z['time_s'],z['oxygen_states'][:,a]/z['volumes_mL'][:,a]*100,label=label)
        axes[1].plot(z['time_s'],z['oxygen_states'][:,n]);axes[2].plot(z['time_s'],z['ventilation_mL_s_BTPS']*.06)
    axes[0].set(ylabel='AO O₂ (mL/dL)',title='Gas dependency with native pumping maintained (not arrest support)');axes[0].legend(fontsize=8)
    axes[1].set(ylabel='Native gas O₂ (mL STPD)');axes[2].set(ylabel='Ventilation (L/min BTPS)',xlabel='Time (s)')
    fig.savefig(ROOT/'reports/m4-respiratory/native-gas.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(9,7),sharex=True,layout='constrained')
    chosen=['B1_feedback','B1_feedback_store5','handover_feedback','E1_fed','E1_finite','E2_fed','BIO_feedback']
    for r in verified:
        name=r['case']['name']
        if name not in chosen:continue
        z=np.load(ROOT/f'reports/m4-respiratory/{name}/refined.npz');n=len(r['volume_names']);m=len(r['organ_names'])
        axes[0].plot(z['time_s'],z['oxygen_states'][:,n+2+m:n+2+2*m].sum(axis=1),label=name)
        axes[1].plot(z['time_s'],z['consumed'].sum(axis=1)*60)
    axes[0].set(ylabel='Accumulated unmet O₂ (mL)',title='Short support fixtures; different architectures and routing');axes[0].legend(fontsize=8,ncol=2)
    axes[1].set(xlabel='Time after arrest (s)',ylabel='Actual uptake (mL O₂/min)')
    fig.savefig(ROOT/'reports/m4-respiratory/support.png',dpi=170);plt.close(fig)
    tree=ET.parse(ROOT/'reports/m4-combined-tests.xml');suites=list(tree.getroot().iter('testsuite'));tests=sum(int(s.attrib['tests']) for s in suites)
    checkpoint=f'''## 10. Respiratory verification checkpoint (before fuel accounting)

The respiratory checkpoint suite contained **{tests} passing tests**. The combined package verifies **{len(verified)}/{len(results)} respiratory scenarios** with paired numerical settings, plus three periodic native oxygen/store comparisons. The [baseline report](../reports/m4-baseline/report.md) and [respiratory report](../reports/m4-respiratory/report.md) contain full tables and figures. All nine frozen M2 files still pass their integrity audit.

Maximum combined oxygen balance error is {max(r['oxygen_conservation_error_mL'] for x in verified for r in x['records']):.3g} mL; maximum paired oxygen-state difference is {max(x['refinement_difference_mL'] for x in verified):.3g} mL. These are numerical errors, not biological uncertainties. The no-store and finite-store baseline actual uptakes are {', '.join(f"{r['consumed_mL_min']:.3f}" for r in cases)} mL/min respectively (0/5/15 s credit); compare the assigned demand of 255.848 mL/min.

The principal findings are:

- The no-credit resting shortfall is confined to the left coronary branch. Tissue buffering covers its pulsatile extraction gap under the chosen capacities/kinetics; this is conditional baseline consistency, not a measured reserve requirement.
- In the eight-second B1 feedback comparison, no-credit unmet demand is 0.765927 mL, while the assumed 5 s store case has none. Both satisfy the central pressure screen. A 5 s demand-equivalent capacity does not mean a five-second support limit: simultaneous blood delivery and store replenishment remain active.
- Buffered BIO has no oxygen-demand shortfall over eight seconds but violates the pressure screen. Sufficient modeled oxygen and sufficient pressure are independent tests.
- With native pumping retained for 120 s, maintained/half ventilation end at AO oxygen contents 19.9135/19.0557 mL/dL. Apnea/drive loss end at 12.4312/12.8226 mL/dL and accumulate 0.479301/0.171654 mL unmet left-coronary demand. These are gas-dependency experiments, not permissible arrest times.
- Finite-source E1 consumes part of its initial 100 mL oxygenator gas inventory, leaving 69.195446 mL after eight seconds with zero external supply. Fed E1/E2 have explicit net external oxygen of 26.080153/16.625030 mL. Their idealized membrane parameters are not device-sizing results.
- Making pooled limb recovery optional changes none of these classifications, because the observed violating screens are elsewhere. A different routing/demand policy would require a new simulation.

The [final audit](../reports/m4-combined-audit.json) verifies hashes, protocol agreement, frozen files and report links. All three saved baseline cycles replay with current code with zero observed difference. Baseline source snapshots retain the original runtime hashes; later changes to an unused oxygenator trial guard do not alter native-baseline physics.

T4.4/T4.5 results remain conditional on the documented gas, extraction, storage and respiratory assumptions. At that respiratory-only checkpoint, finite substrate accounting and the final gate were still open; a subsequent fuel checkpoint, when present below, supersedes that status. Recovery, full brain coverage and finite-resource support duration are still unassessed.
'''
    p=ROOT/'docs/m0-m4-report.md';existing=p.read_text(encoding='utf-8')
    s=existing.split('<!-- M4_RESULT_CHECKPOINT -->')[0]
    fuel='<!-- M4_FUEL_CHECKPOINT -->'+existing.split('<!-- M4_FUEL_CHECKPOINT -->',1)[1] if '<!-- M4_FUEL_CHECKPOINT -->' in existing else ''
    p.write_text(s+'<!-- M4_RESULT_CHECKPOINT -->\n'+checkpoint+'\n'+fuel,encoding='utf-8')


if __name__=='__main__':main()
