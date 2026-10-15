"""Normal-heart comparison and explicitly assumed auxiliary metabolic demand."""
from pathlib import Path
import json,sys,hashlib,os,xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=ROOT/'reports/m5-normal-comparison'
def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))
def main():
    normal=read('reports/m5-normal-comparison/verification.json');assert normal['status']=='verified'
    cycle=normal['runs'][-1]['metrics']['cycle'];matched=normal['runs'][-1]['metrics']['matched_8s']
    baseline=sum(s['work']['chamber_J'] for s in cycle['sides'].values())/.75
    native8=sum(s['work']['chamber_J'] for s in matched['sides'].values())
    scenarios=read('reports/m5-self-demand/scenarios.json');native=scenarios['normal_heart_O2_mL_min'];budget=read('data/m1-demand-budget.json')['scenarios']['unsuppressed']['whole_body_mL_min']
    mechanics=read('reports/m5-biological-energy/verification.json')['results'];hydraulics=read('reports/m5-t51/verification.json')['summaries'];h={x['case']:x for x in hydraulics}
    lines=['# T5.2: auxiliary heart versus the normally beating heart','','2026-10-10. Comparison of project model outputs and declared metabolic assumptions; not a patient-specific or empirically calibrated auxiliary design. [Equations, sources and limits](../../docs/m5-t52-normal-comparison.md).','','## Normal reference', '',f"The normal regional model beats at 80/min with CO {cycle['CO_L_min']:.4f} L/min and mean AO {cycle['mean_AO_mmHg']:.2f} mmHg (maximum/minimum {cycle['maximum_AO_mmHg']:.2f}/{cycle['minimum_AO_mmHg']:.2f}). Combined ventricular work is {baseline:.4f} W averaged over a complete 0.75 s cycle. M1 supplies the independent constructed myocardial demand of {native:.4f} mL O2/min.",'','## Mechanical and demand comparison','','Normal mechanical power is a full-cycle mean; BIO mechanical power is an eight-second transient mean including activation. BIO mean AO and output in this table also use the full eight seconds, unlike T5.1\'s final-two-second table. Oxygen estimates are the new assumed workload/activation/maintenance closure, not the old fixed 10 mL/min.','','| Model | Rate /min | Mean AO mmHg | Left output L/min | Chamber W | % normal W | O2 target mL/min | % normal O2 target |','|---|---:|---:|---:|---:|---:|---:|---:|',f"| Normal resting reference | 80 | {cycle['mean_AO_mmHg']:.2f} | {cycle['CO_L_min']:.3f} | {baseline:.4f} | 100.0 | {native:.3f} | 100.0 |"]
    summaries=[];names=['Normal'];powers=[baseline];demands=[native]
    for case in scenarios['scenarios']:
        if not case['name'].endswith('_reference'):continue
        name=case['trajectory'];pars=read(f'reports/m5-biological/{name}/result.json')['parameters'];a=np.load(ROOT/f'reports/m5-biological-energy/{name}-refined.npz');t=a['time_s'];y=a['states'];p=a['pressures_mmHg']
        power=sum(s['work']['chamber_J'] for s in mechanics[name]['runs'][-1]['sides'].values())/8
        ao=float(np.trapezoid(p[:,2],t)/8);co=float(np.trapezoid(y[23],t)/8*.06);oxygen=case['auxiliary_mL_min']
        lines.append(f"| BIO {name} | {60/pars['period_s']:.0f} | {ao:.2f} | {co:.3f} | {power:.4f} | {100*power/baseline:.1f} | {oxygen:.3f} | {100*oxygen/native:.1f} |")
        summaries.append(dict(case=name,chamber_power_W=power,fraction_normal_power=power/baseline,eight_second_work_J=8*power,fraction_matched_native_8s_work=8*power/native8,mean_AO_mmHg=ao,left_output_L_min=co,auxiliary_O2_mL_min=oxygen,total_two_hearts_O2_mL_min=native+oxygen,total_body_O2_mL_min=budget+oxygen))
        names.append(name);powers.append(power);demands.append(oxygen)
    lines+=['',f"For a strictly matched eight-second normal interval beginning at phase 0.375 s, actual native work is {native8:.4f} J and measured average output {matched['CO_L_min']:.4f} L/min. The cycle-average eight-second energy equivalent is {baseline*8:.4f} J. Their difference reflects noninteger-cycle boundaries; do not use them interchangeably.",'','Lower work occurs with lower output and generally lower pressure. It does not establish superior efficiency or adequate preservation. BIO also uses a different, selectively restricted circulation and auxiliary blood allocation, so this is a scale comparison rather than a matched architecture superiority test.','','## Oxygen demand and the retained native heart','','The native heart remains off in BIO runs but retains its unchanged unsuppressed M1 target in these experiments. Therefore auxiliary demand is an additional cost. For the default BIO reference, total cardiac target is '+f"{summaries[0]['total_two_hearts_O2_mL_min']:.3f} mL/min versus {native:.3f} for the normal heart alone. Whole-body target is {summaries[0]['total_body_O2_mL_min']:.3f} versus {budget:.3f} mL/min. Lower native-preservation demand is a separate explicit scenario, not an automatic consequence of stopping contraction.",'','The split between maintenance, activation and workload is assumed 15/25/60% at the normal reference. Tissue mass ratio is an independent assumed 0.25; it is not derived from cavity scale. Increased rate and contractility increase the activation term; workload cost follows the positive transmural-work proxy. This is not PVA and does not quantify all isovolumic costs. No calibrated chemical efficiency or physical mass has been inferred.','','## Joint supply checks and sensitivities','','The historical 10 mL/min auxiliary target is replaced, and absolute auxiliary oxygen stores are held fixed. Original-demand, oxygen, glucose and stoichiometry ledgers are checked with paired replays. Once a sampled auxiliary shortfall occurs, continued prescribed pumping is a diagnostic counterfactual.','','| Case | Auxiliary target mL/min | Total unmet O2 in 8 s, mL | Auxiliary unmet O2, mL | First auxiliary shortfall bracket, s |','|---|---:|---:|---:|---|']
    replay=[]
    for case in scenarios['scenarios']:
        r=read(f"reports/m5-self-demand/{case['name']}/result.json");assert r['status']=='verified' and r['case']==case
        v=r['records'][-1];aux=sum(x for key,x in v['per_organ_unmet_oxygen_mL'].items() if key.startswith('aux_'))
        bracket=v['first_sampled_auxiliary_shortfall_bracket_s']
        bracket_label='none observed' if bracket is None else f'{bracket[0]:.3f}–{bracket[1]:.3f}'
        lines.append(f"| {case['name']} | {case['auxiliary_mL_min']:.3f} | {v['unmet_oxygen_mL']:.5f} | {aux:.5f} | {bracket_label} |")
        replay.append(dict(case=case['name'],unmet_oxygen_mL=v['unmet_oxygen_mL'],auxiliary_unmet_oxygen_mL=aux,first_auxiliary_shortfall_bracket_s=v['first_sampled_auxiliary_shortfall_bracket_s'],driver_status=v['prescribed_driver_status']))
    lines+=['','Illustrative sensitivity values are neither confidence intervals nor established biological limits. L-M4-F01 remains: no tissue glucose buffer and no validated joint-fuel physiological resting baseline. A shortfall is not a tissue injury time or proven duration of viable support.','','## Verification and scope','','191 tests pass. The normal reference has paired ODE/output-grid work checks; all seven new demand scenarios have paired oxygen/glucose replays. Historical M2/M4/T5.1 and the previous BIO energy results remain unchanged. The archive includes equations, assumptions, intermediate files and provenance hashes.','','T5.2 remains open for a better identified metabolic relation, PVA/septal energy interpretation, maintenance/mass constraints and justified supply-limited pumping. No case is declared physiologically feasible.','','![Normal and auxiliary comparison](comparison.png)','','## Reproduce','','```powershell','.\\.venv\\Scripts\\python.exe scripts/run_m5_normal_comparison.py','.\\.venv\\Scripts\\python.exe scripts/derive_m5_self_demand.py','.\\.venv\\Scripts\\python.exe scripts/run_m5_self_demand.py','.\\.venv\\Scripts\\python.exe scripts/report_m5_normal_comparison.py','```']
    fig,ax=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    colors=['#355c7d']+['#d58b40']*4
    ax[0].bar(names,powers,color=colors);ax[0].set(ylabel='Mean chamber boundary power (W)',title='Normal full cycle / BIO eight-second transient')
    ax[1].bar(names,demands,color=colors);ax[1].set(ylabel='Requested oxygen (mL/min)',title='M1 normal target / assumed BIO demand')
    for a in ax:a.tick_params(axis='x',rotation=25)
    fig.savefig(OUT/'comparison.png',dpi=160);plt.close(fig)
    (OUT/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    (OUT/'comparison.json').write_text(json.dumps(dict(status='verified_at_declared_scope',normal_periodic_power_W=baseline,normal_matched_work_J=native8,normal_heart_demand_mL_min=native,summaries=summaries,replays=replay),indent=2)+'\n')
    suites=list(ET.parse(ROOT/'reports/m5-t52-demand-tests.xml').getroot().iter('testsuite'));assert sum(int(s.get('tests',0)) for s in suites)==191 and all(int(s.get('errors',0))+int(s.get('failures',0))==0 for s in suites)
    for manifest in ('reports/m5-t51/manifest.json','reports/m5-biological-energy/audit.json'):
        for p,h in read(manifest)['hashes'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    files={Path(__file__),ROOT/'scripts/run_m5_normal_comparison.py',ROOT/'scripts/derive_m5_self_demand.py',ROOT/'scripts/run_m5_self_demand.py',ROOT/'configs/m5-myocardial-demand.json',ROOT/'reports/m5-t52-demand-tests.xml',ROOT/'data/m1-demand-budget.json',ROOT/'configs/m4-organ-mapping.json',ROOT/'configs/m4-oxygen-reference.json',ROOT/'configs/m4-fuel-reference.json',ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml',*list((ROOT/'src/cardiac_redundancy').glob('*.py')),*list((ROOT/'tests').glob('test_myocardial_demand.py'))}
    for directory in ('reports/m5-self-demand','reports/m5-biological-energy','reports/m5-normal-comparison'):
        files.update(p for p in (ROOT/directory).rglob('*') if p.is_file() and p.name!='manifest.json')
    for case in scenarios['scenarios']:
        files.update(ROOT/p for p in read(f"reports/m5-self-demand/{case['name']}/result.json")['inputs'])
    (OUT/'manifest.json').write_text(json.dumps(dict(status='verified',tests=191,paired_metabolic_cases=7,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}),indent=2)+'\n')
if __name__=='__main__':main()
