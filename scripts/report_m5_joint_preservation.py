"""Verify archived joint-grid provenance and summarize conditional outcomes."""
from pathlib import Path
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'reports/m5-joint-preservation'
def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))
def main():
    cfg=read('configs/m5-joint-preservation.json');records=[]
    normal=read('reports/m5-normal-comparison/verification.json')['runs'][-1]['metrics']['cycle']
    old=read('reports/m5-biological-energy/verification.json')['results']
    olddemands={c['trajectory']:c for c in read('reports/m5-self-demand/scenarios.json')['scenarios'] if c['name'] in ['BIO_reference_reference','BIO_faster_reference']}
    lines=['# Joint BIO and preservation-conductance exploration','','2026-10-14. Four paired hydraulic/energy and four paired metabolic cases. Fixed active/preservation demand profile; auxiliary demand is recalculated for each workload. [Protocol and limitations](../../docs/m5-joint-preservation.md).','','## Pressure, flow and pump cost','','Pressure/flow use the final two seconds; work is over eight seconds. All fractions are exploratory assumptions. The normal full-cycle chamber-power reference is 1.5486 W and its myocardial oxygen target is 20.7904 mL/min. These are reference denominators, not a phase-matched transient.','','| Case | Mean/min/max AO mmHg | Left pump L/min | Chamber power W | Auxiliary O2 mL/min | Pressure screen |','|---|---|---:|---:|---:|---|']
    for c in cfg['cases']:
        h=read(f"reports/m5-joint-preservation/{c['name']}/hydraulics.json");r=read(f"reports/m5-joint-preservation/{c['name']}/result.json")
        assert h['status']==r['status']=='verified' and h['case']==c
        assert all(r['case'][k]==v for k,v in c.items())
        assert h['volume_names']==r['volume_names']
        assert h['flow_names']==read(f"reports/m5-biological/{c['pump']}/result.json")['flow_names']
        profile=next(x for x in read('configs/m5-organ-preservation.json')['cases'] if x['name']==cfg['demand_case'])
        assert r['case']['profile']==profile['profile'] and r['case']['liver_function_fraction']==profile['liver_function_fraction']
        hr=h['runs'][-1];mr=r['records'][-1]
        work=sum(x['work']['chamber_J'] for x in hr['sides'].values());d=sum(x['total_mL_min'] for x in hr['demand'].values())
        assert abs(d-r['case']['auxiliary_mL_min'])<1e-10 and hr['demand']==r['case']['sides']
        if c['preservation_conductance']==.25:
            prior=sum(x['work']['chamber_J'] for x in old[c['pump']]['runs'][-1]['sides'].values())
            assert abs(work-prior)<.005 and abs(d-olddemands[c['pump']]['auxiliary_mL_min'])<.005
        a=np.load(OUT/c['name']/'fine-refined.npz');t=a['time_s'];mask=t>=cfg['duration_s']-2
        flow_summary={edge:float(np.trapezoid(a['flows_mL_s'][mask,h['flow_names'].index(edge)],t[mask])/(t[mask][-1]-t[mask][0])*60) for edge in ['brain','hepatic_outflow','respiratory_muscle','lungs','kidney_protected','kidney_other','aux_coronary_left','aux_coronary_right']}
        ao=a['pressures_mmHg'][:,h['volume_names'].index('AO')]
        full_pressure=dict(mean=float(np.trapezoid(ao,t)/(t[-1]-t[0])),minimum=float(ao.min()),maximum=float(ao.max()))
        records.append(dict(case=c,hydraulics=hr,metabolic=mr,work_J=work,auxiliary_mL_min=d,regional_flow_mL_min=flow_summary,full_interval_AO_mmHg=full_pressure))
        lines.append('| %s | %.2f / %.2f / %.2f | %.3f | %.4f | %.4f | %s |'%(c['name'],hr['mean_AO_mmHg'],hr['min_AO_mmHg'],hr['max_AO_mmHg'],hr['auxiliary_left_L_min'],work/cfg['duration_s'],d,'pass' if hr['pressure_screen'] else 'fail'))
    lines+=['','The screen uses mean >=65, minimum >=50 and maximum <=90 mmHg over the fixed window. It does not establish beat-defined systolic/diastolic pressures or tissue preservation. Chamber power is mechanical output, not metabolic input.','','## Demand and delivery','','| Case | Total target mL O2/min | Active unmet mL | Preservation unmet mL | Auxiliary unmet mL | First active shortfall s | First preservation shortfall s | First auxiliary shortfall s |','|---|---:|---:|---:|---:|---|---|---|']
    def bracket(b):return 'none observed' if b is None else '%.3f–%.3f'%tuple(b)
    for r in records:
        g=r['metabolic']['groups'];total=sum(v['target_final_mL_min'] for v in g.values())
        lines.append('| %s | %.3f | %.6f | %.6f | %.6f | %s | %s | %s |'%(r['case']['name'],total,g['active']['unmet_O2_mL'],g['preservation']['unmet_O2_mL'],g['pump']['unmet_O2_mL'],bracket(g['active']['first_shortfall_bracket_s']),bracket(g['preservation']['first_shortfall_bracket_s']),bracket(g['pump']['first_shortfall_bracket_s'])))
    lines+=['','Branches with integrated unmet oxygen above 1e-7 mL (a reporting tolerance, not an injury threshold):','']
    for r in records:
        failed={k:v for k,v in r['metabolic']['per_organ_unmet_oxygen_mL'].items() if v>1e-7}
        lines.append('- '+r['case']['name']+': '+(', '.join('%s %.6f mL'%(k,v) for k,v in failed.items()) or 'none'))
    lines+=['','## Regional flow response' ,'','Final-two-second signed mean flows. Hepatic outflow includes actual arterial and portal return. These are blood flows, not tissue oxygen uptake.','','| Case | Brain mL/min | Hepatic outflow mL/min | Respiratory muscle mL/min | Both kidneys mL/min | Auxiliary coronaries mL/min |','|---|---:|---:|---:|---:|---:|']
    for r in records:
        f=r['regional_flow_mL_min']
        lines.append('| %s | %.2f | %.2f | %.2f | %.2f | %.2f |'%(r['case']['name'],f['brain'],f['hepatic_outflow'],f['respiratory_muscle'],f['kidney_protected']+f['kidney_other'],f['aux_coronary_left']+f['aux_coronary_right']))
    lines+=['','## Findings and next decision','','No case clears every modeled supply requirement. All have an auxiliary-left coronary deficit. The 0.10 conductance cases additionally have native-left/right coronary deficits despite higher aortic pressure. Uniform preservation restriction therefore cannot be selected solely from its pressure or active-flow benefit.','','At fixed rate, stronger restriction increases modeled brain flow and reduces renal flow; hepatic outflow falls despite higher pressure because arterial and portal routes interact. At either conductance, the faster pump consumes more oxygen and has greater integrated auxiliary deficit. These are conditional comparisons within this grid, not universal monotonic laws.','','Next test native coronary controls separately from other preservation branches, and characterize auxiliary coronary delivery versus workload. Retain the same demand floors initially so an apparent improvement cannot be manufactured by lowering a failing target. The inherited missing tissue-glucose-buffer limitation also needs to remain visible when interpreting coronary shortages. No configuration is selected as an optimized or viable design.']
    lines+=['','Full-interval pressure extrema and means are also retained in summary.json. A final-window screen is not a statement that every earlier transient sample satisfies the pressure criteria.','','An eight-second replay with finite assumed oxygen stores is not evidence of steady-state adequacy. Any auxiliary-shortfall flag limits the interpretation of later prescribed pumping. Per-organ shortages remain in each result.json; favorable group totals do not waive pressure, recovery or other constraints.','','## Verification','','Paired hydraulic trajectory/event comparisons, conserved blood inventory, septal-root checks, independent PV and inertial-energy checks, and paired metabolic conservation/refinement checks pass. Unchanged 0.25-conductance reference/faster cases reproduce the historical mechanical work and self-demand within 0.005 J and 0.005 mL/min. Existing core regression evidence remains the 206-test record; this new script adds the paired integration checks above.','','All earlier manifests are checked before this report is emitted. The new manifest records code, configurations, profiles, source seeds, trajectories and outputs. No historical reference is retuned.','','## Reproduce','','```powershell','./.venv/Scripts/python.exe scripts/run_m5_joint_preservation.py','./.venv/Scripts/python.exe scripts/report_m5_joint_preservation.py','```','']
    for manifest in ['reports/m5-normal-comparison/manifest.json','reports/m5-biological-energy/audit.json','reports/m5-native-preservation/manifest.json','reports/m5-organ-preservation/manifest.json']:
        for p,h in read(manifest)['hashes'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    (OUT/'report.md').write_text('\n'.join(lines),encoding='utf-8')
    (OUT/'summary.json').write_text(json.dumps(dict(status='verified_conditional_exploration',records=records),indent=2)+'\n',encoding='utf-8')
    files={Path(__file__),ROOT/'scripts/run_m5_joint_preservation.py',ROOT/'scripts/run_m5_organ_preservation.py',ROOT/'configs/m5-joint-preservation.json',ROOT/'configs/m5-organ-preservation.json',ROOT/'configs/m5-myocardial-demand.json',ROOT/'docs/m5-joint-preservation.md',ROOT/'reports/m5-organ-preservation-tests.xml',ROOT/'reports/m5-normal-comparison/verification.json',ROOT/'data/m1-demand-budget.json'}
    files.update((ROOT/'src/cardiac_redundancy').glob('*.py'));files.update((ROOT/'configs').glob('*'))
    for c in cfg['cases']:
        for name in ['hydraulics.json','result.json']:
            files.update(ROOT/p for p in read(f"reports/m5-joint-preservation/{c['name']}/{name}")['inputs'])
    files.update(p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json')
    (OUT/'manifest.json').write_text(json.dumps(dict(status='verified',paired_hydraulic_cases=len(records),paired_metabolic_cases=len(records),hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files) if p.is_file()}),indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
