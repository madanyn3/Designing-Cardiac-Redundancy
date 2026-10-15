"""Assemble an audited fixed-hydraulics organ-demand exploration report."""
from pathlib import Path
import json,hashlib,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reports/m5-organ-preservation'
def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))
def main():
    config=read('configs/m5-organ-preservation.json')
    sources={c['name']:c for c in read('reports/m5-self-demand/scenarios.json')['scenarios']}
    baseline=sum(read('data/m1-demand-budget.json')['baseline_mL_min'].values())
    rows=['# Tunable active and preservation demand sweep','','2026-10-14. Fixed reference BIO hydraulics, eight-second paired metabolic replays. No tissue recovery or survivability duration is established.','','The active set is whole-brain coverage, lung nutritive tissue, respiratory muscle and liver. The pulmonary gas-exchange circuit remains separately represented. All other native tissues retain positive preservation targets. Auxiliary myocardium retains its independent workload-derived target.','','## Assumptions','','Brain 60%, liver 50%, native heart 30%: inherited M1 exploratory fractions. Lung and respiratory-muscle targets remain 100%. Other native tissues vary at 10%, 25% and 50%; these new fractions are design variables, not published viable minima. A five-second ramp tests transition cost. Liver export-capacity fraction is separately assumed at 0.5 in every new case; oxygen demand does not automatically determine export capacity.','','Saved hydraulics retain preservation supply conductances at 0.25 and active supply at their original settings. Conductance fraction is not flow fraction. The unchanged reference has a final-two-second mean aortic pressure about 59.91 mmHg, below the inherited 65 mmHg mean-pressure screen. No new pressure, flow, pump work or size saving can be inferred from changing demand alone. Absolute oxygen stores and initial inventories remain unchanged.','','## Demand targets','','Targets are fractions of normal regional demand, with auxiliary demand added separately. Normal whole-body reference is %.4f mL O2/min.'%baseline,'','| Case | Active target | Preservation target | Auxiliary target | Total target mL/min | Reduction vs normal |','|---|---:|---:|---:|---:|---:|']
    results=[]
    for c in config['cases']:
        r=read('reports/m5-organ-preservation/'+c['name']+'/result.json')
        assert r['status']=='verified' and r['case']==sources[c['source']]|c
        g=r['records'][-1]['groups'];total=sum(v['target_final_mL_min'] for v in g.values())
        rows.append('| %s | %.3f | %.3f | %.3f | %.3f | %.2f%% |'%(c['name'],g['active']['target_final_mL_min'],g['preservation']['target_final_mL_min'],g['pump']['target_final_mL_min'],total,100*(1-total/baseline)))
        results.append(r)
    rows+=['','The ramp has the same final target as its step counterpart but greater integrated demand. These reductions describe requested metabolism under assumptions, not achieved physiological suppression or measured delivery.','','## Delivery screen','','| Case | Active unmet O2 mL | Preservation unmet O2 mL | Auxiliary unmet O2 mL | First active shortfall s | First auxiliary shortfall s |','|---|---:|---:|---:|---|---|']
    def bracket(b):return 'none observed' if b is None else '%.3f–%.3f'%tuple(b)
    for r in results:
        g=r['records'][-1]['groups']
        rows.append('| %s | %.6f | %.6f | %.6f | %s | %s |'%(r['case']['name'],g['active']['unmet_O2_mL'],g['preservation']['unmet_O2_mL'],g['pump']['unmet_O2_mL'],bracket(g['active']['first_shortfall_bracket_s']),bracket(g['pump']['first_shortfall_bracket_s'])))
    rows+=['','An eight-second replay with credited oxygen stores is not a steady-state or long-duration test. Group sums retain every individual organ; individual results are in each result.json. A group passes its demand screen only when every member does. A shortfall is not an injury threshold, and zero shortfall does not validate recovery. Results after an auxiliary shortage are prescribed-driver counterfactuals under L-M5-E03.','','## Design-space interpretation','','This is the demand axis of the exploration, not a completed optimization. Future joint sweeps must rerun hydraulics for each conductance/pump change and recompute auxiliary self-demand from its changed workload. Keep active-organ objectives, preservation assumptions, liver export capacity and transition times explicit. Report per-organ deficits and constraint failures; do not select the lowest preservation target as a proven survival optimum. Preservation floors require evidence or declared sensitivity bounds outside the optimizer.','','## Verification and reproduction','','Paired solver checks cover oxygen/glucose conservation, the original scheduled-demand integral, oxidation stoichiometry, nonnegative inventories, oxygen-store capacity and refinement. Unit tests check exact profile coverage, independent targets, transition integrals and rejection of zero preservation or auxiliary-demand discounts. Historical manifests remain unchanged.','','```powershell','./.venv/Scripts/python.exe scripts/run_m5_organ_preservation.py','./.venv/Scripts/python.exe scripts/report_m5_organ_preservation.py','```','','[Configuration and usage](../../docs/m5-organ-preservation.md).']
    suites=list(ET.parse(ROOT/'reports/m5-organ-preservation-tests.xml').getroot().iter('testsuite'))
    assert all(int(s.get('failures',0))+int(s.get('errors',0))==0 for s in suites)
    count=sum(int(s.get('tests',0)) for s in suites)
    rows.insert(-1,'%d tests pass in the complete suite.'%count)
    for name in ['reports/m5-normal-comparison/manifest.json','reports/m5-biological-energy/audit.json','reports/m5-native-preservation/manifest.json']:
        for p,h in read(name)['hashes'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    (OUT/'report.md').write_text('\n'.join(rows)+'\n',encoding='utf-8')
    files={ROOT/'configs/m5-organ-preservation.json',ROOT/'scripts/run_m5_organ_preservation.py',Path(__file__),ROOT/'src/cardiac_redundancy/organ_preservation.py',ROOT/'tests/test_organ_preservation.py',ROOT/'reports/m5-organ-preservation-tests.xml',ROOT/'reports/m5-self-demand/scenarios.json',ROOT/'configs/m4-fuel-reference.json',ROOT/'configs/m4-oxygen-reference.json',ROOT/'configs/m4-organ-mapping.json',ROOT/'data/m1-demand-budget.json',ROOT/'docs/m5-organ-preservation.md'}
    files.update((ROOT/'src/cardiac_redundancy').glob('*.py'))
    for r in results:files.update(ROOT/p for p in r['inputs'])
    files.update(p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json')
    (OUT/'manifest.json').write_text(json.dumps(dict(status='verified',tests=count,paired_cases=len(results),hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}),indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
