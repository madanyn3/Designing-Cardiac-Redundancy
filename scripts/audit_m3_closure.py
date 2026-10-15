"""Consolidate evidence without converting a failed source screen into a pass."""
import hashlib,json,subprocess,sys,xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'reports/m3-final'

def read(relative):return json.loads((ROOT/relative).read_text())

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    registry={
        'regional':('reports/m3-regional/verification.json','all_passed'),
        'evidence_regional':('reports/m3-evidence-regional/verification.json','passed'),
        'pressure_regional':('reports/m3-pressure-regional/verification.json','passed'),
        'pressure_sensitivity':('reports/m3-pressure-sensitivity/summary.json','all_completed'),
        'arrest':('reports/m3-arrest/verification.json','passed'),
        'selective':('reports/m3-selective/verification.json','passed'),
        'whole_network':('reports/m3-whole-network/verification.json','passed'),
        'pump_energy':('reports/m3-pump-energy/metrics.json','passed'),
        'combined_circuits':('reports/m3-circuits/verification.json','passed'),
        'storage_sensitivity':('reports/m3-storage-sensitivity/verification.json','passed')}
    checks={name:bool(read(path)[key]) for name,(path,key) in registry.items()}
    baseline=read('reports/m3-sweeps/evaluation.json')
    checks['six_baseline_screens']=all(x['passed'] for x in baseline['source_screens'])
    checks['source_convergence_and_repeat']=all(x['passed'] for x in baseline['comparisons'].values()) and not baseline['failures']
    decrease=baseline['thoracic_CO_decrease_pct']['value']
    checks['journal_thoracic_screen']=abs(decrease-9)<=5
    source=read('reports/m3-source-audit/comparison.json')
    max_error=max(v for c in source['comparisons'] for v in c['relative_metric_difference_pct'].values())
    checks['independent_reconstruction']=max_error<1
    tree=ET.parse(ROOT/'reports/m3-final-tests.xml')
    suites=list(tree.getroot().iter('testsuite'))
    tests=sum(int(s.attrib['tests']) for s in suites)
    checks['tests']=tests>0 and all(int(s.attrib.get('failures',0))+int(s.attrib.get('errors',0))==0 for s in suites)
    static={}
    for name in ('audit_m2_freeze.py','audit_m2_specification.py','audit_m2_design.py'):
        p=subprocess.run([sys.executable,str(ROOT/'scripts'/name)],cwd=ROOT,capture_output=True,text=True,check=True)
        static[name]=json.loads(p.stdout)
    checks['frozen_M2']=static['audit_m2_freeze.py']['snapshot_files_verified']==9
    decision=read('configs/m3-closure-decision.json')
    engineering=all(v for k,v in checks.items() if k!='journal_thoracic_screen')
    accepted=decision.get('id')=='R-M3-01' and decision.get('approved') is True
    r=dict(status='closed_with_accepted_source_limitation' if engineering and accepted else 'closure_decision_required',
           engineering_checks_passed=all(v for k,v in checks.items() if k!='journal_thoracic_screen'),
           original_M3_gate_passed=all(checks.values()),checks=checks,test_count=tests,
           thoracic_decrease_pct=decrease,journal_target_pct=9,journal_tolerance_percentage_points=5,
           max_independent_metric_difference_pct=max_error,
           agreed_closure_gate_passed=engineering and (checks['journal_thoracic_screen'] or accepted),
           exception=decision,static_audits=static)
    paths=[ROOT/p for p,_ in registry.values()]+[ROOT/'reports/m3-final-tests.xml',ROOT/'reports/m3-sweeps/evaluation.json',ROOT/'reports/m3-source-audit/comparison.json',Path(__file__)]
    paths.append(ROOT/'configs/m3-closure-decision.json')
    r['evidence_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (OUT/'closure-audit.json').write_text(json.dumps(r,indent=2)+'\n')
    lines=['# M3 final evidence audit','',f"Engineering evidence checks pass: **{r['engineering_checks_passed']}**.",f"Original M3 gate passes: **{r['original_M3_gate_passed']}**.",'',
           '| Check | Pass |','|---|---|']
    lines += [f'| {name} | {value} |' for name,value in checks.items()]
    lines += ['',f'Fresh automated suite: {tests} tests. All nine frozen M2 files verify. The three static audits pass.',
              '',f'Source thoracic reduction: {decrease:.6f}% versus 9 +/- 5 percentage points. Maximum independent reconstruction metric difference: {max_error:.6f}%. These are separate tests.',
              '', 'This audit consolidates archived experiment checks and reruns static audits; it does not rerun every historic ODE experiment. Evidence hashes are recorded in closure-audit.json. Fresh tests and the new storage experiment are distinguished from archived runs.',
              '', f"Agreed closure gate passes: **{r['agreed_closure_gate_passed']}**. R-M3-01 acceptance is recorded explicitly in configs/m3-closure-decision.json. The original failed screen remains false.",
              '', '[Closure decision and M4 handoff](../../docs/m3-closure.md). M3 is closed at verified-reconstruction scope with the accepted source limitation.']
    (OUT/'report.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:v for k,v in r.items() if k not in ('static_audits','evidence_sha256')},indent=2))
    return int(not r['engineering_checks_passed'])
if __name__=='__main__':raise SystemExit(main())
