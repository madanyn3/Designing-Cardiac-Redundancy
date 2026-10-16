"""Verify bounded M7 closure, preserve prior evidence, freeze new artifacts."""
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]


def verify(path,field):
    hashes=json.loads((ROOT/path).read_text())[field]
    for name,digest in hashes.items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise RuntimeError(f'Changed evidence: {name}')
    return len(hashes)


def main():
    prior={}
    for path,field in [('reports/m5-final-comparison/manifest.json','hashes'),
                       ('reports/m7-global/manifest.json','sha256'),
                       ('reports/m7-global/input-manifest.json','sha256'),
                       ('reports/m7-identifiability/recovery/manifest.json','sha256'),
                       ('reports/m7-structural-replay/manifest.json','sha256'),
                       ('reports/m7-regulation/input-manifest.json','hashes'),
                       ('reports/m7-finite-reserve/summary.json','hashes')]:
        prior[path]=verify(path,field)
    reg=json.loads((ROOT/'reports/m7-regulation/summary.json').read_text())
    if len(reg)!=3:raise RuntimeError('Three regulation cases required')
    for row in reg:
        for name in ['result.json','hydraulics.json']:
            if json.loads((ROOT/'reports/m7-regulation'/row['name']/name).read_text())['status']!='verified':raise RuntimeError('Unverified case')
    reserve=json.loads((ROOT/'reports/m7-finite-reserve/summary.json').read_text())['rows']
    if len(reserve)!=9 or any(r['paired_difference_mL']>1e-6 or r['balance_error_mL']>1e-8 or r['residual_mL'] < -1e-7 for r in reserve):
        raise RuntimeError('Incomplete or unverified reserve diagnostic')
    suites=ET.parse(ROOT/'reports/m7-final-tests.xml').getroot().findall('testsuite')
    if not suites or any(int(s.get('errors',0))+int(s.get('failures',0)) for s in suites):raise RuntimeError('Failed regression')
    paths=[]
    for directory in ['reports/m7-regulation','reports/m7-finite-reserve']:
        paths += [p for p in (ROOT/directory).rglob('*') if p.is_file()]
    paths += list((ROOT/'src/cardiac_redundancy').glob('*.py'))
    paths += [ROOT/'scripts'/name for name in ['run_m7_regulation.py','report_m7_regulation.py','report_m7_reserve.py','finalize_m7.py']]
    paths += [ROOT/'tests/test_regulated_circulation.py',ROOT/'tests/test_finite_reserve.py',ROOT/'reports/m7-final-tests.xml']
    paths += [ROOT/'docs'/name for name in ['m7-closure.md','m7-t74-t75-closure.md','m7-evidence-priorities.md','m7-structural-physiology-evidence.md','m7-final-structural-protocol.md']]
    result=dict(status='M7 closed at bounded assessment scope',tests=sum(int(s.get('tests',0)) for s in suites),
                prior_verified=prior,sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(paths))})
    out=ROOT/'reports/m7-final';out.mkdir(parents=True,exist_ok=True)
    (out/'manifest.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='sha256'},indent=2))


if __name__=='__main__':main()
