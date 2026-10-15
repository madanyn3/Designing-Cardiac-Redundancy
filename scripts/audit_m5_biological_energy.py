"""Audit mechanical T5.2 results and preserve implementation/output provenance."""
from pathlib import Path
import hashlib,json,sys,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
def main():
    out=ROOT/'reports/m5-biological-energy';r=json.loads((out/'verification.json').read_text())
    assert r['status']=='verified' and len(r['results'])==4
    for p,h in r['hashes'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    old=json.loads((ROOT/'reports/m5-t51/manifest.json').read_text())
    for p,h in old['hashes'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    tree=ET.parse(ROOT/'reports/m5-t52-first-tests.xml');suites=list(tree.getroot().iter('testsuite'))
    assert sum(int(s.get('tests',0)) for s in suites)==183
    assert all(int(s.get('errors',0))==int(s.get('failures',0))==0 for s in suites)
    maxima={k:max(run['sides'][side][k] for case in r['results'].values() for run in case['runs'] for side in ('left','right')) for k in ('storage_integral_error_J','coarse_grid_difference_J','PV_quadrature_difference_J','maximum_balance_residual_W')}
    files={Path(__file__),ROOT/'scripts/report_m5_biological_energy.py',ROOT/'reports/m5-t52-first-tests.xml',*list((ROOT/'src/cardiac_redundancy').glob('*.py')),*list((ROOT/'tests').glob('test_biological_energy.py'))}
    files.update(ROOT/p for p in r['hashes']);files.update(ROOT/p for p in ('configs/baseline.yaml','configs/m3-evidence-regional.yaml','configs/m2-frozen-overlay.yaml','configs/m3-whole-network-protocol.json'))
    files.update(p for p in out.iterdir() if p.is_file() and p.name!='audit.json')
    audit=dict(status='verified',tests=183,paired_cases=4,scope='mechanical first increment; T5.2 remains open',maximum_errors=maxima,historical_T51_hashes='unchanged',hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)})
    (out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps({k:v for k,v in audit.items() if k!='hashes'},indent=2))
if __name__=='__main__':main()
