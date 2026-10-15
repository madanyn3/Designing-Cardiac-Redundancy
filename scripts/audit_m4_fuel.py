"""Audit current fuel evidence and preserved historical provenance without reruns."""
import hashlib,json,re,sys,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]


def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))
def digest(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()


def main():
    fuel=read('reports/m4-fuel/manifest.json')
    assert fuel['status']=='verified'
    for manifest in (fuel,read('reports/m4-respiratory/manifest.json')):
        for group in ('hashes','output_hashes'):
            for path,expected in manifest[group].items():assert digest(path)==expected,('hash',path)
    baseline=read('reports/m4-baseline/result.json')
    for path,expected in baseline['hashes'].items():
        assert digest(baseline.get('archived_runtime',{}).get(path,path))==expected,('baseline',path)
    for path,expected in read('data/m2-freeze-manifest.json')['sha256'].items():
        assert digest(path)==expected,('frozen M2',path)
    maximum_release_error=maximum_nonox_error=0.
    for r in fuel['results']:
        n=len(r['volume_names']);m=len(r['organ_names']);p=r['parameters']
        name=r['case']['name'] if 'case' in r else r['name']
        for label in ('nominal','refined'):
            z=np.load(ROOT/f'reports/m4-fuel/{name}/{label}.npz');y=z['states'];t=z['time_s']
            if 'case' in r:y=y[:,n+5+3*m:]
            initial=p['liver_volume_L']*p['glycogen_mmol_L_liver']*p['accessible_fraction']
            er=float(np.max(abs(y[:,n]+y[:,-3]-initial)))
            en=float(np.max(abs(y[:,-2]+y[:,-1]-t*p['nonoxidative_mmol_min']/60)))
            assert er<1e-7 and en<1e-7,(name,er,en)
            maximum_release_error=max(maximum_release_error,er);maximum_nonox_error=max(maximum_nonox_error,en)
    suites=list(ET.parse(ROOT/'reports/m4-fuel-tests.xml').getroot().iter('testsuite'))
    tests=sum(int(s.attrib['tests']) for s in suites)
    assert tests>=161 and all(int(s.attrib.get('failures',0))+int(s.attrib.get('errors',0))==0 for s in suites)
    docs=['docs/m4-fuel-accounting.md','reports/m4-fuel/report.md','docs/m4-closure.md']
    for name in docs:
        path=ROOT/name
        for target in re.findall(r'\]\(([^)]+)\)',path.read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'):continue
            if (path.parent/target.split('#')[0]).resolve()==(ROOT/'reports/m4-fuel/audit.json').resolve():continue
            assert (path.parent/target.split('#')[0]).exists(),('broken link',name,target)
    result=dict(status='verified',test_count=tests,fuel_cases=fuel['cases'],paired_runs=fuel['paired_runs'],
        glycogen_release_identity_error_mmol=maximum_release_error,nonoxidative_demand_identity_error_mmol=maximum_nonox_error,
        frozen_M2_files=9,historical_oxygen_provenance='unchanged',local_links='resolved',
        milestone='M4 closed at declared exploratory accounting scope; recovery and physiological duration unassessed')
    (ROOT/'reports/m4-fuel/audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
