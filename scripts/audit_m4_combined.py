"""Audit saved M4 evidence and replay each baseline cycle with current code."""
import hashlib,json,re,sys,xml.etree.ElementTree as ET
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.whole_network import WholeNetworkModel
from cardiac_redundancy.organ_oxygen import map_demands
from cardiac_redundancy.oxygen import OxygenParameters
from cardiac_redundancy.respiratory import RespiratoryOxygenTransport
from cardiac_redundancy.oxygen_numerics import bounded_jacobian


def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))
def digest(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()


def replay(store):
    h=WholeNetworkModel(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml',arterial_fraction=.2,inlet_resistance_fraction=.1)
    organs=map_demands(read('data/m1-demand-budget.json')['baseline_mL_min'],read('configs/m4-organ-mapping.json')['mapping'],h.flow_names,.75)
    m=RespiratoryOxygenTransport(h.volume_names,h.edges,organs,OxygenParameters(**read('configs/m4-oxygen-reference.json')['parameters']),store_seconds=store)
    z=np.load(ROOT/'reports/m4-baseline/native-cycle.npz');saved=np.load(ROOT/f'reports/m4-baseline/store-{store}-refined.npz')
    t=z['time_s'];v=interp1d(t,z['states'][:len(h.volume_names)],axis=1);q=interp1d(t,z['flows'],axis=0)
    rhs=lambda ti,y:m.rhs(v(min(ti,.75)),y,q(min(ti,.75)))
    sol=solve_ivp(rhs,(0,.75),saved['oxygen_states'][0],method='Radau',jac=bounded_jacobian(rhs),rtol=1e-10,atol=1e-12,max_step=.000625,t_eval=t)
    error=float(np.max(abs(sol.y.T-saved['oxygen_states'])))
    assert sol.success and error<1e-7,(store,error)
    return dict(store_seconds=store,current_code_replay_difference_mL=error)


def main():
    baseline=read('reports/m4-baseline/result.json');manifest=read('reports/m4-respiratory/manifest.json')
    for p,expected in baseline['hashes'].items():
        actual=baseline.get('archived_runtime',{}).get(p,p)
        assert digest(actual)==expected,('baseline input/runtime hash',p)
    for group in ('hashes','output_hashes'):
        for p,expected in manifest[group].items():assert digest(p)==expected,('respiratory hash',p)
    protocol=read('configs/m4-respiratory-protocol.json')
    assert manifest['status']=='verified' and len(protocol['cases'])==17
    for c in protocol['cases']:
        r=read(f"reports/m4-respiratory/{c['name']}/result.json")
        assert r['case']==c and r['status']=='verified' and len(r['records'])==2
        assert r['refinement_difference_mL']<.002 and r['refinement_ventilation_difference_mL_s_BTPS']<.002
        for a in r['records']:
            assert a['oxygen_conservation_error_mL']<1e-6 and a['demand_accounting_error_mL']<1e-6
            assert a['assessments']['all_tissues']['biological_recovery']=='unassessed'
    for p,expected in read('data/m2-freeze-manifest.json')['sha256'].items():assert digest(p)==expected,('M2 changed',p)
    tree=ET.parse(ROOT/'reports/m4-combined-tests.xml');suites=list(tree.getroot().iter('testsuite'))
    tests=sum(int(s.attrib['tests']) for s in suites)
    assert tests==148 and all(int(s.attrib.get('failures',0))+int(s.attrib.get('errors',0))==0 for s in suites)
    docs=['docs/m0-m4-report.md','docs/m4-respiratory-verification.md','reports/m4-baseline/report.md','reports/m4-respiratory/report.md']
    for name in docs:
        p=ROOT/name
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'):continue
            assert (p.parent/target.split('#')[0]).exists(),('broken link',name,target)
    with ProcessPoolExecutor(max_workers=3) as pool:replays=list(pool.map(replay,(0,5,15)))
    result=dict(status='verified',test_count=tests,respiratory_scenarios=17,paired_runs=34,
        baseline_replays=replays,frozen_M2_files=9,local_report_links='resolved',
        provenance='Original baseline runtime snapshots match their recorded hashes; current-code cycle replay independently passes.',
        milestone='M4 open: finite glycogen/glucose accounting and final gate review',
        biological_recovery='unassessed')
    (ROOT/'reports/m4-combined-audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
