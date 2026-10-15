"""Identical-boundary energy comparison and activation-timing audit."""
from pathlib import Path
import json,sys,hashlib,xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.architecture_energy import external_input
from cardiac_redundancy.activation_timing import ActivationTiming
from cardiac_redundancy.pump_energy import MMHG_ML_TO_J as K
from cardiac_redundancy.preservation_energy import oxygen_energy
OUT=ROOT/'reports/m5-architecture-foundation'
def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))
def main():
    config=read('configs/m5-architecture-reference.json');results=[]
    rows=['# T5.4 first increment: matched-boundary energy benchmark','','2026-10-14. BIO is primary; mechanical support is an abstract externally powered reference. These calculations reuse the historical joint-grid boundary histories. They are not independent mechanical circulation simulations and do not implement the revised one-lung/one-kidney organ set.','','All comparisons use the exact same two-stage port pressures and inlet/outlet flow histories, hence the same timing and pulsatility. No extra oxygenator, flow, organ-delivery benefit or oxygen credit is introduced. The comparator assumes zero overhead and no recovery of negative port work; efficiencies are assumed sensitivities.','','| Case | BIO chamber output W | BIO requested metabolic equivalent W | Positive port-work rate W | Mechanical input equivalent at eta=0.5 W |','|---|---:|---:|---:|---:|']
    for case in read('configs/m5-joint-preservation.json')['cases']:
        name=case['name'];path=f'reports/m5-joint-preservation/{name}'
        h=read(path+'/hydraulics.json');r=read(path+'/result.json');a=np.load(ROOT/path/'fine-refined.npz')
        t=a['time_s'];p=a['pressures_mmHg'];y=a['states'];vn=h['volume_names']
        left=K*(p[:,vn.index('AO')]*y[23]-p[:,vn.index('PU')]*y[22])
        right=K*(p[:,vn.index('PA')]*y[25]-p[:,vn.index('VC')]*y[24])
        powers=np.c_[left,right];duration=float(t[-1]-t[0])
        variants=[external_input(t,powers,e,config['overhead_W']) for e in config['mechanical_efficiency_sensitivities']]
        signed=sum(v['work']['port_J'] for v in h['runs'][-1]['sides'].values())
        assert abs(signed-variants[0]['signed_port_J'])<1e-8
        # Coarse-output quadrature is separate from already verified ODE pairs.
        idx=np.unique(np.r_[np.arange(0,len(t),2),len(t)-1])
        check=external_input(t[idx],powers[idx],.5)
        selected=next(v for v in variants if v['efficiency']==.5)
        assert abs(check['external_input_J']-selected['external_input_J'])<.005
        bio=float(oxygen_energy(r['case']['auxiliary_mL_min']*duration/60))
        work=sum(v['work']['chamber_J'] for v in h['runs'][-1]['sides'].values())
        rows.append('| %s | %.4f | %.4f | %.4f | %.4f |'%(name,work/duration,bio/duration,selected['positive_port_J']/duration,selected['external_input_J']/duration))
        results.append(dict(case=name,duration_s=duration,BIO_requested_energy_equivalent_J=bio,BIO_chamber_work_J=work,mechanical_variants=variants))
    rows+=['','BIO input is an assumed oxygen-energy equivalent, not measured ATP expenditure. Mechanical input is an assumed external-energy ledger with zero overhead. Different energy sources and uncalibrated efficiencies prevent selecting a winner from this table. Initial elastic/inertial storage and pericardial work also matter in a transient.','','The formula sums max(Pout*Qout-Pin*Qin,0) separately for each stage, integrates over time, divides by eta, then adds overhead times duration. Negative port work is reported separately and never cancels positive cost in another stage. This is an accounting convention at prescribed boundaries, not proof a physical actuator can reproduce them.','','## Remaining T5.4 work','','Implement the amended mandatory-organ roles and separate supported/preserved pulmonary paths; ensure one-lung ventilation, gas exchange, perfusion and storage are explicit. Then run matched closed-circuit BIO and bounded mechanical comparators with documented blood-volume allocation, native failure, organ targets and endpoints. Keep BIO self-demand versus external mechanical energy separate. Materials/device optimization and E1/E2 architecture searches are outside the user-narrowed scope. T5.4 remains open.','','[Protocol and scope](../../docs/m5-t54-comparison.md).']
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'report.md').write_text('\n'.join(rows)+'\n',encoding='utf-8')
    timing={name:dict(parameters=cfg,release_s=ActivationTiming(**cfg).latency_s,ramp_complete_s=ActivationTiming(**cfg).latency_s+cfg['ramp_s'],effective_support_s=None) for name,cfg in read('configs/m5-activation-timing.json')['cases'].items()}
    maintenance=read('configs/m5-myocardial-demand.json')['reference'];native=read('data/m1-demand-budget.json')['baseline_mL_min']['heart']
    dormant=native*maintenance['maintenance_fraction']*maintenance['mass_ratio']
    (OUT/'results.json').write_text(json.dumps(dict(architecture=results,timing=timing,dormant_assumption=dict(O2_mL_min=dormant,energy_equivalent_W=float(oxygen_energy(dormant)/60),classification='T5.2 proxy extrapolation, not empirical dormant-organ cost')),indent=2)+'\n',encoding='utf-8')
    suites=list(ET.parse(ROOT/'reports/m5-t54-t55-foundation-tests.xml').getroot().iter('testsuite'))
    assert all(int(s.get('errors',0))+int(s.get('failures',0))==0 for s in suites)
    tests=sum(int(s.get('tests',0)) for s in suites)
    for name in ['reports/m5-joint-preservation/manifest.json','reports/m5-primer-audit/manifest.json']:
        for p,h in read(name)['hashes'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    files={Path(__file__),ROOT/'src/cardiac_redundancy/architecture_energy.py',ROOT/'src/cardiac_redundancy/activation_timing.py',ROOT/'tests/test_architecture_activation.py',ROOT/'configs/m5-architecture-reference.json',ROOT/'configs/m5-activation-timing.json',ROOT/'reports/m5-t54-t55-foundation-tests.xml',ROOT/'docs/m5-t54-comparison.md',ROOT/'docs/m5-t55-activation.md',OUT/'report.md',OUT/'results.json'}
    files.update((ROOT/'src/cardiac_redundancy').glob('*.py'))
    for c in read('configs/m5-joint-preservation.json')['cases']:
        files.update(ROOT/f"reports/m5-joint-preservation/{c['name']}/{f}" for f in ['fine-refined.npz','hydraulics.json','result.json'])
    files.update(ROOT/p for p in ['configs/m5-joint-preservation.json','configs/m5-myocardial-demand.json','data/m1-demand-budget.json'])
    (OUT/'manifest.json').write_text(json.dumps(dict(status='verified_foundation_only',tests=tests,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}),indent=2)+'\n',encoding='utf-8')
    print('Verified foundation; tests:',tests)
if __name__=='__main__':main()
