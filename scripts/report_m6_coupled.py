"""Consolidate verified M6 cases without changing any frozen M5 artifact."""
from pathlib import Path
import json,hashlib,sys,xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'reports/m6-coupled'
def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))

def main():
    cfg=read('configs/m6-coupled-protocol.json');results=[]
    lines=['# M6 coupled controller, handover and fault experiments','','2026-10-15. Thirteen scenarios, nominal/refined hydraulic runs, and paired metabolic replay where the hydraulic domain permits. These experiments test policy effects, not a biologically implemented controller or safe recovery duration.','','The observed policy uses causal 0.75 s pressure/native-output filtering, 50 ms sampling, 0.25 s confirmation and 0.25 s signal/recruitment delay. One case halves sampling to 25 ms. The sensor prehistory is an explicit constant initial-pressure/5.213513 L/min assumption for every policy, not newly simulated pre-arrest physiology. The fixed schedule knows arrest onset; that privileged information is not evidence of better sensing.','','| Case | Actual release s | Final 2 s mean AO mmHg | Auxiliary left L/min | Unmet metabolism mL O2-equivalent | Outcome |','|---|---:|---:|---:|---:|---|']
    fig,axes=plt.subplots(3,1,figsize=(10,10),sharex=True)
    selected={'fixed_arrest','observed_arrest','graded_handover','abrupt_handover','rearrest'}
    for case in cfg['cases']:
        name=case['name'];h=read(f'reports/m6-coupled/{name}/hydraulics.json');r=read(f'reports/m6-coupled/{name}/result.json')
        assert h['status']=='verified';last=h['runs'][-1]
        item=dict(case=name,release_s=last['release_s'],MAP_mmHg=last['final_MAP_mmHg'],aux_left_L_min=last['aux_left_L_min'],hydraulic_scaled_errors=h['scaled_errors'],termination=last['termination'],controller_events=last['controller_events'])
        if last['termination'] is None:
            assert r['status']=='verified'
            a=np.load(OUT/name/'hydraulic-refined.npz');b=np.load(OUT/name/'metabolic-refined.npz')
            t=b['time_s'];requested=b['requested_O2_mL_s'];ceiling=b['glucose_ceiling_mL_s'];actual=b['oxygen_consumed_mL_s']
            glucose=requested-ceiling;oxygen=ceiling-actual
            assert glucose.min()>-1e-10 and oxygen.min()>-1e-10
            assert np.max(abs(glucose+oxygen-b['oxygen_unmet_mL_s']))<1e-10
            final=r['records'][-1]
            item.update(unmet_O2_equivalent_mL=final['unmet_oxygen_mL'],groups=final['groups'],
                sampled_glucose_limited_mL=float(np.trapezoid(glucose.sum(axis=1),t)),
                sampled_additional_oxygen_limited_mL=float(np.trapezoid(oxygen.sum(axis=1),t)),
                glycogen_remaining_mmol=final['glycogen_remaining_mmol'],released_mmol=final['released_mmol'],
                paired_oxygen_difference_mL=r['paired_max_oxygen_state_difference'],
                paired_glucose_difference_mmol=r['paired_max_glucose_state_difference_mmol'])
            for rec in r['records']:
                assert rec['oxygen_balance_error_mL']<1e-6 and rec['glucose_balance_error_mmol']<1e-8
            outcome='completed horizon; adequacy separately screened'
            if name in selected:
                axes[0].plot(a['time_s'],a['pressures_mmHg'][:,2],label=name)
                n=len(r['volume_names']);m=len(r['organ_names'])
                axes[1].plot(t,b['states'][:,n+2+m:n+2+2*m].sum(axis=1),label=name)
            if name in ('observed_arrest','sensor_loss','stuck_high','inlet_obstruction','low_preload'):
                axes[2].plot(a['time_s'],a['pressures_mmHg'][:,2],label=name)
        else:
            assert r['status']=='hydraulic_domain_limit';outcome='domain stop; no extrapolated metabolism'
        results.append(item)
        (OUT/name/'status.json').write_text(json.dumps(dict(status='verified_domain_stop' if last['termination'] else 'hydraulic_and_metabolic_verified'),indent=2)+'\n')
        number=lambda v:'—' if v is None else f'{v:.4f}'
        lines.append('| '+name+' | '+number(item['release_s'])+' | '+number(item['MAP_mmHg'])+' | '+number(item['aux_left_L_min'])+' | '+number(item.get('unmet_O2_equivalent_mL'))+' | '+outcome+' |')
    byname={r['case']:r for r in results}
    sampling=dict(default_release_s=byname['observed_arrest']['release_s'],fine_release_s=byname['observed_fast_sample']['release_s'],default_MAP_mmHg=byname['observed_arrest']['MAP_mmHg'],fine_MAP_mmHg=byname['observed_fast_sample']['MAP_mmHg'])
    lines+=['','## Reading the results','','SUPPORT denotes actuator command, not physiological adequacy. The coupled observed controller conservatively receives recovery adequacy=false because online tissue/recovery certification is absent; it does not automatically withdraw on pressure alone. Abrupt/graded handovers are explicitly prescribed scenarios. Synthetic policy tests separately exercise qualified recovery and aborted withdrawal. This is a retained limitation, not a claim of autonomous metabolic control.','','Native restart starts at 4 s with a 1 s ramp. Systemic reopening starts at 5 s over 2 s. Native-organ demand restoration starts independently at 5.5 s over 2 s. Preserved-lung ventilation/perfusion allocation remains unchanged; all restored demand must compete for the actual available supply. Scheduled withdrawal starts at 7 s over 0 or 2 s. Re-arrest occurs at 8 s in its own observed-policy case. No repair or biochemical injury process is simulated.','','Metabolic demand uses the measured auxiliary command envelope; maintenance persists when inactive and workload cost is recalculated from each hydraulic trajectory. The offline work normalization is not causal energetic feedback. The glucose-versus-oxygen components in summary.json are saved-sample quadrature; the combined deficit above is the ODE ledger and can differ slightly at narrow events. No numerical or domain failure is labeled a tissue-survival result.','','Sampling sensitivity: '+json.dumps(sampling)+'. This is distinct from paired ODE-tolerance refinement, which retains the same sensor sampling.','','## Verification and limits','','All hydraulic cases check conserved blood, continuous physical states across integration segments, finite septal residuals, nonnegative valve flows and chamber energy identities. Saved PV quadrature is checked separately. Paired state/pressure/flow scaled differences must remain below 0.005; activation time and controller transition traces must match. Metabolic replay checks paired convergence, oxygen/glucose balances, finite stores and oxidation stoichiometry. Recovery ramps do not reset the auxiliary beat phase. Abrupt withdrawal intentionally changes the pressure law without teleporting blood.','','No test establishes a viable auxiliary design. Source/model limitations, assumed reduced demands, missing local glucose buffering, lack of supply-dependent contractility, absent CO2 and absent injury/blood-compatibility models remain. Low preload and obstruction are synthetic stresses; no real fault probability is assigned. The fixed and observed policies share initial physical inventory, but their capabilities differ in knowledge of arrest onset.','','[Failure-mode coverage](../../docs/failure-modes.md) · [Controller specification](../../docs/m6-controller-specification.md) · [M6 closure](../../docs/m6-closure.md)','','![Pressure and metabolic deficits](comparison.png)']
    axes[0].set_ylabel('Aortic pressure (mmHg)');axes[1].set_ylabel('Unmet metabolism (mL O2-equivalent)')
    axes[2].set_ylabel('Fault-case aortic pressure (mmHg)');axes[2].set_xlabel('Time after native arrest (s)')
    for ax in axes:ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(OUT/'comparison.png',dpi=150);plt.close(fig)
    (OUT/'summary.json').write_text(json.dumps(dict(cases=results,sampling=sampling),indent=2)+'\n')
    (OUT/'verification.json').write_text(json.dumps(dict(status='verified',hydraulic_cases=len(results),metabolic_cases=sum(r['termination'] is None for r in results),cases=[r['case'] for r in results]),indent=2)+'\n')
    suites=list(ET.parse(ROOT/'reports/m6-coupled-tests.xml').getroot().iter('testsuite'))
    assert all(int(s.get('errors',0))+int(s.get('failures',0))==0 for s in suites)
    count=sum(int(s.get('tests',0)) for s in suites)
    lines+=['',f'Full regression suite: {count} passed. Reproduce hydraulics with `scripts/run_m6_coupled.py`, metabolic replay with `scripts/run_m6_metabolism.py`, then this report with `scripts/report_m6_coupled.py`.']
    (OUT/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    old=read('reports/m5-final-comparison/manifest.json')
    for p,digest in old['hashes'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==digest,p
    files=set((ROOT/'src/cardiac_redundancy').glob('*.py'))
    files.update(ROOT/p for p in old['hashes'])
    files.add(ROOT/'reports/m5-final-comparison/manifest.json')
    files.update(ROOT/p for p in ['scripts/run_m6_coupled.py','scripts/run_m6_metabolism.py','scripts/report_m6_coupled.py','configs/m6-controller.json','configs/m6-coupled-protocol.json','docs/m6-closure.md','docs/failure-modes.md','reports/m6-coupled-tests.xml'])
    files.update(OUT.glob('*/*.json'));files.update(OUT.glob('*/hydraulic-*.npz'));files.update(OUT.glob('*/metabolic-*.npz'));files.update(OUT/p for p in ['report.md','summary.json','comparison.png','verification.json'])
    (OUT/'manifest.json').write_text(json.dumps(dict(tests=count,status='verified_qualified_M6',hashes={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}),indent=2)+'\n')
    print('M6 report verified:',count,'tests;',len(results),'cases')

if __name__=='__main__':main()
