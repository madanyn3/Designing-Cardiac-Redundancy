"""Energy-equivalent ledgers and native-preservation comparisons for T5.2."""
from pathlib import Path
import json,sys,hashlib,xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.preservation_energy import oxygen_energy,integrated_native_fraction
from cardiac_redundancy.fuel import O2_ML_STPD_PER_MMOL_GLUCOSE
OUT=ROOT/'reports/m5-native-preservation'
def read(p):return json.loads((ROOT/p).read_text(encoding='utf-8'))
def main():
    protocol=read('configs/m5-native-preservation.json');sources={c['name']:c for c in read('reports/m5-self-demand/scenarios.json')['scenarios']}
    budget=read('data/m1-demand-budget.json')['baseline_mL_min'];normal=budget['heart'];body=sum(budget.values());mechanics=read('reports/m5-biological-energy/verification.json')['results']
    summaries=[]
    rows=['# T5.2 native preservation and energy accounting','','2026-10-14. Six paired runs; native preservation varies, all other native-organ demands and hydraulics remain fixed. [Evidence and equations](../../docs/m5-t52-native-preservation.md).','','## Demand and supply','','Normal-heart demand reference: '+f'{normal:.4f} mL O2/min. The 30% and 50% native targets are transferred M1 assumptions, not validated human minima.','', '| Case | Final native + auxiliary target mL/min | % normal heart target | Native unmet O2 mL | Auxiliary unmet O2 mL | First auxiliary shortfall s |','|---|---:|---:|---:|---:|---|']
    for c in protocol['cases']:
        case=sources[c['source']]|c;r=read(f"reports/m5-native-preservation/{c['name']}/result.json");assert r['status']=='verified' and r['case']==case
        record=r['records'][-1];a=np.load(OUT/c['name']/'refined.npz');t=a['time_s'];z=a['states'];n=len(r['volume_names']);m=len(r['organ_names']);k=n+5+3*m
        consumed=z[-1,n+2:n+2+m];debt=z[-1,n+2+m:n+2+2*m];requested=consumed+debt
        groups={'native':[i for i,name in enumerate(r['organ_names']) if name in ('coronary_left','coronary_right')], 'auxiliary':[i for i,name in enumerate(r['organ_names']) if name.startswith('aux_coronary_')], 'whole_body':list(range(m))}
        ledger={}
        for group,idx in groups.items():
            req=float(requested[idx].sum());used=float(consumed[idx].sum());unmet=float(debt[idx].sum())
            request_glucose=req/O2_ML_STPD_PER_MMOL_GLUCOSE
            actual_glucose=float(z[-1,k+n+1+np.array(idx)].sum())
            assert abs(actual_glucose-used/O2_ML_STPD_PER_MMOL_GLUCOSE)<1e-8
            ledger[group]=dict(requested_O2_mL=req,supported_O2_mL=used,unmet_O2_mL=unmet,
                requested_energy_J=float(oxygen_energy(req)),supported_energy_J=float(oxygen_energy(used)),unmet_energy_J=float(oxygen_energy(unmet)),
                requested_mean_equivalent_W=float(oxygen_energy(req)/8),requested_glucose_equivalent_mmol=request_glucose,actual_glucose_used_mmol=actual_glucose,
                caloric_sensitivity_requested_J={str(e):float(oxygen_energy(req,e)) for e in protocol['energy_equivalent_kcal_per_L']})
        expected=normal/60*integrated_native_fraction(8,c['native_fraction'],c['native_ramp_s']);assert abs(ledger['native']['requested_O2_mL']-expected)<1e-6
        final_native=normal*c['native_fraction'];final_pair=final_native+case['auxiliary_mL_min'];final_body=body-normal+final_pair
        bracket=record['first_sampled_auxiliary_shortfall_bracket_s'];prefix=t if bracket is None else t[t<=bracket[0]+1e-12]
        # This prefix audits the auxiliary driver only. It certifies neither
        # adequate pressure nor all-organ supply, and is not a safe duration.
        count=len(prefix);np.savez_compressed(OUT/c['name']/'driver-audit-prefix.npz',time_s=t[:count],states=z[:count],volumes_mL=a['volumes_mL'][:count])
        ws=mechanics[case['trajectory']]['runs'][-1]['sides'];work=sum(v['work']['chamber_J'] for v in ws.values());port=sum(v['work']['port_J'] for v in ws.values())
        s=dict(case=c['name'],final_native_target_mL_min=final_native,auxiliary_target_mL_min=case['auxiliary_mL_min'],final_combined_target_mL_min=final_pair,fraction_normal_heart_target=final_pair/normal,final_body_target_mL_min=final_body,ledger=ledger,
            auxiliary_chamber_work_J=work,auxiliary_port_work_J=port,apparent_chamber_work_over_requested_energy=work/ledger['auxiliary']['requested_energy_J'],
            first_auxiliary_shortfall_bracket_s=bracket,driver_audit_prefix_end_s=None if not count else float(prefix[-1]),recovery='unassessed; prefix is not a safe support interval')
        summaries.append(s);b='none observed' if bracket is None else f'{bracket[0]:.3f}–{bracket[1]:.3f}'
        rows.append(f"| {c['name']} | {final_pair:.3f} | {100*final_pair/normal:.1f} | {ledger['native']['unmet_O2_mL']:.5f} | {ledger['auxiliary']['unmet_O2_mL']:.5f} | {b} |")
    rows+=['','The five-second ramp has the same final target as the instantaneous 30% case but greater integrated demand. All six cases still have auxiliary shortages under the declared glucose/oxygen model. Lowering the native demand is not a cure for the auxiliary local supply shortfall.','','## Requested energy equivalents and mechanical work','','Use the inherited reporting convention 4.85 kcal/L O2, equal to 20.2924 J/mL. These are requested metabolic-energy equivalents, not measured ATP costs or new fuel credits. Input-equivalent sensitivities at 4.7 and 5.05 kcal/L are retained in JSON.','','| Case | Native requested J | Auxiliary requested J | Auxiliary supported J | Auxiliary unmet equivalent J | Auxiliary chamber work J | Apparent work/request % |','|---|---:|---:|---:|---:|---:|---:|']
    for s in summaries:
        a=s['ledger']['auxiliary'];n=s['ledger']['native'];rows.append(f"| {s['case']} | {n['requested_energy_J']:.4f} | {a['requested_energy_J']:.4f} | {a['supported_energy_J']:.4f} | {a['unmet_energy_J']:.5f} | {s['auxiliary_chamber_work_J']:.4f} | {100*s['apparent_chamber_work_over_requested_energy']:.2f} |")
    rows+=['','The apparent ratio is not a calibrated muscle efficiency: transient net work, initial stored energy, pericardial exchange and assumed demand all matter. No ratio is improved by replacing requested energy with the lower actual uptake after supply failure. Full signed port work, valve losses and inertial storage remain in the earlier mechanical ledger.','','Requested and actual glucose are reported separately; actual glucose matches actual supported oxygen through the existing 134.484 mL/mmol stoichiometry. No second glycogen charge is added for the reported joules. The mixed-substrate caloric convention is not asserted to calibrate the all-glucose comparator.','','## Normal-heart and total-body comparison','',f"Normal myocardial budget is {normal:.3f} mL/min, or {float(oxygen_energy(normal)/60):.4f} W by the reporting convention. Reference BIO plus native 30% preservation requests {summaries[0]['final_combined_target_mL_min']:.3f} mL/min, {100*summaries[0]['fraction_normal_heart_target']:.1f}% of normal cardiac demand. Whole-body target is {summaries[0]['final_body_target_mL_min']:.3f} versus {body:.3f} mL/min: other organ demands remain unchanged, so cardiac savings do not imply the same whole-body percentage saving.",'','## Numerical verification and disposition','','198 tests pass. Six paired replays check oxygen/glucose conservation, time-varying original-demand integrals, stoichiometry and numerical refinement. They preserve all earlier mechanical and metabolic evidence. Energy-equivalent balances and glucose oxidation are independently checked during report assembly. Each case exports a driver-audit prefix ending before its first sampled auxiliary shortfall; subsequent prescribed-driver results remain diagnostic counterfactuals. This does not establish an actual collapse mechanism, tissue injury or recovery time.','','The task-level disposition and carried limitations are recorded in the [T5.2 closure](../../docs/m5-t52-closure.md).','','## Reproduce','','```powershell','.\\.venv\\Scripts\\python.exe scripts/run_m5_native_preservation.py','.\\.venv\\Scripts\\python.exe scripts/report_m5_native_preservation.py','```']
    (OUT/'report.md').write_text('\n'.join(rows)+'\n',encoding='utf-8')
    (OUT/'energy-summary.json').write_text(json.dumps(dict(status='verified_conditional_accounting',normal_heart_target_mL_min=normal,summaries=summaries),indent=2)+'\n')
    for manifest in ('reports/m5-normal-comparison/manifest.json','reports/m5-biological-energy/audit.json'):
        for p,h in read(manifest)['hashes'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    suites=list(ET.parse(ROOT/'reports/m5-t52-preservation-tests.xml').getroot().iter('testsuite'));assert sum(int(s.get('tests',0)) for s in suites)==198 and all(int(s.get('errors',0))+int(s.get('failures',0))==0 for s in suites)
    files={Path(__file__),ROOT/'scripts/run_m5_native_preservation.py',ROOT/'configs/m5-native-preservation.json',ROOT/'reports/m5-t52-preservation-tests.xml',ROOT/'src/cardiac_redundancy/preservation_energy.py',ROOT/'tests/test_preservation_energy.py',ROOT/'reports/m5-self-demand/scenarios.json',ROOT/'reports/m5-biological-energy/verification.json',ROOT/'data/m1-demand-budget.json',ROOT/'configs/m4-organ-mapping.json',ROOT/'configs/m4-oxygen-reference.json',ROOT/'configs/m4-fuel-reference.json',ROOT/'configs/m5-myocardial-demand.json',*list((ROOT/'src/cardiac_redundancy').glob('*.py'))}
    for c in protocol['cases']:
        files.update(ROOT/p for p in read(f"reports/m5-native-preservation/{c['name']}/result.json")['inputs'])
    files.update(p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json')
    (OUT/'manifest.json').write_text(json.dumps(dict(status='verified',tests=198,paired_cases=6,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}),indent=2)+'\n')
if __name__=='__main__':main()
