"""Compare the archived thoracic diagnostic runs without changing any target."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main():
    directory=ROOT/"reports/m3-thoracic-diagnostic"
    baseline=json.loads((ROOT/"reports/m3-sweeps/baseline/result.json").read_text())
    nominal=json.loads((ROOT/"reports/m3-sweep-repairs-thoracic/thoracic_zero/result.json").read_text())
    rows=[]
    lines=['# Thoracic-pressure reproduction discrepancy','',
           '2026-10-04. Post-sweep numerical diagnosis; the original 9±5 percentage-point screen is unchanged. The baseline and the zero-thoracic-pressure configuration each warm up independently. These are steady-cycle comparisons, not an acute transient timing experiment.','',
           '| Run | CO (L/min) | Decrease from baseline (%) | CO difference from nominal perturbation (%) |','|---|---:|---:|---:|']
    for name,path in [('nominal',ROOT/'reports/m3-sweep-repairs-thoracic/thoracic_zero/result.json'),
                      ('refined',directory/'thoracic_zero_refined/result.json'),
                      ('BDF',directory/'thoracic_zero_BDF/result.json')]:
        result=json.loads(path.read_text())
        if result['status']!='passed_numerical_checks':
            raise RuntimeError(f'Diagnostic incomplete/failed: {name}')
        co=result['metrics']['mean']['CO_L_min']
        decrease=100*(1-co/baseline['metrics']['mean']['CO_L_min'])
        difference=100*abs(co-nominal['metrics']['mean']['CO_L_min'])/nominal['metrics']['mean']['CO_L_min']
        rows.append({'run':name,'CO_L_min':co,'decrease_pct':decrease,'difference_from_nominal_pct':difference,
                     'source_screen_passed':4<=decrease<=14})
        lines.append(f'| {name} | {co:.6f} | {decrease:.6f} | {difference:.6f} |')
    max_difference=max(r['difference_from_nominal_pct'] for r in rows)
    lines += ['',f'Largest tested change in perturbed CO is {max_difference:.6f}%. '
              + ('This is below the 0.5% numerical screen; the failed source comparison persists under refinement and a second solver.' if max_difference<.5 and not any(r['source_screen_passed'] for r in rows) else 'Investigate the differing diagnostic outcomes; do not claim numerical insensitivity.'),'',
              'The original [Smith manuscript](https://ir.canterbury.ac.nz/server/api/core/bitstreams/bdd17d64-212e-46af-ac38-6eee0f28f06b/content) abstract and Results/Figure 8 report approximately 9%. We rechecked the runtime mapping against pinned PMR equations and CVSX `models.py`: external thoracic pressure enters PA, PU and pericardium. The dedicated unit test confirms a +4 mmHg shift is added once to LV/RV/PA/PU pressures at unchanged volume, leaving AO/VC unchanged.','',
              'The discrepancy is therefore not explained by the tested integration/sampling settings or a missing/double-applied external offset. Its cause remains unresolved. Candidate questions are exact parameter-set identity, waveform/driver convention, perturbation timing and whether the reconstructed source matches the author run. No alternative physiological configuration was tuned to obtain 9%.','',
              'Next: execute the matched independent CVSX comparison, inspect the original parameter/protocol provenance and compare the perturbation under matched conditions. M3 source-reproduction closure remains open. The main baseline screen, integration verification and conditional parameter experiments remain usable with this limitation.','',
              'Reproduce with `scripts/run_m3_sweeps.py --protocol configs/m3-thoracic-diagnostic.json --workers 2 --out reports/m3-thoracic-diagnostic`, followed by `scripts/report_m3_thoracic.py`. Each run retains an implementation snapshot, effective config, warm-up history and measured trajectory.','']
    (directory/'diagnosis.json').write_text(json.dumps(rows,indent=2)+'\n')
    (directory/'diagnosis.md').write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps(rows,indent=2))


if __name__=='__main__': main()
