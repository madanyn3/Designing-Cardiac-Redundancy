"""Summarize coupled one-factor sensitivities without treating them as convergence tests."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reports/m3-pressure-sensitivity'


def main():
    protocol=json.loads((OUT/'manifest.json').read_text())['protocol']
    base=json.loads((ROOT/'reports/m3-pressure-regional/baseline/result.json').read_text())
    results={'nominal':base}
    for case in protocol['cases']:
        results[case['name']]=json.loads((OUT/case['name']/'result.json').read_text())
    complete=all(r['status']=='passed_numerical_checks' for r in results.values())
    rows=['# Coupled cerebral/coronary sensitivity results','',
        '2026-10-06. One-factor experiments about the pressure-aware nominal configuration; not a full factorial sweep or physiological validation. Each case has independent periodic initialization and ten measured cycles. Nominal solver convergence is reported separately.','',
        'The 3x3 M2 ICP/tone grid yields seven distinct sums: 5, 10, 15, 20, 25, 30, 35 mmHg. Because this model depends only on the sum, duplicates are exactly the same constitutive parameterization. Pc=20 is the nominal reference. R multipliers 0.5/1.5 and kappa=0/1 are varied separately about the nominal setting. Cerebral R is calibrated only at nominal Pc=20, never at each tested threshold.','',
        '| Case | Warm-up | MAP mmHg | CO L/min | Brain mL/min | Left coronary mL/min | Right coronary mL/min | Status |',
        '|---|---:|---:|---:|---:|---:|---:|---|']
    summaries={}
    for n,r in results.items():
        if r['status']!='passed_numerical_checks':
            rows.append(f"| {n} | — | — | — | — | — | — | {r['status']}: {r.get('error','not finished')} |")
            continue
        m=r['metrics']; q=r['mean_branch_flows_mL_s']
        summaries[n]={'MAP_mmHg':m['MAP_mmHg'],'CO_L_min':m['CO_L_min'],
                      **{k+'_mL_min':q[k]*60 for k in ('brain','coronary_left','coronary_right')}}
        rows.append(f"| {n} | {r['warmup_cycles']} | {m['MAP_mmHg']:.4f} | {m['CO_L_min']:.5f} | {q['brain']*60:.3f} | {q['coronary_left']*60:.3f} | {q['coronary_right']*60:.3f} | passed numerical checks |")
    rows += ['',f'All declared cases completed numerical checks: **{complete}**.','']
    if complete:
        brain_nominal=summaries['nominal']['brain_mL_min']
        brain_low=summaries['Pc35']['brain_mL_min']
        brain_high=summaries['Pc5']['brain_mL_min']
        r_half=summaries['R0.5']['brain_mL_min']/brain_nominal
        r_high=summaries['R1.5']['brain_mL_min']/brain_nominal
        cor0=summaries['kappa0']['coronary_left_mL_min']
        cor1=summaries['kappa1']['coronary_left_mL_min']
        rows += [f'At fixed cerebral R, raising closing pressure from 5 to 35 mmHg changes mean cerebral flow from {brain_high:.3f} to {brain_low:.3f} mL/min. The half-resistance and 1.5-times-resistance cases produce {r_half:.4f} and {r_high:.4f} times nominal cerebral flow, respectively. Changing kappa from 0 to 1 changes left-coronary flow from {cor0:.3f} to {cor1:.3f} mL/min. These comparisons include the response of the coupled circulation.','']
    rows += [
        'Changes between these rows reflect a different assumed model, not numerical error. Flow ratios are solved outputs; a half-resistance case is not forced to double achieved flow, because aortic/venous pressures and cardiac filling respond. Likewise higher closing pressure can redistribute blood and alter the heart operating point. No row establishes an oxygen minimum or recovery duration.','',
        'The four nominal comparison runs test method, tolerance and initialization. The sensitivity runs use nominal numerical settings; their pass label checks conservation, domain, periodicity and algebraic residuals, not a separate refined/BDF comparison at every sensitivity point. Higher-order interactions require a later joint sweep.','',
        'Run `scripts/run_m3_regional.py --protocol configs/m3-pressure-sensitivity-protocol.json --out reports/m3-pressure-sensitivity --workers 4` and `scripts/report_m3_pressure_sensitivity.py` with the project Python. Each case retains its effective closure parameters, configuration/code snapshot, result, trajectory and valve history.','']
    (OUT/'report.md').write_text('\n'.join(rows),encoding='utf-8')
    (OUT/'summary.json').write_text(json.dumps({'all_completed':complete,'cases':summaries},indent=2)+'\n')
    print(json.dumps({'all_completed':complete,'cases':summaries},indent=2))
    return 0 if complete else 1


if __name__=='__main__':raise SystemExit(main())
