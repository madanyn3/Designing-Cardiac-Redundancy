"""Report all registered recovery cases and enforce declared numerical gates."""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'reports/m7-identifiability/recovery'


def main():
    protocol = json.loads((OUT/'input-manifest.json').read_text())
    for path, expected in protocol['hashes'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != expected:
            raise RuntimeError(f'Changed input: {path}')
    rows = []
    lines = ['# T7.3 nonlinear recovery and compensation tests', '',
             'Synthetic nominal truth; fixed four-parameter subset; conditional local study, not global uniqueness or patient calibration.', '',
             '| Experiment | Iterations | Recovered multipliers: Rsys, Rpul, Cao, contractility | Max refined scaled residual | Gate |',
             '|---|---:|---|---:|---|']
    for spec in protocol['specs']:
        result = json.loads((OUT/spec['name']/'recovery.json').read_text())
        error = float(np.max(np.abs(result['final_log_multipliers'])))
        passed = result['converged_update']
        if not spec.get('noise'):
            passed = passed and result['max_refined_residual'] < .0005
            if 'free' not in spec:
                passed = passed and error < .002
        rows.append(dict(name=spec['name'], passed=bool(passed), max_log_parameter_error=error,
                         multipliers=result['multipliers'], max_scaled_residual=result['max_refined_residual'],
                         iterations=len(result['history'])))
        lines.append(f"| {spec['name']} | {len(result['history'])} | {', '.join(f'{v:.6f}' for v in result['multipliers'])} | {result['max_refined_residual']:.6g} | {'pass' if passed else 'FAIL'} |")
    lines += ['', 'The noisy case is one illustrative perturbation, not an uncertainty distribution or coverage study. A cases hold Rpul and Cao away from truth and fit the remaining two coordinates: successful MAP/CO matching demonstrates finite nonlinear ambiguity. B/C recoveries only test two nearby starts and a nominal truth.',
              '', 'Every iteration solves the nonlinear periodic model. A fixed local Jacobian only proposes updates. Final residuals come from refined integration; no trajectory surrogate supplies accepted residuals. Numerical tolerances, bounds and observation scales are specified in `docs/m7-recovery-protocol.md`.']
    summary = dict(all_gates_pass=all(r['passed'] for r in rows), results=rows)
    (OUT/'summary.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    (OUT/'report.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    print('\n'.join(lines))
    if not summary['all_gates_pass']:
        raise RuntimeError('T7.3 recovery gates not satisfied')


if __name__ == '__main__':
    main()
