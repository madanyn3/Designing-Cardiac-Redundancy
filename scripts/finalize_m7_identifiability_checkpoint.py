"""Verify and freeze the local T7.3 checkpoint, without declaring T7.3 closed."""
import hashlib
import json
from pathlib import Path
import platform
import importlib.metadata
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'reports/m7-identifiability/pulsatile'


def main():
    manifest = json.loads((OUT/'input-manifest.json').read_text())
    for path, expected in manifest['hashes'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != expected:
            raise RuntimeError(f'Changed input: {path}')
    results = [json.loads((OUT/case[0]/'result.json').read_text()) for case in manifest['cases']]
    if any(r['status'] != 'verified' for r in results):
        raise RuntimeError('Unverified cases')
    summary = json.loads((OUT/'summary.json').read_text())
    checks = {k: v['step_relative_frobenius'] < .02 and v['solver_relative_frobenius'] < .02
              for k, v in summary['observations'].items()}
    if not all(checks.values()):
        raise RuntimeError('Jacobian stability requires investigation')
    suites = ET.parse(ROOT/'reports/m7-identifiability-regression-tests.xml').getroot().findall('testsuite')
    if not suites or any(int(s.get('failures', 0))+int(s.get('errors', 0)) for s in suites):
        raise RuntimeError('Regression tests incomplete or failed')
    checkpoint = dict(status='local experiment verified; T7.3 remains open', case_count=len(results),
        tests=sum(int(s.get('tests', 0)) for s in suites), jacobian_stability=checks,
        sum_case_wall_seconds=sum(r['elapsed_s'] for r in results),
        python=platform.python_version(), packages={n: importlib.metadata.version(n) for n in ['numpy','scipy','pytest']},
        remaining=['nonlinear synthetic recovery', 'profile/noise analysis', 'final conditional-identifiability decisions'])
    (OUT/'checkpoint.json').write_text(json.dumps(checkpoint, indent=2)+'\n', encoding='utf-8')
    paths = [p for p in OUT.rglob('*') if p.is_file() and p.name != 'manifest.json']
    paths += [ROOT/'scripts'/name for name in ['report_m7_pulsatile_identifiability.py', 'finalize_m7_identifiability_checkpoint.py']]
    paths += [ROOT/'tests/test_m7_pulsatile_identifiability.py', ROOT/'docs/m7-pulsatile-identifiability-protocol.md',
              ROOT/'reports/m7-identifiability-regression-tests.xml']
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (OUT/'manifest.json').write_text(json.dumps({'sha256': hashes}, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(checkpoint, indent=2))


if __name__ == '__main__':
    main()
