"""Freeze the first structural block, explicitly leaving T7.4/T7.5 open."""
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'reports/m7-structural-replay'


def main():
    summary = json.loads((OUT/'summary.json').read_text())
    if len(summary) != 6:
        raise RuntimeError('Six structural cases required')
    for row in summary:
        result = json.loads((OUT/row['name']/'result.json').read_text())
        if result['status'] != 'verified':
            raise RuntimeError('Unverified structural case')
    verified = {}
    for path, field in [('reports/m7-identifiability/recovery/manifest.json', 'sha256'),
                         ('reports/m7-global/input-manifest.json', 'sha256'),
                         ('reports/m7-structural-replay/input-manifest.json', 'hashes')]:
        hashes = json.loads((ROOT/path).read_text())[field]
        for source, expected in hashes.items():
            if hashlib.sha256((ROOT/source).read_bytes()).hexdigest() != expected:
                raise RuntimeError(f'Changed input: {source}')
        verified[path] = len(hashes)
    suites = ET.parse(ROOT/'reports/m7-structural-regression-tests.xml').getroot().findall('testsuite')
    if not suites or any(int(s.get('errors', 0))+int(s.get('failures', 0)) for s in suites):
        raise RuntimeError('Regression evidence failed')
    paths = [p for p in OUT.rglob('*') if p.is_file() and p.name != 'manifest.json']
    paths += [ROOT/'scripts/report_m7_structural_replay.py', Path(__file__),
              ROOT/'tests/test_structural_closures.py', ROOT/'reports/m7-structural-regression-tests.xml']
    result = dict(status='first structural block verified; T7.4/T7.5 remain open',
                  tests=sum(int(s.get('tests', 0)) for s in suites), prior_verified=verified,
                  sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (OUT/'manifest.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k != 'sha256'}, indent=2))


if __name__ == '__main__':
    main()
