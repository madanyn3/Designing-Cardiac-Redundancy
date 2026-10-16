"""Freeze a qualified T7.3 closure only after all declared recovery gates pass."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'reports/m7-identifiability/recovery'


def main():
    summary = json.loads((OUT/'summary.json').read_text())
    if len(summary['results']) != 7 or not summary['all_gates_pass']:
        raise RuntimeError('Incomplete recovery evidence')
    verified = {}
    for relative, field in [('reports/m7-global/manifest.json', 'sha256'),
                            ('reports/m7-global/input-manifest.json', 'sha256'),
                            ('reports/m7-identifiability/pulsatile/manifest.json', 'sha256'),
                            ('reports/m7-identifiability/recovery/input-manifest.json', 'hashes')]:
        hashes = json.loads((ROOT/relative).read_text())[field]
        for path, expected in hashes.items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != expected:
                raise RuntimeError(f'Changed evidence: {path}')
        verified[relative] = len(hashes)
    paths = [p for p in OUT.rglob('*') if p.is_file() and p.name != 'manifest.json']
    paths += [ROOT/'scripts'/name for name in ['report_m7_recovery.py', 'finalize_m7_t73.py']]
    paths += [ROOT/'docs/m7-recovery-protocol.md', ROOT/'docs/m7-t73-closure.md']
    manifest = dict(status='T7.3 closed at conditional identifiability scope', verified_prior_hash_counts=verified,
                    sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (OUT/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status':manifest['status'], 'files':len(paths), 'verified_prior_hash_counts':verified}, indent=2))


if __name__ == '__main__':
    main()
