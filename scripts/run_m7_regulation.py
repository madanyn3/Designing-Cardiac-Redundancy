"""New paired circulation plus metabolic replays for bounded regulation."""
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
import run_m7_global_hydraulics as hydraulic
import run_m7_global_metabolism as metabolic
from cardiac_redundancy.regulated_circulation import RegulatedCirculation

OUT = ROOT/'reports/m7-regulation'


def run(spec):
    class Circuit(RegulatedCirculation):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, regulation_gain=spec['gain'], **kwargs)
    case = json.loads((ROOT/'reports/m7-factorial/hydraulic_nominal/result.json').read_text())['case']
    case.update(name=spec['name'], hydraulic_source=str((OUT/spec['name']).relative_to(ROOT)))
    hydraulic.SensitivityCirculation = Circuit
    hydraulic.OUT = OUT
    hydraulic.run(case)
    metabolic.OUT = OUT
    result = metabolic.run(case)
    return dict(name=spec['name'], status=result['status'])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    specs = [dict(name='passive', gain=None), dict(name='gain2', gain=2.), dict(name='gain4', gain=4.)]
    paths = [Path(__file__), ROOT/'src/cardiac_redundancy/regulated_circulation.py',
             ROOT/'scripts/run_m7_global_hydraulics.py', ROOT/'scripts/run_m7_global_metabolism.py',
             ROOT/'docs/m7-final-structural-protocol.md', ROOT/'reports/m7-global/input-manifest.json']
    manifest = dict(specs=specs, hashes={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    hydraulic.write(OUT/'input-manifest.json', manifest)
    with ProcessPoolExecutor(max_workers=3) as pool:
        for future in as_completed([pool.submit(run, s) for s in specs]):
            print(future.result(), flush=True)


if __name__ == '__main__':
    main()
