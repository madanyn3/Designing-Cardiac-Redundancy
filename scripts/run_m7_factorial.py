"""Run all preregistered factorial cases; preserve errors and resume verified cases."""
from pathlib import Path
import json,sys,time,hashlib
from concurrent.futures import ProcessPoolExecutor
ROOT=Path(__file__).resolve().parents[1]
import run_m7_factorial_hydraulics as hyd
import run_m7_factorial_metabolism as met
OUT=ROOT/'reports/m7-factorial'

def run(case):
    dest=OUT/case['name'];dest.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter()
    result=dest/'result.json'
    if result.exists() and json.loads(result.read_text(encoding='utf-8'))['status']=='verified':return
    try:
        if case['block']=='hydraulic':hyd.run(case)
        met.run(case)
        met.write(dest/'completion.json',dict(status='complete',elapsed_s=time.perf_counter()-start))
    except Exception as exc:
        met.write(dest/'completion.json',dict(status='error',error=repr(exc),elapsed_s=time.perf_counter()-start))
        return str(exc)

if __name__=='__main__':
    frozen=json.loads((OUT/'input-manifest.json').read_text(encoding='utf-8'))['sha256']
    changed=[p for p,h in frozen.items() if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h]
    if changed:raise RuntimeError('Frozen inputs changed: '+str(changed))
    cases=json.loads((ROOT/'configs/m7-factorial.json').read_text(encoding='utf-8'))['cases']
    if len(sys.argv)>1:cases=[c for c in cases if c['name'] in sys.argv[1:]]
    with ProcessPoolExecutor(max_workers=3) as pool:errors=list(pool.map(run,cases))
    if any(errors):raise SystemExit('Some cases failed; inspect completion.json')
