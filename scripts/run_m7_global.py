"""Staged Morris experiment with exact hydraulic reuse and explicit failures."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import json,sys,time,hashlib,os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[name]='1'
import run_m7_global_hydraulics as hyd
import run_m7_global_metabolism as met
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'reports/m7-global'
def run_hyd(case):
    dest=hyd.OUT/case['name'];p=dest/'hydraulics.json'
    if p.exists() and met.read(p)['status']=='verified':return None
    start=time.perf_counter()
    try:
        hyd.run(case);met.write(dest/'completion.json',dict(status='complete',elapsed_s=time.perf_counter()-start))
    except Exception as exc:
        met.write(dest/'completion.json',dict(status='error',error=repr(exc),elapsed_s=time.perf_counter()-start));return repr(exc)
def run_met(case):
    dest=OUT/case['name'];p=dest/'result.json'
    if p.exists() and met.read(p)['status']=='verified':return None
    start=time.perf_counter()
    try:
        # Hydraulic reuse is valid only for identical settings in this fixed
        # scenario. Metabolic knobs cannot alter the prescribed mechanics.
        h=met.read(ROOT/case['hydraulic_source']/'hydraulics.json')
        for key in ('preservation_conductance','contractility','coronary_conductance'):
            if abs(h['case']['knobs'][key]-case['knobs'][key])>1e-10:raise ValueError('Mismatched hydraulic cache')
        try:
            result=met.run(case)
        except RuntimeError as first:
            # Keep the failed attempt; retry numerically at stricter settings.
            # This never bypasses a conservation/capacity/domain check.
            attempt=met.read(dest/'result.json') if (dest/'result.json').exists() else dict(error=repr(first))
            met.write(dest/'first-attempt.json',attempt)
            result=met.run(dict(case,numerical_refinement=True))
        if result['status']!='verified':raise RuntimeError(result['status'])
        met.write(dest/'completion.json',dict(status='complete',elapsed_s=time.perf_counter()-start))
    except Exception as exc:
        met.write(dest/'completion.json',dict(status='error',error=repr(exc),elapsed_s=time.perf_counter()-start));return repr(exc)
if __name__=='__main__':
    n=int(sys.argv[1]) if len(sys.argv)>1 else 8
    cfg=met.read('configs/m7-global.json')
    if n not in cfg['stages']:raise ValueError('Unregistered stage')
    frozen=met.read(OUT/'input-manifest.json')['sha256']
    changed=[p for p,h in frozen.items() if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h]
    if changed:raise RuntimeError('Inputs changed '+str(changed))
    cases=[c for c in cfg['cases'] if c['trajectory']<n];needed={c['hydraulic_source'] for c in cases}
    hydraulics=[c for c in cfg['hydraulics'] if 'reports/m7-global/hydraulics/'+c['name'] in needed]
    with ProcessPoolExecutor(max_workers=6) as pool:errors=list(pool.map(run_hyd,hydraulics))
    if any(errors):raise SystemExit('Hydraulic errors retained; metabolism not extrapolated')
    with ProcessPoolExecutor(max_workers=6) as pool:errors=list(pool.map(run_met,cases))
    if any(errors):raise SystemExit('Metabolic errors retained; no successful-subset analysis')
    met.write(OUT/f'stage-{n}.json',dict(status='complete',trajectories=n,cases=len(cases),new_hydraulic_cases=len(hydraulics)))
