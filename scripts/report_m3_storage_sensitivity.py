"""Numerical and scenario comparisons for the final storage/compliance audit."""
import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'reports/m3-storage-sensitivity'
def read(p):return json.loads(p.read_text())
def load(folder,r):
    result=[]
    for s in r['stages']:
        with np.load(folder/f"segment-{s['index']}.npz") as z:result.append({k:z[k].copy() for k in z.files})
    return result

def main():
    manifest=read(OUT/'manifest.json');comparisons={};summaries={};passed=True;maxdrift=0.
    for f in manifest['factors']:
        folder=OUT/f'factor-{f}';ref=read(folder/'nominal/result.json')
        if ref['status'] not in ('completed','model_domain_limit'):raise RuntimeError(str(ref))
        aa=load(folder/'nominal',ref)
        for label in ('refined','BDF'):
            r=read(folder/label/'result.json')
            if r['status'] not in ('completed','model_domain_limit'):raise RuntimeError(str(r))
            bb=load(folder/label,r);error=0.
            for a,b in zip(aa,bb):
                ta=a['time']-ref['start_s'];tb=b['time']-r['start_s'];mask=ta<=min(ta[-1],tb[-1])
                for key in ('states','pressures','flows'):
                    xa=a[key] if key=='states' else a[key].T;xb=b[key] if key=='states' else b[key].T
                    interp=np.array([np.interp(ta[mask],tb,row) for row in xb])
                    error=max(error,float(np.max(abs(xa[:,mask]-interp)/np.maximum(1,np.max(abs(xa[:,mask]),axis=1))[:,None])))
            ea=[e for s in ref['stages'] for e in s['events']];eb=[e for s in r['stages'] for e in s['events']]
            same=[(e['valve'],e['action']) for e in ea]==[(e['valve'],e['action']) for e in eb]
            dt=max([abs((a['time_s']-ref['start_s'])-(b['time_s']-r['start_s']))*1000 for a,b in zip(ea,eb)],default=0) if same else None
            term=(ref['termination'] or {}).get('node')==(r['termination'] or {}).get('node')
            endpoint=abs((aa[-1]['time'][-1]-ref['start_s'])-(bb[-1]['time'][-1]-r['start_s']))*1000
            good=bool(error<.005 and same and dt<1 and term and endpoint<1 and len(aa)==len(bb))
            comparisons[f'{f}-{label}']=dict(passed=good,error_pct=100*error,event_ms=dt,endpoint_ms=endpoint);passed &= good
            maxdrift=max(maxdrift,r['relative_blood_drift'])
        maxdrift=max(maxdrift,ref['relative_blood_drift'])
        t=np.hstack([d['time']-ref['start_s'] for d in aa]);q=np.vstack([d['flows'] for d in aa]);p=np.vstack([d['pressures'] for d in aa]);fn=ref['flow_names']
        mask=(t>=2)&(t<=4)
        avg=lambda x:float(np.trapezoid(x[mask],t[mask])/(t[mask][-1]-t[mask][0]))
        summaries[str(f)]=dict(status=ref['status'],warm_cycles=len(read(folder/'warmup.json')),initial_limb_inventory_mL=600*f,
             mean_AO_2_4_mmHg=avg(p[:,2]),brain_2_4_mL_min=60*avg(q[:,fn.index('brain')]),
             left_support_2_4_L_min=.06*avg(q[:,fn.index('support_left')]),limb_blood_at_4_mL=float(aa[1]['states'][6:16,-1].sum()))
    result=dict(passed=passed,comparisons=comparisons,scenarios=summaries,maximum_relative_blood_drift=maxdrift)
    (OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    lines=['# M3 limb storage/compliance sensitivity','',
      'Scale initial limb V, unstressed U and compliance C together by .5, 1 or 1.5; compensate central VC V and U so physical blood remains 5,250 mL. Resistances stay fixed. Each case warms independently. These are correlated engineering sensitivities, not measured anatomical bounds.',
      '', 'After arrest, support heads ramp to 20/100 mmHg in one second. Preservation conductances are .25; reopening occurs during seconds 4–6; native pumping stays off. Compare nominal/refined Radau and BDF through eight seconds or a model-domain stop.',
      '', '| Factor | Initial limb mL | Warm-up cycles | Mean AO 2–4 s mmHg | Brain 2–4 s mL/min | Left support 2–4 s L/min | Limb blood at 4 s mL |',
      '|---|---:|---:|---:|---:|---:|---:|']
    for f,s in summaries.items():lines.append(f"| {f} | {s['initial_limb_inventory_mL']:.1f} | {s['warm_cycles']} | {s['mean_AO_2_4_mmHg']:.4f} | {s['brain_2_4_mL_min']:.3f} | {s['left_support_2_4_L_min']:.5f} | {s['limb_blood_at_4_mL']:.3f} |")
    lines+=['',f'All declared comparisons pass: **{passed}**. Maximum relative blood drift {maxdrift:.3e}.',
            '', '| Comparison | Max scaled error % | Valve timing ms | Pass |','|---|---:|---:|---|']
    for n,c in comparisons.items():lines.append(f"| {n} | {c['error_pct']:.6f} | {c['event_ms']:.6f} | {c['passed']} |")
    lines+=['', 'Acceptance: <0.5% coordinate difference, scaled by max(1 native unit, reference peak); matched event sequence and domain node; valve/end timing <1 ms. Per-run invariants enforce blood drift <1e-8, septal residual <1e-8 mmHg, state continuity and valve directionality.',
            '', 'This completes the prescribed V/U/C scaling check for this declared support scenario. Dependence of retained limb blood and delivered flows on storage parameters must remain in later uncertainty analyses. The experiment does not establish a preservation interval, fully isolated anatomy or an optimal pump.',
            '', '[Anatomical disposition](../../docs/m3-anatomical-assumptions.md). Rerun `scripts/run_m3_storage_sensitivity.py`, then `scripts/report_m3_storage_sensitivity.py`. Effective configurations, warm-up histories, seeds, solver runs and implementation snapshots are saved here.']
    (OUT/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2));return 0 if passed else 1
if __name__=='__main__':raise SystemExit(main())
