"""Quantify F0 transient convergence; physiological survival is not evaluated."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reports/m3-arrest'


def read(path): return json.loads(path.read_text())


def main():
    protocol=read(OUT/'manifest.json')['protocol']
    results={c['name']:read(OUT/c['name']/'result.json') for c in protocol['cases']}
    comparisons={}; passed=True
    for name,r in results.items():
        if r['status'] not in ('completed','model_domain_limit'):
            comparisons[name]={'passed':False,'error':r.get('error',r['status'])}; passed=False; continue
        if name.endswith('-nominal'): continue
        reference=name.rsplit('-',1)[0]+'-nominal'; base=results[reference]
        if base['status'] not in ('completed','model_domain_limit'):
            comparisons[name]={'passed':False,'error':'reference incomplete'}; passed=False; continue
        a=np.load(OUT/reference/f"segment-{base['post_segment']}.npz")
        b=np.load(OUT/name/f"segment-{r['post_segment']}.npz")
        ta=a['time_s']-base['arrest_time_s']; tb=b['time_s']-r['arrest_time_s']
        # Different localized domain endpoints are compared explicitly. Never
        # extrapolate a terminated solution beyond its valid time interval.
        mask=ta<=min(ta[-1],tb[-1]); tc=ta[mask]
        errors={}
        for key in ('states','pressures_mmHg','flows_mL_s'):
            xa=a[key] if key=='states' else a[key].T
            xb=b[key] if key=='states' else b[key].T
            interp=np.array([np.interp(tc,tb,row) for row in xb])
            scale=np.maximum(1,np.max(abs(xa[:,mask]),axis=1))
            errors[key]=float(np.max(abs(interp-xa[:,mask])/scale[:,None]))
        ea=read(OUT/reference/'valves.json'); eb=read(OUT/name/'valves.json')
        same=[(e['valve'],e['action']) for e in ea]==[(e['valve'],e['action']) for e in eb]
        dt=max([abs(x['time_s']-y['time_s'])*1000 for x,y in zip(ea,eb)],default=0) if same else None
        term_match=(base['termination'] is None and r['termination'] is None) or (
            base['termination'] is not None and r['termination'] is not None and
            base['termination']['node']==r['termination']['node'])
        endpoint_ms=abs(ta[-1]-tb[-1])*1000
        ok=bool(max(errors.values())<.005 and same and dt<1 and term_match and endpoint_ms<1)
        comparisons[name]={'scaled_maximum_errors':errors,'valve_sequence_matches':same,
            'valve_time_difference_ms':dt,'termination_matches':term_match,'endpoint_difference_ms':endpoint_ms,'passed':ok}
        passed &= ok
    rows=['# M3 F0 — scheduled native arrest without support','',
        '2026-10-06. The user closed the cerebral/coronary pressure-law task. This next task removes native LV/RV/septal activation at a prescribed time, retaining passive mechanics, all blood and the inertial valve-flow states. No support, oxygen transport, injury or recovery model is present.','',
        '## Protocol and provenance','',
        'Five arrest phases are each tested with nominal Radau, refined Radau and BDF, over ten seconds after arrest. Two phases bracket a recorded aortic-valve closing time by 10 microseconds. All start from the same recorded periodic nominal pressure-aware state; its original model/configuration hashes are checked before reuse. The new intervention wrapper does not change those constitutive laws. Files and the seed are preserved in the implementation snapshot.','',
        'The ten-second horizon is an engineering observation window, not a supported survival duration. Each segment retains its own activation flag; the intervention record saves identical states and separately evaluated pressures/flows on either side. Pressure jumps are allowed by abrupt activation removal. There is no interpolation across that jump.','',
        '## Nominal outcomes','',
        '| Arrest phase | Phase s | Observed post-arrest s | Final aortic pressure mmHg | Final brain flow mL/min | Result |',
        '|---|---:|---:|---:|---:|---|']
    for name,r in results.items():
        if not name.endswith('-nominal'):continue
        if 'final_AO_mmHg' in r:
            rows.append(f"| {name[:-8]} | {r['case']['phase_s']:.8f} | {r['final_elapsed_s']:.4f} | {r['final_AO_mmHg']:.4f} | {r['final_brain_flow_mL_min']:.4f} | {r['status']} |")
        else: rows.append(f"| {name} | — | — | — | — | {r['status']} |")
    rows += ['', 'Aortic pressure here is instantaneous at the final observation; it is not a cardiac-cycle MAP after pumping stops. Nonzero transient flow can come from stored pressure/volume and inertial motion despite zero active contraction. Zero cerebral flow is a result of this model’s closing-pressure law, not an independently validated brain-injury time.','',
        '## Numerical verification','',
        'Compare post-event states, pressures and flows on aligned elapsed times. For each coordinate, maximum absolute difference is divided by max(1 in its native unit, its nominal peak magnitude); the limit is 0.5%. This provides an absolute error scale near zero, without unstable pointwise relative errors. Matched valve sequences must agree within 1 ms; domain endpoint identity and timing are checked separately. These are numerical criteria, not biological thresholds.','',
        '| Comparison | Maximum scaled error % | Valve timing difference ms | Pass |',
        '|---|---:|---:|---|']
    for name,c in comparisons.items():
        value=max(c.get('scaled_maximum_errors',{}).values(),default=float('nan'))*100
        rows.append(f"| {name} | {value:.6g} | {c.get('valve_time_difference_ms')} | {c['passed']} |")
    valid=[r for r in results.values() if 'relative_volume_drift' in r]
    rows += ['',f'Overall declared numerical comparisons pass: **{passed}**.',
        f"Maximum relative blood-volume drift: {max((r['relative_volume_drift'] for r in valid),default=float('nan')):.3e}.",
        f"Maximum continuous-state jump at arrest: {max((r['state_jump'] for r in valid),default=float('nan')):.3e}.",
        '', 'Automated tests additionally check uninterrupted versus unchanged-restart agreement, high-activation arrest continuity, conservation and suppression of scheduled arrest after an initial domain violation. Full results are in [the JUnit record](../m3-arrest-tests.xml). This does not close all simultaneous-event, support-port or domain stress cases.','',
        '## Representative transient','', '![F0 transient](transient.png)','',
        '## Rerun and next task','',
        'Run `.venv/Scripts/python.exe scripts/run_m3_arrest.py`, then `.venv/Scripts/python.exe scripts/report_m3_arrest.py`. The protocol is [m3-arrest-protocol.json](../../configs/m3-arrest-protocol.json); orchestration is [interventions.py](../../src/cardiac_redundancy/interventions.py).','',
        'Next implement B1 abstract right/left support ports, with finite head/flow, conservative transfers, ramp and trigger timing, and reservoir-domain checks. Then selective and oxygenator/unloading topologies, followed by the chamber-based auxiliary fixture. The source-paper perturbation discrepancy remains open; M3 as a whole is not closed.','']
    r=results['peak_activation-nominal']
    if 'post_segment' in r:
        d=np.load(OUT/'peak_activation-nominal'/f"segment-{r['post_segment']}.npz")
        # Recorded branch order is immutable for this configuration.
        import sys
        sys.path.insert(0,str(ROOT/'src'))
        from cardiac_redundancy.pressure_regional import PressureRegionalModel
        m=PressureRegionalModel(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml')
        t=d['time_s']-r['arrest_time_s']; fig,axes=plt.subplots(3,1,figsize=(9,8),sharex=True)
        for n in ('AO','VC','PA','PU'):axes[0].plot(t,d['pressures_mmHg'][:,m.volume_names.index(n)],label=n)
        axes[0].set_ylabel('Pressure (mmHg)'); axes[0].legend(ncol=4)
        for n in ('brain','coronary_left','AV'):axes[1].plot(t,d['flows_mL_s'][:,m.flow_names.index(n)]*60,label=n)
        axes[1].set_ylabel('Flow (mL/min)'); axes[1].legend()
        for n in ('LV','RV'):axes[2].plot(t,d['states'][m.volume_names.index(n)],label=n)
        axes[2].set_ylabel('Cavity volume (mL)'); axes[2].set_xlabel('Time after native arrest (s)'); axes[2].legend()
        fig.suptitle('No support: arrest at peak native activation'); fig.tight_layout()
        fig.savefig(OUT/'transient.png',dpi=150); plt.close(fig)
    (OUT/'report.md').write_text('\n'.join(rows),encoding='utf-8')
    (OUT/'verification.json').write_text(json.dumps({'passed':passed,'comparisons':comparisons},indent=2)+'\n')
    print(json.dumps({'passed':passed,'comparisons':comparisons},indent=2))
    return 0 if passed else 1


if __name__=='__main__':raise SystemExit(main())
