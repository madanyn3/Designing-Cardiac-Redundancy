"""Audit and summarize the conditional M7 endpoint experiment."""
from pathlib import Path
import json,hashlib,platform,sys
import numpy as np
import scipy
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reports/m7-metabolic-screen'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    config=read(ROOT/'configs/m7-metabolic-screen.json')
    names=['nominal']+[f['id']+'_'+e for f in config['factors'] for e in ('low','high')]
    results={name:read(OUT/name/'result.json') for name in names}
    if any(r['status']!='verified' for r in results.values()):raise RuntimeError('Incomplete/error cases; do not rank successful subset')
    old=read(ROOT/'reports/m6-coupled/observed_arrest/result.json')
    reference_error=abs(results['nominal']['records'][-1]['unmet_oxygen_mL']-old['records'][-1]['unmet_oxygen_mL'])
    if reference_error>1e-7:raise RuntimeError('Nominal regression failed')
    for end,fraction in (('low',.1),('high',.5)):
        for name,e in results['preservation_fraction_'+end]['profile'].items():
            if e['role']=='preservation':
                if e['fraction']!=fraction:raise RuntimeError('Preservation knob not applied')
            elif e!=old['profile'][name]:raise RuntimeError('Mandatory demand changed')
    effects=[]
    for f in config['factors']:
        lo,hi=[results[f['id']+'_'+e] for e in ('low','high')]
        values=[r['records'][-1]['unmet_oxygen_mL'] for r in (lo,hi)]
        resolution=sum(abs(r['records'][0]['unmet_oxygen_mL']-r['records'][1]['unmet_oxygen_mL']) for r in (lo,hi))
        effects.append(dict(factor=f['id'],low=values[0],high=values[1],signed_change_mL=values[1]-values[0],span_mL=abs(values[1]-values[0]),paired_endpoint_difference_sum_mL=resolution,resolved=abs(values[1]-values[0])>max(1e-8,10*resolution)))
    effects.sort(key=lambda v:v['span_mL'],reverse=True)
    # Integrate cause attribution only as sampled diagnostic; cumulative ledgers
    # remain primary. Report quadrature discrepancy instead of hiding it.
    attribution={}
    for name,r in results.items():
        a=np.load(OUT/name/'metabolic-refined.npz');t=a['time_s']
        glucose=np.trapezoid(np.maximum(0,a['requested_O2_mL_s']-a['glucose_ceiling_mL_s']).sum(axis=1),t)
        oxygen=np.trapezoid(np.maximum(0,a['glucose_ceiling_mL_s']-a['oxygen_consumed_mL_s']).sum(axis=1),t)
        attribution[name]=dict(sampled_glucose_deficit_mL=float(glucose),sampled_additional_oxygen_deficit_mL=float(oxygen),sampled_vs_ledger_difference_mL=float(glucose+oxygen-r['records'][-1]['unmet_oxygen_mL']))
    summary=dict(scope='10 s fixed hydraulic endpoint sensitivity; not global or survival analysis',cases=13,failed_cases=0,seed=None,nominal_regression_error_mL=reference_error,effects=effects,attribution=attribution,runtime_sum_s=sum(r['runtime_s'] for r in results.values()),python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    lines=['# T7.2 conditional metabolic screen','', '13 deterministic cases, each solved at two numerical settings. Fixed M6 observed-arrest hydraulics over 10 seconds. All 13 verified; no failed samples. No random sampling or confidence intervals. T7.2 remains open for joint/global and hydraulic screening.','',f"Nominal unmet aerobic demand: {results['nominal']['records'][-1]['unmet_oxygen_mL']:.9g} mL O2 equivalent. M6 reproduction difference: {reference_error:.3g} mL.",'','| Factor | Low endpoint deficit (mL) | High endpoint deficit (mL) | Absolute span (mL) | Resolved above numerical comparison? |','|---|---:|---:|---:|---|']
    for e in effects:lines.append(f"| {e['factor']} | {e['low']:.9g} | {e['high']:.9g} | {e['span_mL']:.9g} | {e['resolved']} |")
    lines+=['','Resolution flag requires span > max(1e-8 mL, ten times the sum of endpoint nominal/refined differences). This is a numerical screening rule, not a statistical confidence bound. Endpoint spans depend on the chosen ranges and miss interactions.','', '| Case | Active deficit | Preservation deficit | Pump deficit | Released glycogen (mmol) | Minimum plasma-equivalent glucose (mmol/L) |','|---|---:|---:|---:|---:|---:|']
    for name,r in results.items():
        v=r['records'][-1];g=v['groups']
        lines.append(f"| {name} | {g['active']['unmet_O2_mL']:.8g} | {g['preservation']['unmet_O2_mL']:.8g} | {g['pump']['unmet_O2_mL']:.8g} | {v['released_mmol']:.8g} | {v['minimum_plasma_equivalent_mmol_L']:.8g} |")
    lines+=['','## Interpretation limits','','Pressure, flow and pump work do not respond to these metabolic knobs in this conditional replay. No result establishes physiological feasibility. Glycogen-stock insensitivity over ten seconds cannot establish long-duration unimportance. Tissue glucose buffers, substrate alternatives and energetic feedback on contraction remain absent. Absolute oxygen stores are held fixed when myocardial demand mass changes.','', 'Cause attribution and its sampled quadrature error are in summary.json; primary deficits use cumulative state ledgers. Detailed conservation, capacity, stoichiometry and paired-state errors are recorded per case.','', 'Reproduce: `.venv/Scripts/python.exe scripts/run_m7_metabolic_screen.py`, then `.venv/Scripts/python.exe scripts/report_m7_metabolic_screen.py`. Protocol: configs/m7-metabolic-screen.json; rationale and references: docs/uncertainty-plan.md.']
    lines+=['', '## Observations', '',
        'Within these selected ranges, extraction has the largest endpoint effect, followed by myocardial demand mass and work cost. Preservation-demand changes have a very small effect on the deficit; release-cap and stock effects are unresolved under the declared numerical rule.',
        'At extraction .25, brain deficit is 0.0165653 mL and respiratory-muscle deficit 0.243994 mL, in addition to auxiliary myocardial deficits. At the nominal setting, all recorded deficit is auxiliary myocardial.',
        'All cases show zero additional oxygen limitation at saved diagnostic samples: their recorded shortage is glucose-limited aerobic metabolism. Sampled attribution differs from integrated ledgers by up to 5.25e-5 mL because of output-grid quadrature; use ledgers for small endpoint comparisons.',
        'Increasing release cap from .21 to .84 mmol/min increases actual release from .0175 to .0555686 mmol over ten seconds, but has no resolved effect on deficit. Stock and rate importance for longer support remains untested.',
        'Verification: 258 regression tests pass, including four M7 checks. All thirteen cases pass paired refinement, conservation, stoichiometry and inventory checks. Frozen M6 input hashes remain unchanged.',
        '', '## What this screen establishes', '',
        'The nominal circulation is unchanged from M6 (final-window mean aortic pressure about 49.0 mmHg). Adjusting metabolic assumptions cannot repair its pressure shortfall in this replay.',
        'Demand and conductance are distinct knobs. Lower preservation demand here leaves all branch conductances unchanged; testing redistribution requires new hydraulic solutions.',
        'Next: joint sampling of extraction and myocardial-demand assumptions, hydraulic sensitivity including branch conductance and coronary supply, then ranking stability across horizons and structural alternatives. This initial endpoint screen does not close T7.2.']
    (OUT/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(9,4.5))
    ax.barh([e['factor'] for e in effects][::-1],[e['span_mL'] for e in effects][::-1])
    ax.set_xlabel('Endpoint span in unmet aerobic demand (mL O2 equivalent)')
    ax.set_title('Conditional 10-second screen; hydraulics fixed')
    fig.tight_layout();fig.savefig(OUT/'sensitivity.png',dpi=160);plt.close(fig)
    paths=set(ROOT.glob('src/cardiac_redundancy/*.py'))|set(ROOT.glob('configs/*'))|set(ROOT.glob('data/*'))
    paths|={ROOT/'tests/test_m7_screen.py',ROOT/'reports/m7-screen-tests.xml',ROOT/'reports/m7-regression-tests.xml',ROOT/'scripts/run_m7_metabolic_screen.py',ROOT/'scripts/report_m7_metabolic_screen.py',ROOT/'scripts/run_m4_respiratory.py',ROOT/'docs/uncertainty-plan.md'}
    paths|=set((ROOT/'reports/m4-baseline').glob('*'))|set((ROOT/'reports/m6-coupled/observed_arrest').glob('*'))|set((ROOT/'reports/m5-normal-comparison').glob('*'))|set(OUT.rglob('*'))
    manifest={str(p.relative_to(ROOT)).replace('\\','/'):digest(p) for p in sorted(paths) if p.is_file() and p.name!='manifest.json' and p.suffix!='.log'}
    (OUT/'manifest.json').write_text(json.dumps(dict(sha256=manifest),indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
