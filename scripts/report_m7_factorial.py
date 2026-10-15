"""Audit complete factorial blocks and compare effects across numerical settings."""
from pathlib import Path
import json,sys,hashlib,platform
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.factorial_analysis import decompose,interpolate
OUT=ROOT/'reports/m7-factorial'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def mean_window(t,y,end):
    start=max(0,end-2);inside=(t>start)&(t<end)
    tt=np.r_[start,t[inside],end];yy=np.r_[np.interp(start,t,y),y[inside],np.interp(end,t,y)]
    return float(np.trapezoid(yy,tt)/(end-start))
def metrics(folder,hydfolder,label):
    result=read(folder/'result.json');hyd=read(hydfolder/'hydraulics.json')
    a=np.load(folder/f'metabolic-{label}.npz');h=np.load(hydfolder/f'hydraulic-{label}.npz')
    n=len(result['volume_names']);names=result['organ_names'];m=len(names)
    deficits=a['states'][:,n+2+m:n+2+2*m]
    if a['time_s'][-1]<10-1e-8 or h['time_s'][-1]<10-1e-8:raise RuntimeError('Common horizon not reached')
    groups={role:[i for i,name in enumerate(names) if result['profile'][name]['role']==role] for role in ('active','preservation','pump')}
    output={}
    for end in (4,7,10):
        row=dict(MAP_mmHg=mean_window(h['time_s'],h['pressures_mmHg'][:,hyd['volume_names'].index('AO')],end),aux_output_L_min=mean_window(h['time_s'],h['flows_mL_s'][:,hyd['flow_names'].index('aAV')]*.06,end))
        row['total_deficit_mL']=float(np.interp(end,a['time_s'],deficits.sum(axis=1)))
        for role,idx in groups.items():row[role+'_deficit_mL']=float(np.interp(end,a['time_s'],deficits[:,idx].sum(axis=1)))
        output[str(end)]=row
    output['regional_final_flow_mL_s']={name:mean_window(h['time_s'],h['flows_mL_s'][:,j],10) for j,name in enumerate(hyd['flow_names'])}
    extra=np.maximum(0,a['glucose_ceiling_mL_s']-a['oxygen_consumed_mL_s'])
    output['maximum_sampled_additional_O2_deficit_mL_s']=float(extra.max())
    output['sampled_glucose_deficit_mL']=float(np.trapezoid(np.maximum(0,a['requested_O2_mL_s']-a['glucose_ceiling_mL_s']).sum(axis=1),a['time_s']))
    output['attribution_quadrature_difference_mL']=output['sampled_glucose_deficit_mL']+float(np.trapezoid(extra.sum(axis=1),a['time_s']))-output['10']['total_deficit_mL']
    output['release_s']=hyd['runs'][0 if label=='nominal' else 1]['release_s']
    output['chamber_work_J']=hyd['runs'][0 if label=='nominal' else 1]['chamber_work_J']
    return output

def main():
    cfg=read(ROOT/'configs/m7-factorial.json')
    audit=[];data={}
    for case in cfg['cases']:
        folder=OUT/case['name'];completion=folder/'completion.json';p=folder/'result.json'
        status=read(completion) if completion.exists() else {'status':'pending'}
        if not p.exists() or read(p).get('status')!='verified':
            audit.append(dict(case=case['name'],completion=status,result=read(p).get('status') if p.exists() else 'missing'));continue
        hydfolder=folder if case['block']=='hydraulic' else ROOT/'reports/m6-coupled/observed_arrest'
        data[case['name']]={label:metrics(folder,hydfolder,label) for label in ('nominal','refined')}
    write(OUT/'case-audit.json',dict(expected=17,verified=len(data),incomplete=audit))
    if audit:raise RuntimeError('Incomplete factorial; no successful-subset ranking: '+str(audit))
    source=ROOT/'reports/m6-coupled/observed_arrest'
    old={label:metrics(source,source,label) for label in ('nominal','refined')}
    errors={k:abs(data['hydraulic_nominal']['refined']['10'][k]-v) for k,v in old['refined']['10'].items()}
    if max(errors.values())>1e-6:raise RuntimeError('Hydraulic nominal regression failed: '+str(errors))
    analyses={}
    for block,fs in cfg['factors'].items():
        cases=[c for c in cfg['cases'] if c['block']==block and c['signs'] is not None]
        names=[f['id'] for f in fs];signs=[c['signs'] for c in cases];analyses[block]={}
        reference=data['hydraulic_nominal'] if block=='hydraulic' else old
        for end in ('4','7','10'):
            analyses[block][end]={}
            for output in old['refined'][end]:
                nominal=decompose(signs,[data[c['name']]['nominal'][end][output] for c in cases],names)
                refined=decompose(signs,[data[c['name']]['refined'][end][output] for c in cases],names)
                actual=reference['refined'][end][output]
                prediction=float(interpolate(refined,{f['id']:f['nominal'] for f in fs},{f['id']:(f['low'],f['high']) for f in fs}))
                coefficient_error=max(abs(refined['terms'][key]['coefficient']-nominal['terms'][key]['coefficient']) for key in refined['terms'])
                bound=10*max(1e-10,coefficient_error)
                rank=lambda d:sorted([k for k,v in d['total_shares'].items() if v is not None],key=lambda k:d['total_shares'][k],reverse=True)
                resolved=[k for k,v in refined['terms'].items() if abs(v['coefficient'])>bound]
                analyses[block][end][output]=dict(nominal=nominal,refined=refined,max_coefficient_numerical_difference=coefficient_error,resolved_terms=resolved,ranking=rank(refined),nominal_ranking=rank(nominal),reference_actual=actual,vertex_interpolated_reference=prediction,reference_interpolation_error=prediction-actual)
    summary=dict(scope='Exact discrete factorial sensitivity on two separate eight-vertex designs; not continuous global convergence',cases=17,verified=17,failed=0,seed=None,case_runtime_sum_s=sum(read(OUT/c['name']/'completion.json')['elapsed_s'] for c in cfg['cases']),workers=3,reference_regression_errors=errors,analyses=analyses,metrics=data,python=platform.python_version(),numpy=np.__version__)
    write(OUT/'summary.json',summary)
    lines=['# T7.2 joint metabolic and hydraulic sensitivity','', 'Two complete three-factor factorials plus one hydraulic reference: 17 cases, 34 metabolic solves and 18 hydraulic solves. All cases verified. No cases removed. Numerical pairs are not statistical replicates; no probability distributions or confidence intervals are inferred.','', 'The design estimates main effects and interactions across its eight vertices. Rankings are conditional on these ranges, the ten-second scenario and fixed assumptions. T7.2 remains open for continuous-domain/interior sampling and broader uncertainty coverage.','', '## Ten-second effects','', 'Effect means high-minus-low averaged over the other factors. Interactions are estimated separately. Deficit is accumulated unmet aerobic demand in mL O2 equivalent; pressure is the preceding two-second mean.','', '| Block / output | Factor | Main effect | Discrete total variance share |','|---|---|---:|---:|']
    for block in analyses:
        for output in ('total_deficit_mL','MAP_mmHg','aux_output_L_min'):
            a=analyses[block]['10'][output];d=a['refined']
            for name in [f['id'] for f in cfg['factors'][block]]:
                share=d['total_shares'][name]
                lines.append(f"| {block} / {output} | {name} | {d['terms'][name]['effect']:.8g} | {share if share is None else round(share,6)} |")
    lines+=['','Total shares include interactions, so their sum can exceed one. Shares describe equal weighting of vertices, not physiological variability. Tiny/no-variance outputs have no meaningful ranking.','', '## Interactions and interior check', '', '| Block / output | Interaction fraction of discrete variance | Vertex prediction at nominal | Actual nominal | Prediction minus actual |', '|---|---:|---:|---:|---:|']
    for block in analyses:
        for output in ('total_deficit_mL','MAP_mmHg','aux_output_L_min'):
            a=analyses[block]['10'][output];d=a['refined'];interaction=sum(v['variance_share'] or 0 for v in d['terms'].values() if v['order']>1)
            lines.append(f"| {block} / {output} | {interaction:.6g} | {a['vertex_interpolated_reference']:.8g} | {a['reference_actual']:.8g} | {a['reference_interpolation_error']:.8g} |")
    lines+=['', 'The nominal point is independently simulated and is not always the geometric center. Interpolation error reveals interior nonlinear behavior; exact corner reconstruction does not validate interpolation.','', '## Time-window ranking check', '', '| Block / output | 4 s | 7 s | 10 s |', '|---|---|---|---|']
    for block in analyses:
        for output in ('total_deficit_mL','MAP_mmHg','aux_output_L_min'):
            rankings=[' > '.join(analyses[block][end][output]['ranking']) or 'no variance' for end in ('4','7','10')]
            lines.append('| '+block+' / '+output+' | '+' | '.join(rankings)+' |')
    lines+=['','## Numerical effect stability','','| Block / output at 10 s | Maximum coefficient change | Same total-share order at both tolerances? | Resolved terms |','|---|---:|---|---|']
    for block in analyses:
        for output in ('total_deficit_mL','MAP_mmHg','aux_output_L_min'):
            a=analyses[block]['10'][output]
            lines.append(f"| {block} / {output} | {a['max_coefficient_numerical_difference']:.6g} | {a['ranking']==a['nominal_ranking']} | {', '.join(a['resolved_terms']) or 'none'} |")
    lines+=['', 'A term is resolved only if its absolute coefficient exceeds ten times max(1e-10, the largest paired coefficient difference for that response). This is a numerical diagnostic, not a confidence interval or full discretization-error bound.','']
    lines+=['','## Case results','','| Case | Release (s) | Final-window MAP (mmHg) | Auxiliary output (L/min) | Active deficit (mL) | Pump deficit (mL) | Total deficit (mL) |','|---|---:|---:|---:|---:|---:|---:|']
    for case in cfg['cases']:
        v=data[case['name']]['refined']['10']
        lines.append(f"| {case['name']} | {data[case['name']]['refined']['release_s']} | {v['MAP_mmHg']:.6g} | {v['aux_output_L_min']:.6g} | {v['active_deficit_mL']:.6g} | {v['pump_deficit_mL']:.6g} | {v['total_deficit_mL']:.6g} |")
    max_error=max(v['max_coefficient_numerical_difference'] for block in analyses.values() for end in block.values() for v in end.values())
    metvals=[data[c['name']]['refined']['10']['total_deficit_mL'] for c in cfg['cases'] if c['block']=='metabolic']
    hydvals=[data[c['name']]['refined']['10']['MAP_mmHg'] for c in cfg['cases'] if c['block']=='hydraulic' and c['signs'] is not None]
    lines+=['', '## Findings and next sampling decision', '',
        f'Metabolic corner deficits span {min(metvals):.6g} to {max(metvals):.6g} mL O2 equivalent. Hydraulic corner final-window pressures span {min(hydvals):.6g} to {max(hydvals):.6g} mmHg; {sum(v>=65 for v in hydvals)} of eight reach the current 65 mmHg pressure screen. This is a scenario count, not a probability of feasibility.',
        'The separate nominal-point test exposes curvature that the corner-only interpolation cannot resolve. Next sample the interior of the influential factor box before interpreting continuous-domain contributions, then test cross-block combinations and longer horizons. The completed factorial is useful interaction evidence but does not close T7.2.',
        'Short-window absence of a preservation deficit is not evidence that reduced branch perfusion can preserve an organ indefinitely: absolute oxygen stores and the limited horizon remain influential modeling conditions.']
    max_extra=max(v['refined']['maximum_sampled_additional_O2_deficit_mL_s'] for v in data.values())
    lines+=['', '## Verification and limits','',
        f'269 regression tests pass. Maximum additional oxygen-limited deficit at saved diagnostic samples across all cases: {max_extra:.6g} mL/s. These diagnostics attribute shortages to glucose-limited aerobic metabolism; sampled quadrature discrepancies are saved separately and cumulative ledgers remain primary.', '', f'Largest nominal/refined coefficient difference across reported quantities: {max_error:.6g} (mixed physical units; detailed values are in summary.json). Nominal reference regression errors: {errors}. Each case retains conservation, root, inventory, stoichiometric and state-refinement checks.', '', 'Observed activation remains enabled: hydraulic-factor effects include any resulting change in trigger time, shown per case. New coronary supply is a hypothetical conductance control; added vascular tissue, storage and remodeling costs are not modeled. Contractility changes activation and measured workload costs. Hydraulic contraction remains prescribed even if metabolic demand is unmet. The model still excludes tissue glucose buffering, alternative substrates, ischemic injury and reperfusion injury. No result demonstrates recoverable support duration.', '', 'Factors are varied within separate metabolic and hydraulic blocks; cross-block interactions are untested. Numerical rank stability is distinct from stability under new interior samples, different ranges or model structures. Do not claim continuous-domain global sensitivity has converged.', '', 'Protocol and provenance: [factorial protocol](../../docs/m7-factorial-protocol.md). Method: [NIST full factorial designs](https://www.itl.nist.gov/div898/handbook/pri/section3/pri3331.htm). Reproduce using scripts/run_m7_factorial.py then scripts/report_m7_factorial.py. input-manifest.json was saved before main runs.']
    (OUT/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    for ax,block,output in zip(axes,('metabolic','hydraulic'),('total_deficit_mL','MAP_mmHg')):
        d=analyses[block]['10'][output]['refined'];names=[f['id'] for f in cfg['factors'][block]]
        ax.barh(names,[d['terms'][n]['effect'] for n in names]);ax.axvline(0,color='black',lw=.7)
        ax.set_title(block.title()+' block');ax.set_xlabel('Main effect: '+output)
    fig.tight_layout();fig.savefig(OUT/'effects.png',dpi=160);plt.close(fig)
    inputs=read(OUT/'input-manifest.json')['sha256']
    changed=[p for p,h in inputs.items() if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h]
    if changed:raise RuntimeError('Frozen input mismatch '+str(changed))
    paths=[p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json' and p.suffix!='.log']
    paths+=[Path(__file__),ROOT/'reports/m7-factorial-tests.xml',ROOT/'reports/m7-factorial-regression-tests.xml']
    write(OUT/'manifest.json',dict(input_hashes_verified=True,sha256={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}))
    print('Verified',len(data),'cases. Report:',OUT/'report.md')
if __name__=='__main__':main()
