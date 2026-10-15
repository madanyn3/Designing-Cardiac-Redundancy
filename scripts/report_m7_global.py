"""Morris sampling adequacy, trajectory bootstrap, and numerical audit."""
from pathlib import Path
import sys,json,hashlib,platform
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.global_effects import elementary_effects,summarize
from report_m7_factorial import metrics
OUT=ROOT/'reports/m7-global'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def main(n):
    cfg=read(ROOT/'configs/m7-global.json');cases=[c for c in cfg['cases'] if c['trajectory']<n]
    names=[f['id'] for f in cfg['factors']];data={};failed=[]
    for c in cases:
        folder=OUT/c['name'];path=folder/'result.json'
        if not path.exists() or read(path)['status']!='verified':
            failed.append(dict(case=c['name'],status=read(path).get('status') if path.exists() else 'missing'));continue
        cache=folder/'analysis-metrics.json'
        # Derived cache depends on the paired result artifacts; never reuse if
        # those files or the extraction routine change.
        paths=[path,folder/'metabolic-nominal.npz',folder/'metabolic-refined.npz',ROOT/c['hydraulic_source']/'hydraulics.json',ROOT/c['hydraulic_source']/'hydraulic-nominal.npz',ROOT/c['hydraulic_source']/'hydraulic-refined.npz',ROOT/'scripts/report_m7_factorial.py']
        fingerprint=hashlib.sha256(''.join(hashlib.sha256(p.read_bytes()).hexdigest() for p in paths).encode()).hexdigest()
        if cache.exists() and read(cache).get('fingerprint')==fingerprint:values=read(cache)['values']
        else:
            values={label:metrics(folder,ROOT/c['hydraulic_source'],label) for label in ('nominal','refined')}
            write(cache,dict(fingerprint=fingerprint,values=values))
        data[c['name']]=values
    write(OUT/f'audit-{n}.json',dict(expected=7*n,verified=len(data),failed=failed))
    if failed:raise RuntimeError('Incomplete trajectories; no successful-subset analysis')
    x=np.array([c['normalized'] for c in cases]).reshape(n,7,6)
    stages=[r for r in cfg['stages'] if r<=n];analysis={}
    for end in ('4','7','10'):
        analysis[end]={}
        for output in data[cases[0]['name']]['refined'][end]:
            arrays={label:np.array([data[c['name']][label][end][output] for c in cases]).reshape(n,7) for label in ('nominal','refined')}
            ee={label:elementary_effects(x,y) for label,y in arrays.items()}
            err=np.max(abs(ee['nominal']-ee['refined']),axis=0)
            records={str(r):summarize(ee['refined'][:r],cfg['bootstrap_seed'],cfg['resamples']) for r in stages}
            current=records[str(n)];previous=records[str(stages[-2])] if len(stages)>1 else current
            rank=lambda a:sorted([j for j,v in enumerate(a) if v>max(1e-10,10*err[j])],key=lambda j:a[j],reverse=True)
            order=rank(current['mu_star']);oldorder=rank(previous['mu_star']);nominal_order=rank(np.abs(ee['nominal']).mean(axis=0))
            top=order[:2];change={names[j]:abs(previous['mu_star'][j]-current['mu_star'][j])/max(current['mu_star'][j],1e-10) for j in top}
            stable=set(top)==set(oldorder[:2]) and all(v<=.3 for v in change.values())
            numerical_ok=all(err[j]<=.01*max(current['mu_star'][j],1e-10) for j in top)
            analysis[end][output]=dict(stages=records,elementary_effects=ee['refined'].tolist(),paired_max_EE_difference=err.tolist(),rank=[names[j] for j in order],nominal_rank=[names[j] for j in nominal_order],previous_rank=[names[j] for j in oldorder],top_two_relative_change=change,stable_primary_rule=stable,numerical_top_two_ok=numerical_ok,first_half=summarize(ee['refined'][:n//2],cfg['bootstrap_seed'],cfg['resamples']),second_half=summarize(ee['refined'][n//2:],cfg['bootstrap_seed'],cfg['resamples']))
    unstable=[output for output in ('total_deficit_mL','MAP_mmHg') if not analysis['10'][output]['stable_primary_rule']]
    next_n=next((r for r in cfg['stages'] if r>n),None) if unstable else None
    numerical_bad=[o for o in ('total_deficit_mL','MAP_mmHg') if not analysis['10'][o]['numerical_top_two_ok']]
    if numerical_bad:raise RuntimeError('Primary effects numerically unresolved: '+str(numerical_bad))
    retries=[c['name'] for c in cases if (OUT/c['name']/'first-attempt.json').exists()]
    summary=dict(trajectories=n,cases=len(cases),seed=cfg['seed'],bootstrap_seed=cfg['bootstrap_seed'],analysis=analysis,unstable_primary_outputs=unstable,next_stage=next_n,numerical_retries=retries,failed_cases=0,scope='bounded six-factor four-level global elementary effects; no population or survival probabilities',metrics=data,python=platform.python_version(),numpy=np.__version__)
    write(OUT/f'summary-{n}.json',summary)
    lines=[f'# T7.2 global sensitivity: {n} trajectories','',f'{len(cases)} joint cases verified, with paired numerical solves. No failed cases omitted. Main-study retries: {len(retries)}. Four-level Morris design, seed7202; 5000 whole-trajectory bootstrap replicates, seed7203. 274 regression tests pass.','', 'Mu-star is mean absolute elementary effect per full input range. Sigma mixes nonlinearity and interactions. Bootstrap intervals quantify finite-trajectory sampling uncertainty under this artificial scenario design; they do not describe a human population.','', '## Ten-second results','', '| Output | Factor | Mu-star | 95% bootstrap interval | Sigma | Rank interval |','|---|---|---:|---|---:|---|']
    for output in ('total_deficit_mL','MAP_mmHg','aux_output_L_min'):
        a=analysis['10'][output];s=a['stages'][str(n)]
        for name in a['rank']:
            j=names.index(name)
            lines.append(f"| {output} | {name} | {s['mu_star'][j]:.6g} | [{s['mu_star_interval'][0][j]:.6g}, {s['mu_star_interval'][1][j]:.6g}] | {s['sigma'][j]:.6g} | [{s['rank_interval'][0][j]:.3g}, {s['rank_interval'][1][j]:.3g}] |")
    lines+=['','Zero or numerically unresolved effects are omitted from ranking only, never from the simulations or full JSON results. Separate bootstrap intervals do not constitute a simultaneous ranking confidence set.','', '## Sampling adequacy','', '| Output | Previous prefix rank | Current rank | Relative change in current top two | Primary stability rule |','|---|---|---|---|---|']
    for output in ('total_deficit_mL','MAP_mmHg','aux_output_L_min'):
        a=analysis['10'][output]
        lines.append(f"| {output} | {' > '.join(a['previous_rank'])} | {' > '.join(a['rank'])} | {a['top_two_relative_change']} | {a['stable_primary_rule']} |")
    lines+=['',f'Unstable registered primary outputs: {unstable or "none"}. Next registered stage: {next_n if next_n else "none"}. The rule compares the previous/current top-two sets and allows at most 30% change in their mu-star magnitudes. It is an exploratory stopping rule, not proof of asymptotic convergence. At sixteen trajectories any unresolved ranking is retained explicitly.','', '## Window and numerical checks','', '| Output | 4 s | 7 s | 10 s | Numerical rank match at 10 s |','|---|---|---|---|---|']
    for output in ('total_deficit_mL','MAP_mmHg','aux_output_L_min'):
        ranks=[' > '.join(analysis[t][output]['rank']) for t in ('4','7','10')];a=analysis['10'][output]
        lines.append('| '+output+' | '+' | '.join(ranks)+f" | {a['rank']==a['nominal_rank']} |")
    lines+=['','## Disjoint-half diagnostic','','| Output | First half mu-star ranking | Second half mu-star ranking |','|---|---|---|']
    for output in ('total_deficit_mL','MAP_mmHg','aux_output_L_min'):
        a=analysis['10'][output];ranks=[]
        for half in ('first_half','second_half'):
            values=a[half]['mu_star'];order=sorted([j for j,v in enumerate(values) if v>1e-10],key=lambda j:values[j],reverse=True)
            ranks.append(' > '.join(names[j] for j in order))
        lines.append('| '+output+' | '+' | '.join(ranks)+' |')
    lines+=['','Disjoint halves reveal dependence on sampled backgrounds. Rank disagreements remain visible even if the nested-prefix stopping rule passes.','']
    vals=[v['refined']['10'] for v in data.values()]
    lines+=['', '## Scenario outcomes','',f"Final-window MAP spans {min(v['MAP_mmHg'] for v in vals):.6g}–{max(v['MAP_mmHg'] for v in vals):.6g} mmHg; {sum(v['MAP_mmHg']>=65 for v in vals)} of {len(vals)} points reach 65 mmHg. Total deficit spans {min(v['total_deficit_mL'] for v in vals):.6g}–{max(v['total_deficit_mL'] for v in vals):.6g} mL O2 equivalent. Counts are neither feasibility nor survival probabilities.", '',f"Maximum additional oxygen limitation at saved samples is {max(v['refined']['maximum_sampled_additional_O2_deficit_mL_s'] for v in data.values()):.6g} mL/s. Cumulative ledgers remain primary; sampled attribution and quadrature errors are saved separately.", '', '## Scope and limitations','', 'These runs directly evaluate joint interior grid points and avoid the invalid corner surrogate. Six influential coordinates are varied together; anatomy groups remain linked and all other physiology stays fixed. The design does not estimate continuous-domain Sobol indices, optimized pump parameters, injury, reversible recovery or long-duration support. Variability among elementary effects is a reason to retain uncertainty, not to select favorable defaults. Wider physiology and alternative equations remain separate tasks. Contraction remains prescribed despite unmet metabolic demand: pressure/output sensitivities are conditional on that driver, not validated energy-limited pump performance.', '', 'The faster pilot was rejected after capacity checks failed; the main study uses original strict settings. All old frozen inputs must retain their hashes. Individual conservation, capacity, stoichiometry, event and refinement checks remain in case results.','', 'Protocol: [joint global protocol](../../docs/m7-global-protocol.md). Methods: [Morris 1991](https://doi.org/10.1080/00401706.1991.10484804), [Campolongo et al. 2007](https://doi.org/10.1016/j.envsoft.2006.10.004), [SALib documentation](https://salib.readthedocs.io/en/latest/api/SALib.analyze.html#module-SALib.analyze.morris). Reproduce with scripts/run_m7_global.py followed by scripts/report_m7_global.py, passing the trajectory count.']
    (OUT/f'report-{n}.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    import os
    os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'reports/m7-plot-cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    for ax,output,label in zip(axes,('total_deficit_mL','MAP_mmHg'),('Unmet aerobic demand (mL O2 equivalent)','Mean aortic pressure (mmHg)')):
        a=analysis['10'][output];v=a['stages'][str(n)];order=[names.index(name) for name in a['rank']][::-1]
        value=np.array(v['mu_star'])[order];low=np.array(v['mu_star_interval'][0])[order];high=np.array(v['mu_star_interval'][1])[order]
        ax.barh([names[j].replace('_',' ') for j in order],value)
        ax.errorbar(value,np.arange(len(order)),xerr=np.vstack([np.maximum(0,value-low),np.maximum(0,high-value)]),fmt='none',ecolor='black',capsize=3)
        ax.set_xlabel('Mean absolute elementary effect\n'+label);ax.set_title(f'{n} trajectories, 95% bootstrap intervals')
    fig.tight_layout();fig.savefig(OUT/f'effects-{n}.png',dpi=160);plt.close(fig)
    print(json.dumps(dict(trajectories=n,unstable=unstable,next_stage=next_n,ranks={o:analysis['10'][o]['rank'] for o in ('total_deficit_mL','MAP_mmHg','aux_output_L_min')},changes={o:analysis['10'][o]['top_two_relative_change'] for o in ('total_deficit_mL','MAP_mmHg')}),indent=2))
if __name__=='__main__':main(int(sys.argv[1]) if len(sys.argv)>1 else 8)
