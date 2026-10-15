"""Consolidate qualified T7.2 closure only after the registered sampling decision."""
from pathlib import Path
import json,sys,hashlib,shutil
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'reports/m7-global'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def main(n):
    summary=read(OUT/f'summary-{n}.json');cfg=read(ROOT/'configs/m7-global.json')
    if summary['next_stage'] is not None:raise RuntimeError('Registered expansion remains required')
    if summary['failed_cases']:raise RuntimeError('Unresolved main-study failures')
    frozen=read(OUT/'input-manifest.json')['sha256']
    changed=[p for p,h in frozen.items() if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h]
    if changed:raise RuntimeError('Frozen inputs changed: '+str(changed))
    cases=[c for c in cfg['cases'] if c['trajectory']<n]
    for c in cases:
        if read(OUT/c['name']/'result.json')['status']!='verified':raise RuntimeError('Missing verified case')
    suites=ET.parse(ROOT/'reports/m7-global-regression-tests.xml').getroot().findall('testsuite')
    tests=sum(int(s.attrib['tests']) for s in suites)
    if any(int(s.attrib.get('errors',0))+int(s.attrib.get('failures',0)) for s in suites):raise RuntimeError('Regression failed')
    names=[f['id'] for f in cfg['factors']];rows=[]
    for output in ('total_deficit_mL','MAP_mmHg','aux_output_L_min'):
        a=summary['analysis']['10'][output];stats=a['stages'][str(n)];j=names.index(a['rank'][0])
        rows.append(f"| {output} | {a['rank'][0]} | {stats['mu_star'][j]:.6g} | [{stats['mu_star_interval'][0][j]:.6g}, {stats['mu_star_interval'][1][j]:.6g}] |")
    paths={c['hydraulic_source'] for c in cases};new=[p for p in paths if p.startswith('reports/m7-global/')]
    vals=[v['refined']['10'] for v in summary['metrics'].values()]
    unstable=summary['unstable_primary_outputs']
    status='Registered leading-effect stability rule passed.' if not unstable else 'Maximum registered sample reached; unresolved primary ordering/magnitude retained: '+', '.join(unstable)+'. No convergence claim for these outputs.'
    text=f"""# T7.2 closure: bounded global sensitivity characterization

2026-10-15. **Closed at the stated exploratory scope.** M7 remains open for T7.3 identifiability, T7.4 structural alternatives and T7.5 evidence priorities. Preprint writing remains parked; optimization belongs to M8.

## Evidence completed

- Initial thirteen-case metabolic endpoint screen, followed by seventeen factorial cases with an independent nominal regression.
- Joint Morris study: **{n} trajectories, {7*n} cases**, all verified with paired metabolic solves. {len(paths)} distinct hydraulic configurations were used: {len(new)} newly computed with paired solves, {len(paths)-len(new)} reused from the frozen factorial study.
- **{tests} regression tests pass.** Main-study numerical retries: {len(summary['numerical_retries'])}; no failed main-study points omitted. Four failed accelerated pilot cases are explicitly retained in the pilot audit; both accelerated settings were rejected.
- Direct interior and cross-block evaluations replace the inaccurate corner surrogate. Sampling seed7202; trajectory bootstrap seed7203, 5000 replicates. Compare nested sample counts, disjoint halves, observation windows and numerical settings.

{status}

## Leading effects at ten seconds

| Output | Largest mean absolute elementary effect | Mu-star | 95% trajectory-bootstrap interval |
|---|---|---:|---|
"""+'\n'.join(rows)+f"""

Effects have output units per full chosen parameter range; they are not variance fractions. Intervals quantify sampling variability under this grid design, not patient variability. The detailed report retains overlapping secondary rank intervals, disjoint-half disagreements and horizon-dependent ordering. Passing the leading-effect stopping rule does not establish a precise order for every factor.

Pressure spans **{min(v['MAP_mmHg'] for v in vals):.4g}–{max(v['MAP_mmHg'] for v in vals):.4g} mmHg**; **{sum(v['MAP_mmHg']>=65 for v in vals)} of {len(vals)}** sampled cases reaches the 65 mmHg final-window pressure screen. Deficit spans **{min(v['total_deficit_mL'] for v in vals):.6g}–{max(v['total_deficit_mL'] for v in vals):.6g} mL O2 equivalent**. These results exclude neither designs outside the box nor different model structures. No feasible auxiliary heart or recovery window is established.

## Closure criteria and boundaries

| Requirement | Disposition |
|---|---|
| Screening followed by global analysis of influential inputs | Complete: six-factor joint elementary-effect study over a four-level grid |
| Interior nonlinearity and cross-block interactions | Directly sampled; sigma reflects nonlinearity plus interactions, not a pure interaction index |
| Sample size, seeds, failures and numerical verification | Recorded in configurations, per-case results, stage audits and manifests |
| Sampling adequacy and ranking uncertainty | Explicit nested-prefix rule plus bootstrap, disjoint-half, time-window and numerical checks; unresolved secondary ordering retained |
| Correlations and dependencies | Anatomical groups stay linked; independent scenario controls are not asserted to be independent human physiology |
| Wider physiological uncertainty | Not globally quantified here; registry families held fixed remain limitations, not proven unimportant |

This closure concerns a **bounded six-parameter, ten-second global screening study**. It does not claim continuous-domain Sobol convergence, calibrated distributions, robust long-duration preservation, or a universal ranking across model forms. Glucose buffering, substrate mixture and oxygen-store assumptions remain consequential structural limits. Contraction remains prescribed despite unmet demand. The short-window stock/release results cannot establish long-term glycogen irrelevance.

## Next tasks and artifacts

T7.3 should test what baseline observations can identify and expose compensating parameter combinations. T7.4 should test the instantaneous-glucose closure and other declared structural alternatives before M8 treats a favorable modeled point as credible. T7.5 can prioritize measurements using sensitivity plus these structural/identifiability limits.

- [Global report](../reports/m7-global/report.md), [stage {n} report](../reports/m7-global/report-{n}.md), [quantitative results](../reports/m7-global/summary-{n}.json), [figure](../reports/m7-global/effects-{n}.png).
- [Frozen global protocol](m7-global-protocol.md), [configuration](../configs/m7-global.json), [input manifest](../reports/m7-global/input-manifest.json), [pilot audit](../reports/m7-global/pilot-audit.json).
- [Initial screening](../reports/m7-metabolic-screen/report.md), [factorial study](../reports/m7-factorial/report.md), [uncertainty registry](uncertainty-plan.md).
"""
    (ROOT/'docs/m7-t72-closure.md').write_text(text,encoding='utf-8')
    shutil.copyfile(OUT/f'report-{n}.md',OUT/'report.md')
    (OUT/'closure.json').write_text(json.dumps(dict(task='T7.2',status='closed_at_bounded_global_screening_scope',trajectories=n,cases=len(cases),tests=tests,unstable_primary_outputs=unstable,frozen_inputs_verified=True),indent=2)+'\n',encoding='utf-8')
    artifacts=[p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    artifacts += [ROOT/'docs/m7-t72-closure.md',ROOT/'scripts/report_m7_global.py',ROOT/'scripts/finalize_m7_t72.py',ROOT/'reports/m7-global-tests.xml',ROOT/'reports/m7-global-regression-tests.xml',ROOT/'configs/m7-global-pilot.json',ROOT/'configs/m7-global-pilot-strict.json']
    (OUT/'manifest.json').write_text(json.dumps(dict(sha256={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in artifacts}),indent=2)+'\n',encoding='utf-8')
    print('T7.2 qualified closure:',n,'trajectories;',len(cases),'cases;',tests,'tests')
if __name__=='__main__':main(int(sys.argv[1]))
