"""Reproducible T5.3 geometry, feedback and unit audit; no new design claim."""
from pathlib import Path
from dataclasses import asdict,replace
import sys,json,hashlib,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.primer_sizing import PrimerInputs,fixed_point,geometry,iterate,steady_power_W
OUT=ROOT/'reports/m5-primer-audit'
def main():
    p=PrimerInputs();thin=fixed_point(p);exact=fixed_point(p,True)
    same=geometry(thin['flow_mL_min'],p,True)
    iterations={mode:{str(seed):iterate(p,seed,1e-6,mode=='exact_geometry') for seed in [0,700,1500]} for mode in ['thin','exact_geometry']}
    old=json.loads((ROOT/'data/primer-arithmetic-audit.json').read_text())
    assert abs(thin['flow_mL_min']-old['flow_fixed_point_mL_min'])<1e-4
    assert abs(thin['two_shell_mass_g']-old['two_thin_shells_mass_g'])<1e-4
    power=[dict(left_head_mmHg=l,right_head_mmHg=r,left_W=steady_power_W(thin['flow_mL_min'],l),right_W=steady_power_W(thin['flow_mL_min'],r),total_W=steady_power_W(thin['flow_mL_min'],l+r)) for l in [50,70,90] for r in [5,10,20]]
    sensitivity={str(f):fixed_point(replace(p,self_flow_mL_min_g=f/300)) for f in [200,250]}
    out=dict(classification='Conditional historical arithmetic, not a physical heart specification',inputs=asdict(p),thin=thin,exact_same_dimensions=same,exact_geometry_reclosed=exact,iterations=iterations,self_flow_normalization_sensitivity=sensitivity,illustrative_two_head_power=power,one_path_head70_W=steady_power_W(thin['flow_mL_min'],70))
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    rows=['# T5.3 primer arithmetic verification','','Conditional on historical inputs. Exact-shell variants correct geometry only, not wall mechanics or anatomical completeness.','','| Calculation | Flow mL/min | Radius cm | Wall mm | Two-shell mass g | Feedback gain |','|---|---:|---:|---:|---:|---:|']
    for label,r in [('Thin-shell fixed point',thin),('Exact shell, same dimensions',same),('Exact geometry, reclosed feedback',exact)]:
        rows.append('| %s | %.6f | %.6f | %.6f | %.6f | %.8f |'%(label,r['flow_mL_min'],r['radius_cm'],10*r['wall_cm'],r['two_shell_mass_g'],r['feedback_gain']))
    rows+=['','The same-dimension row is not a fixed point for its changed self-flow coefficient. Reclosing that loop is an arithmetic sensitivity, not a validated thick-wall solution.','','## Convergence','','| Geometry | Seed mL/min | Updates to <=1e-6 mL/min error | Final error bound mL/min |','|---|---:|---:|---:|']
    for mode,group in iterations.items():
        for seed,r in group.items():rows.append('| %s | %s | %d | %.3g |'%(mode,seed,r['updates'],r['error_bound_mL_min']))
    rows+=['','## Power units','',f"At {thin['flow_mL_min']:.6f} mL/min and an assumed single-path 70 mmHg pressure rise: {out['one_path_head70_W']:.8f} W.",'','| Assumed left head mmHg | Assumed right head mmHg | Total steady hydraulic W |','|---:|---:|---:|']
    for r in power:rows.append('| %g | %g | %.8f |'%(r['left_head_mmHg'],r['right_head_mmHg'],r['total_W']))
    rows+=['','The head grid is an arithmetic illustration, not a physiological range. Equal path flows are assumed here only. Actual M5 circuits retain separate flows and port pressures. Neither table supplies muscle input energy or a support duration.','','[Full audit, assumptions and closure](../../docs/m5-t53-closure.md).','','Reproduce: `./.venv/Scripts/python.exe scripts/audit_m5_primer.py`.']
    (OUT/'report.md').write_text('\n'.join(rows)+'\n',encoding='utf-8')
    suites=list(ET.parse(ROOT/'reports/m5-t53-tests.xml').getroot().iter('testsuite'))
    assert sum(int(x.get('tests',0)) for x in suites)==217
    assert all(int(x.get('failures',0))+int(x.get('errors',0))==0 for x in suites)
    files=[ROOT/'reports/m5-t53-tests.xml',ROOT/'docs/m5-t53-closure.md',ROOT/'docs/m5-supported-organ-set.md',Path(__file__),ROOT/'src/cardiac_redundancy/primer_sizing.py',ROOT/'tests/test_primer_sizing.py',ROOT/'docs/archive/primer-pre-m1-2026-10-01.md',ROOT/'data/primer-arithmetic-audit.json',OUT/'results.json',OUT/'report.md']
    (OUT/'manifest.json').write_text(json.dumps(dict(status='arithmetic_verified',hashes={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}),indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
