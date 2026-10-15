"""Build a versioned regional amendment without modifying the M2 freeze.

Published parent fractions and conditional tissue estimates are kept distinct.
Unresolved limb bone/other allocations remain explicitly provisional; all
non-limb complements are subtracted from named whole-body budgets.
"""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=json.loads((ROOT/'configs/project-extension.yaml').read_text())
    evidence=json.loads((ROOT/'data/m3-regional-flow-candidates.json').read_text())
    estimates={r['branch']:r for r in evidence['candidates']}
    ref='https://doi.org/10.1097/00004032-199508000-00003'
    fractions={'brain':.12,'heart':.04,'kidneys':.19,'liver_arterial':.065,
               'stomach_oesophagus':.01,'small_intestine':.10,'large_intestine':.04,
               'pancreas':.01,'spleen':.03,'lungs':.025}
    old={r['name']:r for r in p['branches']}
    targets={n:5600*f for n,f in fractions.items()}
    targets.update({n:r['flow_ml_min'] for n,r in estimates.items()})
    # No matching skeletal volume or isolated other-tissue rate has yet been
    # admitted. Retain these historical construction inputs transparently.
    for n in ('limb_bone','limb_other'):
        targets[n]=old[n]['target_flow_mL_s']*60
    complements={'nonlimb_nonrespiratory_muscle':5600*.17-targets['limb_muscle']-targets['respiratory_muscle'],
                 'nonlimb_skin':5600*.05-targets['limb_skin'],
                 'nonlimb_adipose':5600*.05-targets['limb_adipose'],
                 'nonlimb_skeleton':5600*.05-targets['limb_bone'],
                 'other_source_tissues_less_limb_other':5600*.05-targets['limb_other']}
    assert all(v>0 for v in complements.values())
    targets['remaining']=sum(complements.values())
    assert abs(sum(targets.values())-5600)<1e-8
    portal=['stomach_oesophagus','small_intestine','large_intestine','pancreas','spleen']
    branches=[]
    for n,q in targets.items():
        recipient='H' if n in portal+['liver_arterial'] else ('PU' if n=='lungs' else ('B_'+n if n.startswith('limb_') else 'VC'))
        downstream=12 if recipient=='H' else 4 if recipient=='PU' else 5 if recipient.startswith('B_') else 2
        classification='published reference fraction transferred to project CO'
        source=ref
        if n in estimates:
            classification=estimates[n]['status']; source=estimates[n]['source']
        elif n in ('limb_bone','limb_other'):
            classification='UNRESOLVED historical construction allocation; not evidence-qualified'
            source=None
        elif n=='remaining': classification='derived named parent-budget complements; inherits unresolved partitions'
        branches.append(dict(name=n,**{'from':'AO','to':recipient},baseline_fraction=q/5600,
            target_flow_mL_s=q/60,R_mmHg_s_mL=(100-downstream)/(q/60),
            **{'class':classification},source=source))
    hepatic=sum(targets[n] for n in portal+['liver_arterial'])
    branches.append(dict(name='hepatic_outflow',**{'from':'H','to':'VC'},baseline_fraction=None,
        target_flow_mL_s=hepatic/60,R_mmHg_s_mL=10/(hepatic/60),**{'class':'derived serial return'}))
    for n in targets:
        if n.startswith('limb_'):
            branches.append(dict(name=n+'_return',**{'from':'B_'+n,'to':'VC'},baseline_fraction=None,
                target_flow_mL_s=targets[n]/60,R_mmHg_s_mL=3/(targets[n]/60),
                **{'class':'derived using unchanged assumed 5 mmHg pool pressure'}))
    p.update(schema='M3 evidence amendment v1',status='referenced estimates with explicit unresolved partitions; not physiological validation',
        branches=branches,portal_branches=portal,remaining_flow_components_mL_min=complements,
        evidence_notes=['docs/m3-splanchnic-split.md','docs/m3-regional-flow-candidates.md'],
        unresolved=['limb_bone and limb_other flows','M1 sub-organ demand reconciliation',
                    'pressure/storage anchors','kidney/coronary sub-organ split','anatomical transfers'])
    (ROOT/'configs/m3-evidence-regional.yaml').write_text(json.dumps(p,indent=2)+'\n')
    protocol=json.loads((ROOT/'configs/m3-regional-protocol.json').read_text())
    protocol.update(version='2026-10-05-evidence-regional-v1',project_config='configs/m3-evidence-regional.yaml',
        variant='Referenced regional estimates and five portal beds; fixed-R comparator; unresolved partitions explicit')
    (ROOT/'configs/m3-evidence-regional-protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    print(json.dumps({'arterial_total_mL_min':sum(targets.values()),'hepatic_mL_min':hepatic,'remaining_components':complements},indent=2))


if __name__=='__main__': main()
