"""Arm B coupled replay: fat supply set by the depot's actual perfusion.

Arm A charges every joule of oxidative demand to glucose. Arm B measures how
much fatty acid the visceral depot can actually deliver at the perfusion the
support policy grants it, converts that into per-organ glucose shares, and
replays the metabolism with those shares.

The chain is deliberately one-way and auditable:

    support policy -> depot conductance factor -> measured adipose flow
                   -> fat supply (lipolysis or perfusion, whichever binds)
                   -> per-organ glucose shares -> coupled replay

Nothing is assumed available. The fat supply comes from the saved trajectory's
own adipose flow, so a design that does not perfuse the depot does not get the
fuel. The brain stays wholly on glucose throughout.

This reuses the verified M7 replay unchanged apart from its declared
`glucose_share` input, so Arm A and Arm B differ by the share array and the
network variant, and by nothing else.
"""
from pathlib import Path
import json, sys, time
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import run_m7_global_metabolism as met
import run_m8_metabolism as m8
from cardiac_redundancy.adipose_policy import ADIPOSE_BRANCH, AdiposePolicy
from cardiac_redundancy.design import DesignPoint, ScenarioPoint, mass_ratio
from cardiac_redundancy.fat import LipolysisParameters, fat_supply, glucose_shares
from cardiac_redundancy.fuel import FuelParameters

OUT = 'reports/m8b-coupled'
#: The branch the depot was disaggregated out of; its baseline oxygen store
#: splits between the two in proportion to the demand each now carries.
PARENT_ORGAN = 'remaining'
NETWORK = 'configs/m8-adipose-regional.yaml'
MAPPING = 'configs/m8-adipose-organ-mapping.json'


def read(path): return json.loads((ROOT / path).read_text(encoding='utf-8'))


def measured_adipose_flow(hydraulic_dir, flow_names):
    """Mean depot perfusion over the settled part of the saved trajectory.

    Taken from the trajectory rather than from the branch's reference flow,
    because the reference is defined at the construction pressure of 100 mmHg
    while the supported state runs well below it.
    """
    a = np.load(ROOT / hydraulic_dir / 'hydraulic-refined.npz')
    t = a['time_s']; q = a['flows_mL_s']
    if ADIPOSE_BRANCH not in flow_names:
        raise ValueError('This trajectory has no %s branch; it is an Arm A run' % ADIPOSE_BRANCH)
    column = q[:, list(flow_names).index(ADIPOSE_BRANCH)]
    settled = t >= max(float(t[0]), float(t[-1])-30.)
    window = t[settled]
    if len(window) < 2:
        raise ValueError('Trajectory too short to average depot perfusion')
    return float(np.trapezoid(column[settled], window)/(window[-1]-window[0])*60.)


def evaluate(design, scenario, policy, name, hydraulic_dir, out=OUT,
             lipolysis=LipolysisParameters()):
    hydraulics = read(Path(hydraulic_dir) / 'hydraulics.json')
    if hydraulics['runs'][-1]['termination'] is not None:
        return dict(status='hydraulic_domain_limit', name=name)
    flow_names = hydraulics.get('flow_names') or hydraulics['runs'][-1].get('flow_names')
    if flow_names is None:
        raise ValueError('Saved hydraulics do not record flow names')

    adipose_flow = measured_adipose_flow(hydraulic_dir, flow_names)
    supply = fat_supply(adipose_flow, lipolysis)

    # Organ demands under this design, needed to apportion the fat supply.
    met.OUT = ROOT / out
    assumptions, sides = m8.auxiliary_demand(design, scenario, hydraulics)
    case = dict(name=name, device='BIO', architecture='BIO', sides=sides,
                network=NETWORK,
                liver_function_fraction=scenario.liver_function_fraction,
                hydraulic_source=str(hydraulic_dir),
                knobs=dict(extraction_limit=scenario.extraction_limit,
                           preservation_fraction=scenario.preservation_fraction))
    pars = FuelParameters(**{**read('configs/m4-fuel-reference.json')['parameters'],
        'glycogen_mmol_L_liver': scenario.glycogen_mmol_L_liver,
        'release_cap_mmol_min': scenario.release_cap_mmol_min})
    original = met.read
    store = [sides[s]['total_mL_min']*m8.STORE_SECONDS/60 for s in ('left', 'right')]

    def coupled_read(path):
        data = original(path)
        if str(path).endswith('m5-myocardial-demand.json'):
            data = dict(data, auxiliary_store_capacity_mL=store)
        if str(path).endswith('m4-organ-mapping.json'):
            data = read(MAPPING)
        if str(path).endswith('m5-organ-preservation.json'):
            # The depot needs a profile entry, and it must agree with the
            # conductance the same policy granted it hydraulically.
            data = json.loads(json.dumps(data))
            profile = data['cases'][1]['profile']
            profile[ADIPOSE_BRANCH] = policy.organ_entry(scenario.preservation_fraction)
        return data

    # The organ list must be read under the same patched configuration the
    # replay will use, or the share array is sized for the other arm.
    met.read = coupled_read
    try:
        organs = met.network(case)[1]
    finally:
        met.read = original
    names = [o.name for o in organs]
    demands = [o.demand_mL_min for o in organs]
    shares = glucose_shares(names, demands, supply)
    if len(shares['shares']) != len(organs):
        raise ValueError('Share array does not match the organ list')
    if ADIPOSE_BRANCH not in names:
        raise ValueError('Arm B replay is missing the %s organ' % ADIPOSE_BRANCH)
    case['glucose_share'] = shares['shares']
    # The depot was disaggregated out of the residual branch, so the baseline's
    # stored oxygen for that branch must split with the demand that moved.
    demand_by_name = dict(zip(names, demands))
    parent_total = demand_by_name[PARENT_ORGAN]+demand_by_name[ADIPOSE_BRANCH]
    case['store_split'] = {
        PARENT_ORGAN: (PARENT_ORGAN, demand_by_name[PARENT_ORGAN]/parent_total),
        ADIPOSE_BRANCH: (PARENT_ORGAN, demand_by_name[ADIPOSE_BRANCH]/parent_total)}

    start = time.perf_counter()
    met.read = coupled_read
    try:
        result = met.short(case, pars)
    finally:
        met.read = original
    record = result['records'][-1]
    summary = dict(status='verified', arm='B', name=name, design=vars(design),
        scenario=vars(scenario), policy=policy.describe(),
        runtime_s=time.perf_counter()-start,
        adipose_flow_mL_min=adipose_flow, fat_supply=supply, glucose_shares=shares,
        coupled_mass_ratio=mass_ratio(design, scenario), rate_bpm=design.rate_bpm,
        auxiliary_demand_mL_min={s: v['total_mL_min'] for s, v in sides.items()},
        groups=m8.group_ledgers(result, Path(out) / name),
        unmet_oxygen_mL=record['unmet_oxygen_mL'],
        glycogen_remaining_mmol=record['glycogen_remaining_mmol'],
        glucose_remaining_mmol=record['glucose_remaining_mmol'],
        released_mmol=record['released_mmol'],
        minimum_plasma_equivalent_mmol_L=record['minimum_plasma_equivalent_mmol_L'],
        scope='fat supply measured from the depot perfusion this design actually '
              'delivers; the brain stays wholly on glucose')
    dest = ROOT / out / name; dest.mkdir(parents=True, exist_ok=True)
    (dest / 'm8b-coupled.json').write_text(
        json.dumps(summary, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return summary


if __name__ == '__main__':
    source = Path(sys.argv[1])
    perfused = len(sys.argv) > 2 and sys.argv[2].lower() in ('1', 'true', 'perfused')
    policy = AdiposePolicy(perfused=perfused)
    design = DesignPoint(**read(source / 'hydraulics.json')['design'])
    name = source.name+('-perfused' if perfused else '-preserved')
    r = evaluate(design, ScenarioPoint(), policy, name, source)
    if r['status'] != 'verified':
        print(name, r['status'], flush=True); sys.exit(0)
    s = r['fat_supply']; g = r['glucose_shares']
    print('%s: depot flow %.1f mL/min (%s binds), fat %.1f mL O2-eq/min, covers %.1f%%' % (
        name, r['adipose_flow_mL_min'], s['binding'], s['supply_mL_O2_eq_min'],
        100*g['covered_fraction']), flush=True)
    for role, grp in sorted(r['groups'].items()):
        print('   %-12s requested %8.3f mL  unmet sustained %.6f  %s' % (
            role, grp['requested_mL'], grp['unmet_sustained_mL'], grp['sustained_status']), flush=True)
