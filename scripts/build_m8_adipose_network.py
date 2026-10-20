"""Arm B step 1: disaggregate visceral adipose out of the residual branch.

Arm A charges all oxidative demand to glucose. Arm B asks what changes when
tissue can oxidise fat, which requires adipose tissue to be perfused well
enough to export fatty acid. That makes adipose conductance a design variable,
so the circuit has to represent it explicitly.

The existing network already accounts for the whole cardiac output: its 23
branches have baseline fractions summing to exactly 1 and a resistance
construction R = (P_AO - P_VC)/Q throughout. Visceral adipose is therefore
already inside the budget, carried by the residual `remaining` branch, and the
existing `limb_adipose` branch is classed in the frozen config as an
"abdominal-to-limb transfer candidate, not limb measurement".

So this **disaggregates** rather than adds: `remaining` is split into
`remaining` and `visceral_adipose`, preserving total flow, the fraction sum and
the total oxygen demand. Adding a branch instead would inflate cardiac output
and demand beyond the closed M1/M3 accounting.

The frozen M3 configuration is not modified. This writes an Arm B variant, so
Arm A stays reproducible and the two arms differ by one declared change.

Declared inputs, each an assumption rather than a measurement of this subject:

  * visceral adipose mass, 2.0 kg for the reference young adult male;
  * adipose blood flow at 4 mL per 100 g per minute, the midpoint of the
    3-5 range reported for abdominal adipose after an overnight fast;
  * oxygen demand taken from the project's **own** limb_adipose demand-to-flow
    ratio rather than an external number, so the split stays internally
    consistent with the M1 budget.

Adipose perfusion is expected to bind before lipolytic capacity does: carrying
the fatty acid that whole-body lipolysis can release, at the flow this branch
provides, needs a venous concentration at or above the plausible plasma
ceiling. Establishing that is the next step, not this one.
"""
from pathlib import Path
import json, sys
ROOT = Path(__file__).resolve().parents[1]

SOURCE_NETWORK = 'configs/m3-evidence-regional.yaml'
SOURCE_MAPPING = 'configs/m4-organ-mapping.json'
OUT_NETWORK = 'configs/m8-adipose-regional.yaml'
OUT_MAPPING = 'configs/m8-adipose-organ-mapping.json'

VISCERAL_ADIPOSE_MASS_KG = 2.0
ADIPOSE_FLOW_ML_PER_100G_MIN = 4.0
PARENT_BRANCH = 'remaining'
PARENT_BUDGET_ENTRY = 'remaining'


def read(path):
    return json.loads((ROOT / path).read_text(encoding='utf-8'))


def visceral_flow_mL_s(flow_per_100g_min=ADIPOSE_FLOW_ML_PER_100G_MIN,
                       mass_kg=VISCERAL_ADIPOSE_MASS_KG):
    """Blood flow to the disaggregated depot, in mL/s."""
    return mass_kg*10.*flow_per_100g_min/60.   # 100 g units


def build(flow_per_100g_min=ADIPOSE_FLOW_ML_PER_100G_MIN, mass_kg=VISCERAL_ADIPOSE_MASS_KG,
          suffix=''):
    """Build one Arm B variant. Perfusion rate is a declared axis, not a constant:
    reported abdominal adipose blood flow after an overnight fast spans 3-5 mL
    per 100 g per minute, and the rate is actively regulated, rising with
    lipolytic stimulation. A variant is written per value rather than one
    number being chosen silently."""
    net = read(SOURCE_NETWORK)
    branches = net['branches']
    parent = next(b for b in branches if b['name'] == PARENT_BRANCH)
    limb = next(b for b in branches if b['name'] == 'limb_adipose')
    budget = read('data/m1-demand-budget.json')['baseline_mL_min']

    total_flow = sum(b['target_flow_mL_s'] for b in branches if b['from'] == 'AO')
    total_fraction = sum(b['baseline_fraction'] for b in branches if b['from'] == 'AO')
    drive = parent['R_mmHg_s_mL']*parent['target_flow_mL_s']        # (P_AO - P_VC)

    flow = visceral_flow_mL_s(flow_per_100g_min, mass_kg)
    if flow >= parent['target_flow_mL_s']:
        raise ValueError('Disaggregated depot cannot exceed its parent branch')
    fraction = parent['baseline_fraction']*flow/parent['target_flow_mL_s']

    # Oxygen demand per unit flow, taken from the project's own limb adipose
    # entry so the split does not import an outside specific-metabolic-rate.
    demand_per_flow = budget['limb_adipose']/limb['target_flow_mL_s']
    visceral_demand = demand_per_flow*flow

    visceral = dict(name='visceral_adipose', **{'from': 'AO'}, to='VC',
                    baseline_fraction=fraction, target_flow_mL_s=flow,
                    R_mmHg_s_mL=drive/flow,
                    **{'class': 'disaggregated from the residual branch; declared mass and '
                                'flow-per-mass, not a measurement of this subject'},
                    source='abdominal adipose blood flow 3-5 mL/100 g/min after overnight fast; '
                           'demand per unit flow taken from the project limb_adipose entry')
    reduced = dict(parent)
    reduced['target_flow_mL_s'] = parent['target_flow_mL_s']-flow
    reduced['baseline_fraction'] = parent['baseline_fraction']-fraction
    reduced['R_mmHg_s_mL'] = drive/reduced['target_flow_mL_s']
    reduced['class'] = parent['class']+'; visceral adipose disaggregated out under Arm B'

    out = dict(net)
    out['branches'] = [visceral if False else (reduced if b['name'] == PARENT_BRANCH else b)
                       for b in branches]
    out['branches'].insert(out['branches'].index(reduced)+1, visceral)
    out['status'] = ('Arm B variant of the frozen M3 evidence network: visceral adipose '
                     'disaggregated from the residual branch, totals preserved')
    out['arm_b_disaggregation'] = dict(
        parent=PARENT_BRANCH, mass_kg=mass_kg,
        flow_mL_per_100g_min=flow_per_100g_min,
        visceral_flow_mL_s=flow, visceral_demand_mL_O2_min=visceral_demand,
        demand_per_flow_mL_O2_min_per_mL_s=demand_per_flow,
        preserved=['total AO flow', 'baseline fraction sum', 'total oxygen demand',
                   'resistance construction R = (P_AO - P_VC) / Q'])

    # Checks: the split must conserve what the frozen network already closed.
    new_flow = sum(b['target_flow_mL_s'] for b in out['branches'] if b['from'] == 'AO')
    new_fraction = sum(b['baseline_fraction'] for b in out['branches'] if b['from'] == 'AO')
    assert abs(new_flow-total_flow) < 1e-9, 'total AO flow changed'
    assert abs(new_fraction-total_fraction) < 1e-12, 'baseline fraction sum changed'
    # The construction is R = (P_from - P_to)/Q, so only branches sharing the
    # parent's endpoints share its driving pressure: liver_arterial drains to
    # the hepatic node at a different pressure and has its own drive.
    for b in out['branches']:
        if b['from'] == parent['from'] and b['to'] == parent['to']:
            assert abs(b['R_mmHg_s_mL']*b['target_flow_mL_s']-drive) < 1e-6, b['name']

    # Metabolic side: move the same share of the parent budget entry across.
    mapping = read(SOURCE_MAPPING)
    weights = dict(mapping['mapping'][PARENT_BUDGET_ENTRY])
    share = visceral_demand/budget[PARENT_BUDGET_ENTRY]
    assert share < weights[PARENT_BRANCH], 'visceral share exceeds the residual weight'
    weights[PARENT_BRANCH] -= share
    weights['visceral_adipose'] = share
    assert abs(sum(weights.values())-1.) < 1e-12, 'mapping weights no longer sum to one'
    out_map = dict(mapping)
    out_map['mapping'] = dict(mapping['mapping'], **{PARENT_BUDGET_ENTRY: weights})
    out_map['status'] = ('Arm B variant: visceral adipose disaggregated from the residual '
                         'budget entry, weights still sum to one')

    network_path = OUT_NETWORK.replace('.yaml', suffix+'.yaml')
    mapping_path = OUT_MAPPING.replace('.json', suffix+'.json')
    (ROOT / network_path).write_text(json.dumps(out, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    (ROOT / mapping_path).write_text(json.dumps(out_map, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return out, out_map, visceral, reduced, visceral_demand


if __name__ == '__main__':
    rate = float(sys.argv[1]) if len(sys.argv) > 1 else ADIPOSE_FLOW_ML_PER_100G_MIN
    tag = '' if rate == ADIPOSE_FLOW_ML_PER_100G_MIN else '-atbf%g' % rate
    net, mapping, visceral, reduced, demand = build(rate, suffix=tag)
    print('wrote %s and %s' % (OUT_NETWORK.replace('.yaml', tag+'.yaml'),
                               OUT_MAPPING.replace('.json', tag+'.json')))
    print('  adipose blood flow %.1f mL/100 g/min over %.1f kg' % (rate, VISCERAL_ADIPOSE_MASS_KG))
    print('  visceral_adipose: flow %.4f mL/s (%.1f mL/min), fraction %.6f, R %.3f mmHg s/mL'
          % (visceral['target_flow_mL_s'], visceral['target_flow_mL_s']*60,
             visceral['baseline_fraction'], visceral['R_mmHg_s_mL']))
    print('  oxygen demand %.4f mL O2/min, from the project limb_adipose demand-to-flow ratio'
          % demand)
    print('  residual branch: flow %.4f -> %.4f mL/s, R %.3f -> %.3f'
          % (reduced['target_flow_mL_s']+visceral['target_flow_mL_s'],
             reduced['target_flow_mL_s'],
             (reduced['R_mmHg_s_mL']*reduced['target_flow_mL_s'])
             / (reduced['target_flow_mL_s']+visceral['target_flow_mL_s']),
             reduced['R_mmHg_s_mL']))
    print('  totals preserved: AO flow, fraction sum, oxygen demand, R construction')
