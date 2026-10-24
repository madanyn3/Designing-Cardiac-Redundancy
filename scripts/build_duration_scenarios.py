"""Supported duration across the three declared physiology scenarios.

Table 3 of the preprint had no artifact: the scenario figures were carried in
prose only. This computes them from the declared scenario parameters and the
corrected duration law, which applies the liver-function fraction the earlier
analytic block omitted, and writes one artifact the table reads.
"""
from pathlib import Path
import json, sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
from cardiac_redundancy.duration import SubstrateBudget, supported_duration

OUT = 'reports/m8-duration'
PROTOCOL = json.loads((ROOT / 'configs/m8-protocol.json').read_text(encoding='utf-8'))
#: Total oxidative demand per scenario, from the coupled replays.
DEMAND_ML_O2_MIN = dict(favourable=87.2, reference=108.2, adverse=144.8)
FLOOR_MMOL_L = 3.0
#: Fields of SubstrateBudget a scenario may set; others are scenario factors
#: that act on the circulation rather than on the substrate budget.
BUDGET_FIELDS = ('glycogen_mmol_L_liver', 'release_cap_mmol_min', 'liver_function_fraction')


def budget_for(scenario, floor=0.):
    spec = PROTOCOL['scenarios'][scenario]
    kwargs = {k: v for k, v in spec.items() if k in BUDGET_FIELDS}
    return SubstrateBudget(glucose_floor_mmol_L=floor, **kwargs)


def main():
    rows = []
    for scenario in ('favourable', 'reference', 'adverse'):
        demand = DEMAND_ML_O2_MIN[scenario]
        free = supported_duration(demand, budget_for(scenario))
        floored = supported_duration(demand, budget_for(scenario, FLOOR_MMOL_L))
        rows.append(dict(
            scenario=scenario, demand_mL_O2_min=demand,
            export_mmol_min=free['effective_export_mmol_min'],
            release_cap_mmol_min=free['release_cap_mmol_min'],
            liver_function_fraction=free['liver_function_fraction'],
            glucose_demand_mmol_min=free['glucose_demand_mmol_min'],
            net_drain_mmol_min=free['net_drain_mmol_min'],
            duration_min=free['minutes'], duration_at_floor_min=floored['minutes'],
            limited_by=free['limited_by'],
            glycogen_total_mmol=free['glycogen_total_mmol'],
            store_limited_minutes=free['store_limited_minutes']))
    result = dict(
        scope='supported duration across the three declared physiology scenarios',
        protocol_version=PROTOCOL['version'], floor_mmol_L=FLOOR_MMOL_L,
        demand_source='total oxidative demand from the coupled replays',
        rows=rows,
        note='the liver-function fraction is applied, which the earlier analytic block in '
             'verification-60s.json omitted: that block reports 39.57 min at reference against '
             '33.8 here, and the difference is exactly that factor',
        interpretation='interval over which the declared delivery screens can be met; never '
                       'survival, recovery or injury')
    dest = ROOT / OUT; dest.mkdir(parents=True, exist_ok=True)
    (dest / 'duration-scenarios.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('%-12s %9s %9s %9s %9s' % ('scenario', 'demand', 'export', 'minutes', 'at floor'))
    for r in rows:
        print('%-12s %9.1f %9.3f %9.1f %9.1f' % (
            r['scenario'], r['demand_mL_O2_min'], r['export_mmol_min'],
            r['duration_min'], r['duration_at_floor_min']))
    print('\nstore at reference: %.1f mmol, %.0f min at the effective cap' % (
        rows[1]['glycogen_total_mmol'], rows[1]['store_limited_minutes']))
    return result


if __name__ == '__main__':
    main()
