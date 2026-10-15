"""Reproduce the M1 assumed demand budget; no circulation or survival simulation."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def calculate(c):
    eps = c['oxygen_energy_equivalent_kcal_per_L']
    convert = lambda energy: energy * 1000 / (1440 * eps)
    masses = c['mass_partition_kg']
    assert abs(sum(masses.values()) - c['body_mass_kg']) < 1e-9
    rates = c['energy_rates_kcal_per_kg_day']
    base = {k: convert(masses[k]*rates[k]) for k in ['brain','heart','kidneys','liver','limb_muscle']}
    base['respiratory_muscle'] = c['respiratory_reference_mL_min']
    base['gut'] = c['gut_reference_mL_min']
    base['lungs'] = c['lung_reference_mL_min']
    # All skeletal muscle includes limb, respiratory and axial muscle.
    residual_mass = c['body_mass_kg'] - sum(masses[k] for k in ['brain','heart','kidneys','liver']) - c['total_skeletal_muscle_kg'] - c['adipose_kg']
    axial_mass = c['total_skeletal_muscle_kg'] - masses['limb_muscle'] - masses['respiratory_muscle']
    base['remaining'] = convert(residual_mass*rates['residual'] + c['adipose_kg']*rates['adipose'] + axial_mass*rates['skeletal_muscle']) - base['gut'] - base['lungs']
    # Repartition existing residual/adipose energy; do not add new body mass or energy.
    for tissue, spec in c['limb_nonmuscle_partition'].items():
        base[tissue] = convert(masses[tissue]*spec['energy_rate_kcal_kg_day'])
        base['remaining'] -= base[tissue]
    assert min(base.values()) > 0
    expected = convert(sum(masses[k]*rates[k] for k in ['brain','heart','kidneys','liver']) + c['total_skeletal_muscle_kg']*rates['skeletal_muscle'] + c['adipose_kg']*rates['adipose'] + residual_mass*rates['residual'])
    # Replace, rather than add twice, respiratory muscle's allocated resting energy.
    expected += base['respiratory_muscle'] - convert(masses['respiratory_muscle']*rates['skeletal_muscle'])
    assert abs(sum(base.values())-expected) < 1e-9
    scenarios = {}
    for name, multipliers in c['demand_multipliers'].items():
        assert set(multipliers) == set(base)
        values = {k: base[k]*multipliers[k] for k in base}
        assert all(0 < m <= 1 for m in multipliers.values())
        scenarios[name] = {'organ_mL_min': values, 'whole_body_mL_min': sum(values.values()),
                           'protected_set_mL_min': sum(v for k,v in values.items() if not k.startswith('limb_')),
                           'limb_demand_mL_min': sum(v for k,v in values.items() if k.startswith('limb_'))}
    store = c['muscle_store_assumptions']
    concentration = store['typeI_fraction']*.70 + (1-store['typeI_fraction'])*.49
    sat = lambda p: p/(p+store['P50_mmHg'])
    bound_per_kg = concentration*store['gas_mL_per_mmol']/store['density_kg_per_L']*(sat(store['initial_PO2_mmHg'])-sat(store['floor_PO2_mmHg']))
    return {'classification':'D conditional on A inputs; not pump sizing or viability',
            'baseline_mL_min':base, 'scenarios':scenarios,
            'muscle_bound_store_mL_per_kg':bound_per_kg,
            'limb_bound_store_mL':bound_per_kg*masses['limb_muscle'],
            'respiratory_bound_store_mL':bound_per_kg*masses['respiratory_muscle'],
            'checks':'mass balance; disjoint energy allocation; positive demands; scenario coverage'}

if __name__ == '__main__':
    config = json.loads((ROOT/'data/m1-exploratory-assumptions.json').read_text(encoding='utf-8'))
    result = calculate(config)
    (ROOT/'data/m1-demand-budget.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
