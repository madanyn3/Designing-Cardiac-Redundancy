"""Static inverse-design arithmetic, not a circulation simulation or viability test."""
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def audit():
    cfg = json.loads((ROOT / 'configs/project-extension.yaml').read_text())
    design = json.loads((ROOT / 'configs/auxiliary-heart-design.yaml').read_text())
    source = json.loads((ROOT / 'configs/baseline.yaml').read_text())
    branches = {b['name']: b for b in cfg['branches']}
    resistance = {name: b['R_mmHg_s_mL'] for name, b in branches.items()}
    ref = cfg['regional_reference']
    ao, vc, pu = (ref[k] for k in ['P_AO_mmHg', 'P_VC_mmHg', 'P_PU_mmHg'])
    # Reduce the portal/arterial junction analytically; hepatic return is in series.
    r_hepatic = 1 / (1 / resistance['gut'] + 1 / resistance['liver_arterial'])
    r_hepatic += resistance['hepatic_outflow']
    g_retained = 1 / r_hepatic
    g_retained += sum(1 / b['R_mmHg_s_mL'] for b in branches.values()
                      if b['from'] == 'AO' and b['to'] == 'VC')
    limb_names = [n for n, b in branches.items()
                  if b['from'] == 'AO' and b['to'].startswith('B_limb_')]
    g_limb = sum(1 / (resistance[n] + resistance[n + '_return']) for n in limb_names)
    r_open, r_closed = 1 / (g_retained + g_limb), 1 / g_retained
    q_lung = (ao - pu) / resistance['lungs']
    q_open = (ao - vc) / r_open
    q_closed = (ao - vc) / r_closed
    q_total_open, q_total_closed = .06 * (q_open + q_lung), .06 * (q_closed + q_lung)
    assert math.isclose(q_total_open, ref['Q_systemic_L_min'])
    assert math.isclose(q_total_closed, 4.592)
    assert r_closed > r_open
    assert math.isclose(r_closed / r_open, .99 / .81)
    # Equal kidney subdivision must preserve original parallel conductance.
    r_kidney_each = 2 * resistance['kidneys']
    assert math.isclose(2 / r_kidney_each, 1 / resistance['kidneys'])
    q_brain = (ao - vc) / resistance['brain']
    mean_target = design['pressure']['mean_AO_nominal_mmHg']
    required_brain_r = (mean_target - vc) / q_brain
    fixed_r_ratio = (mean_target - vc) / (ao - vc)
    assert math.isclose(required_brain_r / resistance['brain'], fixed_r_ratio)
    assert design['main_scenario']['brain_demand_multiplier'] == .8
    assert design['main_scenario']['brain_demand_multiplier_cases'] == [1, .8, .6]
    assert design['main_scenario']['brain_minimum_mean_flow_ratio_to_project_baseline'] is None
    assert design['main_scenario']['automatic_flow_to_demand_derating'] is False
    assert design['pressure']['proven_minimum_SBP_mmHg'] is None
    assert design['pressure']['proven_minimum_DBP_mmHg'] is None
    assert design['inverse_problem']['required_CO_L_min'] is None
    # Algebraic scaling at multiple states; not a validation of biological allometry.
    max_error = 0
    pars = source['source_parameters']
    for comp, suffix, coordinate in [('lvf_calculator', 'lvf', 80),
                                      ('rvf_calculator', 'rvf', 80),
                                      ('septum', 'spt', 4)]:
        p = pars[comp]
        for e in [0, .3, 1]:
            for alpha in [.25, .5, 1, 1.5]:
                original = (e * p['E_es_' + suffix] * (coordinate - p['V_d_' + suffix])
                            + (1 - e) * p['P_0_' + suffix]
                            * math.expm1(p['lambda_' + suffix] * (coordinate - p['V_0_' + suffix])))
                scaled = (e * p['E_es_' + suffix] / alpha
                          * (alpha * coordinate - alpha * p['V_d_' + suffix])
                          + (1 - e) * p['P_0_' + suffix]
                          * math.expm1(p['lambda_' + suffix] / alpha
                                       * (alpha * coordinate - alpha * p['V_0_' + suffix])))
                max_error = max(max_error, abs(original - scaled))
                assert math.isclose(original, scaled, abs_tol=1e-10)
    return {
        'status': 'passed static arithmetic only; no ODE, pump sizing or viability result',
        'classification': 'D conditional on A branch allocation and fixed boundary pressures',
        'AO_to_VC_R_mmHg_s_mL': {'open': r_open, 'limbs_closed': r_closed,
                                  'percent_increase': 100 * (r_closed / r_open - 1)},
        'AO_to_VC_Q_L_min': {'open': .06 * q_open, 'limbs_closed': .06 * q_closed},
        'lung_nutritive_Q_L_min': .06 * q_lung,
        'total_AO_outflow_L_min': {'open': q_total_open, 'limbs_closed': q_total_closed},
        'whole_outflow_operating_point_R_percent_increase': 100 * (q_total_open / q_total_closed - 1),
        'brain': {'anchor_Q_L_min': .06 * q_brain, 'anchor_R_mmHg_s_mL': resistance['brain'],
                  'R_to_preserve_baseline_flow_comparator_only_mmHg_s_mL': required_brain_r,
                  'fixed_R_flow_ratio_at_target_MAP': fixed_r_ratio},
        'kidney_each_R_mmHg_s_mL_equal_split': r_kidney_each,
        'scaled_wall_law_max_error_mmHg': max_error,
        'not_established': ['organ functional minima', 'minimum viable SBP/DBP',
                            'required auxiliary output', 'physical septal tissue volume',
                            'biological pump size', 'recoverable support duration']
    }


if __name__ == '__main__':
    result = audit()
    (ROOT / 'data/m2-design-audit.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, indent=2))
