"""Reproduce historical primer arithmetic; not a physiological pump model."""
import json
import math
from pathlib import Path

def calculate():
    # Historical assumptions, deliberately NOT admitted as design inputs.
    external, rate, ef, pressure, stress, density = 700, 60, .45, 90, 22, 1.05
    alpha = pressure * .133322368 / (2 * stress)
    mass_per_flow = 6 * density * alpha / (rate * ef)  # two identical thin shells
    feedback = mass_per_flow * 250 / 300
    q = external / (1 - feedback)
    radius = (3 * q / (rate * ef) / (4 * math.pi)) ** (1 / 3)
    ratio = ((1 + alpha) ** 3 - 1) / (3 * alpha)
    iterates = [external]
    for _ in range(3):
        iterates.append(external + feedback * iterates[-1])
    return {
        "classification": "D conditional on historical assumptions; not validated design",
        "h_over_r": alpha, "feedback_gain": feedback,
        "flow_fixed_point_mL_min": q, "radius_cm": radius,
        "wall_cm": alpha * radius, "two_thin_shells_mass_g": mass_per_flow * q,
        "exact_to_thin_shell_mass_ratio_same_dimensions": ratio,
        "two_exact_shells_mass_g_same_dimensions": mass_per_flow * q * ratio,
        "flow_iterates_mL_min": iterates,
        "error_after_two_updates_mL_min": abs(q - iterates[2]),
        "one_path_hydraulic_power_W_if_head_70_mmHg": q / 1000 * 70 * .0022220394667,
        "beats_70_years_if_100000_per_day": 100000 * 365.25 * 70,
    }

if __name__ == "__main__":
    results = calculate()
    assert 0 < results["feedback_gain"] < 1
    assert abs(results["flow_fixed_point_mL_min"] - (700 + results["feedback_gain"] * results["flow_fixed_point_mL_min"])) < 1e-10
    output = Path(__file__).resolve().parents[1] / "data" / "primer-arithmetic-audit.json"
    output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2))
