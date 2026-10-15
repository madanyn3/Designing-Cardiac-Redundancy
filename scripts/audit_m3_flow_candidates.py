"""Calculate literature-linked flow candidates; never modify the hydraulic model.

Masses remain the M1 assumptions. Published perfusion coefficients are transferred
to those masses with the limitations described in m3-regional-flow-candidates.md.
This audit checks conversions, not physiological validity or oxygen adequacy.
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    masses = json.loads((ROOT / "data/m1-exploratory-assumptions.json").read_text())["mass_partition_kg"]
    rows = []
    for branch, omega in (("limb_muscle", 0.538e-3), ("limb_skin", 1.05e-3)):
        density = 1085.0  # Kalra et al., Table 1: kg/m3
        volume_ml = masses[branch] / density * 1e6
        flow = omega * volume_ml * 60
        rows.append(dict(branch=branch, mass_kg=masses[branch], density_kg_m3=density,
                         perfusion_s_inverse=omega, volume_ml=volume_ml,
                         flow_ml_min=flow, fraction_at_5600_ml_min=flow / 5600,
                         source="https://doi.org/10.1007/s11517-014-1229-8",
                         status="candidate model-coefficient transfer; assumed M1 mass"))
    flow = 3.53 * masses["limb_adipose"] * 10
    rows.append(dict(branch="limb_adipose", mass_kg=masses["limb_adipose"],
                     perfusion_ml_min_100g=3.53, flow_ml_min=flow,
                     fraction_at_5600_ml_min=flow / 5600,
                     source="https://pmc.ncbi.nlm.nih.gov/articles/PMC3287351/",
                     status="abdominal-to-limb transfer candidate, not limb measurement"))
    volume_ml = masses["respiratory_muscle"] / 1085 * 1e6
    flow = 9.9 * volume_ml / 100
    rows.append(dict(branch="respiratory_muscle", mass_kg=masses["respiratory_muscle"],
                     density_kg_m3=1085, volume_ml=volume_ml,
                     perfusion_ml_min_100ml=9.9, flow_ml_min=flow,
                     fraction_at_5600_ml_min=flow / 5600,
                     source="https://doi.org/10.1152/japplphysiol.01160.2007",
                     status="local intercostal-to-respiratory-pool surrogate; not diaphragm measurement"))
    assert all(row["flow_ml_min"] > 0 for row in rows)
    # Independent dimensional check: s^-1 -> mL/min per 100 mL.
    assert abs(rows[0]["flow_ml_min"] - 3.228 * rows[0]["volume_ml"] / 100) < 1e-9
    result = dict(status="research candidates only; incomplete distribution; no ODE run",
                  denominator_note="5600 mL/min is the existing project anchor, not a matched study CO",
                  candidates=rows,
                  limb_bone=dict(perfusion_ml_min_100ml=[5.5, 6.5],
                                 flow_ml_min_per_litre=[55, 65],
                                 source="https://doi.org/10.1186/1749-799X-3-36",
                                 status="requires matching bone-plus-marrow volume; do not substitute 4 kg as 4 L"),
                  unresolved=["limb_other composition/perfusion", "bulk skeletal volume",
                              "whole-body disjoint reconciliation", "respiratory anatomical transfer"])
    output = ROOT / "data/m3-regional-flow-candidates.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    for row in rows:
        print(f"{row['branch']}: {row['flow_ml_min']:.3f} mL/min; {100*row['fraction_at_5600_ml_min']:.3f}% of project anchor")


if __name__ == "__main__":
    main()
