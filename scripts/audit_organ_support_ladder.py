"""Reproduce protected-set demand subtotals, not flow, sizing or survival."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ADDITIONS = [
    ("L0", "Brain + native-lung/respiratory-muscle proxy", ["brain", "lungs", "respiratory_muscle"]),
    ("L1", "+ native myocardium", ["heart"]),
    ("L2", "+ kidney A", ["kidney_A"]),
    ("L3", "+ liver", ["liver"]),
    ("L4", "+ gut", ["gut"]),
    ("L5", "+ kidney B", ["kidney_B"]),
    ("L6", "+ remaining-body pool", ["remaining"]),
    ("L7", "+ limb skin", ["limb_skin"]),
    ("L8", "+ limb bone", ["limb_bone"]),
    ("L9", "+ other limb tissue", ["limb_other"]),
    ("L10", "+ limb adipose", ["limb_adipose"]),
    ("L11", "+ limb skeletal muscle", ["limb_muscle"]),
]


def main():
    budget = json.loads((ROOT / "data/m1-demand-budget.json").read_text())
    output = {
        "classification": "D conditional on M1 A inputs; protected-set demand only",
        "scope": "Native tissues only; auxiliary self-demand and unresolved respiratory dependencies are additional",
        "exclusion": "Unprotected tissue demand is retained, not set to zero",
        "kidney_split": "A: half of the existing paired demand for each kidney",
        "horizons_min": [1, 5, 15, 30, 60],
        "exploratory_horizons_min": [120],
        "scenarios": {},
    }
    for name, scenario in budget["scenarios"].items():
        demand = dict(scenario["organ_mL_min"])
        paired = demand.pop("kidneys")
        demand.update(kidney_A=paired / 2, kidney_B=paired / 2)
        whole = sum(demand.values())
        assert abs(whole - scenario["whole_body_mL_min"]) < 1e-8
        selected = []
        rows = []
        for stage, label, added in ADDITIONS:
            assert not set(added).intersection(selected)
            selected.extend(added)
            protected = sum(demand[key] for key in selected)
            rows.append({"stage": stage, "label": label, "protected_pools": list(selected),
                         "added_demand_mL_O2_min": sum(demand[key] for key in added),
                         "protected_demand_mL_O2_min": protected,
                         "unprotected_demand_mL_O2_min": whole - protected})
        assert set(selected) == set(demand)
        assert abs(rows[-1]["unprotected_demand_mL_O2_min"]) < 1e-8
        assert abs(rows[6]["protected_demand_mL_O2_min"] - scenario["protected_set_mL_min"]) < 1e-8
        output["scenarios"][name] = {"whole_native_demand_mL_O2_min": whole, "stages": rows}
    (ROOT / "data/organ-support-ladder.json").write_text(json.dumps(output, indent=2) + "\n")
    print("PASS: disjoint additions, kidney split, full-body totals and M1 protected-set reconciliation")
    for row in output["scenarios"]["modest_reduction"]["stages"]:
        print(f'{row["stage"]}: {row["protected_demand_mL_O2_min"]:.2f} mL O2/min')


if __name__ == "__main__":
    main()
