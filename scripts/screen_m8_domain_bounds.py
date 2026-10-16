"""T8.1 pre-search screen: which cavity scales the fixed inventory admits.

The auxiliary cavities are filled by a one-time physical transfer from the
pulmonary-vein and vena-cava pools (configs/m2-frozen-overlay.yaml donors).
Total blood volume is fixed, so a larger auxiliary cavity is paid for by a
smaller donor pool, and at some scale a donor reaches its unstressed volume.

This screen reports that boundary on the declared initial state before any
integration. It is an inventory/model-domain bound, not a biological size
limit, and it says nothing about whether an admissible scale performs well.
"""
from pathlib import Path
import json, sys
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / 'src'))
from cardiac_redundancy.sensitivity_circulation import SensitivityCirculation
from cardiac_redundancy.biological_pump import BiologicalParameters

OUT = ROOT / 'reports/m8-protocol'
STATE = 'reports/m5-joint-preservation/reference_g25/fine-refined.npz'


def model(scale, **knobs):
    return SensitivityCirculation(
        ROOT / 'configs/baseline.yaml', ROOT / 'configs/m3-evidence-regional.yaml',
        ROOT / 'configs/m2-frozen-overlay.yaml', arterial_fraction=.2,
        inlet_resistance_fraction=.1,
        biological=BiologicalParameters(left_scale=scale, right_scale=scale,
                                        septal_scale=scale, pericardial_scale=scale,
                                        **knobs))


def rescaled_initial_state(m, saved):
    """Move the saved reference state to this model's cavity allocation.

    The saved trajectory was produced at the reference scale, so its auxiliary
    cavities hold the reference volumes. Reallocating with the same donor rule
    keeps total inventory fixed; only the split changes.
    """
    y = np.asarray(saved, dtype=float).copy()
    cavities = m.initial_state[16:18]
    for j, donor in enumerate(('PU', 'VC')):
        y[m.volume_names.index(donor)] += y[16+j]-cavities[j]
    y[16:18] = cavities
    return y


def main():
    saved = np.load(ROOT / STATE)['states'][:, 0]
    total = float(saved[:18].sum())
    records = []
    for scale in [round(.05*k, 3) for k in range(2, 21)]:
        try:
            m = model(scale)
        except ValueError as exc:
            records.append(dict(scale=scale, admissible=False, reason=repr(exc))); continue
        y = rescaled_initial_state(m, saved)
        margins = m.domain_margins(y)
        records.append(dict(scale=scale, admissible=min(margins.values()) > 0,
            cavity_mL=[float(v) for v in y[16:18]],
            inventory_error_mL=float(abs(y[:18].sum()-total)),
            donor_margins_mL={k: margins[k] for k in ('PU', 'VC')},
            tightest=min(margins, key=margins.get),
            tightest_margin_mL=float(min(margins.values()))))
    admissible = [r['scale'] for r in records if r['admissible']]
    result = dict(
        scope='initial-state inventory and model-domain screen at equal cavity scales',
        reference_state=STATE, total_inventory_mL=total,
        declared_parameter_bounds=dict(scale=[.1, 1.], period_s=[.4, 1.5],
                                       contractility=[.5, 2.]),
        admissible_scale_range=[min(admissible), max(admissible)] if admissible else None,
        interpretation='an inventory/domain bound on the declared initial state; '
                       'admissibility here is not adequacy and not an anatomical limit',
        records=records)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'domain-bounds.json').write_text(
        json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    for r in records:
        print(f"scale {r['scale']:.2f} admissible={r['admissible']}"
              + (f" tightest={r['tightest']} {r['tightest_margin_mL']:.1f} mL"
                 if r['admissible'] else f" reason={r.get('reason','domain')[:60]}"), flush=True)
    return result


if __name__ == '__main__':
    main()
