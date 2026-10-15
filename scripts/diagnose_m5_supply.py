"""Disaggregate saved unmet demand without changing the frozen M5 experiment."""
from pathlib import Path
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
rows = []
for name in ('bio_immediate', 'bio_neural_delay', 'bio_slow_delay'):
    path = ROOT / 'reports/m5-final-comparison' / name
    result = json.loads((path / 'result.json').read_text(encoding='utf-8'))
    a = np.load(path / 'refined.npz')
    t = a['time_s']; actual = a['oxygen_consumed_mL_s']
    unmet = a['oxygen_unmet_mL_s']; ceiling = a['glucose_ceiling_mL_s']
    requested = actual + unmet
    glucose_limited = requested - ceiling
    oxygen_limited = ceiling - actual
    assert np.allclose(glucose_limited + oxygen_limited, unmet, atol=1e-12)
    assert np.max(abs(oxygen_limited)) < 1e-12
    j = result['organ_names'].index('aux_coronary_left')
    n = len(result['volume_names']); m = len(result['organ_names'])
    stores = a['states'][:, n + 5 + 2 * m + j]
    rows.append(dict(case=name,
        integrated_ledger_unmet_O2_equivalent_mL=result['records'][-1]['groups']['pump']['unmet_O2_mL'],
        sampled_glucose_limited_O2_equivalent_mL=float(np.trapezoid(glucose_limited[:, j], t)),
        sampled_additional_oxygen_limited_mL=float(np.trapezoid(oxygen_limited[:, j], t)),
        minimum_left_auxiliary_oxygen_store_mL=float(stores.min())))
out = ROOT / 'reports/m5-supply-diagnosis.json'
out.write_text(json.dumps(dict(scope='Saved-sample attribution; no changed hydraulics or physiological calibration',
    note='Sampled trapezoidal totals differ slightly from the ODE deficit ledger at narrow pulsatile events.', results=rows), indent=2) + '\n', encoding='utf-8')
print(out.read_text(encoding='utf-8'))
