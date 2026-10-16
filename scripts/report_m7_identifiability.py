"""Exact observation audit of the M3 parallel-resistor comparator, not BIO fitting.

Pressure differences are known inputs here. Unknown venous pressure, regional
storage, closing pressures and autoregulation require separate experiments.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def observation_matrix(pressure_drop, branches, observed=()):
    """Stack total-flow observations and optional individual branch flow traces.

    Columns are conductances [mL/s/mmHg]. Repeated pressure samples change
    precision, but cannot distinguish identical parallel paths by total flow.
    """
    dp = np.asarray(pressure_drop, dtype=float)
    if dp.ndim != 1 or not len(dp) or not np.all(np.isfinite(dp)):
        raise ValueError('Pressure drops must be a finite nonempty vector')
    if branches < 1 or len(set(observed)) != len(observed) or any(i < 0 or i >= branches for i in observed):
        raise ValueError('Invalid branch observation set')
    rows = [dp[:, None] * np.ones((1, branches))]
    for i in observed:
        rows.append(dp[:, None] * np.eye(branches)[i])
    return np.vstack(rows)


def main():
    path = ROOT / 'configs/m3-evidence-regional.yaml'
    config = json.loads(path.read_text(encoding='utf-8-sig'))
    rows = [r for r in config['branches'] if r['from'] == 'AO' and r['to'] == 'VC']
    conductance = np.array([1/r['R_mmHg_s_mL'] for r in rows])
    n = len(rows)
    # Synthetic excitation: not a measured waveform or a physiological target.
    dp = np.array([70., 80., 90., 100.])
    result = {'scope': 'Exact linear parallel-resistor comparator only',
              'branch_names': [r['name'] for r in rows], 'branch_count': n,
              'pressure_drop_mmHg': dp.tolist(), 'observation_sets': {}}
    for name, observed in [('total_only', ()), ('total_plus_one_branch', (0,)),
                           ('total_plus_all_but_one', tuple(range(n-1)))]:
        a = observation_matrix(dp, n, observed)
        y = a @ conductance
        recovered = np.linalg.lstsq(a, y, rcond=None)[0]
        rank = int(np.linalg.matrix_rank(a))
        result['observation_sets'][name] = dict(rank=rank, nullity=n-rank,
            max_flow_residual_mL_s=float(np.max(np.abs(a @ recovered-y))),
            max_relative_parameter_error=float(np.max(np.abs(recovered/conductance-1))))
    # Explicit positive alternative: transfer half the smaller conductance.
    alternate = conductance.copy()
    transfer = .5 * min(conductance[:2])
    alternate[0] += transfer
    alternate[1] -= transfer
    a = observation_matrix(dp, n)
    result['indistinguishable_pair'] = dict(original_conductance=conductance.tolist(),
        alternate_conductance=alternate.tolist(),
        max_total_flow_difference_mL_s=float(np.max(np.abs(a @ (alternate-conductance)))),
        branch_flow_change_at_90_mmHg=(90*(alternate-conductance)).tolist())
    result['sha256'] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in [path, Path(__file__)]}
    out = ROOT / 'reports/m7-identifiability'
    out.mkdir(parents=True, exist_ok=True)
    (out/'parallel-audit.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result['observation_sets'], indent=2))


if __name__ == '__main__':
    main()
