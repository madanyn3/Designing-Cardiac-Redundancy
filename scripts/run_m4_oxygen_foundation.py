"""Verify zero-uptake oxygen transport against saved M3 moving blood volumes.

No hydraulic rerun or physiological support claim. Interpolate each saved
segment separately so regional control boundaries are never smoothed across.
"""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from cardiac_redundancy.oxygen import ConservativeOxygenTransport, OxygenParameters
from cardiac_redundancy.whole_network import WholeNetworkModel


def main():
    config_path = ROOT/'configs/m4-oxygen-reference.json'
    config = json.loads(config_path.read_text(encoding='utf-8'))
    source = ROOT/'reports/m3-whole-network/preservation_quarter-nominal'
    metadata = json.loads((source/'result.json').read_text(encoding='utf-8'))
    hydraulics = WholeNetworkModel(ROOT/'configs/baseline.yaml', ROOT/'configs/m3-evidence-regional.yaml',
        ROOT/'configs/m2-frozen-overlay.yaml', arterial_fraction=.2, inlet_resistance_fraction=.1)
    assert metadata['volume_names'] == list(hydraulics.volume_names)
    assert metadata['flow_names'] == list(hydraulics.flow_names)
    assert metadata['status'] == 'completed'
    model = ConservativeOxygenTransport(hydraulics.volume_names, hydraulics.edges, OxygenParameters(**config['parameters']))
    n = len(model.names)
    files = [source/f'segment-{i}.npz' for i in range(4)]
    trajectories = [dict(np.load(f)) for f in files]
    out = ROOT/'reports/m4-oxygen-foundation'
    out.mkdir(parents=True, exist_ok=True)
    results = []
    samples = {}
    for label, ventilation in config['ventilation_cases_mL_min_BTPS'].items():
        for solver, rtol, step in [('nominal', 1e-7, .02), ('refined', 1e-9, .005)]:
            y = model.initial_state(trajectories[0]['states'][:n, 0], config['initial_test_po2_mmHg'], config['initial_alveolar_fraction'])
            initial = y.copy(); all_t = []; all_y = []
            for segment in trajectories:
                t = segment['time']
                volumes = interp1d(t, segment['states'][:n], axis=1, bounds_error=True)
                flows = interp1d(t, segment['flows'], axis=0, bounds_error=True)
                sol = solve_ivp(lambda time, state: model.rhs(volumes(time), state, flows(time), ventilation/60),
                    (t[0], t[-1]), y, method='Radau', rtol=rtol, atol=rtol*.01, max_step=step, t_eval=t)
                if not sol.success:
                    raise RuntimeError(sol.message)
                y = sol.y[:, -1].copy()
                all_t.extend(t); all_y.extend(sol.y.T)
            states = np.array(all_y); times = np.array(all_t)
            drift = float(np.max(np.abs(model.balance(states)-model.balance(initial))))
            if drift > 1e-7 or np.min(states[:, :-1]) < 0:
                raise RuntimeError('Oxygen balance/positivity failure')
            samples[(label, solver)] = states
            np.savez_compressed(out/f'{label}-{solver}.npz', time_s=times, oxygen_states=states)
            results.append(dict(case=label, solver=solver, method='Radau', rtol=rtol, atol=rtol*.01, max_step_s=step,
                duration_s=float(times[-1]-times[0]), maximum_balance_error_mL=drift,
                initial_blood_oxygen_mL=float(initial[:n].sum()), final_blood_oxygen_mL=float(y[:n].sum()),
                initial_lung_oxygen_mL=float(initial[n]), final_lung_oxygen_mL=float(y[n]),
                net_external_oxygen_mL=float(y[-1]), minimum_store_mL=float(states[:, :-1].min())))
    differences = {label: float(np.max(np.abs(samples[(label, 'nominal')]-samples[(label, 'refined')])))
                   for label in config['ventilation_cases_mL_min_BTPS']}
    if max(differences.values()) > .001:
        raise RuntimeError('Oxygen solver refinement exceeds 0.001 mL state tolerance')
    paths = [Path(__file__), config_path, source/'result.json', *files,
             *sorted((ROOT/'src/cardiac_redundancy').glob('*.py')),
             *[ROOT/f'configs/{name}' for name in ('baseline.yaml', 'm3-evidence-regional.yaml', 'm2-frozen-overlay.yaml')]]
    result = dict(status='passed', scope='zero-uptake transport verification only',
        hydraulic_input=str(source.relative_to(ROOT)), volume_names=model.names,
        state_order=[*model.names, 'alveolar_oxygen', 'net_external_oxygen'],
        refinement_maximum_absolute_state_difference_mL=differences, cases=results,
        software=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
        hashes={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (out/'result.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    rows = ['| Gas input | Blood O2 at 8 s (mL) | Lung O2 at 8 s (mL) | Net external O2 (mL) | Balance error (mL) |',
            '|---|---:|---:|---:|---:|']
    for r in results:
        if r['solver'] == 'refined':
            rows.append(f"| {r['case']} | {r['final_blood_oxygen_mL']:.6f} | {r['final_lung_oxygen_mL']:.6f} | {r['net_external_oxygen_mL']:.6f} | {r['maximum_balance_error_mL']:.3g} |")
    report = '\n'.join([
        '# M4 oxygen transport foundation — 2026-10-07', '',
        '**Passed for numerical transport only. Tissue consumption is zero in every run.**', '',
        'Six runs reuse the eight-second M3 quarter-preservation hydraulic trajectory. Blood volumes and signed flows vary; native pumping remains off and branches reopen on the saved schedule. Uniform initial blood PO2 = 40 mmHg is a test fixture, not a physiological initial condition.', '',
        *rows, '',
        f'Maximum nominal/refined oxygen-state difference: {max(differences.values()):.6g} mL (acceptance < 0.001 mL). Balance acceptance: < 1e-7 mL; all saved stores nonnegative. Refinement covers the oxygen solver on the SAME interpolated hydraulic samples, not hydraulic sampling error.', '',
        'With apnea, lung oxygen transferred to blood comes entirely from the finite gas store. Oxygen does not vanish without consumption; this is conservation, not evidence that tissues tolerate apnea. Maintained/reduced labels refer to prescribed ventilation inputs, not demonstrated adequate oxygenation.', '',
        'Missing: organ uptake, extraction limits, normal metabolic baseline, tissue stores, oxygenator gas supply, fuel stock/flux, functional preservation and recovery. No repair-window estimate follows from these runs.', '',
        'Reproduce: `.venv\\Scripts\\python.exe scripts/run_m4_oxygen_foundation.py`. Inputs, versions, solver settings and hashes: [result.json](result.json). Equations and sources: [specification](../../docs/m4-oxygen-specification.md).'])
    (out/'report.md').write_text(report+'\n', encoding='utf-8')
    print(json.dumps(dict(status='passed', runs=len(results), maximum_refinement_difference_mL=max(differences.values()))))


if __name__ == '__main__':
    main()
