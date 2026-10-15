"""Regenerate M5 comparison tables, plots and an immutable-input audit."""
from pathlib import Path
import hashlib
import json
import sys
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from cardiac_redundancy.architecture_energy import external_input
from cardiac_redundancy.activation_timing import ActivationTiming
from cardiac_redundancy.pump_energy import MMHG_ML_TO_J as K
from cardiac_redundancy.preservation_energy import oxygen_energy
OUT = ROOT / 'reports/m5-final-comparison'

def read(path):
    return json.loads((ROOT / path).read_text(encoding='utf-8'))

def main():
    config = read('configs/m5-final-comparison.json')
    timings = read('configs/m5-activation-timing.json')['cases']
    seed = np.load(ROOT / 'reports/m5-joint-preservation/reference_g25/fine-refined.npz')['states'][:, 0]
    results = []
    common_initial = None
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    rows = ['# M5 final architecture and startup comparison', '',
        '2026-10-14. Four paired eight-second experiments; exploratory assumptions, not a validated support duration.', '',
        'All cases start from the identical 5250 mL hydraulic inventory and native failure. Both devices have the same Smith-shaped chamber law, valve parameters, cavity scale and rate. The mechanical reference supplies that actuation externally and has no auxiliary myocardium. It is a finite chamber reference, not a realizable device or unlimited flow source.', '',
        'Brain, native coronary branches, one kidney, hepatic artery and respiratory muscle branches retain full reference conductance. Other systemic branches have conductance fraction 0.25. Supported/preserved lungs have separate pulmonary and nutritive flow paths and gas stores: assumed equal lung shares, pulmonary factors 1/0.1, nutritive factors 1/0.25 and ventilation factors 1/0. These are sensitivity settings, not physiological minima. Both lungs share PA/PU blood reservoirs.', '',
        '| Case | Release / ramp s | Final 2 s mean AO mmHg | AO min–max mmHg | Left output L/min | Chamber work J over 8 s |',
        '|---|---:|---:|---:|---:|---:|']
    energy_rows = ['| Case | Requested auxiliary O2 mL | Unmet auxiliary O2 mL | First auxiliary shortfall s | Requested BIO equivalent J | External input J at eta=0.5 |', '|---|---:|---:|---|---:|---:|']
    for case in config['cases']:
        name = case['name']; base = f'reports/m5-final-comparison/{name}'
        h = read(base + '/hydraulics.json'); r = read(base + '/result.json')
        assert h['status'] == r['status'] == 'verified'
        a = np.load(ROOT / base / 'fine-refined.npz'); b = np.load(ROOT / base / 'refined.npz')
        n = len(r['volume_names']); m = len(r['organ_names']); initial = b['states'][0]
        k = n + 5 + 3 * m
        native_stores = {name: float(initial[n + 5 + 2 * m + j]) for j, name in enumerate(r['organ_names']) if not name.startswith('aux_coronary_')}
        comparable = (initial[:n + 1].tolist(), float(initial[n + 3 + 2 * m]), native_stores, initial[k:k + n + 1].tolist())
        if common_initial is None: common_initial = comparable
        assert comparable == common_initial, 'Native blood/gas/tissue/glucose initialization differs'
        t, y, p = a['time_s'], a['states'], a['pressures_mmHg']
        assert np.array_equal(y[:, 0], seed)
        vn = h['volume_names']; final = h['runs'][-1]; record = r['records'][-1]
        timing = ActivationTiming(**timings[case['timing']])
        left = K * (p[:, vn.index('AO')] * y[23] - p[:, vn.index('PU')] * y[22])
        right = K * (p[:, vn.index('PA')] * y[25] - p[:, vn.index('VC')] * y[24])
        variants = [external_input(t, np.c_[left, right], eta) for eta in (.25, .5, .75, 1.)]
        signed = sum(v['work']['port_J'] for v in final['sides'].values())
        assert abs(signed - variants[0]['signed_port_J']) < 1e-8
        work = sum(v['work']['chamber_J'] for v in final['sides'].values())
        pump = record['groups']['pump']; requested = pump['requested_O2_mL']
        if case['device'] == 'BIO':
            integral_g = max(0., 8 - timing.latency_s - timing.ramp_s / 2)
            expected = sum(d['maintenance_mL_min'] * 8 + d['activation_mL_min'] * integral_g + d['work_mL_min'] * 8 for d in r['case']['sides'].values()) / 60
            assert abs(expected - requested) < 1e-7
        else:
            assert requested == 0 and r['case']['auxiliary_mL_min'] == 0
        # Diagnostic only: a post-ramp two-cycle window satisfying the declared
        # pressure envelope, forward output, and no previous supply deficit.
        # This cannot establish recovery and never equates ramp completion to support.
        bt = b['time_s']; unmet = b['oxygen_unmet_mL_s']
        first_window = None
        for end in np.arange(timing.latency_s + timing.ramp_s + 1.5, 8.0001, .005):
            mask = (t >= end - 1.5) & (t <= end)
            tt = t[mask]; ao = p[mask, vn.index('AO')]
            if len(tt) < 2: continue
            avg = np.trapezoid(ao, tt) / (tt[-1] - tt[0])
            if (avg >= 65 and ao.min() >= 50 and ao.max() <= 90
                    and np.mean(y[23, mask]) > 0 and np.mean(y[25, mask]) > 0
                    and not np.any(unmet[bt <= end] > 1e-8)):
                first_window = [float(end - 1.5), float(end)]; break
        entry = dict(case=name, device=case['device'], release_s=timing.latency_s,
            ramp_s=timing.ramp_s, final_hydraulics={k: final[k] for k in ['mean_AO_mmHg','min_AO_mmHg','max_AO_mmHg','auxiliary_left_L_min','pressure_screen']},
            chamber_work_J=work, groups=record['groups'],
            requested_BIO_energy_equivalent_J=float(oxygen_energy(requested)),
            port_energy_sensitivities=variants, first_passing_post_ramp_window_s=first_window,
            final_lung_gas_fractions=b['lung_fractions'][-1].tolist(),
            paired_hydraulic_scaled_errors=h['scaled_errors'],
            paired_oxygen_max_difference_mL=r['paired_max_oxygen_state_difference'],
            paired_glucose_max_difference_mmol=r['paired_max_glucose_state_difference_mmol'])
        results.append(entry)
        rows.append(f"| {name} | {timing.latency_s:g} / {timing.ramp_s:g} | {final['mean_AO_mmHg']:.3f} | {final['min_AO_mmHg']:.3f}–{final['max_AO_mmHg']:.3f} | {final['auxiliary_left_L_min']:.4f} | {work:.4f} |")
        bracket = pump['first_shortfall_bracket_s']
        bracket = 'none' if bracket is None else f'{bracket[0]:.3f}–{bracket[1]:.3f}'
        energy_rows.append(f"| {name} | {requested:.6f} | {pump['unmet_O2_mL']:.6f} | {bracket} | {entry['requested_BIO_energy_equivalent_J']:.4f} | {variants[1]['external_input_J']:.4f} |")
        axes[0].plot(t, p[:, vn.index('AO')], label=name)
        dt = np.diff(bt)
        debt = np.r_[0., np.cumsum((unmet[:-1].sum(axis=1) + unmet[1:].sum(axis=1)) * dt / 2)]
        axes[1].plot(bt, debt, label=name)
        axes[2].plot(bt, b['lung_fractions'][:, 0], label=name + ' supported')
        axes[2].plot(bt, b['lung_fractions'][:, 1], linestyle='--', label=name + ' preserved')
    rows += ['', 'All final-window pressure screens fail. Flow is auxiliary left output, not a guarantee of organ delivery. The mechanical reference passes the first post-ramp diagnostic window at 1.0–2.5 s, but fails the final pressure screen; no BIO case passes. This is transient screening success, not sustained support. The diagnostic is: a 1.5 s window with mean AO ≥65, minimum ≥50, maximum ≤90 mmHg, positive mean bilateral output and no previous sampled oxygen shortfall in any modeled organ. Windows are scanned every 5 ms after the ramp. These are declared project screens, not validated clinical survival thresholds.', '', *energy_rows, '',
        'External-input columns are hypothetical port-energy equivalents for each trajectory; only the mechanical row draws that external source in the modeled comparison. Efficiencies 0.25/0.5/0.75/1 and zero overhead are assumptions, with full outputs in summary.json. Positive port work is summed separately for both stages before division by efficiency; returned work is not credited. BIO requested energy uses the existing oxygen-energy conversion and is not measured ATP consumption. Neither column establishes physical pump efficiency or a winning architecture.', '',
        'All four cases have zero modeled unmet demand in the mandatory native group (66.9539 mL O2/min) and the other preservation group (31.2158 mL O2/min) over eight seconds. This uses reduced demand assumptions and finite initial tissue oxygen stores. It does not prove delivery alone sustains those organs, nor their eventual recovery. All BIO cases have auxiliary left-myocardial shortfall; subsequent prescribed contraction is counterfactual under the declared demand model. Delayed shortfall is not a longer repair window.', '',
        'BIO maintenance runs before release. Activation follows the configured ramp. Work-related demand is distributed with that ramp and normalized to preserve the integrated eight-second work-related oxygen cost. This is an offline workload proxy, not a causal metabolic controller. Unmeasured one-time activation energy remains unknown and excluded, not measured zero. Blood stasis, thrombosis, dormant viability, anatomy, CO2/pH and injury remain unresolved.', '',
        'The two gas stores partition the existing initial alveolar store; they do not add oxygen. The unventilated preserved lung can initially exchange oxygen from its finite store. Equal lung volume/capacity allocation is an explicit assumption. No extra regional blood storage is introduced. The mechanical case removes auxiliary nutritive branches and oxygen demand while preserving the initial physical blood state.', '',
        'Verification: paired hydraulic tolerances/time steps, valve-event agreement, blood conservation, septal residual, chamber/port/valve/inertial energy identities, paired metabolic integrations, oxygen/glucose conservation, demand identities and oxidation stoichiometry. Separate unit tests verify bilateral gas allocation/conservation and mandatory routing. Full suite: 235 passed. Prior frozen manifests are checked by this generator.', '',
        'Reproduce: `.venv/Scripts/python.exe scripts/run_m5_final_comparison.py`, then `.venv/Scripts/python.exe scripts/report_m5_final_comparison.py`. Run `.venv/Scripts/python.exe -m pytest -q --junitxml=reports/m5-final-tests.xml` for the suite.', '',
        '[Closure, evidence and limitations](../../docs/m5-closure.md) · [Implementation explanation](../../docs/m5-bilateral-implementation.md)', '', '![Pressure, deficits and separate lung stores](comparison.png)']
    print('Post-ramp windows:', [(r['case'], r['first_passing_post_ramp_window_s']) for r in results])
    axes[0].axhline(65, color='black', linestyle=':', label='mean-pressure target (not instantaneous floor)')
    for ax, ylabel in zip(axes, ['Aortic pressure (mmHg)', 'Total unmet oxygen (mL)', 'Alveolar oxygen fraction']):
        ax.set_ylabel(ylabel); ax.grid(alpha=.2); ax.legend(fontsize=7, ncol=2)
    axes[-1].set_xlabel('Seconds after native arrest')
    fig.suptitle('M5: matched initial state; prescribed support, not demonstrated viability')
    fig.tight_layout(); fig.savefig(OUT / 'comparison.png', dpi=160); plt.close(fig)
    (OUT / 'summary.json').write_text(json.dumps(results, indent=2) + '\n', encoding='utf-8')
    (OUT / 'report.md').write_text('\n'.join(rows) + '\n', encoding='utf-8')
    old = ['data/m2-freeze-manifest.json', 'reports/m5-t51/manifest.json', 'reports/m5-biological-energy/audit.json', 'reports/m5-normal-comparison/manifest.json', 'reports/m5-native-preservation/manifest.json', 'reports/m5-organ-preservation/manifest.json', 'reports/m5-joint-preservation/manifest.json', 'reports/m5-primer-audit/manifest.json', 'reports/m5-architecture-foundation/manifest.json']
    checked = 0
    for path in old:
        manifest = read(path)
        hashes = manifest.get('hashes', manifest.get('sha256', {}))
        assert hashes, path
        for file, expected in hashes.items():
            assert hashlib.sha256((ROOT / file).read_bytes()).hexdigest() == expected, file
            checked += 1
    suites = list(ET.parse(ROOT / 'reports/m5-final-tests.xml').getroot().iter('testsuite'))
    assert sum(int(s.get('tests', 0)) for s in suites) == 235
    assert all(int(s.get('errors', 0)) + int(s.get('failures', 0)) == 0 for s in suites)
    files = set((ROOT / 'src/cardiac_redundancy').glob('*.py'))
    files.update((ROOT / 'tests').glob('test_bilateral_support.py'))
    files.update(ROOT / p for p in ['scripts/run_m5_final_comparison.py', 'scripts/run_m5_bilateral_metabolism.py', 'scripts/report_m5_final_comparison.py', 'configs/m5-final-comparison.json', 'configs/m5-activation-timing.json', 'configs/m5-myocardial-demand.json', 'reports/m5-final-tests.xml', 'docs/m5-closure.md', 'docs/m5-bilateral-implementation.md'])
    files.update((ROOT / 'configs').glob('*'))
    files.update(ROOT / p for p in ['scripts/run_m4_respiratory.py', 'data/m1-demand-budget.json'])
    for case in config['cases']:
        base = OUT / case['name']; files.update(base.glob('*'))
        for meta in ['result.json', 'hydraulics.json']:
            files.update(ROOT / p for p in json.loads((base / meta).read_text(encoding='utf-8'))['inputs'])
    files.update(OUT / p for p in ['summary.json', 'report.md', 'comparison.png'])
    (OUT / 'manifest.json').write_text(json.dumps(dict(status='verified_qualified_closure', tests=235, prior_hash_checks=checked, hashes={str(p.relative_to(ROOT)).replace('\\','/'): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}), indent=2) + '\n', encoding='utf-8')
    print(f'Four cases verified; {checked} prior hashes intact; 235 tests; report generated.')

if __name__ == '__main__':
    main()
