"""Finite-difference observation experiment on the frozen native Smith model.

Each case reaches its own periodic state. Old M3 trajectories only initialize
warm-up; they do not substitute for a new equilibrium after a parameter change.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from cardiac_redundancy.source import SourceModel
from cardiac_redundancy.integration import integrate_source
from cardiac_redundancy.experiments import apply_knobs, warm_up, measure, sampled_checks

OUT = ROOT/'reports/m7-identifiability/pulsatile'
PARAMETERS = ('systemic_resistance', 'pulmonary_resistance', 'aortic_compliance', 'common_contractility')


def knobs(x):
    """x are log multipliers; common contractility scales all three active walls."""
    factors = np.exp(x)
    return dict(systemic_resistance_factor=factors[0], pulmonary_resistance_factor=factors[1],
                aortic_compliance_factor=factors[2], lv_elastance_factor=factors[3],
                rv_elastance_factor=factors[3], septal_elastance_factor=factors[3])


def save(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def run(case):
    name, x, refined = case
    folder = OUT/name
    folder.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    try:
        config = apply_knobs(json.loads((ROOT/'configs/baseline.yaml').read_text(encoding='utf-8-sig')), knobs(x))
        settings = config['solver_targets']
        settings.update(relative_tolerance=1e-9 if refined else 1e-7,
                        absolute_volume_tolerance_mL=1e-9 if refined else 1e-7,
                        absolute_flow_tolerance_mL_s=1e-8 if refined else 1e-6,
                        periodic_state_relative_tolerance=1e-7 if refined else 1e-6)
        save(folder/'config.json', config)
        model = SourceModel(folder/'config.json')
        seed = np.load(ROOT/'reports/m3-sweeps/baseline_refined/trajectory.npz')
        model.initial_state = seed['states'][:, -1].copy()
        solver = dict(method='Radau', max_step=.0025 if refined else .005)
        protocol = dict(warmup_phase_intervals=100, consecutive_periodic_comparisons=3)
        y, t, history = warm_up(model, protocol, solver)
        tr = integrate_source(model, config['driver']['period_s'], initial_state=y, start=t,
                              sample_step=.0005 if refined else .001, **solver)
        metrics, pressures = measure(model, tr, 1)
        checks = sampled_checks(model, tr, model.initial_state[:6].sum())
        if metrics['maximum_sampled_septal_residual_mmHg'] >= settings['septal_residual_tolerance_mmHg']:
            raise RuntimeError('Septal residual failed')
        np.savez_compressed(folder/'cycle.npz', phase_s=tr.time-tr.time[0], states=tr.states, pressures=pressures)
        result = dict(status='verified', name=name, log_multipliers=list(x), refined=refined,
                      metrics=metrics, checks=checks, warmup=history, solver=solver)
    except Exception as exc:
        result = dict(status='failed', name=name, error=repr(exc))
    result['elapsed_s'] = time.perf_counter()-started
    save(folder/'result.json', result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=6)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    cases = []
    for refined in (False, True):
        level = 'refined' if refined else 'nominal'
        cases.append((level+'_base', [0.]*4, refined))
        for step in ([.02, .01] if not refined else [.01]):
            for i in range(4):
                for sign in (-1, 1):
                    x = [0.]*4
                    x[i] = sign*step
                    cases.append((f'{level}_h{step:g}_p{i}_{sign:+d}', x, refined))
    paths = [Path(__file__), ROOT/'configs/baseline.yaml',
             ROOT/'reports/m3-sweeps/baseline_refined/trajectory.npz',
             *sorted((ROOT/'src/cardiac_redundancy').glob('*.py'))]
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    manifest = dict(parameters=PARAMETERS, cases=cases, hashes=hashes,
                    scope='Local log-parameter Jacobians; no patient calibration or confidence intervals')
    manifest_path = OUT/'input-manifest.json'
    if manifest_path.exists() and json.loads(manifest_path.read_text())['hashes'] != hashes:
        raise RuntimeError('Input hashes changed; use a separate experiment directory')
    save(manifest_path, manifest)
    pending = [c for c in cases if not (OUT/c[0]/'result.json').exists()]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(run, pending):
            print(result['name'], result['status'], result.get('error', ''), flush=True)


if __name__ == '__main__':
    main()
