"""Bounded local nonlinear recovery using a frozen-Jacobian iteration.

Every residual is a fresh periodic nonlinear solve. The fixed Jacobian proposes
steps; it never substitutes for the forward model in acceptance decisions.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import numpy as np
import run_m7_pulsatile_identifiability as forward
from report_m7_pulsatile_identifiability import observations

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'reports/m7-identifiability/recovery'


def experiment(spec):
    forward.OUT = OUT/spec['name']
    forward.OUT.mkdir(parents=True, exist_ok=True)
    summary = json.loads((ROOT/'reports/m7-identifiability/pulsatile/summary.json').read_text())
    key = spec['observations']
    target = observations(ROOT/'reports/m7-identifiability/pulsatile/refined_base')[key]
    noise = np.zeros_like(target)
    if spec.get('noise'):
        noise = np.random.default_rng(7303).normal(0, .001, len(target))
    target = target + noise
    j = np.array(summary['observations'][key]['jacobian'])
    free = spec.get('free', [0,1,2,3])
    x = np.array(spec['start'], dtype=float)
    history = []
    for iteration in range(6):
        name = f'iteration-{iteration}'
        result = forward.run((name, x.tolist(), False))
        if result['status'] != 'verified':
            raise RuntimeError(result)
        residual = observations(forward.OUT/name)[key]-target
        step = np.linalg.lstsq(j[:, free], residual, rcond=None)[0]
        history.append(dict(x=x.tolist(), residual=residual.tolist(), step=step.tolist(),
                            residual_norm=float(np.linalg.norm(residual))))
        # Stop on a small parameter update, including nonzero noisy residuals.
        if np.max(np.abs(step)) < .0005:
            break
        if iteration == 5:
            break
        proposal = x.copy()
        proposal[free] -= np.clip(step, -.08, .08)
        if np.max(np.abs(proposal)) > .15:
            raise RuntimeError('Proposed parameters exceed registered log bounds')
        x = proposal
    validation = forward.run(('validation-refined', x.tolist(), True))
    if validation['status'] != 'verified':
        raise RuntimeError(validation)
    validation_residual = observations(forward.OUT/'validation-refined')[key]-target
    result = dict(spec=spec, final_log_multipliers=x.tolist(), multipliers=np.exp(x).tolist(),
                  converged_update=bool(np.max(np.abs(step)) < .0005), history=history,
                  refined_residual=validation_residual.tolist(), noise=noise.tolist(),
                  refined_residual_norm=float(np.linalg.norm(validation_residual)),
                  max_refined_residual=float(np.max(np.abs(validation_residual))))
    forward.save(forward.OUT/'recovery.json', result)
    return result


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    specs = []
    for key in ['B_pressure_summaries', 'C_pressure_volume_waveforms']:
        for sign in [-1, 1]:
            specs.append(dict(name=f'{key}_{sign:+d}', observations=key,
                              start=(sign*np.array([.05,-.04,.05,-.03])).tolist()))
    specs.append(dict(name='B_noisy', observations='B_pressure_summaries',
                      start=[.05,-.04,.05,-.03], noise=True))
    for sign in [-1, 1]:
        specs.append(dict(name=f'A_fixed_pulmonary_compliance_{sign:+d}', observations='A_MAP_CO',
                          start=[0.,sign*.08,sign*.08,0.], free=[0,3]))
    paths = [Path(__file__), ROOT/'scripts/run_m7_pulsatile_identifiability.py',
             ROOT/'scripts/report_m7_pulsatile_identifiability.py',
             ROOT/'reports/m7-identifiability/pulsatile/manifest.json']
    manifest = dict(specs=specs, log_bounds=[-.15,.15], max_iterations=6,
                    step_cap=.08, update_stop=.0005, synthetic_noise_scaled_sd=.001,
                    seed=7303, hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    path = OUT/'input-manifest.json'
    if path.exists() and json.loads(path.read_text()) != manifest:
        raise RuntimeError('Recovery protocol changed')
    forward.save(path, manifest)
    pending = [s for s in specs if not (OUT/s['name']/'recovery.json').exists()]
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(experiment, spec) for spec in pending]
        for future in as_completed(futures):
            result = future.result()
            print(result['spec']['name'], result['converged_update'], result['multipliers'],
                  result['refined_residual_norm'], flush=True)


if __name__ == '__main__':
    main()
