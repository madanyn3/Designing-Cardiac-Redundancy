"""Analyze local observation Jacobians; do not interpret SVD as global recovery."""
import json
import os
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'reports/m7-identifiability/pulsatile'
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'reports/m7-plot-cache'))


def observations(folder):
    result = json.loads((folder/'result.json').read_text())
    if result['status'] != 'verified':
        raise ValueError(f'Unverified case: {folder}')
    m = result['metrics']['mean']
    cycle = np.load(folder/'cycle.npz')
    t, p = cycle['phase_s'], cycle['pressures']
    mean_pa = np.trapezoid(p[:, 4], t)/(t[-1]-t[0])
    a = np.array([m['MAP_mmHg']/100, m['CO_L_min']/5])
    b = np.r_[a, m['SBP_mmHg']/100, m['DBP_mmHg']/100,
              mean_pa/20, m['VC_mean_mmHg']/5, m['PU_mean_mmHg']/10]
    # Fixed phase, omit duplicate cycle endpoint. Channel RMS weighting avoids
    # claiming more information merely by oversampling an identical waveform.
    phase = np.arange(50)*t[-1]/50
    channels = [p[:, 2]/100, p[:, 4]/20, cycle['states'][0]/100, cycle['states'][1]/100]
    c = np.concatenate([np.interp(phase, t, channel)/np.sqrt(50) for channel in channels])
    return {'A_MAP_CO': a, 'B_pressure_summaries': b, 'C_pressure_volume_waveforms': c}


def svd_audit(j):
    u, s, vh = np.linalg.svd(j, full_matrices=True)
    rank = int(np.linalg.matrix_rank(j))
    return dict(singular_values=s.tolist(), algebraic_numerical_rank=rank,
                nullity=j.shape[1]-rank,
                right_singular_vectors=vh.tolist(),
                condition_number=float(s[0]/s[-1]) if rank == j.shape[1] else None)


def main():
    bases = {level: observations(OUT/(level+'_base')) for level in ('nominal', 'refined')}
    jacobians = {}
    for level, step in [('nominal', .02), ('nominal', .01), ('refined', .01)]:
        columns = []
        for i in range(4):
            minus = observations(OUT/f'{level}_h{step:g}_p{i}_-1')
            plus = observations(OUT/f'{level}_h{step:g}_p{i}_+1')
            columns.append({key: (plus[key]-minus[key])/(2*step) for key in minus})
        jacobians[f'{level}_{step}'] = {key: np.column_stack([c[key] for c in columns]) for key in columns[0]}
    output = {'scope': 'Local native four-parameter conditional analysis; no nonlinear recovery claim',
              'parameter_order': ['systemic_R', 'pulmonary_R', 'aortic_C', 'common_active_elastance'],
              'scales': 'A: 100 mmHg, 5 L/min; B adds 100,100,20,5,10 mmHg; C: 100,20 mmHg,100,100 mL with channel RMS weighting',
              'observations': {}}
    lines = ['# T7.3 local pulsatile observation audit', '', output['scope'], '',
             '| Observations | Rank / 4 | Singular values | Step-change norm | Solver-change norm |',
             '|---|---:|---|---:|---:|']
    for key in bases['nominal']:
        j = jacobians['refined_0.01'][key]
        nominal = jacobians['nominal_0.01'][key]
        step_error = np.linalg.norm(nominal-jacobians['nominal_0.02'][key])/np.linalg.norm(j)
        solver_error = np.linalg.norm(nominal-j)/np.linalg.norm(j)
        audit = svd_audit(j)
        audit.update(jacobian=j.tolist(), step_relative_frobenius=float(step_error),
                     solver_relative_frobenius=float(solver_error),
                     step_difference_spectral_norm=float(np.linalg.norm(nominal-jacobians['nominal_0.02'][key], 2)),
                     solver_difference_spectral_norm=float(np.linalg.norm(nominal-j, 2)),
                     baseline_scaled_max_change=float(np.max(np.abs(bases['nominal'][key]-bases['refined'][key]))),
                     spectra_by_setting={label: svd_audit(values[key])['singular_values'] for label, values in jacobians.items()})
        # Linearized synthetic inversion tests the observation map, not the
        # nonlinear physiological simulator. Preserve null directions explicitly.
        truth = np.array([.01, -.008, .006, -.004])
        recovered = np.linalg.lstsq(j, j@truth, rcond=None)[0]
        audit['linearized_recovery'] = dict(truth=truth.tolist(), recovered=recovered.tolist(),
            max_parameter_error=float(np.max(np.abs(truth-recovered))),
            max_observation_residual=float(np.max(np.abs(j@(truth-recovered)))))
        output['observations'][key] = audit
        lines.append(f"| {key} | {audit['algebraic_numerical_rank']}/4 | {', '.join(f'{v:.6g}' for v in audit['singular_values'])} | {step_error:.3g} | {solver_error:.3g} |")
    lines += ['', 'Ranks use floating-point matrix rank; practical identifiability requires explicit error assumptions, nonlinear recovery and profiles. Singular vectors use the declared parameter/scaling conventions. C replaces A/B rather than treating summaries of the same waveforms as independent data.',
              '', 'Fixed: heart rate/phase, inventory, passive mechanics, pericardium, valve properties, other vascular compliances and external pressures. Common contractility changes LV/RV/septal active elastances together. No regional or auxiliary parameters are fitted.',
              '', 'Files: `summary.json`, `input-manifest.json`, and per-case configurations, periodic histories, checks and cycle arrays.']
    (OUT/'summary.json').write_text(json.dumps(output, indent=2)+'\n', encoding='utf-8')
    (OUT/'report.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for key, audit in output['observations'].items():
        values = audit['singular_values']
        ax.semilogy(np.arange(1, len(values)+1), values, 'o-', label=key)
    ax.set(xlabel='Singular direction (largest to smallest)', ylabel='Scaled singular value',
           title='Native local identifiability: four conditional parameters', xticks=[1, 2, 3, 4])
    ax.legend(fontsize=8)
    ax.text(.02, .03, 'A has two additional exact null directions.\nScales are not measurement-error estimates.', transform=ax.transAxes, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT/'singular-values.png', dpi=160)
    plt.close(fig)
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
