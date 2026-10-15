"""Organ uptake verification/sensitivity on saved M3 hydraulics, not survival."""
import copy
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
from cardiac_redundancy.oxygen import OxygenParameters
from cardiac_redundancy.organ_oxygen import OrganOxygenTransport, map_demands
from cardiac_redundancy.whole_network import WholeNetworkModel


def read(path):
    return json.loads((ROOT/path).read_text(encoding='utf-8'))


def main():
    config = read('configs/m4-organ-mapping.json'); gas = read('configs/m4-oxygen-reference.json')
    budget = read('data/m1-demand-budget.json')['baseline_mL_min']
    h = WholeNetworkModel(ROOT/'configs/baseline.yaml', ROOT/'configs/m3-evidence-regional.yaml',
                         ROOT/'configs/m2-frozen-overlay.yaml', arterial_fraction=.2, inlet_resistance_fraction=.1)
    cases = [(f'quarter-E{e}', 'preservation_quarter', e, 'reference') for e in config['extraction_limit_cases']]
    cases += [('open-E0.75', 'all_open', .75, 'reference')]
    cases += [(f'quarter-{variant}', 'preservation_quarter', .75, variant)
              for variant in ('gut_equal', 'residual_half', 'residual_double')]
    out = ROOT/'reports/m4-organ-uptake'; out.mkdir(parents=True, exist_ok=True)
    results = []; used = set(); saved = {}; comparisons = {}
    for label, hydraulic_case, extraction, variant in cases:
        mapping = copy.deepcopy(config['mapping'])
        if variant == 'gut_equal':
            mapping['gut'] = dict.fromkeys(mapping['gut'], 1/3)
        elif variant.startswith('residual_'):
            scale = .5 if variant.endswith('half') else 2.
            for organ in ('pancreas', 'spleen'):
                mapping['remaining'][organ] *= scale
            mapping['remaining']['remaining'] = 1-sum(mapping['remaining'][o] for o in ('pancreas', 'spleen'))
        organs = map_demands(budget, mapping, h.flow_names, extraction)
        m = OrganOxygenTransport(h.volume_names, h.edges, organs, OxygenParameters(**gas['parameters']))
        n = len(m.names); count = m.organ_count
        source = ROOT/f'reports/m3-whole-network/{hydraulic_case}-nominal'
        meta = json.loads((source/'result.json').read_text()); used.add(source/'result.json')
        assert meta['volume_names'] == list(h.volume_names) and meta['flow_names'] == list(h.flow_names)
        assert meta['status'] == 'completed'
        segments = []
        for i in range(4):
            path = source/f'segment-{i}.npz'; used.add(path)
            segments.append(dict(np.load(path)))
        for solver, rtol, step in [('nominal', 1e-7, .02), ('refined', 1e-9, .005)]:
            y = m.initial_state(segments[0]['states'][:n, 0], gas['initial_test_po2_mmHg'], gas['initial_alveolar_fraction'])
            initial = y.copy(); states = []; times = []; fluxes = []
            for segment in segments:
                t = segment['time']; volumes = interp1d(t, segment['states'][:n], axis=1, bounds_error=True)
                flows = interp1d(t, segment['flows'], axis=0, bounds_error=True)
                sol = solve_ivp(lambda time, state: m.rhs(volumes(time), state, flows(time), 4690/60),
                    (t[0], t[-1]), y, method='Radau', rtol=rtol, atol=rtol*.01, max_step=step, t_eval=t)
                if not sol.success: raise RuntimeError(sol.message)
                y = sol.y[:, -1].copy(); states.extend(sol.y.T); times.extend(t)
                for ti, yi in zip(t, sol.y.T):
                    _, d = m.evaluate(volumes(ti), yi, flows(ti), 4690/60)
                    fluxes.append(np.array([d[k] for k in ('delivered', 'extractable', 'consumed', 'demand', 'unmet')]))
            states = np.array(states); fluxes = np.array(fluxes); times = np.array(times)
            error = float(np.max(np.abs(m.balance(states)-m.balance(initial))))
            accounting = float(np.max(np.abs(states[:, n+2:n+2+count]+states[:, n+2+count:]
                -(times-times[0])[:, None]*np.array([o.demand_mL_min/60 for o in organs]))))
            if error > 1e-7 or accounting > 1e-7 or states[:, :n+1].min() < 0:
                raise RuntimeError('Store or uptake/deficit conservation failure')
            saved[label, solver] = states
            np.savez_compressed(out/f'{label}-{solver}.npz', time_s=times, oxygen_states=states, organ_fluxes=fluxes)
            result = dict(case=label, solver=solver, hydraulic_input=str(source.relative_to(ROOT)),
                rtol=rtol, atol=rtol*.01, max_step_s=step, method='Radau', duration_s=float(times[-1]-times[0]),
                extraction_limit=extraction, demand_mapping=mapping, organ_names=[o.name for o in organs],
                organ_demands_mL_min=[o.demand_mL_min for o in organs], total_demand_mL_min=sum(budget.values()),
                maximum_balance_error_mL=error, maximum_demand_accounting_error_mL=accounting,
                consumed_by_organ_mL=y[n+2:n+2+count].tolist(), unmet_by_organ_mL=y[n+2+count:].tolist(),
                total_consumed_mL=float(y[n+2:n+2+count].sum()), total_unmet_mL=float(y[n+2+count:].sum()),
                minimum_blood_or_lung_store_mL=float(states[:, :n+1].min()))
            results.append(result); print(label, solver, result['total_unmet_mL'], flush=True)
        comparisons[label] = float(np.max(np.abs(saved[label, 'nominal']-saved[label, 'refined'])))
    if max(comparisons.values()) > .001: raise RuntimeError('Refinement exceeds 0.001 mL')
    used.update([Path(__file__), ROOT/'configs/m4-organ-mapping.json', ROOT/'configs/m4-oxygen-reference.json',
                 ROOT/'data/m1-demand-budget.json'])
    used.update((ROOT/'src/cardiac_redundancy').glob('*.py'))
    used.update(ROOT/f'configs/{n}' for n in ('baseline.yaml', 'm3-evidence-regional.yaml', 'm2-frozen-overlay.yaml'))
    payload = dict(status='passed', scope='uptake verification and assumption sensitivity, not organ viability',
        volume_names=m.names, state_order='blood oxygen, lung oxygen, net external oxygen, organ cumulative consumption, organ cumulative unmet demand',
        flux_order=['delivered', 'extractable', 'consumed', 'demand', 'unmet'], flux_units='mL O2 STPD/s',
        cases=results, refinement_maximum_state_difference_mL=comparisons,
        software=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
        hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(used)})
    (out/'result.json').write_text(json.dumps(payload, indent=2)+'\n', encoding='utf-8')
    lines = ['# M4 organ uptake verification — 2026-10-07', '',
        '**Numerical checks pass. Organ extraction ceilings and fine demand partitions remain assumptions.**', '',
        'All runs retain full M1 resting comparison demand (255.847652 mL O2/min), including preservation tissues. The eight-second M3 trajectory arrests the native heart, applies support and reopens regional branches. Uniform initial blood PO2 40 mmHg is a test condition; this is not a normal metabolic baseline or an accepted arrested-heart demand estimate. Prescribed ventilation is 4.69 L/min BTPS.', '',
        '| Case | Consumed O2 (mL) | Unmet demand (mL) |', '|---|---:|---:|']
    for r in results:
        if r['solver'] == 'refined': lines.append(f"| {r['case']} | {r['total_consumed_mL']:.6f} | {r['total_unmet_mL']:.6f} |")
    lines += ['', f"Maximum oxygen balance error: {max(r['maximum_balance_error_mL'] for r in results):.3g} mL; maximum consumption+unmet-demand accounting error: {max(r['maximum_demand_accounting_error_mL'] for r in results):.3g} mL. Limits: 1e-7 mL each; saved physical stores nonnegative.", '',
        f'Maximum paired oxygen-state refinement difference: {max(comparisons.values()):.6g} mL (limit 0.001 mL). Same saved hydraulic sampling in each pair; sampling-resolution convergence is not established here.', '',
        'Deficits depend on assumed extraction limits, demand partitions, initial oxygenation and hydraulic controls. E=1 is only a conservation upper bound. Total delivered oxygen must not be summed over serial gut/liver pathways as independent external supply. Consumption and unmet-demand ledgers may be summed because demands are disjoint.', '',
        'No zero-deficit result establishes preserved function or recovery. No duration bound, safe branch closure or minimum pump is inferred. Hypothetical demand reductions, finite transitions, normal baseline, tissue stores, oxygenator and fuel coupling remain subsequent work.', '',
        'Reproduce: `.venv\\Scripts\\python.exe scripts/run_m4_organ_uptake.py`. [Machine-readable results and hashes](result.json); [mapping, equations and references](../../docs/m4-organ-uptake.md).']
    (out/'report.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')


if __name__ == '__main__': main()
