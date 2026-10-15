"""Verify finite demand schedules against saved M3 hydraulics and exact integrals."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/'src'))
from cardiac_redundancy.demand_schedule import DemandCommand, FiniteDemandSchedule, scenario_vector
from cardiac_redundancy.oxygen import OxygenParameters
from cardiac_redundancy.organ_oxygen import OrganOxygenTransport, map_demands
from cardiac_redundancy.whole_network import WholeNetworkModel


def read(p): return json.loads((ROOT/p).read_text(encoding='utf-8'))


def main():
    cfg = read('configs/m4-demand-transitions.json'); gas = read('configs/m4-oxygen-reference.json')
    budget = read('data/m1-demand-budget.json'); mapping = read('configs/m4-organ-mapping.json')['mapping']
    h = WholeNetworkModel(ROOT/'configs/baseline.yaml', ROOT/'configs/m3-evidence-regional.yaml',
        ROOT/'configs/m2-frozen-overlay.yaml', arterial_fraction=.2, inlet_resistance_fraction=.1)
    organs = map_demands(budget['baseline_mL_min'], mapping, h.flow_names, .75)
    m = OrganOxygenTransport(h.volume_names, h.edges, organs, OxygenParameters(**gas['parameters']))
    n = len(h.volume_names); count = len(organs)
    source = ROOT/'reports/m3-whole-network/preservation_quarter-nominal'
    meta = json.loads((source/'result.json').read_text()); assert meta['status'] == 'completed'
    assert meta['volume_names'] == list(h.volume_names) and meta['flow_names'] == list(h.flow_names)
    files = [source/f'segment-{i}.npz' for i in range(4)]
    segments = [dict(np.load(p)) for p in files]; arrest = float(segments[0]['time'][0])
    out = ROOT/'reports/m4-demand-transitions'; out.mkdir(parents=True, exist_ok=True)
    results = []; differences = {}; baseline = np.array([o.demand_mL_min for o in organs])
    for case in cfg['cases']:
        target = scenario_vector(organs, budget, case['scenario'])
        schedule = FiniteDemandSchedule(baseline, [DemandCommand(cfg['induction_command_s'], target, case['duration_s']),
            DemandCommand(case['return_command_s'], baseline, case['duration_s'])])
        curve_time = np.unique(np.r_[np.linspace(0, cfg['schedule_horizon_s'], 1261),
                                     schedule.breakpoints(0, cfg['schedule_horizon_s'])])
        np.savez_compressed(out/f"{case['name']}-schedule.npz", time_s=curve_time,
            demand_mL_min=np.array([schedule.at(t) for t in curve_time]),
            integrated_demand_mL=schedule.integral_mL(0, cfg['schedule_horizon_s']))
        saved = {}
        for label, tol, step in [('nominal', 1e-7, .02), ('refined', 1e-9, .005)]:
            y = m.initial_state(segments[0]['states'][:n, 0]); initial = y.copy(); times = []; states = []
            for segment in segments:
                ts = segment['time']; vf = interp1d(ts, segment['states'][:n], axis=1, bounds_error=True)
                qf = interp1d(ts, segment['flows'], axis=0, bounds_error=True)
                # Split at command/ramp boundaries without modifying the saved hydraulics.
                bounds = schedule.breakpoints(float(ts[0]-arrest), float(ts[-1]-arrest))
                for a, b in zip(bounds[:-1], bounds[1:]):
                    left, right = a+arrest, b+arrest
                    samples = np.unique(np.r_[left, ts[(ts > left) & (ts < right)], right])
                    sol = solve_ivp(lambda t, x: m.rhs(vf(t), x, qf(t), 4690/60, schedule.at(t-arrest)),
                        (left, right), y, method='Radau', rtol=tol, atol=tol*.01, max_step=step, t_eval=samples)
                    if not sol.success: raise RuntimeError(sol.message)
                    y = sol.y[:, -1].copy(); times.extend(sol.t-arrest); states.extend(sol.y.T)
            times = np.asarray(times); states = np.asarray(states)
            balance_error = float(np.max(np.abs(m.balance(states)-m.balance(initial))))
            actual = states[:, n+2:n+2+count]+states[:, n+2+count:]
            expected = np.array([schedule.integral_mL(0, t) for t in times])
            demand_error = float(np.max(np.abs(actual-expected)))
            if max(balance_error, demand_error) > 1e-7 or np.min(states[:, :n+1]) < 0:
                raise RuntimeError('Conservation/domain failure')
            np.savez_compressed(out/f"{case['name']}-{label}.npz", time_s=times, oxygen_states=states,
                                demand_mL_min=np.array([schedule.at(t) for t in times]))
            saved[label] = states
            r = dict(case=case, solver=label, method='Radau', rtol=tol, atol=tol*.01, max_step_s=step,
                maximum_balance_error_mL=balance_error, maximum_demand_integral_error_mL=demand_error,
                consumed_mL=float(y[n+2:n+2+count].sum()), unmet_mL=float(y[n+2+count:].sum()),
                demanded_mL=float(expected[-1].sum()), demand_at_8s_mL_min=float(schedule.at(8).sum()),
                target_native_demand_mL_min=float(target.sum()),
                consumed_by_organ_mL=y[n+2:n+2+count].tolist(), unmet_by_organ_mL=y[n+2+count:].tolist())
            results.append(r); print(case['name'], label, r['unmet_mL'], flush=True)
        differences[case['name']] = float(np.max(np.abs(saved['nominal']-saved['refined'])))
    if max(differences.values()) > .001: raise RuntimeError('Refinement failure')
    # No-suppression trajectory must also reproduce the previous constant-demand run.
    previous = np.load(ROOT/'reports/m4-organ-uptake/quarter-E0.75-refined.npz')
    current = np.load(out/'unsuppressed-refined.npz')
    unique, idx = np.unique(current['time_s'], return_index=True)
    old_time = previous['time_s']-previous['time_s'][0]
    interpolated = interp1d(unique, current['oxygen_states'][idx], axis=0)(old_time)
    regression = float(np.max(np.abs(interpolated-previous['oxygen_states'])))
    if regression > .001: raise RuntimeError('Constant-demand regression failure')
    inputs = [Path(__file__), *files, source/'result.json', ROOT/'reports/m4-organ-uptake/quarter-E0.75-refined.npz',
        *sorted((ROOT/'src/cardiac_redundancy').glob('*.py')),
        *[ROOT/p for p in ('configs/m4-demand-transitions.json','configs/m4-organ-mapping.json',
            'configs/m4-oxygen-reference.json','data/m1-demand-budget.json','configs/baseline.yaml',
            'configs/m3-evidence-regional.yaml','configs/m2-frozen-overlay.yaml')]]
    result = dict(status='passed', scope='prescribed demand transition verification; not functional recovery',
        organ_names=[o.name for o in organs], volume_names=h.volume_names, cases=results,
        refinement_maximum_state_difference_mL=differences, constant_demand_regression_mL=regression,
        software=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
        hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs})
    (out/'result.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    lines = ['# Finite demand-transition verification — 2026-10-07', '',
        '**Prescribed numerical scenarios, not demonstrated human suppression or recovery.**', '',
        'Fourteen runs use the same eight-second M3 quarter-preservation trajectory, extraction cap 0.75, uniform initial blood PO2 40 mmHg and prescribed ventilation 4.69 L/min BTPS. Induction starts at 0.5 s. Standard return command is 320 s, outside these hydraulic runs; the interruption case commands return at 3 s. Full schedule curves extend to 630 s but are NOT full-circulation simulations.', '',
        '| Case | Demand at 8 s (mL/min) | Integrated demand (mL) | Unmet demand (mL) |', '|---|---:|---:|---:|']
    for r in results:
        if r['solver'] == 'refined': lines.append(f"| {r['case']['name']} | {r['demand_at_8s_mL_min']:.4f} | {r['demanded_mL']:.6f} | {r['unmet_mL']:.6f} |")
    lines += ['', f"Maximum oxygen balance error: {max(r['maximum_balance_error_mL'] for r in results):.3g} mL; maximum error against exact time-integrated demand: {max(r['maximum_demand_integral_error_mL'] for r in results):.3g} mL (both limits 1e-7). Maximum refinement difference: {max(differences.values()):.6g} mL; constant-demand regression difference: {regression:.6g} mL (both limits 0.001).", '',
        'All saved physical oxygen stores remain nonnegative. No ramp reaches its suppression target within the eight-second hydraulic run. In particular, the 5% target is NOT achieved here. The unit suite separately checks complete 10/60/300 s induction/return schedules and a 620-second two-reservoir oxygen fixture; that fixture is not a prolonged M3 circulation.', '',
        'Deficits reflect these assumed schedules and test initial conditions. Demand restoration does not repair injury, erase prior deficits, restart native pumping or change vascular controls. Working auxiliary-heart demand remains independent and is tested by the API/unit suite, not a BIO run in this experiment.', '',
        'Reproduce: `.venv\\Scripts\\python.exe scripts/run_m4_demand_transitions.py`. [Results and provenance](result.json); [equations and evidence](../../docs/m4-demand-transitions.md).']
    (out/'report.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')


if __name__ == '__main__': main()
