"""M8 design-aware hydraulic evaluation under the frozen T8.1 protocol.

Derived from scripts/run_m7_global_hydraulics.py, which is left unchanged. The
differences are the ones T8.1 requires: a DesignPoint supplies the cavity
scales, period and contractilities, the saved reference state is reallocated to
the chosen cavity size before integration, and the verdict comes from the
frozen pressure evaluator rather than a reported mean.

This stage decides the pressure screen on its own because metabolism is
replayed over saved hydraulics and never feeds back. It does not decide
oxygen or substrate adequacy, and a pressure pass here is not feasibility.
"""
from pathlib import Path
import json, sys
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
from cardiac_redundancy.sensitivity_circulation import SensitivityCirculation
from cardiac_redundancy.circuit_integration import integrate_circuit
from cardiac_redundancy.controller_sensor import TrailingObservation
from cardiac_redundancy.support_controller import SupportController, ControllerParameters
from cardiac_redundancy.biological_energy import biological_power
from cardiac_redundancy.pump_energy import MMHG_ML_TO_J
from cardiac_redundancy.design import DesignPoint
from cardiac_redundancy.feasibility import FeasibilityPolicy, numerical_status, pressure_endpoint
from screen_m8_domain_bounds import rescaled_initial_state

PROTOCOL = json.loads((ROOT / 'configs/m8-protocol.json').read_text(encoding='utf-8'))
STATE = ROOT / 'reports/m5-joint-preservation/reference_g25/fine-refined.npz'


def read(path): return json.loads((ROOT / path).read_text(encoding='utf-8'))


def build(design):
    return SensitivityCirculation(
        ROOT / 'configs/baseline.yaml', ROOT / 'configs/m3-evidence-regional.yaml',
        ROOT / 'configs/m2-frozen-overlay.yaml', arterial_fraction=.2,
        inlet_resistance_fraction=.1, coronary_conductance=design.coronary_conductance,
        biological=design.biological())


def integrate(design, label, rtol, max_step, duration_s, policy='observed', kappa=None):
    """One hydraulic trajectory at one numerical setting.

    `kappa` overrides the coronary compression factor declared in
    configs/m2-frozen-overlay.yaml (`compression_kappa`, 0.5, with declared
    cases 0/0.5/1 that no runner has ever exercised). It is a model-form
    scenario, not a design control: it changes how much of the chamber's
    transmural pressure closes its own coronary bed.
    """
    cfg = read('configs/m6-coupled-protocol.json'); cc = read('configs/m6-controller.json')
    m = build(design)
    if kappa is not None:
        m.native.base.closure_settings['kappa'] = float(kappa)
    y = rescaled_initial_state(m, np.load(STATE)['states'][:, 0])
    m.apply_support_controls(0, design.preservation_conductance)
    if min(m.domain_margins(y).values()) <= 0:
        raise RuntimeError('Initial state outside the model domain')
    controller = SupportController(ControllerParameters(**cc['parameters']))
    sensor = TrailingObservation(cfg['filter_window_s'])
    initial_pressure = m.pressures(0, y)[0]['AO']
    for t in np.round(np.arange(-1, 0, cfg['sample_period_s']), 8):
        sensor.update(float(t), initial_pressure, 5.213513)
    ts = []; ys = []; ps = []; qs = []; amps = []; powers = []
    termination = None; root = 0.; jumps = 0.
    release = .5 if policy == 'fixed' else None
    if release is not None: m.aux_envelope.set(release, 1, 1)
    ticks = np.round(np.arange(0, duration_s+.0001, cfg['sample_period_s']), 8)
    for tick, end in zip(ticks[:-1], ticks[1:]):
        tick = float(tick); end = float(end)
        p, _ = m.pressures(tick, y)
        obs = sensor.update(tick, p['AO'], y[19]*.06)
        if policy == 'observed':
            before = len(controller.events)
            controller.update(tick, *obs, sensor_valid=True,
                              model_valid=min(m.domain_margins(y).values()) > 0,
                              adequate_for_recovery=False)
            for event in controller.events[before:]:
                if event['after'] == 'ACTIVATE':
                    release = controller.release_s
                    m.aux_envelope.set(controller.ramp_start, 1, controller.p.ramp_s)
        bounds = sorted(set([tick, end]+([release] if release is not None and tick <= release <= end else [])))
        for start, stop in zip(bounds[:-1], bounds[1:]):
            if release is not None and m.trigger is None and abs(start-release) < 1e-9:
                m.release_auxiliary(start)
            tr = integrate_circuit(m, stop-start, start=start, initial_state=y,
                                   rtol=rtol, max_step=max_step, sample_step=.005)
            jumps = max(jumps, float(np.max(abs(tr.states[:, 0]-y))))
            for ti, yi in zip(tr.time, tr.states.T):
                if ts and ti <= ts[-1]+1e-12: continue
                p, meta = m.pressures(ti, yi); q = m.flows(ti, yi, [True]*8, pressures=p)
                power = biological_power(m, ti, yi)
                ts.append(ti); ys.append(yi); ps.append([p[n] for n in m.volume_names])
                qs.append([q[n] for n in m.flow_names]); amps.append(m.aux_envelope.value(ti))
                powers.append([[power[s][k] for k in ('chamber_W', 'transmural_W', 'balance_residual_W')]
                               for s in ('left', 'right')])
                root = max(root, abs(meta['septal_residual_mmHg']), abs(meta['aux_septal_residual_mmHg']))
            y = tr.states[:, -1].copy()
            if tr.termination: termination = tr.termination; break
        if termination: break
    t = np.array(ts); v = np.array(ys).T; p = np.array(ps); q = np.array(qs); power = np.array(powers)
    limits = PROTOCOL['numerical']['hydraulic']['conservation_limits']
    drift = float(np.max(abs(v[:18].sum(axis=0)-5250)))
    checks = dict(blood_error_mL=drift, state_jump_mL=jumps, root_residual_mmHg=root)
    work = np.trapezoid(power[:, :, 0], t, axis=0).tolist()
    pv = [float(-MMHG_ML_TO_J*np.sum(.5*(p[1:, j]+p[:-1, j])*np.diff(v[j]))) for j in (16, 17)]
    checks['PV_quadrature_error_J'] = max(abs(a-b) for a, b in zip(work, pv))
    failed = {k: val for k, val in checks.items() if val > limits[k]}
    if abs(power[:, :, 2]).max() > 1e-10: failed['balance_residual_W'] = float(abs(power[:, :, 2]).max())
    mask = t >= max(0, t[-1]-2); win = t[mask]
    mean = lambda z: float(np.trapezoid(z[mask], win)/(win[-1]-win[0])) if len(win) > 1 else None
    return dict(label=label, rtol=rtol, max_step_s=max_step, termination=termination,
                release_s=m.trigger, controller_events=controller.events, checks=checks,
                failed_checks=failed, final_MAP_mmHg=mean(p[:, 2]),
                aux_left_L_min=mean(q[:, m.flow_names.index('aAV')]*.06),
                chamber_work_J=work, volume_names=m.volume_names,
                positive_transmural_mean_W=[float(w/(t[-1]-t[0])) for w in
                    np.trapezoid(np.maximum(power[:, :, 1], 0), t, axis=0)]), (t, p[:, 2], v, p, q, np.array(amps))


def evaluate(design, name, out, duration_s=None, policy='observed', save=True, kappa=None):
    """Paired hydraulic evaluation plus the frozen pressure verdict."""
    duration_s = duration_s or PROTOCOL['evaluation_policy']['horizon_s']
    feasibility = FeasibilityPolicy(horizon_s=duration_s)
    runs = []; arrays = []
    for pair in PROTOCOL['numerical']['hydraulic']['pairs']:
        record, data = integrate(design, pair['label'], pair['rtol'], pair['max_step_s'],
                                 duration_s, policy, kappa)
        runs.append(record); arrays.append(data)
    t, ao, v, p, q, amps = arrays[-1]
    tb = arrays[0][0]
    errors = []
    for x, z in ((arrays[0][2].T, v.T), (arrays[0][3], p), (arrays[0][4], q)):
        interpolated = np.array([np.interp(tb, t, col) for col in z.T]).T
        errors.append(float(np.max(abs(x-interpolated)/np.maximum(1, np.max(abs(x), axis=0)))))
    numerical = numerical_status(runs[-1]['termination'], errors,
                                 PROTOCOL['numerical']['hydraulic']['paired_scaled_error_limit'])
    if any(r['failed_checks'] for r in runs): numerical = 'numerical_failure'
    pressure = pressure_endpoint(feasibility, t, ao, release_s=runs[-1]['release_s'])
    result = dict(status='screened', name=name, design=vars(design), policy=policy,
                  compression_kappa=0.5 if kappa is None else float(kappa),
                  duration_s=duration_s, numerical_status=numerical, scaled_errors=errors,
                  pressure=pressure, runs=runs,
                  scope='hydraulic stage only; decides the pressure screen, not feasibility')
    if save:
        dest = ROOT / out / name; dest.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(dest/'hydraulic-refined.npz', time_s=t, states=v,
                            pressures_mmHg=p, flows_mL_s=q, aux_amplitude=amps)
        (dest/'hydraulics.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n',
                                            encoding='utf-8')
    return result


if __name__ == '__main__':
    # usage: run_m8_hydraulics.py [fixed|observed] [out=DIR] scale[:period]...
    args = sys.argv[1:]
    policy = 'observed'
    if args and args[0] in ('fixed', 'observed'): policy, args = args[0], args[1:]
    out = 'reports/m8-capacity' if policy == 'observed' else 'reports/m8-capacity-oracle'
    if args and args[0].startswith('out='): out, args = 'reports/'+args[0][4:], args[1:]
    # spec = scale[:period[:preservation_conductance[:contractility[:coronary]]]]
    defaults = (.25, .75, .25, 1., 1., .5)
    specs = [tuple([*(float(v) for v in spec.split(':')), *defaults[len(spec.split(':')):]])
             for spec in (args or ['0.25', '0.5', '0.75', '1.0'])]
    summary = []
    for scale, period, preservation, contract, coronary, kappa in specs:
        design = DesignPoint(left_scale=scale, right_scale=scale, septal_scale=scale,
                             pericardial_scale=scale, period_s=period,
                             preservation_conductance=preservation,
                             left_contractility=contract, right_contractility=contract,
                             septal_contractility=contract, coronary_conductance=coronary)
        name = 'scale%g-period%g-cond%g-contract%g-coron%g-kappa%g' % (
            scale, period, preservation, contract, coronary, kappa)
        r = evaluate(design, name, out, policy=policy, kappa=kappa)
        w = r['pressure']['sustained_exposure']
        summary.append(dict(scale=scale, period_s=period, rate_bpm=design.rate_bpm,
                            preservation_conductance=preservation, contractility=contract,
                            coronary_conductance=coronary, compression_kappa=kappa,
                            status=r['pressure']['status'],
                            violating_sustained_windows=r['pressure']['violating_sustained_windows'],
                            first_satisfied_sustained_window_s=r['pressure']['first_satisfied_sustained_window_s'],
                            numerical_status=r['numerical_status'],
                            MAP_mmHg=r['runs'][-1]['final_MAP_mmHg'],
                            aux_left_L_min=r['runs'][-1]['aux_left_L_min'],
                            release_s=r['runs'][-1]['release_s'],
                            chamber_work_J=r['runs'][-1]['chamber_work_J'],
                            violating_window_s={k: w[k]['violating_window_duration_s'] for k in w},
                            shortfall_mmHg_s={k: w[k]['shortfall_mmHg_s'] for k in w}))
        print(f"scale {scale:g} period {period:g} cond {preservation:g} "
              f"contract {contract:g} coron {coronary:g}: "
              f"{summary[-1]['status']} "
              f"({summary[-1]['violating_sustained_windows']} bad windows) "
              f"MAP {summary[-1]['MAP_mmHg']:.4f} "
              f"aux {summary[-1]['aux_left_L_min']:.4f} L/min", flush=True)
    dest = ROOT/out; dest.mkdir(parents=True, exist_ok=True)
    (dest/'summary.json').write_text(json.dumps(dict(
        scope='equal-cavity-scale pressure screen under the frozen T8.1 protocol',
        activation_policy=policy,
        protocol_version=PROTOCOL['version'], thresholds_mmHg=PROTOCOL['evaluation_policy']['thresholds_mmHg'],
        records=summary), indent=2, allow_nan=False)+'\n', encoding='utf-8')
