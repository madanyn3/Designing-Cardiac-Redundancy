"""Synthetic causal signal benchmark; no circulation/clinical outcome claims."""
from pathlib import Path
import sys
import json
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from cardiac_redundancy.controller_sensor import TrailingObservation
from cardiac_redundancy.support_controller import ControllerParameters, SupportController


def main():
    cfg = json.loads((ROOT/'configs/m6-controller.json').read_text())
    out = ROOT/'reports/m6-policy-benchmark'; out.mkdir(parents=True, exist_ok=True)
    results = []
    scenarios = ['normal', 'arrest', 'brief_dip', 'noise_normal', 'noise_arrest',
                 'stuck_high_pressure', 'sensor_loss', 'recovery', 'rearrest']
    for name in scenarios:
        rng = np.random.default_rng(6201)
        sensor = TrailingObservation(cfg['causal_filter_window_s'])
        controller = SupportController(ControllerParameters(**cfg['parameters']))
        rows = []
        for t in np.round(np.arange(-1, 10.0001, cfg['sample_period_s']), 8):
            arrest = name in ('arrest','noise_arrest','stuck_high_pressure','sensor_loss','recovery','rearrest') and t >= 0
            if name in ('recovery','rearrest') and t >= 3: arrest = False
            if name == 'rearrest' and t >= 6: arrest = True
            if name == 'brief_dip' and 0 <= t < .1: arrest = True
            pressure, flow = (50., 0.) if arrest else (85., 5.)
            # Pulsatility and additive noise are hypothetical test inputs.
            pressure += 10*np.sin(2*np.pi*t/.75)
            if name.startswith('noise_'):
                pressure += rng.normal(0, 3)
                flow += rng.normal(0, .1)
            if name == 'stuck_high_pressure' and t >= 0: pressure = 85.
            observation = sensor.update(float(t), pressure, flow)
            if t < 0: continue
            assert observation is not None
            valid = not (name == 'sensor_loss' and t >= 2)
            command = controller.update(float(t), *observation, sensor_valid=valid,
                adequate_for_recovery=name in ('recovery','rearrest') and not arrest)
            rows.append((float(t), *observation, command))
        release = next((e['time_s']+cfg['parameters']['signal_s']+cfg['parameters']['recruitment_s']
                        for e in controller.events if e['reason']=='confirmed'), None)
        # The fixed-delay comparator is an oracle schedule: it knows scenario
        # arrest onset and cannot be credited with detecting false alarms.
        oracle = .5 if name in ('arrest','noise_arrest','stuck_high_pressure','sensor_loss','recovery','rearrest') else None
        result = dict(case=name, first_scheduled_release_s=release,
            oracle_fixed_delay_release_s=oracle, final_state=controller.state,
            fault=controller.fault_reason, events=controller.events)
        results.append(result)
        np.savez_compressed(out/f'{name}.npz', columns=np.array(['time_s','observed_MAP_mmHg','observed_native_L_min','command']), values=np.array(rows))
    assert next(r for r in results if r['case']=='normal')['first_scheduled_release_s'] is None
    assert next(r for r in results if r['case']=='stuck_high_pressure')['first_scheduled_release_s'] is None
    assert next(r for r in results if r['case']=='sensor_loss')['final_state']=='FAULT'
    (out/'results.json').write_text(json.dumps(results, indent=2)+'\n')
    text = ['# M6 causal observation and policy benchmark', '',
        '2026-10-15. Synthetic input signals, not coupled circulation experiments. Seed 6201; 10 ms sampling, 0.75 s trailing integration window with explicit prehistory. Noise is assumed Gaussian SD 3 mmHg / 0.1 L/min. These are policy tests, not measured sensor performance.', '',
        '| Scenario | Observed-policy scheduled release s | Oracle release s | Final state |', '|---|---:|---:|---|']
    for r in results:
        text.append(f"| {r['case']} | {r['first_scheduled_release_s']} | {r['oracle_fixed_delay_release_s']} | {r['final_state']} |")
    text += ['', 'The fixed-delay reference knows arrest onset; its timing is not evidence of a better detector. Observed-policy latency includes the causal filter and confirmation. A stuck-high pressure signal causes a missed trigger because the policy requires both low pressure and low native flow. Signal invalidity latches FAULT and holds the previous command; it does not prove safe perfusion. Recovery adequacy is a synthetic input here, not a physiological result.', '',
        'Next: integrate the policy and continuous actuator envelope with the M5 bilateral circulation; compare pressure, substrate delivery and inventories under identical initial conditions. T6.2 remains open until those experiments are verified.']
    (out/'report.md').write_text('\n'.join(text)+'\n', encoding='utf-8')
    print('\n'.join(text))

if __name__ == '__main__': main()
