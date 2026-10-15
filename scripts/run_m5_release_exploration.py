"""Replay frozen M5 hydraulics with independently configurable glycogen capacity.

Outputs go to a new directory; source experiments and manifests are untouched.
This is metabolic sensitivity only: it does not add supply/contractility feedback.
"""
import argparse
from dataclasses import replace
import hashlib
import math
from pathlib import Path
import run_m5_bilateral_metabolism as replay

ROOT = Path(__file__).resolve().parents[1]

def settings(config):
    cfg = replay.read(config)
    cap = cfg['release_cap_mmol_min']
    fraction = cfg['liver_function_fraction']
    for value in (cap, fraction):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError('Release capacity and liver fraction must be finite numbers')
    if cap < 0 or not 0 <= fraction <= 1:
        raise ValueError('Nonnegative capacity and liver fraction in [0,1] required')
    source = cfg['source_case']
    allowed = {c['name'] for c in replay.read('configs/m5-final-comparison.json')['cases']}
    if source not in allowed:
        raise ValueError('Unknown frozen hydraulic source case')
    archived = replay.read(f'reports/m5-final-comparison/{source}/result.json')
    case = dict(archived['case'], liver_function_fraction=fraction)
    parameters = replace(replay.FuelParameters(**archived['parameters']), release_cap_mmol_min=cap)
    return cfg, case, parameters

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default='configs/m5-glycogen-release-exploration.json')
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    cfg, case, parameters = settings(args.config)
    print('Effective mobilisation cap:', parameters.release_cap_mmol_min * case['liver_function_fraction'], 'mmol/min')
    if args.validate_only:
        return
    # Content-address the settings so different knob values cannot overwrite
    # each other's results. Re-running identical settings is reproducible.
    config_path = ROOT / args.config
    tag = hashlib.sha256(config_path.read_bytes()).hexdigest()[:12]
    replay.OUT = ROOT / 'reports/m5-release-exploration' / tag
    result = replay.short(case, parameters)
    result['inputs'] += [str(config_path.relative_to(ROOT)).replace('\\', '/'),
                        f"reports/m5-final-comparison/{case['name']}/result.json"]
    result['release_exploration'] = cfg
    replay.write(replay.OUT / case['name'] / 'result.json', result)
    replay.write(replay.OUT / 'configuration.json', cfg)
    print('Verified paired metabolic replay:', replay.OUT)

if __name__ == '__main__':
    main()
