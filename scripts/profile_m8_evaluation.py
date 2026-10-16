"""T8.1 cost measurement: wall time of one complete design evaluation.

Reuses the frozen M7 global runners unchanged so the measured cost is the cost
of the pipeline M8 will actually call, not of a simplified stand-in. Both
paired numerical settings are included, because the paired refinement check is
part of an accepted evaluation.

The measurement is machine- and load-specific. It is recorded as the basis for
the declared search budget, not as a portable performance claim.
"""
from pathlib import Path
import json, platform, sys, time
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import run_m7_global_hydraulics as hyd
import run_m7_global_metabolism as met

OUT = ROOT / 'reports/m8-protocol/timing'
NOMINAL = dict(extraction_limit=.5, mass_ratio=.25, work_cost_multiplier=1.,
               preservation_fraction=.25, preservation_conductance=.25,
               contractility=1., coronary_conductance=1.,
               glycogen_mmol_L_liver=229., release_cap_mmol_min=.21)


def main(name='reference'):
    hyd.OUT = OUT / 'hydraulics'; met.OUT = OUT
    case = dict(name=name, policy='observed', knobs=dict(NOMINAL),
                hydraulic_source=str((OUT / 'hydraulics' / name).relative_to(ROOT)))
    start = time.perf_counter(); hyd.run(case); hydraulic_s = time.perf_counter() - start
    start = time.perf_counter(); met.run(case); metabolic_s = time.perf_counter() - start
    record = dict(
        scope='single-evaluation wall time on the frozen M7 global pipeline',
        case=case, hydraulic_s=hydraulic_s, metabolic_s=metabolic_s,
        evaluation_s=hydraulic_s + metabolic_s,
        includes='paired nominal and refined settings for both stages',
        machine=dict(platform=platform.platform(), processor=platform.processor(),
                     python=platform.python_version()))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'single-evaluation.json').write_text(
        json.dumps(record, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(f'hydraulic {hydraulic_s:.1f} s, metabolic {metabolic_s:.1f} s, '
          f'evaluation {record["evaluation_s"]:.1f} s', flush=True)
    return record


if __name__ == '__main__':
    main(*sys.argv[1:])
