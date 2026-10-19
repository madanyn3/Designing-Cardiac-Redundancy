"""Long-horizon verification of the analytic supported-duration bound.

The bound in duration.py assumes constant demand, hepatic export already at
its cap and a fully accessible glucose pool. This runs the optimum design over
a longer horizon and measures the pool's actual depletion rate, so the bound is
verified against a simulation rather than asserted.

A direct longer integration is used rather than tiling a settled cycle: the
circulation is periodic only to about 1 mL per cycle in compartment volume, and
tiling would inject that discontinuity into the substrate bookkeeping the check
depends on.

Verifying the drain *rate* over a minute is what licenses extrapolation to tens
of minutes. Simulating the full interval is not attempted, and the extrapolation
remains conditional on the rate staying constant as the pool falls.
"""
from pathlib import Path
import json, sys, time
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import run_m8_hydraulics as hyd
import run_m8_metabolism as met
from cardiac_redundancy.design import DesignPoint, ScenarioPoint
from cardiac_redundancy.duration import SubstrateBudget, supported_duration
from cardiac_redundancy.fuel import O2_ML_STPD_PER_MMOL_GLUCOSE

OUT = 'reports/m8-duration'
OPTIMUM = dict(left_scale=.651, right_scale=.651, septal_scale=.651, pericardial_scale=.651,
               period_s=1., coronary_conductance=.969)


def main(duration_s=60., kappa=.5):
    design = DesignPoint(**OPTIMUM)
    name = 'optimum-%gs' % duration_s
    folder = ROOT / OUT / name
    start = time.perf_counter()
    if (folder / 'hydraulics.json').is_file():
        record = json.loads((folder / 'hydraulics.json').read_text(encoding='utf-8'))
    else:
        record = hyd.evaluate(design, name, OUT, duration_s=duration_s, policy='fixed',
                              kappa=kappa, release_s=.5)
    hydraulic_s = time.perf_counter()-start
    if record['numerical_status'] != 'valid':
        print('hydraulic stage unusable:', record['numerical_status'], flush=True)
        return record
    start = time.perf_counter()
    if (folder / 'm8-coupled.json').is_file():
        coupled = json.loads((folder / 'm8-coupled.json').read_text(encoding='utf-8'))
    else:
        coupled = met.evaluate(design, ScenarioPoint(), name, Path(OUT) / name, out=OUT)
    metabolic_s = time.perf_counter()-start
    if coupled.get('status') != 'verified':
        print('coupled stage unusable:', coupled.get('status'), flush=True)
        return coupled
    a = np.load(ROOT / OUT / name / 'metabolic-refined.npz')
    t = a['time_s']
    demand = a['requested_O2_mL_s'].sum(axis=1)
    mean_demand = float(np.trapezoid(demand, t)/(t[-1]-t[0])*60)
    # Measured hepatic export and pool trajectory
    released = coupled['released_mmol']
    settled = t >= max(t[0]+5., t[-1]-30.)
    predicted = supported_duration(mean_demand)
    glucose_rate = mean_demand/O2_ML_STPD_PER_MMOL_GLUCOSE
    export_rate = released/((t[-1]-t[0])/60)
    result = dict(
        scope='long-horizon check of the analytic substrate bound at the T8.0 optimum',
        design=vars(design), duration_s=duration_s, compression_kappa=kappa,
        hydraulic_s=hydraulic_s, metabolic_s=metabolic_s,
        pressure_status=record['pressure']['status'],
        status_by_criterion=record['pressure']['status_by_criterion'],
        sustained_windows=record['pressure']['sustained_windows'],
        violating_sustained_windows=record['pressure']['violating_sustained_windows'],
        groups={k: dict(sustained_status=g['sustained_status'],
                        unmet_sustained_mL=g['unmet_sustained_mL']) for k, g in coupled['groups'].items()},
        mean_demand_mL_O2_min=mean_demand,
        measured_glucose_demand_mmol_min=glucose_rate,
        measured_export_mmol_min=export_rate,
        declared_export_cap_mmol_min=SubstrateBudget().release_cap_mmol_min,
        measured_net_drain_mmol_min=glucose_rate-export_rate,
        analytic=predicted,
        glucose_remaining_mmol=coupled.get('glucose_remaining_mmol'),
        glycogen_remaining_mmol=coupled['glycogen_remaining_mmol'],
        minimum_plasma_equivalent_mmol_L=coupled['minimum_plasma_equivalent_mmol_L'],
        interpretation='verifies the drain rate over the simulated interval; the extrapolated '
                       'duration remains conditional on that rate holding as the pool falls')
    if predicted['minutes'] is not None and export_rate > 0:
        measured_minutes = SubstrateBudget().circulating_pool_mmol/(glucose_rate-export_rate)
        result['extrapolated_minutes_from_measured_rate'] = measured_minutes
        result['analytic_error_fraction'] = (measured_minutes-predicted['minutes'])/predicted['minutes']
    dest = ROOT / OUT; dest.mkdir(parents=True, exist_ok=True)
    (dest / f'verification-{duration_s:g}s.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('duration %gs: pressure %s (%s), %d/%d sustained windows violated' % (
        duration_s, record['pressure']['status'], record['pressure']['status_by_criterion'],
        record['pressure']['violating_sustained_windows'], record['pressure']['sustained_windows']), flush=True)
    for role, g in sorted(result['groups'].items()):
        print('   %-12s %s  sustained unmet %.6f mL' % (role, g['sustained_status'],
                                                        g['unmet_sustained_mL']), flush=True)
    print('   demand %.2f mL O2/min -> glucose %.4f mmol/min; export measured %.4f (cap %.2f)' % (
        mean_demand, glucose_rate, export_rate, SubstrateBudget().release_cap_mmol_min), flush=True)
    print('   analytic %s min, extrapolated from measured rate %s min' % (
        'none' if predicted['minutes'] is None else '%.1f' % predicted['minutes'],
        '%.1f' % result['extrapolated_minutes_from_measured_rate']
        if 'extrapolated_minutes_from_measured_rate' in result else 'n/a'), flush=True)
    print('   wall: hydraulic %.1f min, metabolic %.1f min' % (hydraulic_s/60, metabolic_s/60), flush=True)
    return result


if __name__ == '__main__':
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 60.)
