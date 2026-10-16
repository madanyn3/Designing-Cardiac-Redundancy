"""M8 Stage 2: the glucose and oxygen cost of a design that passed Stage 1.

Wraps the verified M7 coupled replay (`run_m7_global_metabolism.short`) without
changing it, and supplies the two T8.1 corrections it cannot make for itself:

  * auxiliary myocardial mass is charged from the chosen cavity scale through
    `design.demand_assumptions`, so a larger pump pays for its muscle;
  * the activation cost is charged at the chosen rate, `60/period_s`, instead
    of the literal 80 beats/min that was correct only while the period was
    fixed;
  * the auxiliary myocardial oxygen store is sized by the same five-second rule
    every other organ uses (`store_seconds=5` in `BilateralRespiration`) applied
    to this design's own demand, instead of the absolute pair frozen in
    configs/m5-myocardial-demand.json. That pair is exactly five seconds of a
    7.0/3.0 mL/min demand, so it is the same rule evaluated once at the
    reference fixture; held fixed, a larger pump carries a store sized for a
    smaller one, and the run aborts on "Oxygen capacity exceeded" as soon as
    supply exceeds demand.

Reported per organ role: requested oxygen, unmet ledger, and the unmet split
at the declared startup allowance. Substrate use remains instantaneous
(L-M4-F01) and contraction remains prescribed, so an unmet auxiliary ledger
ends the supported-performance claim rather than slowing the pump.
"""
from pathlib import Path
import json, sys, time
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import run_m7_global_metabolism as met
from cardiac_redundancy.design import DesignPoint, ScenarioPoint, contractility, demand_assumptions, mass_ratio
from cardiac_redundancy.fuel import FuelParameters
from cardiac_redundancy.myocardial_demand import side_demand

PROTOCOL = json.loads((ROOT / 'configs/m8-protocol.json').read_text(encoding='utf-8'))
ALLOWANCE = PROTOCOL['evaluation_policy']['startup_allowance_s']
#: Store horizon used for every organ in BilateralRespiration; reused here so
#: the auxiliary store follows the same rule rather than a frozen constant.
STORE_SECONDS = 5.


def read(path): return json.loads((ROOT / path).read_text(encoding='utf-8'))


def auxiliary_demand(design, scenario, hydraulics):
    """Per-side auxiliary oxygen demand at this design's size and rate."""
    reference = read('configs/m5-myocardial-demand.json')['reference']
    assumptions = demand_assumptions(design, scenario, reference)
    normal = read('reports/m5-normal-comparison/verification.json')['runs'][-1]['metrics']['cycle']
    native = read('data/m1-demand-budget.json')['baseline_mL_min']['heart']
    weights = read('configs/m4-organ-mapping.json')['mapping']['heart']
    work = hydraulics['runs'][-1]['positive_transmural_mean_W']
    return assumptions, {side: side_demand(
        native*weights['coronary_'+side], work[j],
        normal['sides'][side]['positive_transmural_mean_W'],
        design.rate_bpm, contractility(design, side), assumptions)
        for j, side in enumerate(('left', 'right'))}


def group_ledgers(result, source):
    """Per-role unmet oxygen, split at the declared startup allowance."""
    a = np.load(ROOT / source / 'metabolic-refined.npz')
    t = a['time_s']; unmet = a['oxygen_unmet_mL_s']; requested = a['requested_O2_mL_s']
    roles = {}
    for j, name in enumerate(result['organ_names']):
        roles.setdefault(result['profile'][name]['role'], []).append(j)
    early = t <= ALLOWANCE+1e-12; late = t >= ALLOWANCE
    out = {}
    for role, idx in roles.items():
        u = unmet[:, idx].sum(axis=1)
        bad = np.flatnonzero(u[late] > 1e-8)
        out[role] = dict(organs=[result['organ_names'][j] for j in idx],
            requested_mL=float(np.trapezoid(requested[:, idx].sum(axis=1), t)),
            unmet_total_mL=float(np.trapezoid(u, t)),
            unmet_startup_mL=float(np.trapezoid(u[early], t[early])) if early.sum() > 1 else 0.,
            unmet_sustained_mL=float(np.trapezoid(u[late], t[late])) if late.sum() > 1 else 0.,
            sustained_status='violated' if len(bad) else 'satisfied',
            first_sustained_bracket_s=None if not len(bad) else
                [float(t[late][max(0, bad[0]-1)]), float(t[late][bad[0]])],
            peak_unmet_mL_s=float(u.max()))
    return out


def evaluate(design, scenario, name, hydraulic_dir, out='reports/m8-coupled'):
    met.OUT = ROOT / out
    hydraulics = read(Path(hydraulic_dir) / 'hydraulics.json')
    if hydraulics['runs'][-1]['termination'] is not None:
        return dict(status='hydraulic_domain_limit', name=name)
    assumptions, sides = auxiliary_demand(design, scenario, hydraulics)
    case = dict(name=name, device='BIO', architecture='BIO', sides=sides,
                liver_function_fraction=scenario.liver_function_fraction,
                hydraulic_source=str(hydraulic_dir),
                knobs=dict(extraction_limit=scenario.extraction_limit,
                           preservation_fraction=scenario.preservation_fraction))
    pars = FuelParameters(**{**read('configs/m4-fuel-reference.json')['parameters'],
        'glycogen_mmol_L_liver': scenario.glycogen_mmol_L_liver,
        'release_cap_mmol_min': scenario.release_cap_mmol_min})
    # Size the auxiliary oxygen store by the shared five-second rule at this
    # design's demand. Patching the config read keeps the verified M7 replay
    # untouched while removing a constant frozen at the reference fixture.
    store = [sides[side]['total_mL_min']*STORE_SECONDS/60 for side in ('left', 'right')]
    original = met.read

    def coupled_read(path):
        data = original(path)
        if str(path).endswith('m5-myocardial-demand.json'):
            data = dict(data, auxiliary_store_capacity_mL=store)
        return data

    start = time.perf_counter()
    met.read = coupled_read
    try:
        result = met.short(case, pars)
    finally:
        met.read = original
    record = result['records'][-1]
    summary = dict(status='verified', name=name, design=vars(design), scenario=vars(scenario),
        runtime_s=time.perf_counter()-start,
        coupled_mass_ratio=mass_ratio(design, scenario), rate_bpm=design.rate_bpm,
        coupled_auxiliary_store_mL=store, store_seconds=STORE_SECONDS,
        auxiliary_demand_mL_min={s: v['total_mL_min'] for s, v in sides.items()},
        auxiliary_demand_components=sides,
        chamber_work_J=hydraulics['runs'][-1]['chamber_work_J'],
        groups=group_ledgers(result, Path(out) / name),
        unmet_oxygen_mL=record['unmet_oxygen_mL'],
        glycogen_remaining_mmol=record['glycogen_remaining_mmol'],
        released_mmol=record['released_mmol'],
        minimum_plasma_equivalent_mmol_L=record['minimum_plasma_equivalent_mmol_L'],
        first_auxiliary_shortfall_bracket_s=record['first_sampled_auxiliary_shortfall_bracket_s'],
        prescribed_driver_status=record['prescribed_driver_status'],
        scope='coupled replay over saved hydraulics; instantaneous substrate use and '
              'prescribed contraction remain limitations')
    dest = ROOT / out / name; dest.mkdir(parents=True, exist_ok=True)
    (dest / 'm8-coupled.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n',
                                          encoding='utf-8')
    return summary


if __name__ == '__main__':
    # usage: run_m8_metabolism.py <hydraulic-dir> [scenario]
    source = Path(sys.argv[1])
    scenario = ScenarioPoint(**PROTOCOL['scenarios'][sys.argv[2] if len(sys.argv) > 2 else 'reference'])
    design = DesignPoint(**read(source / 'hydraulics.json')['design'])
    name = source.name + ('' if len(sys.argv) <= 2 else '-'+sys.argv[2])
    r = evaluate(design, scenario, name, source)
    if r['status'] != 'verified':
        print(name, r['status'], flush=True); sys.exit(0)
    print(f"{name}: aux demand {sum(r['auxiliary_demand_mL_min'].values()):.3f} mL O2/min "
          f"(mass ratio {r['coupled_mass_ratio']:.3f}, {r['rate_bpm']:.1f} bpm), "
          f"unmet {r['unmet_oxygen_mL']:.6f} mL, glycogen left {r['glycogen_remaining_mmol']:.2f} mmol",
          flush=True)
    for role, g in sorted(r['groups'].items()):
        print(f"   {role:12s} requested {g['requested_mL']:8.3f} mL  unmet {g['unmet_total_mL']:.6f} "
              f"(startup {g['unmet_startup_mL']:.6f}, sustained {g['unmet_sustained_mL']:.6f}) "
              f"{g['sustained_status']}", flush=True)
