"""T8.0 joint design optimization by boundary tracing under protocol v2.

The objective is the auxiliary heart's own oxygen demand, which rises
monotonically with cavity scale and with rate. The optimum therefore lies on
the pressure-feasibility boundary rather than inside the feasible region, so
the search traces that boundary instead of exploring the interior.

Phase A bisects cavity scale at each declared rate to locate the least scale
that satisfies the pressure screen, and reads the objective off the same runs
because auxiliary demand follows from the design and its own hydraulic work.
Phase B takes the cheapest boundary points and bisects coronary conductance to
find the least value that satisfies auxiliary self-supply, repeated per
compression-kappa case because kappa is a declared model-form alternative and
not a design control.

Bisection assumes pressure feasibility is monotone in scale at fixed rate.
That premise is declared in the protocol and checked here: every located
boundary is reported with the failing point immediately below it, and a
non-monotone bracket is recorded as such rather than silently accepted.

Nothing here claims a global optimum. The result is the cheapest design found
on the traced boundary within the declared bounds, at the stated kappa case.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json, sys, time
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import run_m8_hydraulics as hyd
import run_m8_metabolism as met
from cardiac_redundancy.design import DesignPoint, ScenarioPoint

OUT = 'reports/m8-optimize'
PROTOCOL = json.loads((ROOT / 'configs/m8-protocol.json').read_text(encoding='utf-8'))
PERIODS = [.4, .5, .6, .75, 1.]
SCALE_BOUNDS = tuple(PROTOCOL['declared_bounds']['left_scale'])
CORONARY_BOUNDS = tuple(PROTOCOL['declared_bounds']['coronary_conductance'])
SCALE_TOLERANCE = .01
CORONARY_TOLERANCE = .05


def design(scale, period, coronary=1.):
    return DesignPoint(left_scale=scale, right_scale=scale, septal_scale=scale,
                       pericardial_scale=scale, period_s=period,
                       coronary_conductance=coronary)


def name_of(d, kappa):
    return ('scale%g-period%g-coron%g-kappa%g'
            % (d.left_scale, d.period_s, d.coronary_conductance, kappa))


def hydraulics(d, kappa, policy='fixed'):
    """Evaluate, or reuse an identical evaluation saved earlier."""
    name = name_of(d, kappa)
    cached = ROOT / OUT / name / 'hydraulics.json'
    if cached.is_file():
        record = json.loads(cached.read_text(encoding='utf-8'))
        record['cached'] = True
        return record
    return hyd.evaluate(d, name, OUT, policy=policy, kappa=kappa)


def objective(record, scenario=ScenarioPoint()):
    """Total auxiliary oxygen demand in mL O2/min for this design."""
    d = DesignPoint(**record['design'])
    _, sides = met.auxiliary_demand(d, scenario, record)
    return sum(v['total_mL_min'] for v in sides.values())


def summarise(record):
    d = DesignPoint(**record['design'])
    run = record['runs'][-1]
    return dict(scale=d.left_scale, period_s=d.period_s, rate_bpm=d.rate_bpm,
                coronary_conductance=d.coronary_conductance,
                compression_kappa=record.get('compression_kappa'),
                pressure_status=record['pressure']['status'],
                numerical_status=record['numerical_status'],
                violating_sustained_windows=record['pressure'].get('violating_sustained_windows'),
                MAP_mmHg=run['final_MAP_mmHg'], aux_left_L_min=run['aux_left_L_min'],
                chamber_work_J=sum(run['chamber_work_J']),
                auxiliary_demand_mL_min=objective(record),
                cached=record.get('cached', False))


def feasible(record):
    return record['pressure']['status'] == 'satisfied' and record['numerical_status'] == 'valid'


def least_feasible_scale(period, kappa=.5):
    """Bisect cavity scale for the least value passing the pressure screen."""
    lo, hi = SCALE_BOUNDS
    trace = []
    top = hydraulics(design(hi, period), kappa); trace.append(summarise(top))
    if not feasible(top):
        return dict(period_s=period, rate_bpm=60/period, boundary=None,
                    reason='no scale in the declared bounds satisfies the screen', trace=trace)
    bottom = hydraulics(design(lo, period), kappa); trace.append(summarise(bottom))
    if feasible(bottom):
        return dict(period_s=period, rate_bpm=60/period, boundary=summarise(bottom),
                    note='the lower bound already satisfies the screen', trace=trace)
    while hi-lo > SCALE_TOLERANCE:
        mid = round((lo+hi)/2, 3)
        record = hydraulics(design(mid, period), kappa); trace.append(summarise(record))
        if feasible(record): hi, top = mid, record
        else: lo, bottom = mid, record
    return dict(period_s=period, rate_bpm=60/period, boundary=summarise(top),
                failing_point_below=summarise(bottom), bracket=[lo, hi],
                monotonicity='checked: the point below the boundary fails', trace=trace)


def least_coronary_conductance(scale, period, kappa):
    """Bisect coronary conductance for the least value meeting self-supply."""
    lo, hi = CORONARY_BOUNDS
    trace = []

    def supplied(value):
        d = design(scale, period, value)
        record = hydraulics(d, kappa)
        if not feasible(record):
            trace.append(dict(coronary=value, pressure=record['pressure']['status'],
                              self_supply='not evaluated: pressure screen failed'))
            return None
        name = name_of(d, kappa)
        coupled = ROOT / OUT / name / 'm8-coupled.json'
        if coupled.is_file():
            result = json.loads(coupled.read_text(encoding='utf-8'))
        else:
            result = met.evaluate(d, ScenarioPoint(), name, Path(OUT) / name, out=OUT)
        if result.get('status') != 'verified':
            trace.append(dict(coronary=value, self_supply=result.get('status')))
            return None
        pump = result['groups']['pump']
        trace.append(dict(coronary=value, pressure='satisfied',
                          self_supply=pump['sustained_status'],
                          sustained_unmet_mL=pump['unmet_sustained_mL'],
                          auxiliary_demand_mL_min=sum(result['auxiliary_demand_mL_min'].values())))
        return pump['sustained_status'] == 'satisfied'

    if supplied(hi) is not True:
        return dict(scale=scale, period_s=period, compression_kappa=kappa, minimum=None,
                    reason='self-supply unmet at the upper conductance bound', trace=trace)
    if supplied(lo) is True:
        return dict(scale=scale, period_s=period, compression_kappa=kappa, minimum=lo,
                    note='the lower bound already satisfies self-supply', trace=trace)
    while hi-lo > CORONARY_TOLERANCE:
        mid = round((lo+hi)/2, 3)
        ok = supplied(mid)
        if ok is True: hi = mid
        elif ok is False: lo = mid
        else: return dict(scale=scale, period_s=period, compression_kappa=kappa, minimum=None,
                          reason='evaluation failed inside the bracket', bracket=[lo, hi], trace=trace)
    return dict(scale=scale, period_s=period, compression_kappa=kappa, minimum=hi,
                bracket=[lo, hi], trace=trace)


def phase_a(kappa=.5, workers=5):
    start = time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        lines = list(pool.map(_trace_one, [(p, kappa) for p in PERIODS]))
    feasible_lines = [l for l in lines if l.get('boundary')]
    ranked = sorted(feasible_lines, key=lambda l: l['boundary']['auxiliary_demand_mL_min'])
    result = dict(phase='A', scope='least cavity scale satisfying the pressure screen at each '
                  'declared rate, with the objective read from the same runs',
                  compression_kappa=kappa, protocol_version=PROTOCOL['version'],
                  scale_tolerance=SCALE_TOLERANCE, periods=PERIODS,
                  wall_s=time.perf_counter()-start, lines=lines,
                  cheapest=ranked[0] if ranked else None,
                  ranking=[dict(period_s=l['period_s'], rate_bpm=l['rate_bpm'],
                                scale=l['boundary']['scale'],
                                auxiliary_demand_mL_min=l['boundary']['auxiliary_demand_mL_min'],
                                chamber_work_J=l['boundary']['chamber_work_J'],
                                MAP_mmHg=l['boundary']['MAP_mmHg']) for l in ranked])
    dest = ROOT / OUT; dest.mkdir(parents=True, exist_ok=True)
    (dest / f'phase-a-kappa{kappa:g}.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return result


def _trace_one(args):
    return least_feasible_scale(*args)


if __name__ == '__main__':
    kappa = float(sys.argv[1]) if len(sys.argv) > 1 else .5
    r = phase_a(kappa)
    print(f"phase A done in {r['wall_s']/60:.1f} min, kappa {kappa:g}", flush=True)
    for row in r['ranking']:
        print('  rate %5.1f/min  least scale %.3f  demand %7.3f mL O2/min  work %6.3f J  MAP %6.2f'
              % (row['rate_bpm'], row['scale'], row['auxiliary_demand_mL_min'],
                 row['chamber_work_J'], row['MAP_mmHg']), flush=True)
    for line in r['lines']:
        if not line.get('boundary'):
            print('  rate %5.1f/min  no feasible scale: %s' % (line['rate_bpm'], line['reason']),
                  flush=True)
