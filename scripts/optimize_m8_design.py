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
from scipy.optimize import brentq
from cardiac_redundancy.design import DesignPoint, ScenarioPoint
from cardiac_redundancy.feasibility import pressure_margin

OUT = 'reports/m8-optimize'
PROTOCOL = json.loads((ROOT / 'configs/m8-protocol.json').read_text(encoding='utf-8'))
PERIODS = [.4, .5, .6, .75, 1., 1.25, 1.5]
#: Steps added above the pressure boundary in phase B. Raising coronary
#: conductance diverts flow from the aorta and costs 1-2 mmHg of pressure
#: (87.1 to 85.1 at scale 0.6), which a design sitting exactly on the boundary
#: cannot afford. The ladder buys that margin back with cavity scale.
SCALE_LADDER = [0., .03, .06, .10, .15]
EPSILONS = [1., 1.25, 1.5, 1.75, 2.]
SCALE_BOUNDS = tuple(PROTOCOL['declared_bounds']['left_scale'])
CORONARY_BOUNDS = tuple(PROTOCOL['declared_bounds']['coronary_conductance'])
SCALE_TOLERANCE = .01
CORONARY_TOLERANCE = .05


def design(scale, period, coronary=1.):
    return DesignPoint(left_scale=scale, right_scale=scale, septal_scale=scale,
                       pericardial_scale=scale, period_s=period,
                       coronary_conductance=coronary)


def name_of(d, kappa, release_s=.5, duration_s=None):
    base = ('scale%g-period%g-coron%g-kappa%g'
            % (d.left_scale, d.period_s, d.coronary_conductance, kappa))
    if release_s != .5: base += '-rel%g' % release_s
    # The horizon changes the trajectory, so it belongs in the identity. Runs
    # made at the original ten-second horizon keep their unsuffixed names.
    if duration_s is not None and abs(duration_s-10.) > 1e-9: base += '-h%g' % duration_s
    return base


def hydraulics(d, kappa, policy='fixed', release_s=.5, duration_s=None):
    """Evaluate, or reuse an identical evaluation saved earlier.

    The cache key must carry every setting that changes the trajectory, the
    release time included: keying on the design alone would silently return a
    0.5 s-release run for a 2.5 s-release request.
    """
    duration_s = duration_s or PROTOCOL['evaluation_policy']['horizon_s']
    name = name_of(d, kappa, release_s, duration_s)
    cached = ROOT / OUT / name / 'hydraulics.json'
    if cached.is_file():
        record = json.loads(cached.read_text(encoding='utf-8'))
        record['cached'] = True
        return record
    return hyd.evaluate(d, name, OUT, duration_s=duration_s, policy=policy, kappa=kappa,
                        release_s=release_s)


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


def margin_of(record):
    """Signed mmHg margin, or a large negative value for an unusable run.

    A numerical failure or domain exit is not a pressure result, so it is
    pushed below the boundary rather than being allowed to look like a near
    miss that a root-find could step toward.
    """
    if record['numerical_status'] != 'valid':
        return -1e3
    m = pressure_margin(record['pressure'])
    return -1e3 if m is None else m


def boundary_by_root_find(period, kappa, bracket, tolerance=SCALE_TOLERANCE):
    """Locate the pressure boundary in cavity scale by safeguarded root-find.

    Solves margin(scale) = 0 with Brent's method, which uses the magnitude of
    the margin instead of only its sign and so needs roughly half the
    evaluations a bisection on the boolean verdict does. Brent keeps the
    bracket, which matters because the margin is only piecewise smooth: the
    window and metric attaining the worst margin switch along the boundary.
    """
    trace = {}

    def f(scale):
        scale = round(float(scale), 3)
        if scale not in trace:
            record = hydraulics(design(scale, period), kappa)
            trace[scale] = (margin_of(record), summarise(record))
        return trace[scale][0]

    lo, hi = bracket
    if f(lo) > 0 or f(hi) < 0:
        return dict(period_s=period, rate_bpm=60/period, boundary=None,
                    reason='bracket does not straddle the boundary: margins %.3f at %.3f and '
                           '%.3f at %.3f' % (f(lo), lo, f(hi), hi),
                    trace=[v[1] for v in trace.values()])
    root = brentq(f, lo, hi, xtol=tolerance)
    # Confirm with the boolean verdict, and check the point below still fails.
    at = round(float(root)+tolerance, 3); below = round(float(root)-tolerance, 3)
    f(at); f(below)
    boundary = trace[at][1]; failing = trace[below][1]
    return dict(period_s=period, rate_bpm=60/period,
                boundary=boundary if boundary['pressure_status'] == 'satisfied' else None,
                failing_point_below=failing, root_scale=float(root),
                evaluations=len(trace),
                monotonicity=('checked: the point below the boundary fails'
                              if failing['pressure_status'] != 'satisfied'
                              else 'FAILED: the point below the boundary also passes'),
                margin_at_boundary_mmHg=trace[at][0],
                trace=[v[1] for v in trace.values()])


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


def pair_status(scale, period, coronary, kappa):
    """Both constraints at one (scale, conductance) point, plus the objective.

    Pressure and self-supply are reported separately because conductance moves
    them in opposite directions: it feeds the pump's own myocardium while
    diverting flow from the aorta that drives it.
    """
    d = design(scale, period, coronary)
    record = hydraulics(d, kappa)
    row = dict(scale=scale, period_s=period, rate_bpm=60/period, coronary_conductance=coronary,
               compression_kappa=kappa, pressure=record['pressure']['status'],
               numerical=record['numerical_status'],
               MAP_mmHg=record['runs'][-1]['final_MAP_mmHg'],
               chamber_work_J=sum(record['runs'][-1]['chamber_work_J']),
               auxiliary_demand_mL_min=objective(record))
    if not feasible(record):
        row.update(self_supply='not evaluated: pressure screen failed', feasible=False)
        return row
    name = name_of(d, kappa)
    cached = ROOT / OUT / name / 'm8-coupled.json'
    if cached.is_file():
        result = json.loads(cached.read_text(encoding='utf-8'))
    else:
        result = met.evaluate(d, ScenarioPoint(), name, Path(OUT) / name, out=OUT)
    if result.get('status') != 'verified':
        row.update(self_supply=result.get('status'), feasible=False)
        return row
    pump = result['groups']['pump']
    row.update(self_supply=pump['sustained_status'],
               sustained_unmet_mL=pump['unmet_sustained_mL'],
               startup_unmet_mL=pump['unmet_startup_mL'],
               feasible=pump['sustained_status'] == 'satisfied')
    return row


def least_coronary_conductance(scale, period, kappa):
    """Bisect conductance for the least value meeting self-supply at this scale.

    Self-supply is assumed monotone improving in conductance, which the
    endpoints test. A point that gains self-supply but loses the pressure
    screen is recorded as infeasible, not as a solution.
    """
    lo, hi = CORONARY_BOUNDS
    trace = []

    def check(value):
        row = pair_status(scale, period, value, kappa)
        trace.append(row)
        return row

    top = check(hi)
    if not top['feasible']:
        return dict(scale=scale, period_s=period, rate_bpm=60/period, compression_kappa=kappa,
                    minimum=None, best=None, trace=trace,
                    reason=('self-supply unmet at the upper conductance bound'
                            if top['self_supply'] != 'satisfied'
                            else 'pressure screen lost before self-supply was met'))
    bottom = check(lo)
    if bottom['feasible']:
        return dict(scale=scale, period_s=period, rate_bpm=60/period, compression_kappa=kappa,
                    minimum=lo, best=bottom, trace=trace,
                    note='the lower conductance bound already satisfies both constraints')
    best = top
    while hi-lo > CORONARY_TOLERANCE:
        mid = round((lo+hi)/2, 3)
        row = check(mid)
        if row['feasible']: hi, best = mid, row
        else: lo = mid
    return dict(scale=scale, period_s=period, rate_bpm=60/period, compression_kappa=kappa,
                minimum=hi, best=best, bracket=[lo, hi], trace=trace)


def _ladder_one(args):
    scale, period, kappa = args
    return least_coronary_conductance(scale, period, kappa)


def phase_b(kappa=.5, rates=None, workers=8, boundary_source=.5):
    """Joint (scale, conductance) search along the pressure boundary.

    For each rate, phase A's boundary scale is the cheapest scale that can
    satisfy pressure at nominal conductance. Self-supply may need more
    conductance than that scale can afford, so the ladder steps scale upward
    until a conductance exists that satisfies both, and the cheapest feasible
    pair over the whole grid is reported.
    """
    start = time.perf_counter()
    # Kappa moves MAP by under 1 mmHg (86.7/87.1/87.6 at kappa 0/0.5/1, scale
    # 0.6), so the pressure boundary traced at one kappa is reused as the
    # ladder's starting point at another. The ladder re-evaluates pressure at
    # every rung, so a reused boundary cannot smuggle in a pressure pass.
    source = ROOT / OUT / f'phase-a-kappa{kappa:g}.json'
    if not source.is_file():
        source = ROOT / OUT / f'phase-a-kappa{boundary_source:g}.json'
    phase = json.loads(source.read_text(encoding='utf-8'))
    lines = {l['rate_bpm']: l for l in phase['lines'] if l.get('boundary')}
    chosen = sorted(lines) if rates is None else [r for r in rates if r in lines]
    jobs = [(round(lines[r]['boundary']['scale']+step, 3), 60/r, kappa)
            for r in chosen for step in SCALE_LADDER
            if lines[r]['boundary']['scale']+step <= SCALE_BOUNDS[1]]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(_ladder_one, jobs))
    feasible_pairs = [r['best'] for r in results if r.get('best')]
    ranked = sorted(feasible_pairs, key=lambda r: r['auxiliary_demand_mL_min'])
    frontier = []
    for eps in EPSILONS:
        within = [r for r in feasible_pairs if r['coronary_conductance'] <= eps+1e-9]
        frontier.append(dict(max_conductance=eps,
            cheapest=min(within, key=lambda r: r['auxiliary_demand_mL_min']) if within else None))
    result = dict(phase='B', compression_kappa=kappa, protocol_version=PROTOCOL['version'],
        scope='joint cavity-scale and coronary-conductance search along the pressure boundary; '
              'reference physiology scenario, oracle activation, ten-second horizon',
        objective='total auxiliary oxygen demand (mL O2/min)',
        conductance_caveat='the 0.5-2.0 bracket is a declared assumption with no measured basis; '
                           'docs/uncertainty-plan.md states that using it as a design control '
                           'requires a verified interface, which does not exist',
        boundary_source=str(source.name),
        scale_ladder=SCALE_LADDER, conductance_tolerance=CORONARY_TOLERANCE,
        wall_s=time.perf_counter()-start, rates=chosen,
        best=ranked[0] if ranked else None, feasible_pairs=ranked,
        epsilon_frontier=frontier, lines=results)
    (ROOT / OUT / f'phase-b-kappa{kappa:g}.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return result


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
    phase = sys.argv[1] if len(sys.argv) > 1 else 'a'
    kappa = float(sys.argv[2]) if len(sys.argv) > 2 else .5
    if phase == 'b':
        rates = [float(x) for x in sys.argv[3:]] or None
        r = phase_b(kappa, rates)
        print(f"phase B done in {r['wall_s']/60:.1f} min, kappa {kappa:g}", flush=True)
        for row in r['feasible_pairs']:
            print('  rate %5.1f/min scale %.3f conductance %.3f  demand %7.3f mL O2/min  '
                  'MAP %6.2f  unmet %.6f' % (row['rate_bpm'], row['scale'],
                  row['coronary_conductance'], row['auxiliary_demand_mL_min'],
                  row['MAP_mmHg'], row.get('sustained_unmet_mL', 0.)), flush=True)
        for line in r['lines']:
            if not line.get('best'):
                print('  rate %5.1f/min scale %.3f infeasible: %s'
                      % (line['rate_bpm'], line['scale'], line.get('reason', '')), flush=True)
        for e in r['epsilon_frontier']:
            c = e['cheapest']
            print('  conductance <= %.2f : %s' % (e['max_conductance'],
                  'none feasible' if not c else 'scale %.3f at %.0f/min, %.3f mL O2/min'
                  % (c['scale'], c['rate_bpm'], c['auxiliary_demand_mL_min'])), flush=True)
        sys.exit(0)
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
