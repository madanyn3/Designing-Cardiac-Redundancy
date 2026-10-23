"""Collect the re-measured latency relation for margin-holding designs.

The original relation was measured at three designs, two of which hold +1.33
and +0.69 mmHg of pressure margin, below the 3.0 mmHg protocol v5 declares. A
slope read off a design sitting on the threshold describes the threshold
crossing rather than the design's latency tolerance, so the sweep was re-run
with the reportable designs added.

This reads the saved hydraulic solves and writes one artifact, so the reported
numbers cite a committed file rather than an uncommitted run cache. It solves
nothing itself: a release time with no saved solve is reported as missing.
"""
from pathlib import Path
import json, sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
import optimize_m8_design as opt
from cardiac_redundancy.feasibility import PRESSURE_SCREEN_MMHG, pressure_margin

OUT = 'reports/m8-latency'
PROTOCOL = json.loads((ROOT / 'configs/m8-protocol.json').read_text(encoding='utf-8'))
FLOOR = PROTOCOL['minimum_pressure_margin_mmHg']
RELEASES = (.5, 1., 1.5, 2., 2.5, 3.)
#: (scale, period_s, coronary_conductance, label)
DESIGNS = ((.536, 1.043, 1.748, 'margin-constrained optimum'),
           (.749, 1.304, .953, 'nominal-conductance design'),
           (.65, 1., 1., 'larger comparator'),
           (.475, 1., 1., 'original, below the declared margin'),
           (.387, .75, 1., 'original, below the declared margin'))


def crossing(points, level):
    """Release time at which the margin falls to `level`, by interpolation.

    Linear between the two solved points that bracket it. Returns None when no
    solved pair brackets the level, so an extrapolated figure is never reported
    as if it had been measured.
    """
    ordered = sorted(points)
    for (r0, m0), (r1, m1) in zip(ordered, ordered[1:]):
        if (m0-level) >= 0 > (m1-level):
            if m0 == m1:
                return r0
            return r0+(r1-r0)*(m0-level)/(m0-m1)
    return None


def main():
    horizon = PROTOCOL['evaluation_policy']['horizon_s']
    rows, designs = [], []
    for scale, period, coronary, label in DESIGNS:
        design = opt.design(scale, period, coronary)
        points, missing = [], []
        for release in RELEASES:
            folder = ROOT / opt.OUT / opt.name_of(design, .5, release, horizon)
            f = folder / 'hydraulics.json'
            if not f.is_file():
                missing.append(release); continue
            record = json.loads(f.read_text(encoding='utf-8'))
            margin = pressure_margin(record['pressure'])
            rows.append(dict(scale=design.left_scale, rate_bpm=design.rate_bpm,
                             coronary_conductance=design.coronary_conductance,
                             release_s=release, margin_mmHg=margin,
                             pressure_status=record['pressure']['status'],
                             MAP_mmHg=record['runs'][-1]['final_MAP_mmHg'],
                             source=str(f.relative_to(ROOT)).replace('\\', '/')))
            if margin is not None:
                points.append((release, margin))
        designs.append(dict(
            scale=design.left_scale, rate_bpm=design.rate_bpm,
            coronary_conductance=design.coronary_conductance, label=label,
            margin_at_reference_mmHg=dict(points).get(.5),
            holds_declared_margin=(dict(points).get(.5) or -1) >= FLOOR,
            tolerable_release_screen_s=crossing(points, 0.),
            tolerable_release_margin_s=crossing(points, FLOOR),
            solved_releases=[r for r, _ in sorted(points)], missing_releases=missing))
    result = dict(
        scope='activation-latency relation at designs holding the declared pressure margin',
        protocol_version=PROTOCOL['version'], minimum_pressure_margin_mmHg=FLOOR,
        enforced_threshold_mmHg=PRESSURE_SCREEN_MMHG['mean'], horizon_s=horizon,
        releases_s=list(RELEASES), designs=designs, rows=rows,
        interpretation='prescribed release stands for detection plus recruitment; crossings are '
                       'linear interpolations between solved points, never extrapolations, and '
                       'a design whose solved points do not bracket a level reports None',
        supersedes='reports/m8-optimize/latency.json, measured before a minimum margin was '
                   'declared, at designs holding +1.33 and +0.69 mmHg')
    dest = ROOT / OUT; dest.mkdir(parents=True, exist_ok=True)
    (dest / 'latency-margin.json').write_text(
        json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('latency relation at margin-holding designs (protocol v%d, floor %.1f mmHg)' % (
        PROTOCOL['version'], FLOOR))
    print('  %7s %7s %8s %9s %9s %9s' % ('scale', 'rate', 'coron', 'margin@.5', 'screen_s', 'margin_s'))
    for d in designs:
        print('  %7.3f %7.1f %8.3f %9s %9s %9s  %s' % (
            d['scale'], d['rate_bpm'], d['coronary_conductance'],
            '%+.2f' % d['margin_at_reference_mmHg'] if d['margin_at_reference_mmHg'] is not None else '-',
            '%.2f' % d['tolerable_release_screen_s'] if d['tolerable_release_screen_s'] else '-',
            '%.2f' % d['tolerable_release_margin_s'] if d['tolerable_release_margin_s'] else '-',
            d['label']))
        if d['missing_releases']:
            print('            missing solves at %s s' % d['missing_releases'])
    return result


if __name__ == '__main__':
    main()
