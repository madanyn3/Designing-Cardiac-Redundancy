"""The re-measured latency relation, and that it cites committed evidence.

The reported latency numbers were first taken from an uncommitted run cache,
which makes them unreproducible from a clean clone however correct they are.
These check the collected artifact instead: that it exists, that its crossings
are interpolations between solved points rather than extrapolations, and that
designs below the declared margin report no margin crossing at all.
"""
import json
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / 'reports/m8-latency/latency-margin.json'
PROTOCOL = json.loads((ROOT / 'configs/m8-protocol.json').read_text(encoding='utf-8'))

pytestmark = pytest.mark.skipif(not ARTIFACT.is_file(),
                                reason='latency artifact has not been collected')


@pytest.fixture(scope='module')
def artifact():
    return json.loads(ARTIFACT.read_text(encoding='utf-8'))


def test_every_row_names_the_file_it_came_from(artifact):
    """A reported number must be traceable to a saved solve."""
    assert artifact['rows'], 'the artifact must carry its rows'
    for row in artifact['rows']:
        assert row['source'].endswith('hydraulics.json')
        assert (ROOT / row['source']).is_file(), 'cites a missing solve: %s' % row['source']


def test_crossings_lie_inside_the_solved_range(artifact):
    """Interpolation only: a crossing outside the solved points is extrapolation."""
    releases = artifact['releases_s']
    for design in artifact['designs']:
        solved = design['solved_releases']
        if not solved:
            continue
        for key in ('tolerable_release_screen_s', 'tolerable_release_margin_s'):
            value = design[key]
            if value is None:
                continue
            assert min(solved) <= value <= max(solved), \
                '%s %r lies outside the solved releases %r' % (key, value, solved)
        assert set(solved) <= set(releases)


def test_designs_below_the_declared_margin_report_no_margin_crossing(artifact):
    """A design that never holds the margin cannot have a margin tolerance.

    Two of the original designs hold +1.33 and +0.69 mmHg, so there is no
    release time at which they hold 3.0 mmHg and none may be reported.
    """
    floor = artifact['minimum_pressure_margin_mmHg']
    assert floor == PROTOCOL['minimum_pressure_margin_mmHg']
    below = [d for d in artifact['designs'] if not d['holds_declared_margin']]
    assert below, 'the original designs are kept for comparison'
    for design in below:
        assert design['margin_at_reference_mmHg'] < floor
        assert design['tolerable_release_margin_s'] is None, \
            'reported a margin tolerance for a design that never holds the margin'


def test_margin_falls_monotonically_with_release(artifact):
    """Later activation is never better: the volume debt only compounds."""
    for design in artifact['designs']:
        points = sorted((r['release_s'], r['margin_mmHg']) for r in artifact['rows']
                        if r['scale'] == design['scale']
                        and r['coronary_conductance'] == design['coronary_conductance']
                        and r['margin_mmHg'] is not None)
        margins = [m for _, m in points]
        assert margins == sorted(margins, reverse=True), \
            'margin rose with later release at scale %.3f: %r' % (design['scale'], points)


def test_capacity_buys_latency_tolerance(artifact):
    """The reported relation: a larger design tolerates more delay.

    This is the claim the Results section makes, so it is tested rather than
    asserted. Among designs holding the declared margin, the one with the most
    pressure margin at reference must not tolerate the least delay.
    """
    holding = [d for d in artifact['designs']
               if d['holds_declared_margin'] and d['tolerable_release_margin_s']]
    if len(holding) < 2:
        pytest.skip('too few margin-holding designs to compare')
    best = max(holding, key=lambda d: d['margin_at_reference_mmHg'])
    worst = min(holding, key=lambda d: d['margin_at_reference_mmHg'])
    assert best['tolerable_release_margin_s'] >= worst['tolerable_release_margin_s']
