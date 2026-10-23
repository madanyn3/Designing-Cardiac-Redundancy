"""Protocol v5: the declared pressure margin, and how designs get selected.

These guard the two defects v5 was written for. A design was reported that held
+0.005 mmHg of pressure margin, and a 0.04 per cent change in its activation
period flipped it to violated; the disconfirmation then tested the four
cheapest designs, all of them in that fragile group, and generalized from them.
Both are selection behaviours rather than numerical bugs, so they need tests
that fail if the selection reverts.
"""
import json
import sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))
from cardiac_redundancy.feasibility import PRESSURE_SCREEN_MMHG, pressure_margin

PROTOCOL = json.loads((ROOT / 'configs/m8-protocol.json').read_text(encoding='utf-8'))


def test_protocol_declares_a_positive_margin():
    """The margin must be declared, and declared above zero to mean anything."""
    assert PROTOCOL['version'] >= 5
    floor = PROTOCOL['minimum_pressure_margin_mmHg']
    assert floor > 0, 'a zero margin is the specification that selected +0.005 mmHg'
    assert floor < PRESSURE_SCREEN_MMHG['mean'], 'a margin at the threshold is not a margin'


def test_margin_amendment_records_reason_and_cost():
    """An amendment without its cost cannot be audited against alternatives."""
    entries = PROTOCOL['amendments_v5']
    assert entries, 'v5 must record why the margin was added'
    margin = [e for e in entries if 'margin' in e['change']]
    assert margin, 'the margin amendment must be recorded'
    for field in ('reason', 'reruns'):
        assert all(e.get(field) for e in entries), 'every amendment states %s' % field
    assert any('per cent' in e.get('cost', '') for e in margin), \
        'the margin amendment must state what robustness costs'


def test_pressure_margin_is_signed_against_the_threshold():
    """A satisfied screen has a positive margin and a violated one negative."""
    def endpoint(*means):
        return dict(thresholds_mmHg=dict(PRESSURE_SCREEN_MMHG), enforced_metrics=['mean'],
                    windows=[dict(interval='sustained', metrics_mmHg=dict(mean=m))
                             for m in means])
    threshold = PRESSURE_SCREEN_MMHG['mean']
    assert pressure_margin(dict(windows=[])) is None, 'no sustained window is not a verdict'
    assert pressure_margin(endpoint(70., 68.)) == pytest.approx(68.-threshold)
    assert pressure_margin(endpoint(70., 61.)) < 0, \
        'the worst window sets the margin, not the best'
    # Startup windows are excluded: a design is judged on what it sustains.
    mixed = endpoint(68.)
    mixed['windows'].append(dict(interval='startup', metrics_mmHg=dict(mean=20.)))
    assert pressure_margin(mixed) == pytest.approx(68.-threshold)


@pytest.mark.skipif(not (ROOT / 'reports/m8-margin/margin-candidates.json').is_file(),
                    reason='margin verification has not been run')
def test_leading_designs_hold_the_declared_margin():
    """Disconfirmation must not revert to selecting the cheapest designs."""
    import disconfirm_m8 as dc
    floor = PROTOCOL['minimum_pressure_margin_mmHg']
    designs = dc.leading_designs(4)
    assert designs, 'the selection must return something to stress'
    for row in designs:
        margin = row.get('margin_mmHg')
        assert margin is not None and margin >= floor, \
            'selected a design holding %r mmHg, below the declared %.1f' % (margin, floor)


@pytest.mark.skipif(not (ROOT / 'reports/m8-margin/margin-candidates.json').is_file(),
                    reason='margin verification has not been run')
def test_leading_designs_span_the_conductance_axis():
    """The cheapest corner is not a sample of a flat valley.

    Two conclusions were corrected because the disconfirmation tested four
    designs that shared a corner. The selection must spread over coronary
    conductance, which is the axis the cheap and the assumption-free corners
    differ on.
    """
    import disconfirm_m8 as dc
    designs = dc.leading_designs(4)
    if len(designs) < 2:
        pytest.skip('too few verified designs to span an axis')
    spread = (max(r['coronary_conductance'] for r in designs)
              - min(r['coronary_conductance'] for r in designs))
    assert spread > .2, 'selected designs share a conductance corner (spread %.3f)' % spread
