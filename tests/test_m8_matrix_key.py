"""The matrix cache key must carry the design it solved.

Re-pointing the matrix at a different design returned the previous design's
verdicts byte-identical, because the key was built from the case name and the
compression case alone. A cache that answers a different question is
indistinguishable from a result, so the key is tested rather than trusted.
"""
import json
import sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'scripts'))


def test_matrix_key_distinguishes_designs():
    """Two designs in the same case must not share a cache folder."""
    import run_m8_matrix as matrix
    from cardiac_redundancy.design import DesignPoint
    spec = dict(meaning='test', policy='fixed')
    seen = set()
    for optimum in (matrix.OPTIMUM, matrix.SUPERSEDED_OPTIMUM):
        design = DesignPoint(**optimum)
        name = (f'case-scale{design.left_scale:g}-period{design.period_s:g}'
                f'-coron{design.coronary_conductance:g}-kappa0.5')
        seen.add(name)
    assert len(seen) == 2, 'the two declared designs collide in the cache key'


def test_matrix_records_which_design_it_ran():
    """A saved matrix must say which design produced it."""
    out = ROOT / 'reports/m8-matrix/matrix-kappa0.5.json'
    if not out.is_file():
        pytest.skip('matrix has not been run')
    record = json.loads(out.read_text(encoding='utf-8'))
    assert 'design' in record, 'the matrix must record its design'
    import run_m8_matrix as matrix
    assert record['design']['left_scale'] == pytest.approx(matrix.OPTIMUM['left_scale']), \
        'the saved matrix was produced by a different design than the module declares'


def test_matrix_design_holds_the_declared_margin():
    """The reported matrix must stand on a reportable design.

    The earlier design held +2.26 mmHg, below the declared minimum, so the
    matrix described a design the protocol does not report.
    """
    import run_m8_matrix as matrix
    protocol = json.loads((ROOT / 'configs/m8-protocol.json').read_text(encoding='utf-8'))
    floor = protocol['minimum_pressure_margin_mmHg']
    from cardiac_redundancy.feasibility import pressure_margin
    import optimize_m8_design as opt
    design = opt.design(matrix.OPTIMUM['left_scale'], matrix.OPTIMUM['period_s'],
                        matrix.OPTIMUM['coronary_conductance'])
    folder = ROOT / opt.OUT / opt.name_of(design, .5, .5,
                                          protocol['evaluation_policy']['horizon_s'])
    f = folder / 'hydraulics.json'
    if not f.is_file():
        pytest.skip('the matrix design has no saved reference solve')
    margin = pressure_margin(json.loads(f.read_text(encoding='utf-8'))['pressure'])
    assert margin is not None and margin >= floor, \
        'matrix design holds %r mmHg, below the declared %.1f' % (margin, floor)


def test_matrix_rows_came_from_the_design_it_names():
    """The label is not evidence: check a row against its own saved solve.

    The first re-run wrote a matrix recording design 0.536 whose rows were the
    cached 0.496 solves, because the cache key omitted the design. Recording the
    design is therefore not sufficient — the rows must be shown to come from it.
    """
    out = ROOT / 'reports/m8-matrix/matrix-kappa0.5.json'
    if not out.is_file():
        pytest.skip('matrix has not been run')
    record = json.loads(out.read_text(encoding='utf-8'))
    design = record.get('design')
    if not design:
        pytest.skip('matrix records no design')
    checked = 0
    for row in record['rows']:
        name = ('%s-scale%g-period%g-coron%g-kappa%g' % (
            row['case'], design['left_scale'], design['period_s'],
            design['coronary_conductance'], row['compression_kappa']))
        f = ROOT / 'reports/m8-matrix' / name / 'hydraulics.json'
        if not f.is_file():
            continue
        solve = json.loads(f.read_text(encoding='utf-8'))
        assert solve['design']['left_scale'] == pytest.approx(design['left_scale']), \
            '%s: solve folder holds a different design than the matrix names' % row['case']
        assert row['MAP_mmHg'] == pytest.approx(solve['runs'][-1]['final_MAP_mmHg'], abs=1e-6), \
            ('%s: matrix reports MAP %.4f but its solve gives %.4f, so the row did not come '
             'from the design the matrix names' % (row['case'], row['MAP_mmHg'],
                                                   solve['runs'][-1]['final_MAP_mmHg']))
        checked += 1
    assert checked, 'no row could be checked against a saved solve'
