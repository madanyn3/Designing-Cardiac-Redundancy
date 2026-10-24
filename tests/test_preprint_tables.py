"""The preprint's tables must still match the artifacts they came from.

Tables are generated rather than transcribed, which only helps if the committed
file is actually what the generator produces. This regenerates them and
compares, so a hand edit or a stale table after an artifact changes fails here
rather than reaching the paper.
"""
import re
import subprocess
import sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / 'docs/preprint/tables.md'
RESULTS = ROOT / 'docs/preprint/results.md'

pytestmark = pytest.mark.skipif(not TABLES.is_file(), reason='tables have not been generated')


def test_tables_match_their_artifacts(tmp_path):
    """Regenerating must reproduce the committed file exactly."""
    before = TABLES.read_text(encoding='utf-8')
    result = subprocess.run([sys.executable, str(ROOT / 'scripts/build_preprint_tables.py')],
                            cwd=str(ROOT), capture_output=True, text=True)
    after = TABLES.read_text(encoding='utf-8')
    if before != after:
        TABLES.write_text(before, encoding='utf-8')  # leave the tree as found
        pytest.fail('tables.md differs from what the generator produces; regenerate it '
                    '(stdout: %s stderr: %s)' % (result.stdout[-300:], result.stderr[-300:]))


def test_every_table_results_cites_exists():
    """A cross-reference to a table that was never written is a dangling claim."""
    cited = {int(n) for n in re.findall(r'Table (\d+)', RESULTS.read_text(encoding='utf-8'))}
    present = {int(n) for n in re.findall(r'### Table (\d+)', TABLES.read_text(encoding='utf-8'))}
    missing = sorted(cited-present)
    assert not missing, 'Results cites tables that do not exist: %r' % missing


def test_no_table_cell_is_an_unfilled_gap():
    """The generator leaves gaps rather than inventing values; none may remain."""
    text = TABLES.read_text(encoding='utf-8')
    for marker in ('*artifact absent*', '*none found*', '*scenario curve artifact absent*'):
        assert marker not in text, 'tables carry an unfilled gap: %s' % marker


def test_tables_name_their_sources():
    """Each table states the artifact behind it."""
    for block in TABLES.read_text(encoding='utf-8').split('### Table ')[1:]:
        heading = block.splitlines()[0]
        assert 'Source' in block or 'source' in block, \
            'Table %s names no source' % heading
        for path in re.findall(r'`(reports/[^`]+\.json)`', block):
            assert (ROOT / path).is_file(), 'Table %s cites missing %s' % (heading, path)
