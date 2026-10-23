"""The claim audit, as a test rather than a one-off pass.

The preprint gate requires every numerical claim in the Results section to map
to a saved experiment or a citation. Three of the first forty-eight cited a
path that did not exist, so the audit is kept executable: it fails when a claim
loses its evidence, which is the failure mode that matters after text gets
edited.
"""
import csv
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / 'docs/preprint/results-claims.csv'
KINDS = {'model_output', 'derived', 'cited', 'assumption'}
#: Sources under the run cache are not committed, so a claim resting on one
#: must name the script that regenerates it or it is unreproducible.
CACHE_PREFIX = 'reports/m8-optimize/scale'

pytestmark = pytest.mark.skipif(not LEDGER.is_file(), reason='claim ledger absent')


@pytest.fixture(scope='module')
def claims():
    with LEDGER.open(encoding='utf-8') as f:
        return list(csv.DictReader(f))


def sources_of(row):
    return [p.strip().split(' ')[0] for p in row['source'].split(';') if p.strip()]


def test_ledger_is_well_formed(claims):
    assert claims, 'the ledger must carry claims'
    ids = [c['claim_id'] for c in claims]
    assert len(set(ids)) == len(ids), 'duplicate claim ids'
    for row in claims:
        assert row['kind'] in KINDS, '%s has kind %r' % (row['claim_id'], row['kind'])
        for field in ('claim_text', 'source', 'verification_status', 'qualification'):
            assert row[field].strip(), '%s has an empty %s' % (row['claim_id'], field)


def test_every_claim_cites_something_that_exists(claims):
    """A claim whose evidence has moved or been deleted is not supported."""
    missing = [(c['claim_id'], s) for c in claims for s in sources_of(c)
               if not (ROOT / s).exists()]
    assert not missing, 'claims citing paths that do not exist: %r' % missing


def test_cache_backed_claims_name_their_regenerating_script(claims):
    """Evidence that is not committed must be reproducible from the repo."""
    for row in claims:
        cached = [s for s in sources_of(row) if s.startswith(CACHE_PREFIX)]
        if not cached:
            continue
        named = [w.strip('.,;') for w in row['qualification'].split()
                 if w.strip('.,;').startswith('scripts/')]
        assert named, '%s rests on the run cache but names no script' % row['claim_id']
        for script in named:
            assert (ROOT / script).is_file(), \
                '%s names a missing script %s' % (row['claim_id'], script)


def test_model_outputs_are_not_presented_as_measurements(claims):
    """The study's hardest discipline: a model output is never a measurement.

    Any claim whose text reads as a measurement of a subject must be classified
    as model output or derived, and must carry a qualification saying so.
    """
    for row in claims:
        if 'measured' not in row['claim_text'].lower():
            continue
        assert row['kind'] in {'model_output', 'derived'}, \
            '%s reads as a measurement but is classified %r' % (row['claim_id'], row['kind'])
        assert row['qualification'].strip(), '%s needs an explicit qualification' % row['claim_id']


def test_fractions_are_never_called_probabilities(claims):
    """Scenario fractions are fractions of tested cases, by protocol."""
    for row in claims:
        text = (row['claim_text']+' '+row['qualification']).lower()
        if 'probabilit' in text:
            assert 'never' in text or 'not a probabilit' in text, \
                '%s uses probability language' % row['claim_id']
