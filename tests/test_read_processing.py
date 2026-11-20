"""Read-processing reports distinguish execution evidence from requested settings."""
import json

import pytest

from common import write_tsv
from read_processing import collect, summary
from sample_processing import native_processing


def read_evidence(output, native, run, source, counts=None, **extra):
    directory = output / 'amalgkit_getfastq' / native
    directory.mkdir(parents=True, exist_ok=True)
    manifest = directory / 'getfastq_completion.json'
    manifest.write_text(json.dumps({'schema_version': 4, 'status': 'complete',
                                   'read_source': source, 'runs': [{'run': run, 'status': 'complete'}],
                                   **extra}))
    if counts is not None:
        row = {'run': run, **counts}
        write_tsv(directory / run / 'getfastq_stats.tsv', list(row), [row])
    return manifest


@pytest.mark.parametrize('source,counts,expected', [
    ('amalgkit', {'num_fastp_in': 100, 'num_written': 100}, ('amalgkit', 'yes')),
    ('amalgkit', {'num_fastp_in': 0, 'num_written': 100}, ('amalgkit', 'no')),
    ('amalgkit', None, ('amalgkit', 'unknown')),
    ('amalgkit', {'num_fastp_in': 'nan', 'num_written': 100}, ('amalgkit', 'unknown')),
    ('public-original', {'num_written': 100, 'bp_fastp_in': 10000, 'bp_fastp_out': 10000},
     ('genegalleon_fallback', 'no')),
    ('public-original', {'num_fastp_in': 100, 'num_written': 100}, ('genegalleon_fallback', 'yes')),
])
def test_collect_uses_fastp_read_counts_not_fallback_base_placeholders(tmp_path, source, counts, expected):
    read_evidence(tmp_path, 'Plant', 'SRR1', source, counts)
    assert collect(tmp_path, 'Plant', 'SRR1') == dict(zip(
        ('read_acquisition_method', 'fastp_applied'), expected))


@pytest.mark.parametrize('applied,expected', [(True, 'yes'), (False, 'no'), ('unknown', 'unknown'), ('partial', 'unknown')])
def test_explicit_preprocessing_status_is_independent_of_fallback(tmp_path, applied, expected):
    read_evidence(tmp_path, 'Plant', 'SRR1', 'public-original', fastp_applied=applied)
    assert collect(tmp_path, 'Plant', 'SRR1')['fastp_applied'] == expected


@pytest.mark.parametrize('change', ['missing', 'corrupt', 'unfinished', 'wrong_run', 'duplicate_run'])
def test_missing_or_ineligible_evidence_stays_unknown(tmp_path, change):
    path = read_evidence(tmp_path, 'Plant', 'SRR1', 'amalgkit', {'num_fastp_in': 100})
    if change == 'missing':
        path.unlink()
    elif change == 'corrupt':
        path.write_text('not json')
    else:
        data = json.loads(path.read_text())
        if change == 'unfinished':
            data['status'] = 'running'
        elif change == 'wrong_run':
            data['runs'][0]['run'] = 'SRR2'
        else:
            data['runs'] *= 2
        path.write_text(json.dumps(data))
    assert set(collect(tmp_path, 'Plant', 'SRR1').values()) == {'unknown'}


def test_unrecognized_fallback_schema_does_not_assume_no_fastp(tmp_path):
    read_evidence(tmp_path, 'Plant', 'SRR1', 'public-original', schema_version=5)
    assert collect(tmp_path, 'Plant', 'SRR1')['fastp_applied'] == 'unknown'


def test_stage_reports_can_differ_without_mixing_generating_releases():
    common = {'genegalleon_version': '1.0.0'}
    assembly = {'read_acquisition_method': 'amalgkit', 'fastp_applied': 'yes'}
    quant = {'read_acquisition_method': 'genegalleon_fallback', 'fastp_applied': 'no'}
    # A reused reference can carry an old quant record; the new quant product
    # supplies the authoritative record for that stage.
    original = {**common, 'read_processing': {'assembly': assembly, 'quant': assembly}}
    products = {'reference': {'provenance': {'processing': original}},
                'busco': {'provenance': {'processing': original}},
                'quant': {'provenance': {'processing': {**common, 'read_processing': {'quant': quant}}}}}
    processing = native_processing(products)
    assert processing['read_processing'] == {'assembly': assembly, 'quant': quant}
    assert summary(processing) == {'read_acquisition_method': 'mixed', 'fastp_applied': 'mixed'}
    assert original['read_processing']['quant'] == assembly
    assert summary({'read_processing': {'assembly': assembly}}) == {
        'read_acquisition_method': 'unknown', 'fastp_applied': 'unknown'}
