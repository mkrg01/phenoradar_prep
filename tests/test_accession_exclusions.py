"""Run exclusions select new builds without changing old snapshots or caches."""
import json
import subprocess
from pathlib import Path

import pytest
import yaml

from accession_exclusions import read_exclusions
from common import read_tsv, write_json, write_tsv
from dataset import load, materialize, plan, prepare, status, submit, worker
from dataset_assets import import_existing
from test_datasets import dataset_project, fake_genegalleon, imported, new_dataset


@pytest.mark.parametrize('content,message', [
    ('', 'including accession'),
    ('run\treason\nSRR1\tfailed\n', 'including accession'),
    ('accession\taccession\nSRR1\tSRR1\n', 'unique columns'),
    ('accession\treason\nSRR1\tfirst\nSRR1\tsecond\n', 'duplicate'),
    ('accession\treason\n SRR1\tfailed\n', 'invalid'),
    ('accession\treason\nSRR1\textra\tcolumn\n', 'malformed'),
])
def test_invalid_exclusion_registry_fails_clearly(tmp_path, content, message):
    path = tmp_path/'excluded.tsv'; path.write_text(content)
    with pytest.raises(ValueError, match=message): read_exclusions(path)


def test_registry_allows_empty_policy_local_ids_and_future_notes(tmp_path):
    path = tmp_path/'excluded.tsv'; path.write_text('accession\treason\n')
    assert read_exclusions(path) == {}
    path.write_text('accession\treason\tnotes\nLOCAL1\tassembly_failed\tmanually checked\nSRR123\t\t\n')
    policy = read_exclusions(path)
    assert set(policy) == {'LOCAL1', 'SRR123'}
    assert policy['LOCAL1']['notes'] == 'manually checked'
    assert read_exclusions(None) == {}


def test_new_build_filters_before_inspecting_cached_or_private_inputs(dataset_project):
    root = dataset_project
    store = imported(root)
    # This excluded row would otherwise require a missing private FASTQ and conflict with bad CDS.
    metadata = root/'input/metadata.tsv'
    rows = read_tsv(metadata)
    fields = list(rows[0]) + ['private_file', 'lib_layout', 'read1_path']
    for row in rows: row.update(private_file='', lib_layout='', read1_path='')
    rows[0].update(private_file='yes', lib_layout='single', read1_path='missing.fastq')
    write_tsv(metadata, fields, rows)
    (root/'input/cds/Alpha_plant_longestCDS.fa.gz').write_bytes(b'invalid cached bytes')
    policy = root/'config/excluded_accessions.tsv'
    policy.write_text('accession\treason\nA1\tunusable_run\nSRR999\tprevious_failure\n')
    report = plan(root,root/'config/build.yaml')[-1]
    assert [r['species'] for r in report] == ['Beta_sp-X','Gamma_plant','Alpha_plant']
    assert report[-1]['assembly'] == report[-1]['mapping'] == 'excluded'
    assert report[-1]['reason'] == 'excluded_accession: unusable_run'
    build = prepare(root,'excluded',root/'config/build.yaml')
    assert [i['row']['run'] for i in load(build)['items']] == ['B1','G1']
    assert load(build)['software_lock'] is None
    assert (build/'source_metadata.tsv').read_bytes() == metadata.read_bytes()
    assert (build/'excluded_accessions.tsv').read_bytes() == policy.read_bytes()
    assert read_tsv(build/'excluded_runs.tsv') == [{'species':'Alpha_plant','run':'A1','reason':'unusable_run'}]
    assert submit(build,until='quant',dry_run=True) == []
    inputs = materialize(build)
    assert [r['run'] for r in read_tsv(inputs/'metadata.tsv')] == ['B1','G1']
    assert not (inputs/'cds/Alpha_plant_longestCDS.fa.gz').exists()
    assert (store/'Alpha_plant').is_dir()  # Policy does not delete any source product.
    # Edits to the source policy do not alter the frozen build.
    policy.write_text('accession\treason\nB1\tnew_decision\n')
    assert status(build)[-1]['run'] == 'A1'
    assert status(build)[-1]['assembly'] == 'excluded'
    (build/'excluded_accessions.tsv').write_text('accession\treason\n')
    with pytest.raises(ValueError, match='registered file changed'): load(build)


def test_exclusion_changes_require_a_new_build_and_removal_restores_reuse(dataset_project):
    root = dataset_project; imported(root)
    cfg = root/'config/build.yaml'; policy = root/'config/excluded_accessions.tsv'
    old = prepare(root,'before',cfg)
    policy.write_text('accession\treason\nA1\tfailed_repeatedly\n')
    updated = prepare(root,'after',cfg)
    assert len(load(old)['items']) == 3
    assert len(load(updated)['items']) == 2
    assert len(read_tsv(materialize(old)/'metadata.tsv')) == 3
    policy.write_text('accession\treason\n')
    restored = prepare(root,'restored',cfg)
    assert len(load(restored)['items']) == 3
    assert submit(restored,until='quant',dry_run=True) == []


def test_all_excluded_can_be_planned_but_cannot_start_empty_build(dataset_project):
    root = dataset_project
    (root/'config/excluded_accessions.tsv').write_text('accession\treason\nA1\tfail\nB1\tfail\nG1\tfail\n')
    assert all(r['assembly'] == 'excluded' for r in plan(root,root/'config/build.yaml')[-1])
    with pytest.raises(ValueError,match='all metadata runs are excluded'):
        prepare(root,'empty',root/'config/build.yaml')
    assert not (root/'builds/empty').exists()
    assert not (root/'resources/software').exists()


def test_excluded_species_never_receive_array_indices(dataset_project):
    root = dataset_project; fake_genegalleon(root)
    (root/'config/excluded_accessions.tsv').write_text('accession\treason\nA1\tdownload_failed\n')
    build = prepare(root,'two_species',root/'config/build.yaml')
    commands = submit(build,until='assembly',dry_run=True)
    assert len(commands) == 1 and '--array=1,2%5' in commands[0]
    for index in (1,2): worker(build,'assembly',index)
    events = [json.loads(line)['species'] for line in (build/'genegalleon/events.jsonl').read_text().splitlines()]
    assert events == ['Beta_sp-X','Gamma_plant']
    assert not (build/'genegalleon/input/amalgkit_metadata/Alpha_plant_metadata.tsv').exists()
    assert status(build)[-1]['assembly'] == 'excluded'


def test_legacy_registration_respects_run_exclusions(dataset_project):
    root = dataset_project
    result = import_existing(root/'resources/dataset_assets',root/'input',root/'input/metadata.tsv',excluded_runs={'A1'})
    assert result[0] == {'species':'Alpha_plant','run':'A1','status':'excluded'}
    assert not (root/'resources/dataset_assets/Alpha_plant').exists()
    assert result[1]['status'] == 'registered'


def test_download_failure_is_retried_without_rerunning_completed_species(dataset_project, monkeypatch):
    root = dataset_project
    build = new_dataset(root,('New plant','Other plant'))
    submit(build,until='assembly',dry_run=True)
    worker(build,'assembly',2)
    raw = build/'genegalleon/downloads/SRR1.partial'; raw.write_text('retained resumable download')
    monkeypatch.setenv('FAKE_GG_FAIL_DOWNLOAD','1')
    with pytest.raises(subprocess.CalledProcessError): worker(build,'assembly',1)
    receipt = json.loads((build/'jobs/status/New_plant.assembly.json').read_text())
    assert receipt['state'] == 'failed' and receipt['run'] == 'SRR1' and receipt['stage'] == 'assembly'
    assert [r['assembly'] for r in status(build)] == ['pending','reuse']
    assert '--array=1%5' in submit(build,until='assembly',dry_run=True)[0]
    monkeypatch.delenv('FAKE_GG_FAIL_DOWNLOAD')
    worker(build,'assembly',1)
    assert raw.read_text() == 'retained resumable download'
    assert [r['assembly'] for r in status(build)] == ['reuse','reuse']
    assert submit(build,until='assembly',dry_run=True) == []


def test_interrupted_worker_running_receipt_does_not_block_retry(dataset_project):
    build = new_dataset(dataset_project)
    submit(build,until='assembly',dry_run=True)
    partial = build/'genegalleon/output/transcriptome_assembly/longest_cds/New_plant_longestCDS.fa.gz'
    partial.parent.mkdir(parents=True); partial.write_bytes(b'incomplete output before timeout')
    write_json(build/'jobs/status/New_plant.assembly.json', {'state':'running','run':'SRR1'})
    assert status(build)[0]['assembly'] == 'pending'
    worker(build,'assembly',1)
    assert status(build)[0]['assembly'] == 'reuse'
    assert list((build/'jobs/incomplete/New_plant/assembly').rglob('*longestCDS.fa.gz'))
