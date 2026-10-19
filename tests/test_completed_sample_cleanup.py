"""Completed native work retires without breaking mapping inputs or retries."""
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from common import write_json
from dataset import cleanup_genegalleon, item_products, load, materialize, status, submit, worker
from dataset_assets import record
from genegalleon_cleanup import DIRECTORIES, POLICY_VERSION
from storage_management import run_storage
from test_datasets import dataset_project, native_events, new_dataset


def seed_work(work):
    for name in DIRECTORIES:
        directory = work / name
        directory.mkdir(parents=True, exist_ok=True)
        (directory / 'disposable').write_bytes(b'x' * 8192)


def test_success_retires_work_and_remains_reusable_before_mapping(dataset_project):
    build = new_dataset(dataset_project)
    submit(build, until='quant', dry_run=True)
    work = build / 'work/genegalleon/New_plant_SRR1'
    seed_work(work)
    worker(build, 1)
    assert all(not (work / name).exists() for name in DIRECTORIES)
    manifest = load(build)
    products = item_products(manifest, manifest['items'][0])
    assert Path(products['busco']['short']['path']).exists()
    assert status(build)[0]['quant'] == 'reuse'
    assert (materialize(build) / 'metadata.tsv').is_file()
    events = native_events(build)
    worker(build, 1)
    assert native_events(build) == events
    receipt = json.loads((build / 'jobs/cleanup/New_plant_SRR1.sample.json').read_text())
    assert receipt['completed_sample_cleanup'] == POLICY_VERSION


def test_legacy_complete_receipt_is_upgraded_without_running_tools(dataset_project):
    build = new_dataset(dataset_project)
    submit(build, until='quant', dry_run=True)
    worker(build, 1)
    work = build / 'work/genegalleon/New_plant_SRR1'
    seed_work(work)
    receipt = build / 'jobs/cleanup/New_plant_SRR1.sample.json'
    write_json(receipt, {'state': 'complete', 'removed': [], 'errors': []})
    events = native_events(build)
    preview = run_storage(build, 'build', inspect=True)
    assert preview['counts'] == {'available': 1}
    assert preview['jobs'][0]['inspection']['reclaimable_bytes'] >= 8192 * len(DIRECTORIES)
    assert all((work / name).exists() for name in DIRECTORIES)
    applied = run_storage(build, 'build', inspect=True, apply=True)
    assert applied['counts'] == {'complete': 1}
    assert all(not (work / name).exists() for name in DIRECTORIES)
    assert native_events(build) == events


def test_partial_cleanup_failure_is_retryable_without_recomputation(dataset_project, monkeypatch):
    build = new_dataset(dataset_project)
    submit(build, until='quant', dry_run=True)
    work = build / 'work/genegalleon/New_plant_SRR1'
    seed_work(work)
    blocked = work / 'output/transcriptome_assembly/amalgkit_quant'
    remove = shutil.rmtree

    def fail_one(path, *args, **kwargs):
        if Path(path) == blocked:
            raise PermissionError('temporary failure')
        return remove(path, *args, **kwargs)

    monkeypatch.setattr(shutil, 'rmtree', fail_one)
    with pytest.warns(RuntimeWarning, match='cleanup incomplete'):
        worker(build, 1)
    assert status(build)[0]['quant'] == 'reuse' and blocked.exists()
    assert run_storage(build, 'build')['counts'] == {'pending': 1}
    events = native_events(build)
    monkeypatch.setattr(shutil, 'rmtree', remove)
    assert run_storage(build, 'build', inspect=True, apply=True)['counts'] == {'complete': 1}
    assert not blocked.exists() and native_events(build) == events


@pytest.fixture
def completed(tmp_path):
    work = tmp_path / 'work/genegalleon/Plant_R1'
    seed_work(work)
    published = work / 'products'
    published.mkdir()
    entries = []
    for name in ('cds.fa.gz', 'full.tsv', 'short.txt', 'abundance.tsv.gz'):
        path = published / name
        path.write_text('registered ' + name)
        entries.append(record(path))
    products = {'reference': {'cds': entries[0]}, 'busco': {'full': entries[1], 'short': entries[2]},
                'quant': {'abundance': entries[3]}}
    item = {'species': 'Plant_R1', 'row': {'species_id': 'Plant'}}
    manifest = {'analysis': {'output_root': 'output'}, 'config': {'storage': {'keep_intermediates': False}},
                'raw_inputs': {}}
    receipt = tmp_path / 'jobs/status/Plant_R1.sample.json'
    write_json(receipt, {'state': 'complete'})
    return tmp_path, work, manifest, item, products


@pytest.mark.parametrize('state', ['failed', 'running', 'interrupted'])
def test_failed_or_active_samples_retain_native_work_even_with_products(completed, state):
    build, work, manifest, item, products = completed
    write_json(build / 'jobs/status/Plant_R1.sample.json', {'state': state})
    cleanup_genegalleon(build, manifest, item, 'sample', products)
    assert all((work / name).exists() for name in DIRECTORIES)


def private_read(completed, *, hardlink=False):
    build, work, manifest, _, _ = completed
    original = build / 'original.fastq'
    original.write_bytes(b'original reads')
    staged = work / 'input/reads/Plant_R1/read1_path.fastq'
    staged.parent.mkdir(parents=True)
    if hardlink:
        os.link(original, staged)
    else:
        shutil.copyfile(original, staged)
    manifest['raw_inputs'] = {'Plant_R1': {'read1_path': {'path': str(original)}}}
    return original, staged


@pytest.mark.parametrize('hardlink', [False, True])
def test_private_read_cleanup_preserves_original_and_shared_cache(completed, hardlink):
    build, work, manifest, item, products = completed
    original, staged = private_read(completed, hardlink=hardlink)
    shared = build / 'resources/taxonomy/ete/dataset/taxa.sqlite'
    shared.parent.mkdir(parents=True)
    shared.write_text('shared database')
    directory = work / 'downloads/ete_taxonomy'
    shutil.rmtree(directory)
    directory.symlink_to(shared.parent, target_is_directory=True)
    preview = cleanup_genegalleon(build, manifest, item, 'sample', products, apply=False)
    assert staged.is_file() and original.is_file()
    assert preview['files'] > 0
    result = cleanup_genegalleon(build, manifest, item, 'sample', products)
    assert result['state'] == 'complete' and not staged.exists()
    assert original.read_bytes() == b'original reads'
    assert shared.read_text() == 'shared database'


@pytest.mark.parametrize('change', ['missing', 'changed', 'original_inside', 'untracked'])
def test_cleanup_retains_unique_private_data(completed, change):
    build, work, manifest, item, products = completed
    original, staged = private_read(completed)
    if change == 'missing':
        original.unlink()
    elif change == 'changed':
        original.write_bytes(b'different data')  # Same size, so content verification matters.
    elif change == 'original_inside':
        original = staged.parent / 'original.fastq'
        original.write_bytes(staged.read_bytes())
        manifest['raw_inputs']['Plant_R1']['read1_path']['path'] = str(original)
    else:
        (staged.parent / 'notes.txt').write_text('user file')
    result = cleanup_genegalleon(build, manifest, item, 'sample', products)
    if change == 'untracked':
        assert not staged.exists() and (staged.parent / 'notes.txt').read_text() == 'user file'
    else:
        assert staged.read_bytes() == b'original reads' and result['retained']


def test_legacy_registered_native_files_are_protected(completed):
    build, work, manifest, item, products = completed
    native = work / 'output/transcriptome_assembly/longest_cds/native.fa.gz'
    native.write_bytes(b'registered native CDS')
    products['reference']['cds'] = record(native)
    cleanup_genegalleon(build, manifest, item, 'sample', products)
    assert native.read_bytes() == b'registered native CDS'


def test_corrupt_product_blocks_retirement(completed):
    build, work, manifest, item, products = completed
    Path(products['reference']['cds']['path']).write_text('corrupted')
    with pytest.warns(RuntimeWarning, match='cleanup incomplete'):
        result = cleanup_genegalleon(build, manifest, item, 'sample', products)
    assert result['state'] == 'pending'
    assert all((work / name).exists() for name in DIRECTORIES)


def test_debug_retention_keeps_completed_work_and_private_reads(completed):
    build, work, manifest, item, products = completed
    _, staged = private_read(completed)
    manifest['config']['storage']['keep_intermediates'] = True
    cleanup_genegalleon(build, manifest, item, 'sample', products)
    assert staged.exists() and all((work / name).exists() for name in DIRECTORIES)


def test_shared_cache_is_protected_even_inside_workspace(completed):
    build, work, manifest, item, products = completed
    shared = work / 'downloads/ete_taxonomy/shared/dataset'
    shared.mkdir(parents=True)
    manifest['taxonomy_reference'] = {'dataset': str(shared), 'receipt': {'path': str(shared.parent / 'reference.json')}}
    cleanup_genegalleon(build, manifest, item, 'sample', products)
    assert (work / 'downloads/ete_taxonomy/disposable').exists()


def test_incomplete_sample_retains_native_work(completed):
    build, work, manifest, item, products = completed
    products['quant'] = None
    cleanup_genegalleon(build, manifest, item, 'sample', products)
    assert all((work / name).exists() for name in DIRECTORIES)
