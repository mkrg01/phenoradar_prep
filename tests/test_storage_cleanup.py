"""Large scratch can be retired without losing retry or published reuse."""
import fcntl
import json
import os
import subprocess
from pathlib import Path

import pytest
import yaml

from cleanup_work import cleanup, inspect_job, scratch_usage
from common import file_record
from configuration import validate_keys
from dataset import load, materialize, status, submit, worker
from run_odb_chunk import run as run_odb
from run_kofam import run as run_kofam
from test_datasets import dataset_project, new_dataset, native_events
from test_kofam import kofam_job


@pytest.fixture
def odb_job(fake_odb, frozen_reference, tmp_path, monkeypatch):
    protein = tmp_path / 'Alpha_plant_protein.fa'
    protein.write_text('>Alpha_plant_g1\nMK*\n')
    manifest = tmp_path / 'chunk_000.fs'
    manifest.write_text(str(protein) + '\n')
    events = tmp_path / 'events.txt'
    monkeypatch.setenv('FAKE_ODB_LOG', str(events))
    monkeypatch.delenv('SLURM_CPUS_PER_TASK', raising=False)
    return dict(manifest=manifest, reference=frozen_reference / 'reference.json',
                output_dir=tmp_path / 'out', work_dir=tmp_path / 'work', label='chunk_000',
                command=str(fake_odb), jobs=1, batch_size=1), events


def test_odb_cleanup_preserves_failed_resume_and_published_reuse(odb_job, tmp_path, monkeypatch):
    args, events = odb_job
    reference = file_record(args['reference'])
    fail = tmp_path / 'fail_once'
    fail.touch()
    monkeypatch.setenv('FAKE_ODB_FAIL_ONCE', str(fail))
    base = args['work_dir'] / args['label']
    with pytest.raises(subprocess.CalledProcessError):
        run_odb(**args)
    failed = json.loads((base / 'status.json').read_text())
    work = Path(failed['work'])
    assert work.is_dir()
    assert not (base / 'cleanup.json').exists()
    run_odb(**args)
    assert len(set(events.read_text().splitlines())) == 1
    assert not work.exists()
    assert json.loads((base / 'cleanup.json').read_text())['state'] == 'complete'
    assert not (args['output_dir'] / 'native_results').exists()
    assert file_record(args['reference']) == reference
    before = events.read_text()
    run_odb(**args)
    assert events.read_text() == before
    assert json.loads((base / 'status.json').read_text())['reused']
    # Corrupt results cannot be used as proof that scratch can be deleted.
    output = args['output_dir'] / 'chunk_000.og.annotations'
    output.write_text(output.read_text().replace('OG1', 'OG9'))
    run_odb(**args)
    assert events.read_text() != before
    assert not work.exists()


def test_legacy_cleanup_preview_apply_and_active_job_lock(odb_job):
    args, _ = odb_job
    run_odb(**args, keep_intermediates=True)
    base = args['work_dir'] / args['label']
    work = Path(json.loads((base / 'status.json').read_text())['work'])
    reference_data = Path(json.loads(args['reference'].read_text())['data_dir'])
    original = {p: file_record(p) for p in reference_data.rglob('*') if p.is_file()}
    before = {p: p.stat().st_mtime_ns for p in work.rglob('*') if p.is_file()}
    report = inspect_job('odb', base, args['output_dir'])
    assert report['files'] > 0 and report['reclaimable_bytes'] > 0
    assert work.is_dir() and not (base / 'cleanup.json').exists()
    assert all(p.stat().st_mtime_ns == stamp for p, stamp in before.items())
    with (base / '.lock').open('r') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(BlockingIOError):
            inspect_job('odb', base, args['output_dir'], apply=True)
    inspect_job('odb', base, args['output_dir'], apply=True)
    assert not work.exists()
    assert all(file_record(p) == record for p, record in original.items())
    assert inspect_job('odb', base, args['output_dir'])['reclaimable_bytes'] == 0


def test_cleanup_refuses_changed_output_and_failed_status(odb_job):
    args, _ = odb_job
    run_odb(**args, keep_intermediates=True)
    base = args['work_dir'] / args['label']
    work = Path(json.loads((base / 'status.json').read_text())['work'])
    (args['output_dir'] / 'chunk_000.og.hits').unlink()
    with pytest.raises(ValueError, match='missing or changed'):
        inspect_job('odb', base, args['output_dir'], apply=True)
    assert work.is_dir()
    current = json.loads((base / 'status.json').read_text())
    current['state'] = 'failed'
    (base / 'status.json').write_text(json.dumps(current))
    with pytest.raises(ValueError, match='successful, inactive'):
        inspect_job('odb', base, args['output_dir'], apply=True)
    assert work.is_dir()


@pytest.mark.parametrize('overlap', ['same', 'output_in_work', 'work_in_output'])
def test_odb_rejects_overlapping_directories(odb_job, overlap):
    args, _ = odb_job
    if overlap == 'same': args['output_dir'] = args['work_dir']
    elif overlap == 'output_in_work': args['output_dir'] = args['work_dir'] / 'out'
    else: args['work_dir'] = args['output_dir'] / 'work'
    with pytest.raises(ValueError, match='overlap'):
        run_odb(**args)


def test_kofam_cleanup_and_retry_without_recomputing(kofam_job):
    args = kofam_job['args']
    result = run_kofam(**args)
    work = Path(json.loads((args['work_dir'] / 'status.json').read_text())['work'])
    assert not work.exists()
    assert (args['output_dir'] / 'stdout.log').exists()
    assert run_kofam(**args)['fingerprint'] == result['fingerprint']
    assert len(kofam_job['events'].read_text().splitlines()) == 1
    assert inspect_job('kofam', args['work_dir'], args['output_dir'])['files'] == 0


def test_kofam_publication_failure_keeps_scratch(kofam_job, monkeypatch):
    import run_kofam as module
    def fail(*args): raise OSError('publication failed')
    monkeypatch.setattr(module, '_publish', fail)
    with pytest.raises(OSError, match='publication failed'):
        run_kofam(**kofam_job['args'])
    base = kofam_job['args']['work_dir']
    record = json.loads((base / 'status.json').read_text())
    assert record['state'] == 'failed'
    assert (Path(record['work']) / 'protein.faa').is_file()
    assert not (base / 'cleanup.json').exists()


@pytest.mark.parametrize('keep', [False, True])
def test_genegalleon_keeps_reads_until_verified_quant(dataset_project, monkeypatch, keep):
    config = dataset_project / 'config/build.yaml'
    cfg = yaml.safe_load(config.read_text())
    cfg['storage'] = {'keep_intermediates': keep}
    config.write_text(yaml.safe_dump(cfg))
    build = new_dataset(dataset_project)
    submit(build, until='quant', dry_run=True)
    assert load(build)['analysis']['storage']['keep_intermediates'] is keep
    work = build / 'work/genegalleon/New_plant_SRR1'
    out = work / 'output/transcriptome_assembly'
    fastq = out / 'amalgkit_getfastq/New_plant/SRR1/SRR1.amalgkit.fastq.gz'
    fastq.parent.mkdir(parents=True)
    download_log = fastq.parent / 'download.log'
    download_log.write_text('download QC')
    # A linked original is never altered, even when the staged name is retired.
    original = dataset_project / 'original.fastq.gz'
    original.write_bytes(b'read data')
    os.link(original, fastq)
    temporary = out / 'tmp/1_New_plant/large_intermediate'
    temporary.parent.mkdir(parents=True)
    temporary.write_bytes(b'assembly scratch')
    monkeypatch.setenv('FAKE_GG_INCOMPLETE_MERGE', '1')
    with pytest.raises(subprocess.CalledProcessError):
        worker(build, 1)
    assert fastq.exists()
    assert temporary.exists() is keep
    assert status(build)[0]['assembly'] == 'reuse'
    assert status(build)[0]['busco'] == 'reuse'
    monkeypatch.delenv('FAKE_GG_INCOMPLETE_MERGE')
    worker(build, 1)
    assert temporary.exists() is keep
    assert fastq.exists() is keep
    assert original.read_bytes() == b'read data'
    assert download_log.read_text() == 'download QC'
    before = native_events(build)
    worker(build, 1)
    assert native_events(build) == before
    assert status(build)[0]['quant'] == 'reuse'
    assert (materialize(build) / 'cds/New_plant_SRR1_longestCDS.fa.gz').is_file()


def test_cleanup_does_not_follow_symlinks_or_count_external_hardlinks(tmp_path):
    root, external = tmp_path / 'owned', tmp_path / 'external'
    root.mkdir(); external.mkdir()
    original = external / 'original'
    original.write_bytes(b'x' * 8192)
    work = root / 'scratch'
    work.mkdir()
    (work / 'reference').symlink_to(external, target_is_directory=True)
    os.link(original, work / 'copy')
    (work / 'unique').write_bytes(b'y' * 8192)
    usage = scratch_usage(work)
    assert usage['reclaimable_bytes'] == (work / 'unique').stat().st_blocks * 512
    cleanup(root, ['scratch'], root / 'cleanup.json')
    assert original.read_bytes() == b'x' * 8192
    (root / 'alias').symlink_to(external, target_is_directory=True)
    assert scratch_usage(root / 'alias')['reclaimable_bytes'] == 0
    with pytest.warns(RuntimeWarning, match='incomplete'):
        report = cleanup(root, ['alias/original', '../external', '.'], root / 'cleanup.json')
    assert report['state'] == 'pending' and original.exists()


@pytest.mark.parametrize('value', [None, 'false', 0, 1, [], {}])
def test_storage_flag_is_strict_boolean(value):
    with pytest.raises(ValueError, match='storage.keep_intermediates'):
        validate_keys({'storage': {'keep_intermediates': value}})


def test_odb_defers_cleanup_until_reuse_cache_publication(odb_job):
    from incremental_odb import publish
    from common import write_tsv
    args, events = odb_job
    # The rule must leave scratch if its later cache publication fails.
    run_odb(**args, defer_cleanup=True)
    base = args['work_dir'] / args['label']
    work = Path(json.loads((base / 'status.json').read_text())['work'])
    assert work.is_dir()
    assert not (args['output_dir'] / 'native_results').exists()
    samples = args['manifest'].parent / 'samples.tsv'
    write_tsv(samples, ['species', 'odb_species'], [{'species': 'Alpha_plant', 'odb_species': 'Alpha_plant'}])
    with pytest.raises(FileNotFoundError):
        publish(samples, args['output_dir'], args['manifest'].parent / 'missing-proteins',
                args['manifest'].parent / 'cache', args['label'])
    assert work.exists()
    # The verified chunk can retry publication without rerunning mapping.
    before = events.read_text()
    run_odb(**args, defer_cleanup=True)
    assert events.read_text() == before


def test_cleanup_failure_keeps_published_job_successful(odb_job, monkeypatch):
    import cleanup_work
    args, events = odb_job
    original = cleanup_work.shutil.rmtree
    def fail_work(path, *positional, **kwargs):
        if Path(path).parent == args['work_dir'] / args['label']:
            raise PermissionError('scratch temporarily unavailable')
        return original(path, *positional, **kwargs)
    monkeypatch.setattr(cleanup_work.shutil, 'rmtree', fail_work)
    with pytest.warns(RuntimeWarning, match='incomplete'):
        run_odb(**args)
    base = args['work_dir'] / args['label']
    assert json.loads((base / 'status.json').read_text())['state'] == 'success'
    assert json.loads((base / 'cleanup.json').read_text())['state'] == 'pending'
    monkeypatch.setattr(cleanup_work.shutil, 'rmtree', original)
    before = events.read_text()
    run_odb(**args)
    assert events.read_text() == before
    assert json.loads((base / 'cleanup.json').read_text())['state'] == 'complete'


def test_legacy_kofam_cleanup_recovers_work_from_published_command(kofam_job):
    args = kofam_job['args']
    run_kofam(**args, keep_intermediates=True)
    status_path = args['work_dir'] / 'status.json'
    saved = json.loads(status_path.read_text())
    work = Path(saved.pop('work'))
    saved['reused'] = True  # Earlier releases dropped the work path on reuse.
    status_path.write_text(json.dumps(saved))
    preview = inspect_job('kofam', args['work_dir'], args['output_dir'])
    assert preview['work'] == str(work) and preview['files'] > 0
    before = kofam_job['events'].read_text()
    run_kofam(**args)
    assert not work.exists()
    assert kofam_job['events'].read_text() == before


def test_partial_cleanup_can_be_retried_after_identity_was_removed(odb_job, monkeypatch):
    import cleanup_work
    args, _ = odb_job
    base = args['work_dir'] / args['label']
    original = cleanup_work.shutil.rmtree
    def partially_remove(path, *positional, **kwargs):
        if Path(path).parent == base:
            (Path(path) / 'identity.json').unlink()
            raise PermissionError('remaining scratch is temporarily unavailable')
        return original(path, *positional, **kwargs)
    monkeypatch.setattr(cleanup_work.shutil, 'rmtree', partially_remove)
    with pytest.warns(RuntimeWarning, match='incomplete'):
        run_odb(**args)
    monkeypatch.setattr(cleanup_work.shutil, 'rmtree', original)
    report = inspect_job('odb', base, args['output_dir'], apply=True)
    assert report['cleanup']['state'] == 'complete'
    assert not Path(report['work']).exists()


def test_genegalleon_scan_refuses_linked_parent_without_failing_results(tmp_path):
    from dataset import cleanup_genegalleon
    build = tmp_path / 'build'
    out = build / 'work/genegalleon/Plant_R1/output/transcriptome_assembly'
    out.mkdir(parents=True)
    external = tmp_path / 'private/Plant'
    external.mkdir(parents=True)
    reads = external / 'original.fastq'
    reads.write_bytes(b'private reads')
    (out / 'amalgkit_getfastq').symlink_to(external.parent, target_is_directory=True)
    with pytest.warns(RuntimeWarning, match='incomplete'):
        cleanup_genegalleon(build, {'config': {}, 'analysis': {'output_root': 'test'}},
                           {'species': 'Plant_R1', 'row': {'species_id': 'Plant'}},
                           'quant', {'reference': True, 'quant': True})
    assert reads.read_bytes() == b'private reads'
    report = json.loads((build / 'jobs/cleanup/Plant_R1.quant.json').read_text())
    assert report['state'] == 'pending'
