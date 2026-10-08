"""Reviewed run exclusions free failed reads without sweeping new failures or live jobs."""
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

import dataset
import metadata_catalog as catalog
import sample_review
from common import read_tsv, write_json, write_tsv
from test_sample_curation import sample_project, finish
from test_datasets import dataset_project


def fail(build, index, monkeypatch, stage='BUSCO'):
    monkeypatch.setenv('FAKE_GG_FAIL_' + stage, '1')
    with pytest.raises(subprocess.CalledProcessError):
        finish(build, index)
    monkeypatch.delenv('FAKE_GG_FAIL_' + stage)


def leftovers(build, index):
    item = dataset.load(build)['items'][index - 1]
    work = dataset.workspace(build, dataset.load(build), item)
    out = work / 'output/transcriptome_assembly'
    reads = out / 'amalgkit_getfastq' / item['row']['species_id'] / item['row']['run']
    reads.mkdir(parents=True, exist_ok=True)
    fastq = reads / 'download.amalgkit.fastq.gz'
    fastq.write_bytes(b'failed sample reads' * 1024)
    partial = reads / 'download.sra.part'
    partial.write_bytes(b'partial download')
    log = reads / 'download.log'
    log.write_text('download diagnostics')
    temporary = out / 'tmp/1_native/scratch'
    temporary.parent.mkdir(parents=True, exist_ok=True)
    temporary.write_bytes(b'assembly scratch')
    return work, fastq, partial, log, temporary


def edit_plan(plan, actions):
    rows = read_tsv(plan)
    for row in rows:
        if row['run'] in actions:
            row['action'], row['reason'] = actions[row['run']]
    write_tsv(plan, sample_review.FIELDS, rows)


def test_ordinary_record_holds_failures_and_preserves_failed_reads(sample_project, monkeypatch):
    root, cfg, build = sample_project
    fail(build, 1, monkeypatch)
    _, fastq, partial, log, temporary = leftovers(build, 1)
    result = catalog.record_successes(root, cfg, build)
    row = result['samples'][0]
    assert row['action'] == 'hold' and row['job_state'] == 'failed'
    assert row['stopped_at'] == 'busco'
    assert result['excluded_runs'] == result['cleanup'] == []
    assert all(p.exists() for p in (fastq, partial, log, temporary))
    assert result['samples'][2]['action'] == 'hold'


def test_direct_command_excludes_all_incomplete_samples_with_optional_read_only_preview(sample_project, monkeypatch):
    root, cfg, build = sample_project
    finish(build, 1)
    fail(build, 2, monkeypatch, 'QUANT')
    _, failed_fastq, _, failed_log, failed_scratch = leftovers(build, 2)
    _, unstarted_fastq, _, unstarted_log, unstarted_scratch = leftovers(build, 3)
    cli = [sys.executable, str(root / 'workflow/scripts/metadata_catalog.py'), 'record',
           '--config', str(cfg['_path']), '--build', str(build), '--exclude-failed']
    preview = subprocess.run([*cli, '--dry-run'], cwd=root, capture_output=True, text=True, check=True)
    report = json.loads(preview.stdout)
    assert report['excluded_runs'] == ['B1', 'G1']
    assert report['samples'][1]['stopped_at'] == 'quant'
    assert report['samples'][2]['job_state'] == 'not_started'
    assert not (root / cfg['accepted_samples']).exists()
    assert read_tsv(root / cfg['excluded_accessions']) == []
    assert all(p.exists() for p in (failed_fastq, unstarted_fastq, failed_scratch, unstarted_scratch))
    applied = subprocess.run(cli, cwd=root, capture_output=True, text=True, check=True)
    result = json.loads(applied.stdout)
    assert result['recorded_runs'] == ['A1'] and result['excluded_runs'] == ['B1', 'G1']
    assert all(not p.exists() for p in (failed_fastq, unstarted_fastq, failed_scratch, unstarted_scratch))
    assert failed_log.exists() and unstarted_log.exists()
    assert [row['run'] for row in read_tsv(root / cfg['accepted_samples'])] == ['A1']
    assert dataset.status(build)[1]['assembly'] == dataset.status(build)[1]['busco'] == 'reuse'
    repeated = subprocess.run(cli, cwd=root, capture_output=True, text=True, check=True)
    assert json.loads(repeated.stdout)['excluded_runs'] == []


def test_saved_plan_adopts_success_and_excludes_all_inactive_incomplete_samples(sample_project, monkeypatch):
    root, cfg, build = sample_project
    finish(build, 1)
    fail(build, 2, monkeypatch)
    work, fastq, partial, log, temporary = leftovers(build, 2)
    original = root / 'caller.fastq.gz'
    original.write_bytes(b'caller reads')
    input_reads = work / 'input/reads/original.fastq.gz'
    input_reads.parent.mkdir(parents=True, exist_ok=True)
    input_reads.symlink_to(original)
    os.link(original, fastq.parent / 'linked.fastq.gz')
    frozen = (build / 'build.json').read_bytes()
    plan = root / 'work/review.tsv'
    preview = catalog.record_successes(root, cfg, build, exclude_failed=True, plan=plan)
    assert preview['dry_run'] and preview['excluded_runs'] == ['B1', 'G1']
    assert [r['action'] for r in read_tsv(plan)] == ['accept', 'exclude', 'exclude']
    assert all(p.exists() for p in (fastq, partial, temporary))
    assert not (root / cfg['accepted_samples']).exists()
    assert read_tsv(root / cfg['excluded_accessions']) == []
    assert next(r for r in preview['samples'] if r['run'] == 'B1')['reclaimable_bytes'] > 0
    preview_apply = sample_review.apply_plan(root, cfg, plan, dry_run=True)
    assert preview_apply['excluded_runs'] == ['B1', 'G1'] and fastq.exists()
    result = sample_review.apply_plan(root, cfg, plan)
    assert result['recorded_runs'] == ['A1'] and result['excluded_runs'] == ['B1', 'G1']
    assert [r['run'] for r in read_tsv(root / cfg['accepted_samples'])] == ['A1']
    assert read_tsv(root / cfg['excluded_accessions'])[0]['reason'] == 'processing_abandoned_busco'
    assert all(not p.exists() for p in (fastq, partial, temporary))
    assert log.exists() and input_reads.exists() and original.read_bytes() == b'caller reads'
    assert dataset.status(build)[1]['assembly'] == 'reuse'
    assert (build / 'build.json').read_bytes() == frozen
    repeated = sample_review.apply_plan(root, cfg, plan)
    assert repeated['cleanup_only'] and repeated['cleanup'][0]['state'] == 'complete'
    assert len(read_tsv(root / cfg['excluded_accessions'])) == 2


def test_plan_does_not_sweep_failures_that_happen_after_review(sample_project, monkeypatch):
    root, cfg, build = sample_project
    finish(build, 1)
    plan = root / 'work/review.tsv'
    catalog.record_successes(root, cfg, build, exclude_failed=True, plan=plan)
    edit_plan(plan, {'B1': ('hold', 'retry_later'), 'G1': ('hold', 'retry_later')})
    fail(build, 3, monkeypatch, 'ASSEMBLY')
    _, fastq, _, _, _ = leftovers(build, 3)
    result = sample_review.apply_plan(root, cfg, plan)
    assert result['recorded_runs'] == ['A1'] and result['excluded_runs'] == []
    assert fastq.exists()


def test_plan_can_hold_a_failure_and_exclude_an_already_accepted_low_quality_run(sample_project, monkeypatch):
    root, cfg, build = sample_project
    monkeypatch.setenv('FAKE_GG_HALF_BUSCO', '1')
    finish(build, 1)
    monkeypatch.delenv('FAKE_GG_HALF_BUSCO')
    catalog.record_successes(root, cfg, build, runs=['A1'])
    definition = yaml.safe_load(cfg['_path'].read_text())
    definition['busco_threshold'] = 0.75
    cfg['_path'].write_text(yaml.safe_dump(definition))
    cfg = catalog.configuration(root, cfg['_path'])
    fail(build, 2, monkeypatch)
    _, fastq, _, _, _ = leftovers(build, 2)
    plan = root / 'work/review.tsv'
    catalog.record_successes(root, cfg, build, runs=['A1', 'B1'], exclude_failed=True, plan=plan)
    assert read_tsv(plan)[0]['action'] == 'hold'
    edit_plan(plan, {'A1': ('exclude', 'reviewed_low_quality'), 'B1': ('hold', 'retry_with_more_memory')})
    result = sample_review.apply_plan(root, cfg, plan)
    assert result['excluded_runs'] == ['A1']
    assert read_tsv(root / cfg['accepted_samples']) == []
    assert fastq.exists()


def test_retry_after_review_invalidates_the_planned_exclusion(sample_project, monkeypatch):
    root, cfg, build = sample_project
    fail(build, 1, monkeypatch)
    plan = root / 'work/review.tsv'
    catalog.record_successes(root, cfg, build, runs=['A1'], exclude_failed=True, plan=plan)
    finish(build, 1)
    with pytest.raises(ValueError, match='state changed'):
        sample_review.apply_plan(root, cfg, plan)
    assert read_tsv(root / cfg['excluded_accessions']) == []
    assert not (root / cfg['accepted_samples']).exists()


@pytest.mark.parametrize('lock_name', ['worker', 'native', 'scratch'])
def test_locked_failed_sample_is_held_without_deletion(sample_project, monkeypatch, lock_name):
    root, cfg, build = sample_project
    fail(build, 1, monkeypatch)
    work, fastq, _, _, _ = leftovers(build, 1)
    item = dataset.load(build)['items'][0]
    lock = (Path(dataset.load(build)['config']['store']) / item['species'] / '.worker.lock'
            if lock_name == 'worker' else work / '.native.lock')
    if lock_name == 'scratch':
        lock = work / 'output/transcriptome_assembly/tmp/1_native/.gg_active.lock'
        lock.touch()
    with lock.open('r') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        result = catalog.record_successes(root, cfg, build, runs=['A1'], exclude_failed=True)
    assert result['samples'][0]['status'] == 'active'
    assert result['samples'][0]['allowed_actions'] == ['hold']
    assert result['excluded_runs'] == [] and fastq.exists()


def scheduler(monkeypatch, output, accounting=None):
    run = subprocess.run
    def fake(args, **kwargs):
        if args[0] == 'squeue':
            return subprocess.CompletedProcess(args, 0, output, '')
        return run(args, **kwargs)
    monkeypatch.setattr(subprocess, 'run', fake)
    monkeypatch.setattr('genegalleon_run.scheduler_states', lambda jobs: accounting or {})


def test_pending_retry_is_held_even_when_previous_attempt_failed(sample_project, monkeypatch):
    root, cfg, build = sample_project
    fail(build, 1, monkeypatch)
    write_json(build / 'jobs/submission_9999.json', {'jobs': [
        {'stage': 'sample', 'indices': [1], 'array_offset': 0, 'job_id': '800', 'state': 'submitted'}]})
    scheduler(monkeypatch, '800_1\n')
    result = catalog.record_successes(root, cfg, build, runs=['A1'], exclude_failed=True)
    assert result['samples'][0]['status'] == 'active'
    assert result['excluded_runs'] == []


def test_accounting_proves_oom_after_wrapper_was_killed(sample_project, monkeypatch):
    root, cfg, build = sample_project
    item = dataset.load(build)['items'][0]
    write_json(build / 'jobs/status' / f"{item['species']}.sample.json",
               {'state': 'running', 'slurm_id': '800_1', 'stage': 'assembly'})
    scheduler(monkeypatch, '', {'800_1': {'state': 'OUT_OF_MEMORY', 'exit_code': '0:9'}})
    result = catalog.record_successes(root, cfg, build, runs=['A1'], exclude_failed=True)
    assert result['samples'][0]['job_state'] == 'failed'
    assert result['excluded_runs'] == ['A1']


def test_accounting_lists_failure_before_worker_could_write_a_receipt(sample_project, monkeypatch):
    root, cfg, build = sample_project
    write_json(build / 'jobs/submission_0001.json', {'jobs': [
        {'stage': 'sample', 'indices': [1], 'array_offset': 0, 'job_id': '800', 'state': 'submitted'}]})
    scheduler(monkeypatch, '', {'800_1': {'state': 'OUT_OF_MEMORY', 'exit_code': '0:9'}})
    result = catalog.record_successes(root, cfg, build, runs=['A1'], exclude_failed=True)
    assert result['excluded_runs'] == ['A1']
    assert result['samples'][0]['stopped_at'] == 'assembly'
    assert result['samples'][0]['job_state'] == 'failed'


def test_unavailable_queue_does_not_authorize_exclusion(sample_project, monkeypatch):
    root, cfg, build = sample_project
    fail(build, 1, monkeypatch)
    item = dataset.load(build)['items'][0]
    receipt = build / 'jobs/status' / f"{item['species']}.sample.json"
    job = json.loads(receipt.read_text())
    job['slurm_id'] = '800_1'
    write_json(receipt, job)
    original = subprocess.run
    def unavailable(args, **kwargs):
        if args[0] == 'squeue':
            raise FileNotFoundError('squeue unavailable')
        return original(args, **kwargs)
    monkeypatch.setattr(subprocess, 'run', unavailable)
    result = catalog.record_successes(root, cfg, build, runs=['A1'], exclude_failed=True)
    assert result['samples'][0]['status'] == 'active'
    assert result['excluded_runs'] == []


def test_cleanup_failure_keeps_decision_and_same_plan_retries_cleanup(sample_project, monkeypatch):
    root, cfg, build = sample_project
    fail(build, 1, monkeypatch)
    _, fastq, _, _, temporary = leftovers(build, 1)
    plan = root / 'work/review.tsv'
    catalog.record_successes(root, cfg, build, runs=['A1'], exclude_failed=True, plan=plan)
    import cleanup_work
    remove = cleanup_work.shutil.rmtree
    def broken(path, *args, **kwargs):
        if Path(path).name == 'tmp':
            raise OSError('temporary filesystem error')
        return remove(path, *args, **kwargs)
    monkeypatch.setattr(cleanup_work.shutil, 'rmtree', broken)
    with pytest.warns(RuntimeWarning, match='cleanup incomplete'):
        result = sample_review.apply_plan(root, cfg, plan)
    assert result['excluded_runs'] == ['A1'] and result['cleanup'][0]['state'] == 'pending'
    assert read_tsv(root / cfg['excluded_accessions'])[0]['accession'] == 'A1'
    assert temporary.exists() and not fastq.exists()
    monkeypatch.setattr(cleanup_work.shutil, 'rmtree', remove)
    repeated = sample_review.apply_plan(root, cfg, plan)
    assert repeated['cleanup_only'] and repeated['cleanup'][0]['state'] == 'complete'
    assert not temporary.exists()


def test_review_plan_and_inputs_cannot_be_changed_silently(sample_project, monkeypatch):
    root, cfg, build = sample_project
    fail(build, 1, monkeypatch)
    plan = root / 'work/review.tsv'
    catalog.record_successes(root, cfg, build, exclude_failed=True, plan=plan)
    with pytest.raises(ValueError, match='already exists'):
        catalog.record_successes(root, cfg, build, plan=plan)
    edit_plan(plan, {'G1': ('accept', 'not_finished')})
    with pytest.raises(ValueError, match='invalid review action'):
        sample_review.apply_plan(root, cfg, plan)
    edit_plan(plan, {'G1': ('hold', 'incomplete')})
    write_tsv(root / cfg['excluded_accessions'], ['accession', 'reason'], [{'accession': 'OTHER', 'reason': 'manual'}])
    with pytest.raises(ValueError, match='inputs changed'):
        sample_review.apply_plan(root, cfg, plan)
    assert len(read_tsv(root / cfg['excluded_accessions'])) == 1


def test_cli_saved_review_applies_only_selected_rows(sample_project, monkeypatch):
    root, cfg, build = sample_project
    fail(build, 1, monkeypatch)
    cli = [sys.executable, str(root / 'workflow/scripts/metadata_catalog.py'), 'record',
           '--config', str(cfg['_path'])]
    preview = subprocess.run([*cli, '--build', str(build), '--exclude-failed', '--plan', 'work/review.tsv'],
                             cwd=root, capture_output=True, text=True, check=True)
    assert json.loads(preview.stdout)['excluded_runs'] == ['A1', 'B1', 'G1']
    plan = root / 'work/review.tsv'
    rows = read_tsv(plan)
    write_tsv(plan, sample_review.FIELDS, rows[:1])
    applied = subprocess.run([*cli, '--apply', str(plan)], cwd=root, capture_output=True, text=True, check=True)
    assert json.loads(applied.stdout)['excluded_runs'] == ['A1']
    edit_plan(plan, {'A1': ('hold', 'changed_after_apply')})
    with pytest.raises(ValueError, match='applied review plan changed'):
        sample_review.apply_plan(root, cfg, plan)


def test_busco_rejection_cleans_failed_quant_reads_without_abandoning_other_failures(sample_project, monkeypatch):
    root, cfg, build = sample_project
    monkeypatch.setenv('FAKE_GG_LOW_BUSCO', '1')
    fail(build, 1, monkeypatch, 'QUANT')
    monkeypatch.delenv('FAKE_GG_LOW_BUSCO')
    _, fastq, _, log, _ = leftovers(build, 1)
    result = catalog.record_successes(root, cfg, build, runs=['A1'])
    assert result['excluded_runs'] == ['A1']
    assert result['samples'][0]['status'] == 'busco_below_threshold'
    assert not fastq.exists() and log.exists()
    assert dataset.status(build)[0]['busco'] == 'reuse'


def test_explicit_abandonment_cleans_a_debug_build_and_preserves_linked_reads(sample_project, monkeypatch):
    root, cfg, _ = sample_project
    build_config = root / cfg['build_config']
    settings = yaml.safe_load(build_config.read_text())
    settings['storage'] = {'keep_intermediates': True}
    build_config.write_text(yaml.safe_dump(settings))
    build = dataset.prepare(root, 'debug', build_config)
    dataset.submit(build, until='quant', dry_run=True)
    fail(build, 1, monkeypatch)
    work, fastq, _, _, temporary = leftovers(build, 1)
    runtime = work / f'.genegalleon-runtime-{os.getuid()}'
    assert runtime.is_dir()
    original = root / 'user_reads'
    original.mkdir()
    read = original / 'original.fastq.gz'
    read.write_bytes(b'caller supplied input')
    read_directory = fastq.parent.parent
    import shutil
    shutil.rmtree(read_directory)
    read_directory.symlink_to(original, target_is_directory=True)
    result = catalog.record_successes(root, cfg, build, runs=['A1'], exclude_failed=True)
    assert result['excluded_runs'] == ['A1']
    assert not runtime.exists() and not temporary.exists()
    assert read.read_bytes() == b'caller supplied input'


def test_changed_plan_identity_and_wrong_exclusion_table_are_rejected(sample_project, monkeypatch):
    root, cfg, build = sample_project
    fail(build, 1, monkeypatch)
    plan = root / 'work/review.tsv'
    catalog.record_successes(root, cfg, build, exclude_failed=True, plan=plan)
    rows = read_tsv(plan)
    rows[0]['species'] = 'another_sample'
    write_tsv(plan, sample_review.FIELDS, rows)
    with pytest.raises(ValueError, match='only action and reason'):
        sample_review.apply_plan(root, cfg, plan)
    changed = dict(cfg, excluded_accessions='datasets/another_exclusions.tsv')
    with pytest.raises(ValueError, match='same excluded_accessions'):
        catalog.record_successes(root, changed, build, exclude_failed=True)
    assert read_tsv(root / cfg['excluded_accessions']) == []
