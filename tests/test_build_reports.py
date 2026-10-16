"""Exercise verified plans, Slurm preparation, repair and progress reporting."""
import fcntl
import json
import shutil
import subprocess
from pathlib import Path

import pytest

import dataset
import dataset_assets
from build_dispatch import dispatch, runtime_data, update_verification, verify_and_stage
from build_reports import progress, read_context
from common import file_record, read_tsv, write_json, write_tsv
from database_reuse import PRODUCTS
from test_database_reuse import configure
from test_datasets import dataset_project, fake_genegalleon, imported, new_dataset
from test_submit_cli import invoke, scheduler
from test_portable_build import completed_project, execute


@pytest.fixture
def reuse_project(dataset_project):
    root = dataset_project
    source = imported(root)
    moved = root / 'resources/transferred/database'
    shutil.copytree(source, moved)
    configure(root, str(moved))
    return root, moved


def artifact_paths(source):
    from build_products import load_complete
    data = load_complete(source, verify_files=False)
    return {Path(product[key]['path']) for product in data['products'].values() for key in PRODUCTS}


def prohibit_artifact_hashes(monkeypatch, paths):
    hash_file = dataset_assets.sha256
    def checked(path):
        assert Path(path) not in paths, f'Unexpected artifact checksum: {path}'
        return hash_file(path)
    monkeypatch.setattr(dataset_assets, 'sha256', checked)


@pytest.mark.parametrize('options', [(), ('--verify',)])
def test_plan_verifies_products_by_default_and_saves_details(reuse_project, monkeypatch, capsys, options):
    root, source = reuse_project
    paths = artifact_paths(source)
    calls = []
    hash_file = dataset_assets.sha256
    def counted(path):
        calls.append(Path(path))
        return hash_file(path)
    monkeypatch.setattr(dataset_assets, 'sha256', counted)
    assert invoke(monkeypatch, dataset, 'plan', '--root', root, '--split', '--outdir', 'preview', *options) == 0
    output = capsys.readouterr().out
    assert paths <= set(calls)
    assert 'Reuse: 3' in output
    assert 'pending verification' not in output
    assert '\tspecies\t' not in output
    rows = read_tsv(root / 'preview/plan.tsv')
    assert {row['verification'] for row in rows} == {'verified'}
    assert len(read_tsv(root / 'preview/reuse.tsv')) == 3
    assert read_tsv(root / 'preview/needs_attention.tsv') == []
    assert not (root / 'results').exists()


@pytest.mark.parametrize('missing', [False, True])
def test_plan_verified_damage_becomes_rebuild(reuse_project, monkeypatch, capsys, missing):
    root, source = reuse_project
    damaged = next((source / 'cds').glob('*.gz'))
    if missing:
        damaged.unlink()
    else:
        damaged.write_bytes(damaged.read_bytes() + b'changed')
    assert invoke(monkeypatch, dataset, 'plan', '--root', root) == 0
    assert 'rebuild: 1' in capsys.readouterr().out
    rows = read_tsv(root / 'reports/plan/needs_attention.tsv')
    assert len(rows) == 1
    assert rows[0]['action'] == 'rebuild'
    assert rows[0]['verification'] == 'failed'


def test_plan_scientific_conflicts_still_block(reuse_project, monkeypatch, capsys):
    root, source = reuse_project
    configure(root, str(source), translation={'table': 2})
    assert invoke(monkeypatch, dataset, 'plan', '--root', root) == 1
    assert 'Needs attention: 3' in capsys.readouterr().out
    assert {row['action'] for row in read_tsv(root / 'reports/plan/needs_attention.tsv')} == {'blocked'}


def test_plan_rebuilds_damaged_json_before_interpreting_its_scientific_identity(reuse_project, monkeypatch):
    from build_products import load_complete
    root, source = reuse_project
    products = load_complete(source, verify_files=False)['products']
    product = next(iter(products.values()))
    qc = Path(product['expression_qc']['path'])
    changed = json.loads(qc.read_text())
    changed['expression']['sha256'] = '0' * 64
    write_json(qc, changed)
    assert invoke(monkeypatch, dataset, 'plan', '--root', root) == 0
    attention = read_tsv(root / 'reports/plan/needs_attention.tsv')
    assert len(attention) == 1
    assert attention[0]['action'] == 'rebuild'
    assert 'registered file changed' in attention[0]['reason']


def test_submit_only_schedules_verification_locally(reuse_project, scheduler, monkeypatch, capsys):
    root, source = reuse_project
    prohibit_artifact_hashes(monkeypatch, artifact_paths(source))
    monkeypatch.setattr(dataset, 'resolve_software', lambda *_: pytest.fail('Local software acquisition'))
    assert invoke(monkeypatch, dataset, 'submit', '--root', root, '--name', 'deferred') == 0
    build = root / 'results/deferred'
    assert not (build / 'work/cache').exists()
    assert (build / 'reuse_candidates.json').is_file()
    assert not (build / 'runtime.json').exists()
    assert len(scheduler['commands']) == 1
    assert '_verify' in next(arg for arg in scheduler['commands'][0] if arg.startswith('--job-name='))
    assert 'Verification job: 1001' in capsys.readouterr().out


def test_submit_prepares_plan_cache_hits_locally_without_verification_job(reuse_project, scheduler, monkeypatch, capsys):
    root, source = reuse_project
    assert invoke(monkeypatch, dataset, 'plan', '--root', root) == 0
    prohibit_artifact_hashes(monkeypatch, artifact_paths(source))
    assert invoke(monkeypatch, dataset, 'submit', '--root', root, '--name', 'cached', '--prepare-locally') == 0
    build = root / 'results/cached'
    assert len(scheduler['commands']) == 1
    assert '_mapping' in next(arg for arg in scheduler['commands'][0] if arg.startswith('--job-name='))
    assert not (build / 'jobs/0001_verify.sh').exists()
    assert json.loads((build / 'jobs/verification.json').read_text())['location'] == 'local'
    assert 'preparation completed locally' in capsys.readouterr().out
    assert {row['state'] for row in progress(build, read_context(build))['samples']} == {'reused'}


def test_cached_submit_with_new_samples_schedules_workers_directly(reuse_project, scheduler, monkeypatch):
    root, source = reuse_project
    fake_genegalleon(root)
    metadata = root / 'input/metadata.tsv'
    rows = read_tsv(metadata)
    rows.append(dict.fromkeys(rows[0], '') | {'scientific_name': 'Delta new', 'run': 'NEW1', 'taxid': '42'})
    write_tsv(metadata, list(rows[0]), rows)
    assert invoke(monkeypatch, dataset, 'plan', '--root', root) == 0
    prohibit_artifact_hashes(monkeypatch, artifact_paths(source))
    monkeypatch.setattr('dataset_software.resolve', lambda *_: pytest.fail('Software acquisition on cached submit'))
    assert invoke(monkeypatch, dataset, 'submit', '--root', root, '--name', 'new_sample', '--prepare-locally') == 0
    commands = scheduler['commands']
    assert len(commands) == 2
    assert '_sample' in next(arg for arg in commands[0] if arg.startswith('--job-name='))
    assert next(arg for arg in commands[0] if arg.startswith('--array=')).split('%')[0] == '--array=3'
    assert '--dependency=afterok:1001' in commands[1]
    assert not (root / 'results/new_sample/jobs/0001_verify.sh').exists()
    assert {p.name for p in (root / 'results/new_sample/work/genegalleon').iterdir()} == {'Delta_new_NEW1'}


def test_cached_products_with_uncached_software_still_use_slurm(reuse_project, scheduler, monkeypatch):
    root, source = reuse_project
    metadata = root / 'input/metadata.tsv'
    rows = read_tsv(metadata)
    rows.append(dict.fromkeys(rows[0], '') | {'scientific_name': 'Delta new', 'run': 'NEW1', 'taxid': '42'})
    write_tsv(metadata, list(rows[0]), rows)
    assert invoke(monkeypatch, dataset, 'plan', '--root', root) == 0
    prohibit_artifact_hashes(monkeypatch, artifact_paths(source))
    monkeypatch.setattr('dataset_software.download', lambda *_: pytest.fail('Local software download'))
    assert invoke(monkeypatch, dataset, 'submit', '--root', root, '--name', 'uncached_software', '--prepare-locally') == 0
    assert len(scheduler['commands']) == 1
    assert '_verify' in next(arg for arg in scheduler['commands'][0] if arg.startswith('--job-name='))
    assert not (root / 'results/uncached_software/runtime.json').exists()


def test_file_changed_after_cache_probe_defers_reads_to_slurm(reuse_project, scheduler, monkeypatch):
    import build_dispatch
    root, source = reuse_project
    assert invoke(monkeypatch, dataset, 'plan', '--root', root) == 0
    probe = build_dispatch.cached_preparation
    def changed(path, manifest, context):
        result = probe(path, manifest, context)
        assert result[0]
        cds = next((source / 'cds').glob('*.gz'))
        cds.write_bytes(cds.read_bytes() + b'changed after probe')
        return result
    monkeypatch.setattr(build_dispatch, 'cached_preparation', changed)
    prohibit_artifact_hashes(monkeypatch, artifact_paths(source))
    assert invoke(monkeypatch, dataset, 'submit', '--root', root, '--name', 'race', '--prepare-locally') == 0
    assert len(scheduler['commands']) == 1
    assert '_verify' in next(arg for arg in scheduler['commands'][0] if arg.startswith('--job-name='))


def test_cached_submit_defers_file_preparation_and_dispatch_keeps_cache_hits(reuse_project, scheduler, monkeypatch):
    import build_dispatch
    root, source = reuse_project
    fake_genegalleon(root)
    metadata = root / 'input/metadata.tsv'
    rows = read_tsv(metadata)
    rows.append(dict.fromkeys(rows[0], '') | {'scientific_name': 'Delta new', 'run': 'NEW1', 'taxid': '42'})
    write_tsv(metadata, list(rows[0]), rows)
    assert invoke(monkeypatch, dataset, 'plan', '--root', root) == 0
    prohibit_artifact_hashes(monkeypatch, artifact_paths(source))
    with monkeypatch.context() as local:
        local.setattr(build_dispatch, 'cached_preparation', lambda *_: pytest.fail('Local cache scan before sbatch'))
        local.setattr(build_dispatch, 'verify_and_stage', lambda *_, **__: pytest.fail('Local import before sbatch'))
        local.setattr(dataset, 'stage_workspace', lambda *_: pytest.fail('Local workspace preparation before sbatch'))
        assert invoke(local, dataset, 'submit', '--root', root, '--name', 'fast') == 0
    build = root / 'results/fast'
    assert len(scheduler['commands']) == 1
    assert not (build / 'work').exists()
    assert not (build / 'runtime.json').exists()
    assert json.loads((build / 'jobs/verification.json').read_text())['state'] == 'queued'
    # A second submit uses the saved snapshot and rejects active preparation.
    scheduler['running'] = True
    with pytest.raises(ValueError, match='queued/running'):
        invoke(monkeypatch, dataset, 'submit', '--root', root, '--build', build)
    scheduler['running'] = False
    monkeypatch.setenv('SLURM_JOB_ID', '1001')
    commands = dispatch(build)
    assert len(commands) == 2
    assert '--array=3' in commands[0]
    assert '--dependency=afterok:1001' in commands[0]
    assert '--dependency=afterok:1002' in commands[1]
    assert {p.name for p in (build / 'work/genegalleon').iterdir()} == {'Delta_new_NEW1'}
    assert sum(row['state'] == 'reused' for row in progress(build, read_context(build))['samples']) == 3


def test_completed_build_with_verified_cache_needs_no_slurm_job(completed_project, scheduler, monkeypatch):
    from build_products import load_complete
    from verification_cache import verification_cache
    root, build, _, _ = completed_project
    with verification_cache(root / '.cache/verification.sqlite'):
        completed = load_complete(build)
    prohibit_artifact_hashes(monkeypatch, {Path(entry['path']) for entry in completed['files']})
    assert invoke(monkeypatch, dataset, 'submit', '--root', root, '--build', build) == 0
    assert scheduler['commands'] == []


def test_dispatch_verifies_stages_and_schedules_dependents(reuse_project, scheduler, monkeypatch):
    root, source = reuse_project
    build = dataset.prepare(root, 'deferred', root / 'config/build.yaml', defer=True)
    before = (build / 'build.json').read_bytes()
    monkeypatch.setenv('SLURM_JOB_ID', '900')
    commands = dispatch(build)
    assert len(commands) == 1
    assert '--dependency=afterok:900' in commands[0]
    assert '--kill-on-invalid-dep=yes' in commands[0]
    assert all(row['assembly'] == 'reuse' for row in dataset.status(build))
    assert len(runtime_data(build)['verified_files']) == 3
    assert (build / 'build.json').read_bytes() == before
    # Once imported, resume never needs the source database.
    shutil.rmtree(source)
    verify_and_stage(build, dataset.load(build))
    report = progress(build, read_context(build))
    assert {row['state'] for row in report['samples']} == {'reused'}


def test_dispatch_cli_uses_frozen_submission_resources(reuse_project, scheduler, monkeypatch):
    root, _ = reuse_project
    resources = root / 'resources.yaml'
    resources.write_text('slurm:\n  per_job_resources:\n    controller:\n      cpus: 2\n      mem_gb: 16\n')
    invoke(monkeypatch, dataset, 'submit', '--root', root, '--name', 'resources', '--resources', resources)
    build = root / 'results/resources'
    execution = build / 'jobs/submission_0001.resources.json'
    monkeypatch.setenv('SLURM_JOB_ID', '1001')
    invoke(monkeypatch, dataset, 'dispatch', '--build', build, '--execution', execution)
    assert len(scheduler['commands']) == 2
    assert '--cpus-per-task=2' in scheduler['commands'][1]
    assert '--mem=16000M' in scheduler['commands'][1]
    assert '--dependency=afterok:1001' in scheduler['commands'][1]


def test_dispatch_rejects_changed_submission_resources(reuse_project, scheduler, monkeypatch):
    root, _ = reuse_project
    invoke(monkeypatch, dataset, 'submit', '--root', root, '--name', 'resources')
    build = root / 'results/resources'
    execution = build / 'jobs/submission_0001.resources.json'
    wrapper = json.loads(execution.read_text())
    wrapper['slurm']['stages']['controller']['cpus'] = 12
    write_json(execution, wrapper)
    with pytest.raises(ValueError, match='submission resources changed'):
        invoke(monkeypatch, dataset, 'dispatch', '--build', build, '--execution', execution)
    assert len(scheduler['commands']) == 1


def test_dispatch_repairs_only_damaged_candidate(reuse_project, scheduler, monkeypatch):
    root, source = reuse_project
    fake_genegalleon(root)
    assert invoke(monkeypatch, dataset, 'plan', '--root', root) == 0
    build = dataset.prepare(root, 'repair', root / 'config/build.yaml', defer=True)
    source_cds = source / 'cds/Alpha_plant_A1_longestCDS.fa.gz'
    source_cds.write_bytes(source_cds.read_bytes() + b'changed')
    monkeypatch.setenv('SLURM_JOB_ID', '900')
    commands = dispatch(build)
    assert len(commands) == 2
    assert next(arg for arg in commands[0] if arg.startswith('--array=')).split('%')[0] == '--array=1'
    assert '--dependency=afterok:900' in commands[0]
    assert '--dependency=afterok:1001' in commands[1]
    assert set(runtime_data(build)['rebuilds']) == {'Alpha_plant_A1'}
    assert [row['assembly'] for row in dataset.status(build)] == ['pending', 'reuse', 'reuse']
    dataset.worker(build, 1)
    assert {row['assembly'] for row in dataset.status(build)} == {'reuse'}


def test_retry_invalidates_only_bad_imports_and_downstream_work(reuse_project):
    root, _ = reuse_project
    fake_genegalleon(root)
    build = dataset.prepare(root, 'repair', root / 'config/build.yaml', defer=True)
    verify_and_stage(build, dataset.load(build))
    saved = runtime_data(build)
    files = saved['verified_files']['Alpha_plant_A1']
    missing = next(Path(entry['path']) for entry in files if Path(entry['path']).name == 'abundance.tsv.gz')
    missing.unlink()
    (build / 'work/input').mkdir()
    (build / 'work/input/old.txt').write_text('stale')
    rebuilds = verify_and_stage(build, dataset.load(build))
    assert set(rebuilds) == {'Alpha_plant_A1'}
    assert not (build / 'work/input').exists()
    assert (build / 'work/quarantine/Alpha_plant_A1').is_dir()
    assert set(runtime_data(build)['verified_files']) == {'Beta_sp-X_B1', 'Gamma_plant_G1'}
    snapshot = json.loads((build / 'work/cache/odb/v12_3193/reused/snapshot.json').read_text())
    assert set(snapshot['tables']) == {'Beta_sp-X_B1', 'Gamma_plant_G1'}
    assert len(list((build / 'work/cache/expression').glob('*/receipt.json'))) == 2


def test_status_default_reads_only_progress_and_verify_detects_damage(reuse_project, monkeypatch, capsys):
    root, source = reuse_project
    build = dataset.prepare(root, 'progress', root / 'config/build.yaml', defer=True)
    verify_and_stage(build, dataset.load(build))
    data = runtime_data(build)
    Path(data['verified_files']['Alpha_plant_A1'][0]['path']).unlink()
    prohibit_artifact_hashes(monkeypatch, artifact_paths(source))
    assert invoke(monkeypatch, dataset, 'status', '--build', build, '--outdir', root / 'progress-report') == 0
    assert 'Samples: 3/3 finished' in capsys.readouterr().out
    assert len(read_tsv(root / 'progress-report/status.tsv')) == 3
    assert invoke(monkeypatch, dataset, 'status', '--build', build, '--verify', '--json') == 1
    report = json.loads(capsys.readouterr().out)
    assert report['coverage'] == 'content-verified'
    assert sum(row['verification'] == 'failed' for row in report['samples']) == 1


def test_status_identifies_scheduler_failure_without_completion(reuse_project, monkeypatch):
    root, _ = reuse_project
    build = dataset.prepare(root, 'progress', root / 'config/build.yaml', defer=True)
    write_json(build / 'jobs/verification.json', {'state': 'queued', 'job_id': '123'})
    monkeypatch.setattr('genegalleon_run.scheduler_states', lambda jobs: {'123': {'state': 'OUT_OF_MEMORY'}})
    report = progress(build, read_context(build))
    assert report['verification']['state'] == 'failed'
    assert {row['state'] for row in report['samples']} == {'blocked'}
    assert 'OUT_OF_MEMORY' in report['samples'][0]['reason']
    assert report['verification']['log'] == ''
    assert all(row['log'] == '' for row in report['samples'])


def test_status_reports_each_samples_latest_array_log_files(dataset_project, scheduler, monkeypatch, capsys):
    build = new_dataset(dataset_project, ('Alpha new', 'Beta new', 'Gamma new'))
    monkeypatch.setattr('phase_config.resolve_array_size', lambda: 2)
    dataset.submit(build, until='quant')
    items = dataset.load(build)['items']
    pilot = dataset_project / 'retry.txt'
    pilot.write_text(items[1]['species'] + '\n')
    dataset.submit(build, until='quant', species=pilot)
    # Only Beta was retried. Gamma retains task 1 of the second array batch.
    expected = [('0001_sample_1001_1', '1001_1'),
                ('0002_sample_1003_2', '1003_2'),
                ('0001_sample_3_4_1002_1', '1002_1')]
    for item, (_, identifier), stage in zip(items, expected, ('assembly', 'busco', 'quant')):
        write_json(build / 'jobs/status' / f'{item["species"]}.sample.json',
                   {'state': 'failed', 'stage': stage, 'error': 'native command failed',
                    'slurm_id': identifier})
    monkeypatch.setattr('genegalleon_run.scheduler_states', lambda jobs: {})
    monkeypatch.setattr(subprocess, 'run', lambda cmd, **kwargs: subprocess.CompletedProcess(cmd, 0, stdout=''))
    capsys.readouterr()
    output_dir = dataset_project / 'progress-report'
    assert invoke(monkeypatch, dataset, 'status', '--build', build, '--outdir', output_dir) == 1
    output = capsys.readouterr().out
    exported = read_tsv(output_dir / 'failed.tsv')
    assert len(exported) == 3
    for row, (stem, identifier) in zip(exported, expected):
        stderr = str(build / 'jobs/logs' / f'{stem}.err')
        stdout = str(build / 'jobs/logs' / f'{stem}.out')
        assert row['job_id'] == identifier
        assert row['log'] == row['stderr'] == stderr
        assert row['stdout'] == stdout
        assert f'    Stderr: {stderr}\n    Stdout: {stdout}\n' in output
        assert not Path(stderr).exists()  # Reporting never depends on log creation/content.
    assert '    Logs: ' not in output
    assert invoke(monkeypatch, dataset, 'status', '--build', build, '--json') == 1
    report = json.loads(capsys.readouterr().out)
    assert [row['stderr'] for row in report['samples']] == [row['stderr'] for row in exported]
    assert [row['stdout'] for row in report['samples']] == [row['stdout'] for row in exported]


@pytest.mark.parametrize('controller,stem', [('verification', '0001_verify_1001'),
                                             ('database', '0001_mapping_1001_4294967294')])
def test_status_reports_cancelled_controller_log_files(reuse_project, scheduler, monkeypatch, capsys,
                                                     controller, stem):
    root, _ = reuse_project
    build = dataset.prepare(root, 'cancelled', root / 'config/build.yaml', defer=controller == 'verification')
    dataset.submit(build, defer=controller == 'verification')
    monkeypatch.setattr('genegalleon_run.scheduler_states', lambda jobs: {'1001': {'state': 'CANCELLED'}})
    capsys.readouterr()
    assert invoke(monkeypatch, dataset, 'status', '--build', build) == 1
    output = capsys.readouterr().out
    stderr = str(build / 'jobs/logs' / f'{stem}.err')
    stdout = str(build / 'jobs/logs' / f'{stem}.out')
    assert f'  {controller}: Slurm CANCELLED; completion was not recorded\n' in output
    assert f'    Stderr: {stderr}\n    Stdout: {stdout}\n' in output
    report = progress(build, read_context(build))
    assert report[controller]['stderr'] == stderr
    assert report[controller]['stdout'] == stdout
    if controller == 'verification':
        assert {row['state'] for row in report['samples']} == {'blocked'}
        assert all(row['log'] == row['stderr'] == stderr for row in report['samples'])


def test_status_explains_when_submission_log_paths_are_unavailable(reuse_project, monkeypatch, capsys):
    from build_reports import print_progress
    root, _ = reuse_project
    build = dataset.prepare(root, 'no-submission', root / 'config/build.yaml', defer=True)
    write_json(build / 'jobs/verification.json', {'state': 'failed', 'error': 'preparation failed locally'})
    print_progress(progress(build, read_context(build)))
    output = capsys.readouterr().out
    assert 'Logs: unavailable (no recorded log file paths)' in output
    assert str(build / 'jobs/logs') not in output


def test_changed_shared_manifest_stops_dispatch(reuse_project, scheduler):
    root, source = reuse_project
    build = dataset.prepare(root, 'immutable', root / 'config/build.yaml', defer=True)
    (source / 'manifest.json').write_text('{}')
    with pytest.raises(ValueError, match='registered file changed'):
        dispatch(build)
    assert scheduler['commands'] == []
    assert json.loads((build / 'jobs/verification.json').read_text())['state'] == 'failed'


def test_interrupted_initial_staging_can_be_retried(reuse_project):
    root, _ = reuse_project
    build = dataset.prepare(root, 'interrupted', root / 'config/build.yaml', defer=True)
    (build / 'work/cache').mkdir(parents=True)
    (build / 'work/cache/interrupted.txt').write_text('unfinished')
    verify_and_stage(build, dataset.load(build))
    assert not (build / 'work/cache/interrupted.txt').exists()
    assert (build / 'work/quarantine/unfinished_import/cache/interrupted.txt').is_file()
    assert {row['assembly'] for row in dataset.status(build)} == {'reuse'}


def test_old_builds_have_light_status_fallback(reuse_project, monkeypatch, capsys):
    root, _ = reuse_project
    build = dataset.prepare(root, 'legacy', root / 'config/build.yaml')
    (build / 'status_context.json').unlink()
    assert invoke(monkeypatch, dataset, 'status', '--build', build) == 0
    assert 'Samples: 3/3 finished' in capsys.readouterr().out


def test_help_lists_public_workflow_commands(monkeypatch, capsys):
    with pytest.raises(SystemExit) as result:
        invoke(monkeypatch, dataset, '--help')
    assert result.value.code == 0
    output = capsys.readouterr().out
    assert 'plan' in output and 'submit' in output and 'status' in output
    assert all(name not in output for name in ('dispatch', 'materialize', 'worker', 'prepare'))


def test_deferred_build_publishes_a_portable_database(completed_project):
    from build_products import complete, load_complete
    root, source, env, events = completed_project
    configure(root, str(source / 'database'))
    build = dataset.prepare(root, 'deferred-publication', root / 'config/build.yaml', defer=True)
    prepared = (build / 'build.json').read_bytes()
    verify_and_stage(build, dataset.load(build))
    dataset.materialize(build)
    execute(root, build, 'database', env, mapping=True)
    complete(build)
    assert load_complete(build)['build_id'] == 'deferred-publication'
    assert (build / 'build.json').read_bytes() == prepared
    provenance = build / 'database/provenance'
    assert (provenance / 'prepared_build.json').read_bytes() == prepared
    hashes = json.loads((provenance / 'checksums.json').read_text())
    assert hashes['build.json'] == file_record(provenance / 'build.json')['sha256']
    report = progress(build, read_context(build))
    from build_reports import verify_progress
    assert verify_progress(build, report)['database']['state'] == 'complete'
    assert {row['verification'] for row in report['samples']} == {'verified'}
    assert not events.exists()  # All mapping tables came from the source database.


def test_new_private_sample_does_not_hash_reads(dataset_project, scheduler, monkeypatch):
    root = dataset_project
    fake_genegalleon(root)
    reads = root / 'input/private.fastq'
    reads.write_text('@r1\nACGT\n+\nIIII\n')
    metadata = root / 'input/private.tsv'
    row = {'scientific_name': 'Private plant', 'run': 'LOCAL1', 'taxid': '42',
           'private_file': 'yes', 'lib_layout': 'single', 'read1_path': reads.name}
    write_tsv(metadata, list(row), [row])
    prohibit_artifact_hashes(monkeypatch, {reads})
    assert invoke(monkeypatch, dataset, 'plan', '--root', root, '--metadata', metadata) == 0
    invoke(monkeypatch, dataset, 'submit', '--root', root, '--name', 'private', '--metadata', metadata)
    build = root / 'results/private'
    dispatch(build, until='quant')
    assert dataset.load(build)['raw_inputs']['Private_plant_LOCAL1']['read1_path'] == {'path': str(reads)}
    dataset.worker(build, 1)
    assert dataset.status(build)[0]['quant'] == 'reuse'


def test_private_reads_required_for_known_rebuild_are_reported(dataset_project, monkeypatch, capsys):
    root = dataset_project
    metadata = root / 'input/private.tsv'
    row = {'scientific_name': 'Private plant', 'run': 'LOCAL1', 'taxid': '42',
           'private_file': 'yes', 'lib_layout': 'single', 'read1_path': 'missing.fastq'}
    write_tsv(metadata, list(row), [row])
    assert invoke(monkeypatch, dataset, 'plan', '--root', root, '--metadata', metadata) == 1
    capsys.readouterr()
    assert 'FASTQ' in read_tsv(root / 'reports/plan/needs_attention.tsv')[0]['reason']


def test_every_sample_batch_waits_for_successful_preparation(dataset_project, scheduler, monkeypatch):
    root = dataset_project
    fake_genegalleon(root)
    metadata = root / 'input/new.tsv'
    rows = [{'scientific_name': name, 'run': f'SRR{i}', 'taxid': '42'}
            for i, name in enumerate(('Alpha new', 'Beta new', 'Gamma new'), 1)]
    write_tsv(metadata, list(rows[0]), rows)
    build = dataset.prepare(root, 'batches', root / 'config/build.yaml', metadata, defer=True)
    monkeypatch.setenv('SLURM_JOB_ID', '900')
    monkeypatch.setattr('phase_config.resolve_array_size', lambda: 2)
    commands = dispatch(build)
    assert len(commands) == 3
    assert '--dependency=afterok:900' in commands[0]
    assert '--dependency=afterany:1001,afterok:900' in commands[1]
    assert '--kill-on-invalid-dep=yes' in commands[1]
    assert '--dependency=afterok:1001:1002' in commands[2]


def test_repair_preserves_work_held_by_an_active_native_lock(reuse_project):
    root, _ = reuse_project
    fake_genegalleon(root)
    build = dataset.prepare(root, 'guarded', root / 'config/build.yaml', defer=True)
    verify_and_stage(build, dataset.load(build))
    files = runtime_data(build)['verified_files']['Alpha_plant_A1']
    Path(files[0]['path']).unlink()
    work = build / 'work/genegalleon/Alpha_plant_A1'
    work.mkdir(parents=True)
    (work / 'retained.txt').write_text('active work')
    with (work / '.native.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(BlockingIOError):
            verify_and_stage(build, dataset.load(build))
        assert (work / 'retained.txt').read_text() == 'active work'
        assert not (build / 'work/quarantine').exists()


def test_progress_uses_new_verified_products_after_a_previous_job_failure(reuse_project, monkeypatch):
    import subprocess
    root, _ = reuse_project
    build = dataset.prepare(root, 'recovered', root / 'config/build.yaml', defer=True)
    write_json(build / 'jobs/submission_0001.json', {'jobs': [
        {'stage': 'sample', 'indices': [1], 'job_id': '111', 'state': 'submitted'}]})
    monkeypatch.setattr('genegalleon_run.scheduler_states', lambda jobs: {'111_1': {'state': 'FAILED'}})
    monkeypatch.setattr(subprocess, 'run', lambda command, **kwargs: subprocess.CompletedProcess(command, 0, '', ''))
    verify_and_stage(build, dataset.load(build))
    report = progress(build, read_context(build))
    assert {row['state'] for row in report['samples']} == {'reused'}


def test_verification_progress_preserves_job_identity_and_terminal_state(reuse_project):
    root, _ = reuse_project
    build = dataset.prepare(root, 'identity', root / 'config/build.yaml', defer=True)
    update_verification(build, 'running', job_id='123', started_at='2026-10-09T00:00:00+00:00')
    update_verification(build, 'running', checked=1, total=3, current_sample='Alpha_plant_A1')
    current = json.loads((build / 'jobs/verification.json').read_text())
    assert current['job_id'] == '123'
    update_verification(build, 'complete', finished_at='2026-10-09T00:01:00+00:00')
    update_verification(build, 'queued', job_id='123')
    assert json.loads((build / 'jobs/verification.json').read_text())['state'] == 'complete'
