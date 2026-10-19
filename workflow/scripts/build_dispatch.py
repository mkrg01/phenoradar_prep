"""Prepare cached reuse locally, or verify on Slurm, then schedule sample work."""
import copy
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import time

from common import now, write_json
from odb_settings import cache_namespace
from dataset_assets import (ProductIntegrityError, UncachedVerification, cached_verification,
                            digest, locked, record, verification_cached, verify, verify_many)


def save_candidates(directory, selection):
    from database_reuse import SAMPLE_FIELDS
    selected, locks = {}, {}
    for name, entry in selection['reusable'].items():
        chosen = dict(entry)
        source = chosen['source']['path']
        locks[source] = chosen.pop('software_lock', None)
        if chosen.get('kind') == 'stages':
            chosen['row'] = {key: chosen['row'].get(key, '') for key in SAMPLE_FIELDS}
        else:
            chosen['product'] = dict(chosen['product'])
            chosen['product']['row'] = {key: chosen['product']['row'].get(key, '') for key in SAMPLE_FIELDS}
        selected[name] = chosen
    write_json(Path(directory) / 'reuse_candidates.json',
               {'samples': selected, 'sources': selection['reuse_sources'],
                'software_locks': locks, 'rebuilds': selection.get('rebuilds', {})}, compact=True)


def runtime_data(path):
    filename = Path(path) / 'runtime.json'
    if not filename.is_file():
        return None
    data = json.loads(filename.read_text())
    if data.get('schema_version') != 1 or digest({k: v for k, v in data.items() if k != 'sha256'}) != data.get('sha256'):
        raise ValueError('build runtime record changed')
    return data


def apply_runtime(path, manifest):
    data = runtime_data(path)
    if data is None:
        return
    expected = json.loads((Path(path) / 'checksums.json').read_text())['build.json']
    if data['prepared_sha256'] != expected:
        raise ValueError('runtime record belongs to a different build snapshot')
    manifest['config']['genegalleon'] = data['genegalleon_config']
    manifest['genegalleon'] = data['genegalleon']
    manifest['software_lock'] = data['software_lock']
    manifest['raw_inputs'] = data['raw_inputs']
    manifest['busco_reference'] = data.get('busco_reference', manifest.get('busco_reference'))
    manifest['taxonomy_reference'] = data.get('taxonomy_reference', manifest.get('taxonomy_reference'))
    for item in manifest['items']:
        if item['species'] in data['references']:
            item['reference_id'] = data['references'][item['species']]


def save_runtime(path, manifest, references, raw_inputs, files, rebuilds):
    snapshot = Path(path) / 'work/cache/odb' / cache_namespace(manifest["config"]["odb"]) / 'reused/snapshot.json'
    data = {'schema_version': 1,
            'prepared_sha256': json.loads((Path(path) / 'checksums.json').read_text())['build.json'],
            'genegalleon_config': manifest['config']['genegalleon'], 'genegalleon': manifest['genegalleon'],
            'software_lock': manifest['software_lock'], 'references': references,
            'busco_reference': manifest.get('busco_reference'),
            'taxonomy_reference': manifest.get('taxonomy_reference'),
            'raw_inputs': raw_inputs, 'verified_files': files, 'rebuilds': rebuilds,
            'shared_files': [record(snapshot)] if snapshot.is_file() else []}
    data['sha256'] = digest(data)
    write_json(Path(path) / 'runtime.json', data)


def update_verification(path, state, **values):
    filename = Path(path) / 'jobs/verification.json'
    previous = json.loads(filename.read_text()) if filename.is_file() else {}
    # A short preparation job may finish before sbatch returns to its caller.
    if state == 'queued' and previous.get('job_id') == values.get('job_id') and previous.get('state') in {'running', 'complete', 'failed'}:
        return
    if state == 'queued' or 'started_at' in values:
        previous = {}
    write_json(filename, dict(previous, state=state, updated_at=now(), **values))


def quarantine(path, candidate, label):
    from cleanup_work import owned_path
    root = Path(path)
    checked = owned_path(root, Path(candidate).relative_to(root))
    if not checked.exists() and not checked.is_symlink():
        return
    destination = owned_path(root, Path('work/quarantine') / label)
    if destination.is_symlink():
        raise ValueError(f'unsafe quarantine directory: {destination}')
    destination.mkdir(parents=True, exist_ok=True)
    suffix = 1
    target = destination / checked.name
    while target.exists() or target.is_symlink():
        target = destination / f'{checked.name}_{suffix}'
        suffix += 1
    checked.rename(target)
    return target


def quarantine_sample(path, manifest, item, entries=()):
    from dataset import workspace
    name = item['species']
    root = Path(path)
    candidates = [Path(manifest['config']['store']) / name, workspace(path, manifest, item)]
    reference = item.get('reference_id')
    if reference:
        candidates.append(root / 'work/cache/proteins' / reference)
    expression = root / 'work/cache/expression'
    candidates.extend({Path(entry['path']).parent for entry in entries
                       if Path(entry['path']).is_relative_to(expression)})
    for candidate in candidates:
        quarantine(path, candidate, name)
    for stage in ('assembly', 'busco', 'quant', 'sample'):
        receipt = root / 'jobs/status' / f'{name}.{stage}.json'
        quarantine(path, receipt, name)


def invalidate_mapping(path, manifest, rebuilds):
    """Keep healthy imported tables, and discard downstream work using bad inputs."""
    from mapping_tables import load_tables, relative_file, write_tables
    path = Path(path)
    directory = path / 'work/cache/odb' / cache_namespace(manifest["config"]["odb"]) / 'reused'
    snapshot = directory / 'snapshot.json'
    if snapshot.is_file():
        data = load_tables(snapshot, verify_files=False)
        old = quarantine(path, directory, 'mapping')
        entries = {name: dict(entry, table=relative_file(old, entry['table']))
                   for name, entry in data['tables'].items() if name not in rebuilds}
        if entries:
            write_tables(directory, entries, data['version'], data['node'],
                         references=data.get('reference_sha256s', []))
    for relative in ('work/input', 'input_receipt.json', 'work/database', 'work/mapping'):
        quarantine(path, path / relative, 'mapping')


def candidate_files(chosen):
    if chosen.get('kind') == 'stages':
        products = chosen['products']
        result = [products['reference']['cds']]
        for key in ('full', 'short'):
            if (products.get('busco') or {}).get(key):
                result.append(products['busco'][key])
        if products.get('quant'):
            result.append(products['quant']['abundance'])
        return result
    from database_reuse import PRODUCTS, UPSTREAM
    if chosen.get('reuse_mapping') is False:
        return [chosen['product'][key] for key in UPSTREAM]
    return [chosen['product'][key] for key in PRODUCTS] + [chosen['mapping']['table']]


def verify_and_stage(path, manifest, workers=4, software=None, cached_only=False):
    """Validate frozen candidates before importing any of their products."""
    from database_reuse import stage
    from dataset import item_products, private_fastq_path
    from dataset_software import resolve as resolve_software
    from build_reports import read_context
    path = Path(path)
    previous = runtime_data(path)
    planned = {row['species']: row for row in read_context(path)['items']}
    references, rebuilds = {}, {}
    local_files = (previous or {}).get('verified_files', {})
    last_update = 0

    def progress(checked, total, phase='checksums'):
        nonlocal last_update
        current = time.monotonic()
        if current - last_update >= 1 or checked == total:
            update_verification(path, 'running', phase=phase, checked=checked,
                                total=total, current_sample=None)
            last_update = current

    def failures(samples):
        entries = [(name, entry) for name, files in samples.items() for entry in files]
        result = {}
        for (name, _), error in zip(entries, verify_many((entry for _, entry in entries), workers, progress)):
            if error is not None:
                result.setdefault(name, str(error))
        return result

    for entry in (previous or {}).get('shared_files', []):
        verify(entry)
    if previous is None and manifest.get('deferred_verification'):
        saved = json.loads(verify(manifest['auxiliary']['reuse_candidates.json']).read_text())
        for source in saved['sources']:
            verify(source)  # A changed shared manifest is not a per-sample repair.
        selected = {}
        rebuilds.update(saved['rebuilds'])
        damaged = failures({name: candidate_files(chosen) for name, chosen in saved['samples'].items()})
        rebuilds.update(damaged)
        for name, chosen in saved['samples'].items():
            if name not in damaged:
                chosen['software_lock'] = saved['software_locks'].get(chosen['source']['path'])
                selected[name] = chosen
        update_verification(path, 'running', phase='staging', checked=0,
                            total=len(selected), current_sample=None)
        last_update = 0
        staging = Path(tempfile.mkdtemp(prefix='.verify-', dir=path))
        try:
            local_files = stage(staging, path, selected, manifest['config'],
                                progress=lambda checked, total: progress(checked, total, 'staging'))
            incoming = staging / 'work/cache'
            from cleanup_work import owned_path
            destination = owned_path(path, 'work/cache')
            if incoming.exists():
                if destination.exists():
                    if any(destination.iterdir()):
                        quarantine(path, destination, 'unfinished_import')
                    else:
                        destination.rmdir()
                destination.parent.mkdir(parents=True, exist_ok=True)
                incoming.rename(destination)
        finally:
            shutil.rmtree(staging)
    else:
        rebuilds.update(failures(local_files))
    pending = False
    needs_busco = False
    raw_inputs = copy.deepcopy(manifest['raw_inputs'])
    for item in manifest['items']:
        name = item['species']
        if name not in rebuilds:
            try:
                products = item_products(manifest, item)
            except ProductIntegrityError as error:
                rebuilds[name] = str(error)
        if name in rebuilds:
            if item['row'].get('reference_id'):
                raise ValueError(f'{name}: damaged reference is explicitly pinned in metadata; remove reference_id for a new build')
            from sample_review import sample_guard
            with sample_guard(path, manifest, item, native=True):
                quarantine_sample(path, manifest, item, local_files.get(name, []))
            item['reference_id'] = ''
            references[name] = ''
            local_files.pop(name, None)
            products = {'reference': None, 'busco': None, 'quant': None}
        else:
            reference = products['reference']
            references[name] = reference['reference_id'] if reference else ''
            item['reference_id'] = references[name]
        needs_reads = not products['reference'] or not products['quant']
        if needs_reads and item['row'].get('private_file', '').lower() == 'yes':
            raw_inputs[name] = {key: {'path': str(private_fastq_path(item['row'][key]))}
                                for key in ('read1_path', 'read2_path') if item['row'].get(key)}
        if any(not products[key] for key in ('reference', 'busco', 'quant')):
            pending = True
        else:
            receipt = path / 'jobs/status' / f'{name}.sample.json'
            completed_here = receipt.is_file() and json.loads(receipt.read_text()).get('state') == 'complete'
            reused = planned[name]['action'] == 'reuse' and not completed_here and name not in (previous or {}).get('rebuilds', {})
            write_json(receipt, {'state': 'reused' if reused else 'complete',
                                 'verification': 'verified', 'finished_at': now()})
        if not products['busco']:
            needs_busco = True
    if pending:
        if cached_only and software is None:
            from dataset_software import resolve_cached
            software = resolve_cached(manifest['config']['genegalleon'])
            if software is None:
                raise UncachedVerification('uncached software preparation')
        if software is not None:
            for entry in software[1]:
                verify(entry)
        manifest['config']['genegalleon'], manifest['genegalleon'], manifest['software_lock'] = (
            software if software is not None else resolve_software(manifest['config']['genegalleon']))
        from taxonomy_reference import prepare as prepare_taxonomy
        manifest['taxonomy_reference'] = prepare_taxonomy(manifest['config'].get('taxonomy', {}), manifest['config']['genegalleon'],
                                                         expected=manifest.get('taxonomy_reference'), cached_only=cached_only)
        if needs_busco:
            from busco_reference import prepare as prepare_busco
            manifest['busco_reference'] = prepare_busco(manifest['config']['busco'], manifest['config']['genegalleon'],
                                                       expected=manifest.get('busco_reference'), cached_only=cached_only)
    if rebuilds and previous is not None:
        invalidate_mapping(path, manifest, rebuilds)
    save_runtime(path, manifest, references, raw_inputs, local_files, rebuilds)
    return rebuilds


def dispatch(path, until='database', species=None, execution=None):
    from build_products import load_complete
    from dataset import load, load_execution, submit
    from verification_cache import verification_cache
    path = Path(path).resolve()
    with verification_cache(path.parent.parent / '.cache/verification.sqlite'), locked(path / '.verification.lock'):
        update_verification(path, 'running', started_at=now(), job_id=os.environ.get('SLURM_JOB_ID'))
        try:
            manifest = load(path, check_code=True)
            from phase_config import execution_settings
            slurm = load_execution(execution) if execution else manifest['config']['slurm']
            slurm = execution_settings(slurm)
            if (path / 'completed.json').exists():
                load_complete(path)
                update_verification(path, 'complete', finished_at=now())
                print('Database already complete; recorded products verified.')
                return []
            rebuilds = verify_and_stage(path, manifest, workers=min(4, slurm['stages']['controller']['cpus']))
            update_verification(path, 'running', phase='scheduling')
            commands = submit(path, until, species, resources=execution,
                              dispatch_dependency=os.environ.get('SLURM_JOB_ID'))
            update_verification(path, 'complete', phase='finished', finished_at=now(), rebuilds=rebuilds,
                                job_id=os.environ.get('SLURM_JOB_ID'))
            print(f'Verification complete; rebuilding {len(rebuilds)} samples.')
            return commands
        except BaseException as error:
            update_verification(path, 'failed', error=str(error), job_id=os.environ.get('SLURM_JOB_ID'))
            raise


def save_database_submission(path, batch):
    for job in batch['jobs']:
        if job['stage'] == 'mapping':
            write_json(Path(path) / 'jobs/database.json',
                       {'state': 'queued', 'job_id': job['job_id'], 'updated_at': now()})


def check_submissions(jobs):
    """Reject active or ambiguous submissions before changing imported products."""
    ids = []
    for receipt in jobs.glob('submission_*.json'):
        if receipt.name.endswith('.resources.json'):
            continue
        for job in json.loads(receipt.read_text())['jobs']:
            if job['state'] in {'submitting', 'unknown'}:
                raise ValueError(f'unresolved submission in {receipt}; inspect Slurm before retrying')
            if job.get('job_id'):
                ids.append(job['job_id'])
    if ids:
        import pwd
        queued = subprocess.check_output(['squeue', '--noheader', '--user', pwd.getpwuid(os.getuid()).pw_name,
                                          '--format=%i'], text=True).splitlines()
        if any(identifier.strip().split('_', 1)[0] in ids for identifier in queued):
            raise ValueError('dataset still has queued/running jobs; wait before resubmitting')


def cached_preparation(path, manifest, context):
    """Stat-only eligibility; do not start importing before all candidates agree."""
    from dataset_software import resolve_cached
    previous = runtime_data(path)
    if (path / 'completed.json').exists():
        from build_products import load_complete
        completed = load_complete(path, verify_files=False)
        return all(verification_cached(entry) for entry in completed['files']), None
    if previous is not None:
        entries = [*previous.get('shared_files', []),
                   *(entry for files in previous.get('verified_files', {}).values() for entry in files)]
    elif manifest.get('deferred_verification'):
        saved = json.loads(verify(manifest['auxiliary']['reuse_candidates.json']).read_text())
        products = [entry for chosen in saved['samples'].values() for entry in candidate_files(chosen)]
        # Avoid a large cross-filesystem copy on the submit host.
        device = path.stat().st_dev
        if any(not Path(entry['path']).is_file() or Path(entry['path']).stat().st_dev != device for entry in products):
            return False, None
        entries = [*saved['sources'], *products]
    else:
        entries = []
    if not all(verification_cached(entry) for entry in entries):
        return False, None
    software = None
    if any(row['action'] in {'new', 'partial', 'rebuild'} for row in context['items']):
        software = resolve_cached(manifest['config']['genegalleon'])
        if software is None:
            return False, None
    return True, software


def submit_verification(path, until, species=None, dry_run=False, resources=None, prepare_locally=False):
    """Submit preparation to Slurm, optionally preparing verified cache hits locally."""
    from build_reports import print_plan, read_context
    from dataset import load, submit
    from phase_config import array_concurrency, execution_settings, worker_budgets
    path = Path(path).resolve()
    context = read_context(path)
    # Context contains no data files; loading the frozen configuration does not
    # verify reusable artifacts.
    manifest = load(path, check_code=True)
    slurm = execution_settings(manifest['config']['slurm'], resources)
    if any(row['action'] in {'new', 'partial', 'rebuild'} for row in context['items']):
        array_concurrency(slurm, 'sample')
    if until in {'mapping', 'database'} and species is None:
        worker_budgets(slurm)
    if species is not None:
        from sample_identity import select_samples
        requested = set(Path(species).read_text().splitlines())
        if not requested:
            raise ValueError('pilot species must be present in frozen metadata')
        select_samples([dict(item['row'], species=item['species']) for item in manifest['items']], requested)
    jobs = path / 'jobs'
    jobs.mkdir(exist_ok=True)
    (jobs / 'logs').mkdir(exist_ok=True)
    with locked(path / '.verification.lock'):
        # A cache hit still involves thousands of links and receipts. Keep that
        # work off the submit host unless local preparation was requested.
        if not prepare_locally and not (path / 'completed.json').exists():
            return submit_verification_job(path, manifest, context, slurm, until, species, dry_run)
        with locked(jobs / '.submit.lock'):
            if not dry_run:
                check_submissions(jobs)
        with cached_verification():
            try:
                ready, software = cached_preparation(path, manifest, context)
            except UncachedVerification:
                ready, software = False, None
        if ready:
            update_verification(path, 'running', started_at=now(), job_id=None, location='local')
            try:
                with cached_verification():
                    if (path / 'completed.json').exists():
                        from build_products import load_complete
                        load_complete(path)
                        rebuilds = {}
                    else:
                        rebuilds = verify_and_stage(path, manifest, software=software, cached_only=True)
                    update_verification(path, 'running', phase='scheduling', location='local')
                    # Status and workspace checks also precede sbatch and must
                    # fall back to Slurm if a file changes after preparation.
                    commands = submit(path, until, species, dry_run, resources)
            except UncachedVerification:
                # A file changed after the stat-only check. Let Slurm reverify it.
                pass
            except BaseException as error:
                update_verification(path, 'failed', error=str(error), job_id=None)
                raise
            else:
                update_verification(path, 'complete', phase='finished', finished_at=now(), rebuilds=rebuilds,
                                    job_id=None, location='local')
                print_plan([dict(row, verification='verified') if row['verification'] == 'pending' else row
                            for row in context['items']])
                print('Product checksums matched verified receipts/cache; preparation completed locally.')
                print(f'Status: ./run_build.sh status --build {path}')
                return commands
        return submit_verification_job(path, manifest, context, slurm, until, species, dry_run)


def submit_verification_job(path, manifest, context, slurm, until, species, dry_run):
    from build_reports import print_plan
    jobs = path / 'jobs'
    with locked(jobs / '.submit.lock'):
        if not dry_run:
            check_submissions(jobs)
        number = 1 + max([int(p.name.split('_')[1].split('.')[0])
                          for p in jobs.glob('submission_*.resources.json')] + [0])
        execution = jobs / f'submission_{number:04d}.resources.json'
        write_json(execution, {'slurm': slurm, 'sha256': digest(slurm)})
        script = jobs / f'{number:04d}_verify.sh'
        command = [sys.executable, str(Path(manifest['root']) / 'workflow/scripts/dataset.py'),
                   'dispatch', '--build', str(path), '--until', until, '--execution', str(execution)]
        if species is not None:
            frozen_list = jobs / f'{number:04d}_samples.txt'
            shutil.copy2(species, frozen_list)
            command += ['--species-list', str(frozen_list)]
        script.write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec ' + shlex.join(command) + '\n')
        script.chmod(0o755)
        job = slurm['stages']['controller']
        cmd = ['sbatch', '--parsable', '--nodes=1', '--ntasks=1', f'--cpus-per-task={job["cpus"]}',
               f'--mem={job["mem_gb"] * 1000}M', f'--time={job["time"]}', '--chdir=' + manifest['root'],
               '--job-name=' + manifest['name'] + '_verify',
               f'--output={jobs}/logs/{number:04d}_verify_%j.out',
               f'--error={jobs}/logs/{number:04d}_verify_%j.err']
        if slurm.get('partition'):
            cmd.append('--partition=' + slurm['partition'])
        cmd.append(str(script))
        receipt = jobs / f'submission_{number:04d}.json'
        batch = {'created_at': now(), 'until': until, 'resources': record(execution),
                 'jobs': [{'stage': 'verify', 'indices': None, 'command': cmd,
                           'job_id': None, 'state': 'submitting'}]}
        print_plan(context['items'])
        if dry_run:
            print(shlex.join(cmd))
            print('Worker jobs will be generated after verification on Slurm.')
        else:
            write_json(receipt, batch)
            try:
                identifier = subprocess.check_output(cmd, text=True).strip().split(';')[0]
            except (OSError, subprocess.CalledProcessError):
                batch['jobs'][0]['state'] = 'rejected'
                write_json(receipt, batch)
                raise
            if not identifier.isdigit():
                batch['jobs'][0]['state'] = 'unknown'
                write_json(receipt, batch)
                raise ValueError('unrecognized sbatch result; inspect queue before retrying')
            batch['jobs'][0].update(state='submitted', job_id=identifier)
            write_json(receipt, batch)
            update_verification(path, 'queued', job_id=identifier)
            print(f'Verification job: {identifier}')
        print(f'Status: ./run_build.sh status --build {path}')
        return [cmd]
