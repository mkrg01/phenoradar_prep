"""One reviewed decision per run: adopt, exclude and clean up, or retain for retry."""
from contextlib import contextmanager, ExitStack
import fcntl
import json
import os
from pathlib import Path
import pwd
import re
import subprocess

from accession_exclusions import read_exclusions
from cleanup_work import owned_path
from common import file_record, now, write_json, write_tsv
from dataset_assets import counts, locked
import metadata_catalog as catalog


FIELDS = ['run', 'action', 'reason', 'status', 'species', 'taxid', 'job_state',
          'stopped_at', 'scheduler_state', 'failure_reason',
          'busco_complete', 'busco_total', 'busco_complete_fraction',
          'reclaimable_bytes']
TERMINAL = {'COMPLETED', 'OUT_OF_MEMORY', 'TIMEOUT', 'CANCELLED', 'FAILED',
            'NODE_FAIL', 'PREEMPTED', 'BOOT_FAIL', 'DEADLINE', 'REVOKED'}


def sidecar(path, kind):
    return Path(str(path) + f'.{kind}.json')


@contextmanager
def sample_guard(path, build, item, *, native=False):
    """Probe existing locks without creating any files during a preview."""
    from dataset import workspace
    with ExitStack() as stack:
        worker = Path(build['config']['store']) / item['species'] / '.worker.lock'
        worker = owned_path(build['root'], worker.relative_to(Path(build['root'])))
        paths = [worker]
        if native:
            work = workspace(path, build, item)
            paths.append(owned_path(path, work.relative_to(path) / '.native.lock'))
            scratch = owned_path(path, work.relative_to(path) / 'output/transcriptome_assembly/tmp')
            paths.extend(owned_path(path, lock.relative_to(path)) for lock in scratch.glob('*/.gg_active.lock'))
        for lock in paths:
            if lock.is_symlink():
                raise ValueError(f'unsafe sample lock: {lock}')
            if lock.exists():
                handle = stack.enter_context(lock.open('r'))
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


def jobs(path, build, runs=None):
    """Combine receipts, accounting and live array tasks, including pending retries."""
    from dataset import STAGES, workspace
    from genegalleon_run import observation, scheduler_states
    result, receipts, identifiers, unknown = {}, [], {}, set()
    items = build['items']
    for item in items:
        if runs is not None and item['row']['run'] not in runs:
            continue
        name = item['species']
        rows, evidence = [], []
        for stage in (*STAGES, 'sample'):
            receipt = path / 'jobs/status' / f'{name}.{stage}.json'
            if receipt.is_file():
                row = json.loads(receipt.read_text())
                rows.append(dict(row, receipt_stage=stage))
                evidence.append(file_record(receipt))
                identifier = str(row.get('slurm_id') or row.get('job_id') or '')
                if re.fullmatch(r'[0-9]+(?:_[0-9]+)?', identifier):
                    identifiers.setdefault(identifier, set()).add(name)
        primary = next((r for r in rows if r['receipt_stage'] == 'sample'), None)
        if primary is None:
            primary = next((r for r in reversed(rows) if r.get('state') in {'failed', 'interrupted', 'running'}), {})
        result[name] = {'primary': primary, 'receipts': evidence,
                        'native': observation(workspace(path, build, item)), 'activity': None,
                        'submissions': [], 'submitted_id': None}
        receipts.extend(rows)
    for receipt in sorted((path / 'jobs').glob('submission_*.json')):
        if receipt.name.endswith('.resources.json'):
            continue
        for job in json.loads(receipt.read_text()).get('jobs', []):
            if job.get('stage') not in {*STAGES, 'sample'}:
                continue
            indices = job.get('indices') or []
            names = {items[i - 1]['species'] for i in indices} & result.keys()
            if names:
                evidence = file_record(receipt)
                for name in names:
                    result[name]['submissions'].append(evidence)
            if job.get('state') in {'submitting', 'unknown'}:
                unknown.update(names)
            if job.get('job_id'):
                base = str(job['job_id'])
                offset = job.get('array_offset', 0)
                for index in indices:
                    name = items[index - 1]['species']
                    if name in result:
                        identifier = f'{base}_{index - offset}'
                        identifiers.setdefault(identifier, set()).add(name)
                        result[name]['submitted_id'] = identifier
    accounting = scheduler_states([*receipts, *({'slurm_id': identifier} for identifier in identifiers)])
    if identifiers:
        try:
            output = subprocess.run(['squeue', '--noheader', '--array', '--user',
                                     pwd.getpwuid(os.getuid()).pw_name, '--format=%i'],
                                    capture_output=True, text=True, timeout=10, check=True).stdout
            active = set(output.split())
            for identifier, names in identifiers.items():
                base = identifier.split('_', 1)[0]
                if (identifier in active or base in active
                        or any(j.startswith(base + '_[') for j in active)
                        or ('_' not in identifier and any(j.startswith(base + '_') for j in active))):
                    for name in names:
                        result[name]['activity'] = 'queued_or_running'
        except (OSError, subprocess.SubprocessError):
            unknown.update(name for names in identifiers.values() for name in names)
    for name, state in result.items():
        primary = state['primary']
        identifier = state['submitted_id'] or str(primary.get('slurm_id') or primary.get('job_id') or '')
        scheduler = accounting.get(identifier)
        job_state = primary.get('state', 'not_started')
        if scheduler:
            if scheduler['state'] in TERMINAL - {'COMPLETED'}:
                job_state = 'failed'
            elif scheduler['state'] == 'COMPLETED' and job_state in {'running', 'not_started'}:
                job_state = 'interrupted'
            elif scheduler['state'] not in TERMINAL:
                state['activity'] = 'queued_or_running'
        if name in unknown:
            state['activity'] = state['activity'] or 'unresolved_submission'
        if job_state == 'running' and not state['activity']:
            state['activity'] = 'running_or_unresolved'
        state.update(job_state=job_state, scheduler=scheduler,
                     stopped_at=primary.get('stage') or state['native'].get('current_stage'))
    return result


def review(root, cfg, build_path, runs=None, exclude_failed=False):
    from dataset import load, item_products, cleanup_genegalleon
    path = catalog.path_at(root, build_path)
    if path.name == 'build.json':
        path = path.parent
    build = load(path)
    if Path(build['root']).resolve() != Path(root).resolve():
        raise ValueError('review build belongs to a different project')
    catalog.matching_container(root, cfg, build)
    # A frozen build points to its private copy, so compare the source configuration.
    configured = catalog.yaml.safe_load(catalog.path_at(root, cfg['build_config']).read_text()).get('excluded_accessions')
    if not configured or catalog.path_at(root, configured) != catalog.path_at(root, cfg['excluded_accessions']):
        raise ValueError('selection and build must use the same excluded_accessions table')
    before = catalog.selection_inputs(root, cfg)
    _, previous, accepted = catalog.previous_state(root, cfg)
    exclusions = read_exclusions(catalog.path_at(root, cfg['excluded_accessions']))
    items = {item['row']['run']: item for item in build['items']}
    requested = set(runs) if runs is not None else set(items)
    if requested - items.keys():
        raise ValueError('runs absent from build: ' + ', '.join(sorted(requested - items.keys())))
    states = jobs(path, build, requested)
    report, snapshots = [], {}
    for run in sorted(requested):
        item = items[run]
        row, state = item['row'], states[item['species']]
        tid = catalog.taxid(row['taxid'])
        scheduler = state['scheduler'] or {}
        failure = state['primary'].get('error') or (
            f"{scheduler.get('state')} (exit {scheduler.get('exit_code')})"
            if state['job_state'] in {'failed', 'interrupted'} and scheduler else '')
        entry = dict(run=run, species=item['species'], taxid=tid, action='hold', reason='',
                     job_state=state['job_state'], stopped_at=state['stopped_at'], reclaimable_bytes=0,
                     scheduler_state=scheduler.get('state', ''), failure_reason=failure)
        products, conflict = {}, None
        try:
            with sample_guard(path, build, item, native=True):
                products = item_products(build, item)
        except BlockingIOError:
            state['activity'] = 'locked'
        except (OSError, ValueError) as error:
            conflict = str(error)
        snapshots[run] = dict(receipts=state['receipts'], products=products, conflict=conflict,
                              native=state['native'], job_state=state['job_state'], submissions=state['submissions'])
        if not entry['stopped_at'] and state['job_state'] in {'failed', 'interrupted'}:
            entry['stopped_at'] = next((stage for stage, key in
                                       zip(('assembly', 'busco', 'quant'), ('reference', 'busco', 'quant'))
                                       if not products.get(key)), None)
        allowed = ['hold']
        current = previous.get(tid)
        retained = tid in accepted and accepted[tid]['run'] not in exclusions
        complete = all(products.get(key) for key in ('reference', 'busco', 'quant'))
        busco = products.get('busco') or products.get('assessment')
        if busco:
            qc = counts(busco['counts'])
            entry.update(busco_complete=qc['busco_cds_single'] + qc['busco_cds_duplicated'],
                         busco_total=qc['busco_cds_total'])
            entry['busco_complete_fraction'] = entry['busco_complete'] / entry['busco_total']
        if state['activity']:
            entry.update(status='active', reason=state['activity'])
        elif run in exclusions:
            entry.update(status='excluded', action='exclude', reason=exclusions[run].get('reason') or 'manually_excluded')
            allowed.append('exclude')
        elif conflict:
            entry.update(status='conflict', reason=conflict)
        elif not current or any(current[k] != row[k] for k in catalog.ACCEPTED_FIELDS[:4]):
            entry.update(status='not_current_representative', reason='not_current_representative')
        elif products.get('reference') and busco and entry['busco_complete_fraction'] < cfg['busco_threshold']:
            entry.update(status='busco_below_threshold', retained_accepted=retained,
                         reason=f"busco_completeness_below_{cfg['busco_threshold']}")
            allowed.append('exclude')
            if not retained or (exclude_failed and not complete):
                entry['action'] = 'exclude'
        elif complete:
            entry.update(status='reviewed_complete', action='accept', reason='complete_and_busco_passed')
            allowed.extend(['accept', 'exclude'])
        else:
            entry.update(status='incomplete', reason=failure or 'incomplete')
            allowed.append('exclude')
            if exclude_failed:
                pending = next(stage for stage, key in
                               zip(('assembly', 'busco', 'quant'), ('reference', 'busco', 'quant'))
                               if not products.get(key))
                suffix = 'before_start' if state['job_state'] == 'not_started' else entry['stopped_at'] or pending
                entry.update(action='exclude', reason='processing_abandoned_' + suffix)
        if 'exclude' in allowed:
            with sample_guard(path, build, item):
                cleanup = cleanup_genegalleon(path, build, item, 'excluded', products,
                                             apply=False, discard_reads=True)
            entry['reclaimable_bytes'] = cleanup['reclaimable_bytes']
            entry['cleanup_errors'] = cleanup['errors']
        entry['allowed_actions'] = allowed
        report.append(entry)
    return dict(schema_version=1, kind='sample_review', root=str(Path(root).resolve()),
                build_manifest=file_record(path / 'build.json'), selection_config=file_record(cfg['_path']),
                inputs=before, samples=report, snapshots=snapshots, busco_threshold=cfg['busco_threshold'])


def result_for(data, cfg, *, dry_run):
    exclusions = read_exclusions(catalog.path_at(data['root'], cfg['excluded_accessions']))
    _, _, accepted = catalog.previous_state(data['root'], cfg)
    accepted = {tid: row for tid, row in accepted.items() if row['run'] not in exclusions}
    recorded, excluded = [], []
    for entry in data['samples']:
        run, tid = entry['run'], entry['taxid']
        if entry['action'] == 'accept':
            recorded.append(run)
            accepted[tid] = entry
        elif entry['action'] == 'exclude':
            if run not in exclusions:
                excluded.append(run)
            if accepted.get(tid, {}).get('run') == run:
                accepted.pop(tid)
    return dict(recorded_runs=recorded, excluded_runs=excluded, accepted_samples=len(accepted),
                samples=data['samples'], dry_run=dry_run, busco_threshold=data['busco_threshold'])


def clean_excluded(path, build, samples):
    from dataset import cleanup_genegalleon
    items = {item['row']['run']: item for item in build['items']}
    report = []
    states = jobs(path, build, {row['run'] for row in samples if row['action'] == 'exclude'})
    for entry in samples:
        if entry['action'] != 'exclude':
            continue
        item = items[entry['run']]
        try:
            if states[item['species']]['activity']:
                raise ValueError('sample has queued/running or unresolved jobs')
            with sample_guard(path, build, item):
                cleanup = cleanup_genegalleon(path, build, item, 'excluded', {}, discard_reads=True)
        except (OSError, ValueError) as error:
            cleanup = dict(state='pending', errors=[{'error': str(error)}])
        report.append(dict(run=entry['run'], **cleanup))
    return report


def apply_review(root, cfg, data, *, plan_path=None, plan_record=None):
    from dataset import load
    path = Path(data['build_manifest']['path']).parent
    build = load(path)
    result = result_for(data, cfg, dry_run=False)
    selected = [row for row in data['samples'] if row['action'] != 'hold']
    # Serialize with submit, then check the reviewed evidence while metadata is locked.
    with locked(path / 'jobs/.submit.lock'), locked(catalog.path_at(root, cfg['previous_metadata']).parent / '.metadata.lock'):
        if (file_record(cfg['_path']) != data['selection_config']
                or catalog.selection_inputs(root, cfg) != data['inputs']
                or file_record(path / 'build.json') != data['build_manifest']
                or (plan_path and file_record(plan_path) != plan_record)):
            raise ValueError('sample review inputs changed; repeat the review')
        current = review(root, cfg, path, [row['run'] for row in selected])
        checked = {row['run']: row for row in current['samples']}
        for entry in selected:
            run = entry['run']
            if (current['snapshots'][run] != data['snapshots'][run]
                    or entry['action'] not in checked[run]['allowed_actions']):
                raise ValueError(f'sample review state changed for {run}; repeat the review')
        _, _, accepted = catalog.previous_state(root, cfg)
        exclusion_path = catalog.path_at(root, cfg['excluded_accessions'])
        exclusions = read_exclusions(exclusion_path)
        accepted = {tid: row for tid, row in accepted.items() if row['run'] not in exclusions}
        items = {item['row']['run']: item for item in build['items']}
        decisions = []
        for entry in selected:
            run, tid = entry['run'], entry['taxid']
            row = items[run]['row']
            if entry['action'] == 'exclude':
                decisions.append(dict(accession=run, taxid=tid, bioproject=row['bioproject'],
                                      reason=entry['reason'], source_build=build['name']))
                if accepted.get(tid, {}).get('run') == run:
                    accepted.pop(tid)
            else:
                accepted[tid] = dict(**{k: row[k] for k in catalog.ACCEPTED_FIELDS[:4]},
                                     busco_complete=entry['busco_complete'], busco_total=entry['busco_total'],
                                     source_build=build['name'])
        catalog.append_exclusions(exclusion_path, decisions)
        target = catalog.path_at(root, cfg['accepted_samples'])
        write_tsv(target, catalog.ACCEPTED_FIELDS, [accepted[tid] for tid in sorted(accepted, key=int)])
        write_json(target.with_suffix('.provenance.json'), dict(kind='reviewed_samples', created_at=now(),
                   build_manifest=data['build_manifest'], selection_config=data['selection_config'],
                   inputs=data['inputs'], recorded_runs=result['recorded_runs'], excluded_runs=result['excluded_runs'],
                   busco_threshold=data['busco_threshold'], samples=data['samples'],
                   excluded_accessions=file_record(exclusion_path), accepted_samples=file_record(target)))
        # Commit the decision before deleting reads. A repeated application only retries cleanup.
        if plan_path:
            write_json(sidecar(plan_path, 'applied'), dict(plan=plan_record, samples=data['samples'], result=result))
        result['cleanup'] = clean_excluded(path, build, data['samples'])
    return result


def record(root, cfg, build_path, runs=None, dry_run=False, *, exclude_failed=False, plan=None):
    data = review(root, cfg, build_path, runs, exclude_failed)
    result = result_for(data, cfg, dry_run=dry_run or bool(plan))
    if plan:
        target = catalog.path_at(root, plan)
        if target.exists() or sidecar(target, 'review').exists() or sidecar(target, 'applied').exists():
            raise ValueError('review plan already exists; choose a new path')
        write_tsv(target, FIELDS, [{field: row.get(field) if row.get(field) is not None else ''
                                   for field in FIELDS} for row in data['samples']])
        write_json(sidecar(target, 'review'), data)
        result['plan'] = str(target)
        return result
    if dry_run:
        return result
    return apply_review(root, cfg, data)


def apply_plan(root, cfg, plan, *, dry_run=False):
    from dataset import load
    target = catalog.path_at(root, plan)
    data = json.loads(sidecar(target, 'review').read_text())
    if (data.get('schema_version') != 1 or data.get('kind') != 'sample_review'
            or data.get('root') != str(Path(root).resolve())
            or data['selection_config'] != file_record(cfg['_path'])):
        raise ValueError('review plan does not match this project and selection configuration')
    _, rows = catalog.table(target, FIELDS)
    originals = {row['run']: row for row in data['samples']}
    if len({row['run'] for row in rows}) != len(rows):
        raise ValueError('duplicate run in review plan')
    chosen = []
    for row in rows:
        original = originals.get(row['run'])
        if not original or any(row[key] != (str(original.get(key)) if original.get(key) is not None else '')
                               for key in FIELDS if key not in {'action', 'reason'}):
            raise ValueError('only action and reason may be edited in a review plan; rows may be removed')
        if row['action'] not in original['allowed_actions']:
            raise ValueError(f"invalid review action for {row['run']}: {row['action']}")
        if row['action'] == 'exclude' and not row['reason'].strip():
            raise ValueError('excluded runs require a reason')
        chosen.append(dict(original, action=row['action'], reason=row['reason']))
    data['samples'] = chosen
    if dry_run:
        return result_for(data, cfg, dry_run=True)
    record = file_record(target)
    applied = sidecar(target, 'applied')
    if applied.exists():
        receipt = json.loads(applied.read_text())
        if receipt['plan'] != record:
            raise ValueError('applied review plan changed; prepare a new review')
        path = Path(data['build_manifest']['path']).parent
        if file_record(path / 'build.json') != data['build_manifest']:
            raise ValueError('review build manifest changed')
        with locked(path / 'jobs/.submit.lock'), locked(catalog.path_at(root, cfg['previous_metadata']).parent / '.metadata.lock'):
            exclusions = read_exclusions(catalog.path_at(root, cfg['excluded_accessions']))
            if any(row['run'] not in exclusions for row in chosen if row['action'] == 'exclude'):
                raise ValueError('an applied exclusion was removed; prepare a new review')
            result = dict(receipt['result'], cleanup=clean_excluded(path, load(path), chosen), cleanup_only=True)
        return result
    return apply_review(root, cfg, data, plan_path=target, plan_record=record)
