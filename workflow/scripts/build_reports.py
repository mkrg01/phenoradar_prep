"""Small plan and progress reports, independent of artifact verification."""
from collections import Counter
import json
from pathlib import Path

from common import now, read_tsv, write_json, write_tsv
from dataset_assets import digest


STAGES = ('assembly', 'busco', 'quant')
PLAN_FIELDS = ['species', 'run', 'action', *STAGES, 'mapping', 'verification', 'reason', 'reference_id']
STATUS_FIELDS = ['species', 'run', 'state', 'action', 'verification', 'stage', 'reason', 'job_id', 'log']


def action(row):
    if row.get('action') in {'excluded', 'blocked', 'rebuild', 'reuse', 'partial', 'new'}:
        return row['action']
    values = [row.get(stage) for stage in STAGES]
    if 'excluded' in values:
        return 'excluded'
    if 'conflict' in values:
        return 'blocked'
    if row.get('rebuild'):
        return 'rebuild'
    if all(value == 'reuse' for value in values):
        return 'reuse'
    if 'reuse' in values:
        return 'partial'
    return 'new'


def plan_rows(report, verified=False):
    rows = []
    for entry in report:
        row = dict(entry)
        row['action'] = action(row)
        row['verification'] = ('failed' if row['action'] == 'rebuild' else
                               'verified' if verified and 'reuse' in [row.get(s) for s in STAGES] else
                               'pending' if 'reuse' in [row.get(s) for s in STAGES] else 'not_required')
        rows.append(row)
    return rows


def write_plan(report, directory, verified=False, split=False):
    directory = Path(directory)
    rows = plan_rows(report, verified)
    write_tsv(directory / 'plan.tsv', PLAN_FIELDS, [{k: r.get(k, '') for k in PLAN_FIELDS} for r in rows])
    attention = [r for r in rows if r['action'] in {'blocked', 'rebuild'}]
    write_tsv(directory / 'needs_attention.tsv', PLAN_FIELDS, [{k: r.get(k, '') for k in PLAN_FIELDS} for r in attention])
    if split:
        for filename, actions in (('reuse', {'reuse'}), ('pending', {'new', 'partial', 'rebuild'}),
                                  ('conflicts', {'blocked'}), ('excluded', {'excluded'})):
            write_tsv(directory / f'{filename}.tsv', PLAN_FIELDS,
                      [{k: r.get(k, '') for k in PLAN_FIELDS} for r in rows if r['action'] in actions])
    return rows


def print_plan(rows, directory=None):
    counts = Counter(row['action'] for row in rows)
    included = len(rows) - counts['excluded']
    processing = counts['new'] + counts['partial'] + counts['rebuild']
    unverified = any(row['verification'] == 'pending' for row in rows)
    print(f'Samples: {included} included; {counts["excluded"]} excluded')
    print(f'  Reuse{" candidates" if unverified else ""}: {counts["reuse"]}')
    print(f'  Process: {processing} (new: {counts["new"]}; partial reuse: {counts["partial"]}; rebuild: {counts["rebuild"]})')
    print(f'  Needs attention: {counts["blocked"]}')
    if unverified:
        print('Product checksums are pending verification on Slurm.')
    if directory is not None:
        print(f'Details: {Path(directory) / "plan.tsv"}')
        if counts['blocked'] or counts['rebuild']:
            print(f'Reasons: {Path(directory) / "needs_attention.tsv"}')


def save_context(directory, manifest, report, deferred=False):
    data = {'schema_version': 1, 'name': manifest['name'], 'root': manifest['root'],
            'prepared_at': manifest.get('prepared_at', manifest['created_at']),
            'store': manifest['config']['store'], 'lineage': manifest['config']['busco']['lineage'],
            'deferred': deferred, 'items': plan_rows(report, verified=not deferred)}
    data['sha256'] = digest(data)
    write_json(Path(directory) / 'status_context.json', data)


def read_context(path):
    filename = Path(path) / 'status_context.json'
    if not filename.is_file():
        # Existing builds predate the compact progress index. Read their frozen
        # configuration once, without opening the registered product files.
        from dataset import excluded_report, load
        manifest = load(path)
        planned = {row['species']: row for row in read_tsv(Path(path) / 'plan.tsv')}
        rows = [planned[item['species']] for item in manifest['items']]
        rows += excluded_report(manifest.get('excluded', []))
        deferred = manifest.get('deferred_verification', False)
        return {'name': manifest['name'], 'root': manifest['root'], 'deferred': deferred,
                'items': plan_rows(rows, verified=not deferred)}
    data = json.loads(filename.read_text())
    if data.get('schema_version') != 1 or digest({k: v for k, v in data.items() if k != 'sha256'}) != data.get('sha256'):
        raise ValueError('saved status context changed')
    return data


def controller_state(data):
    from genegalleon_run import scheduler_states
    from sample_review import TERMINAL
    result = dict(data)
    if result.get('job_id') and result.get('state') in {'queued', 'running'}:
        accounting = scheduler_states([result]).get(str(result['job_id']))
        if accounting:
            result['scheduler'] = accounting
            state = accounting['state']
            if state in TERMINAL:
                result.update(state='failed', error=f'Slurm {state}; completion was not recorded')
            elif state in {'RUNNING', 'COMPLETING', 'SUSPENDED'}:
                result['state'] = 'running'
    return result


def progress(path, context):
    from sample_review import jobs
    path = Path(path).resolve()
    active = [row for row in context['items'] if row['action'] != 'excluded']
    build = {'root': context['root'], 'analysis': {'output_root': True},
             'items': [{'species': row['species'], 'row': {'run': row['run']}} for row in active]}
    evidence = jobs(path, build)
    verification_path = path / 'jobs/verification.json'
    verification = controller_state(json.loads(verification_path.read_text())) if verification_path.is_file() else {}
    for receipt in sorted((path / 'jobs').glob('submission_*.json')):
        if receipt.name.endswith('.resources.json'):
            continue
        for job in json.loads(receipt.read_text()).get('jobs', []):
            if job['stage'] == 'verify' and job['state'] in {'submitting', 'unknown'}:
                verification = {'state': 'unknown', 'error': f'Unresolved submission: {receipt}'}
    rows = []
    for planned in active:
        row = dict(planned)
        info = evidence[row['species']]
        primary = info['primary']
        scheduler = info.get('scheduler') or {}
        scheduling = scheduler.get('state')
        if info.get('activity') == 'unresolved_submission':
            state = 'unknown'
        elif scheduling in {'PENDING', 'CONFIGURING', 'REQUEUED', 'RESIZING'}:
            state = 'queued'
        elif scheduling in {'RUNNING', 'COMPLETING', 'SUSPENDED'}:
            state = 'running'
        elif primary.get('verification') == 'verified' and primary.get('state') in {'complete', 'reused'}:
            state = primary['state']
        elif info['job_state'] in {'failed', 'interrupted'}:
            state = 'failed'
        elif primary.get('state') in {'complete', 'reused'}:
            state = 'reused' if primary['state'] == 'reused' else 'complete'
        elif info.get('activity') == 'queued_or_running':
            state = 'active'
        elif info['job_state'] == 'running':
            state = 'unknown'
        elif verification.get('state') == 'running':
            state = 'verifying'
        elif verification.get('state') == 'queued':
            state = 'queued'
        elif verification.get('state') == 'unknown':
            state = 'unknown'
        elif not context['deferred'] and planned['action'] == 'reuse':
            state = 'reused'
        else:
            state = 'not_submitted'
        if verification.get('state') == 'failed' and state in {'verifying', 'not_submitted'}:
            state = 'blocked'
        failure = primary.get('error') or (f'Slurm {scheduling}' if scheduling else '')
        row.update(state=state, stage=info.get('stopped_at') or '',
                   verification='verified' if state in {'complete', 'reused'} else row['verification'],
                   reason=(failure or verification.get('error', '')) if state in {'failed', 'blocked'} else '',
                   job_id=info.get('submitted_id') or primary.get('slurm_id') or primary.get('job_id') or '',
                   log=str(path / 'jobs/logs'))
        rows.append(row)
    for row in context['items']:
        if row['action'] == 'excluded':
            rows.append(dict(row, state='excluded', stage='', job_id='', log=''))
    controller = path / 'jobs/database.json'
    database = controller_state(json.loads(controller.read_text())) if controller.is_file() else {}
    if (path / 'completed.json').is_file():
        database = dict(database, state='complete')
    elif not database:
        database = {'state': 'not_published'}
    return {'build': str(path), 'name': context['name'], 'checked_at': now(),
            'coverage': 'recorded-progress', 'verification': verification,
            'database': database, 'samples': rows}


def verify_progress(path, report):
    """Optional content audit; never stage, repair or run sample work."""
    from build_dispatch import candidate_files, runtime_data
    from dataset import item_products, load
    from dataset_assets import verify
    from build_products import load_complete
    from database_reuse import validate_product
    from mapping_tables import load_tables
    path = Path(path).resolve()
    manifest = load(path)
    runtime = runtime_data(path)
    candidates = None
    if manifest.get('deferred_verification') and runtime is None:
        candidates = json.loads(verify(manifest['auxiliary']['reuse_candidates.json']).read_text())
        for entry in candidates['sources']:
            verify(entry)
    for entry in (runtime or {}).get('shared_files', []):
        verify(entry)
    database = load_complete(path, verify_files=False) if (path / 'completed.json').is_file() else None
    tables = load_tables(verify(database['mapping']), verify_files=False) if database else None
    items = {item['species']: item for item in manifest['items']}
    for row in report['samples']:
        name = row['species']
        if name not in items:
            continue
        try:
            if database:
                validate_product(items[name], database['products'][name], database, tables, manifest['config'])
            elif candidates is not None:
                chosen = candidates['samples'].get(name)
                if chosen:
                    for entry in candidate_files(chosen):
                        verify(entry)
                elif name in candidates['rebuilds']:
                    raise ValueError(candidates['rebuilds'][name])
                else:
                    row['verification'] = 'not_required'
                    continue
            else:
                for entry in (runtime or {}).get('verified_files', {}).get(name, []):
                    verify(entry)
                products = item_products(manifest, items[name])
                if not any(products.values()):
                    row['verification'] = 'not_required'
                    continue
            row['verification'] = 'verified'
        except (ValueError, OSError, KeyError) as error:
            row.update(verification='failed', reason=str(error))
    report['coverage'] = 'content-verified'
    return report


def print_progress(report):
    counts = Counter(row['state'] for row in report['samples'])
    included = len(report['samples']) - counts['excluded']
    completed = counts['reused'] + counts['complete']
    print(f'Build: {report["build"]}')
    print(f'Samples: {completed}/{included} finished (reused: {counts["reused"]}; completed here: {counts["complete"]})')
    print('  ' + '; '.join(f'{state}: {counts[state]}' for state in
                          ('verifying', 'running', 'queued', 'active', 'failed', 'blocked', 'not_submitted', 'unknown')))
    print(f'Database: {report["database"]["state"]}')
    if report.get('verification'):
        preparation = report['verification']
        detail = (f' ({preparation.get("phase", "checksums")}; {preparation["checked"]}/{preparation["total"]} samples checked)'
                  if 'checked' in preparation else '')
        print(f'Preparation: {preparation["state"]}{detail}')
    if report['coverage'] == 'content-verified':
        integrity = Counter(row['verification'] for row in report['samples'])
        print(f'Product verification: {integrity["verified"]} verified; {integrity["failed"]} failed')
    for label in ('verification', 'database'):
        if report[label].get('error'):
            print(f'  {label}: {report[label]["error"]}')
    for row in report['samples']:
        if row['state'] in {'failed', 'blocked'} or row['verification'] == 'failed':
            print(f'  {row["species"]}: {row["stage"] or "verification"}: {row["reason"]}')
            print(f'    Logs: {row["log"]}')


def write_progress(report, directory):
    directory = Path(directory)
    write_tsv(directory / 'status.tsv', STATUS_FIELDS,
              [{k: row.get(k, '') for k in STATUS_FIELDS} for row in report['samples']])
    write_tsv(directory / 'failed.tsv', STATUS_FIELDS,
              [{k: row.get(k, '') for k in STATUS_FIELDS} for row in report['samples']
              if row['state'] in {'failed', 'blocked'} or row.get('verification') == 'failed'])
