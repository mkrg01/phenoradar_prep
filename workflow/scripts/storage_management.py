"""Inspect cleanup receipts and retry successful jobs without rerunning tools."""
import fcntl
import json
from collections import Counter
from pathlib import Path

from cleanup_work import inspect_job, owned_path
from layout import ORTHOGROUP_MAPPING, run_layout


def read(path):
    return json.loads(Path(path).read_text())


def managed(root, path):
    path = Path(path).absolute()
    return owned_path(root, path.relative_to(Path(root).absolute()))


def job_entries(path, manifest, kind):
    root = Path(manifest['root'])
    pipeline = manifest['analysis' if kind == 'build' else 'pipeline']
    output, work, _ = [root / p for p in run_layout(pipeline)]
    for tool, directory, destination, pattern in (
            ('odb', work / ORTHOGROUP_MAPPING, output / ORTHOGROUP_MAPPING / 'chunks', 'chunk_*'),
            ('kofam', work / 'kegg', output / 'kegg/annotation', 'batch_*')):
        directory = managed(root, directory)
        if directory.is_symlink():
            raise ValueError(f'linked work directory: {directory}')
        for status in sorted(directory.glob(f'{pattern}/status.json')):
            base = managed(root, status.parent)
            if base.is_symlink():
                raise ValueError(f'linked job directory: {base}')
            yield {'kind': tool, 'job': base.name, 'work_dir': str(base),
                   'output_dir': str(managed(root, destination / base.name)),
                   'status': str(status), 'receipt': str(base / 'cleanup.json')}
    if kind == 'build':
        from dataset import STAGES, workspace
        for item in manifest['items']:
            work = managed(root, workspace(path, manifest, item))
            if not work.exists():
                continue
            for stage in STAGES:
                status = Path(path) / 'jobs/status' / f"{item['species']}.{stage}.json"
                if status.is_file():
                    yield {'kind': 'genegalleon', 'job': item['species'], 'stage': stage,
                           'work_dir': str(work), 'status': str(status),
                           'receipt': str(Path(path) / 'jobs/cleanup' / f"{item['species']}.{stage}.json")}


def entry_state(entry, keep):
    job = read(entry['status'])
    if job.get('state') not in {'success', 'complete', 'reused'}:
        return job.get('state', 'unknown')
    if keep:
        return 'retained'
    receipt = Path(entry['receipt'])
    if receipt.is_file():
        result = read(receipt)
        if result.get('state') == 'pending':
            entry['errors'] = result.get('errors', [])
            return 'pending'
        if entry['kind'] == 'genegalleon':
            return result.get('state', 'unknown')
    if entry['kind'] != 'genegalleon':
        base = Path(entry['work_dir'])
        scratch = Path(job.get('work', base / job['fingerprint']))
        if scratch.parent != base or scratch.is_symlink():
            raise ValueError('scratch must be a direct managed child')
        if scratch.exists():
            return 'available'
        # Legacy Kofam status may omit its temporary-directory suffix.
        if entry['kind'] == 'kofam' and 'work' not in job:
            from run_kofam import completed_work
            scratch = completed_work(read(Path(entry['output_dir']) / 'provenance.json'), base)
            if scratch.exists():
                return 'available'
        return 'complete'
    return 'available'


def retry_genegalleon(path, manifest, entry, apply):
    from dataset import STAGES, cleanup_genegalleon, item_products
    item = next(item for item in manifest['items'] if item['species'] == entry['job'])
    lock = managed(manifest['root'], Path(manifest['config']['store']) / item['species'] / '.worker.lock')
    # Existing successful jobs have a worker lock. Preview never creates files.
    with lock.open('r') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for stage in STAGES:
            status = Path(path) / 'jobs/status' / f"{item['species']}.{stage}.json"
            if status.is_file() and read(status).get('state') in {'running', 'failed'}:
                raise ValueError(f'{stage} is running or failed; sample scratch retained')
        job = read(entry['status'])
        if job.get('state') not in {'success', 'complete', 'reused'}:
            raise ValueError('cleanup requires a successful, inactive stage')
        products = item_products(manifest, item)
        key = dict(zip(STAGES, ('reference', 'busco', 'quant')))[entry['stage']]
        if not products.get(key):
            raise ValueError('registered stage products are missing; scratch retained')
        result = cleanup_genegalleon(path, manifest, item, entry['stage'], products, apply=False)
        if result['errors']:
            raise ValueError(str(result['errors']))
        if apply:
            result['cleanup'] = cleanup_genegalleon(path, manifest, item, entry['stage'], products)
        return result


def storage_report(path, manifest, kind, *, inspect=False, apply=False):
    path = Path(path).resolve()
    keep = manifest['config'].get('storage', {}).get('keep_intermediates', False)
    jobs = []
    for entry in job_entries(path, manifest, kind):
        try:
            entry['state'] = entry_state(entry, keep)
            if inspect and entry['state'] in {'available', 'pending'}:
                if entry['kind'] == 'genegalleon':
                    result = retry_genegalleon(path, manifest, entry, apply)
                else:
                    result = inspect_job(entry['kind'], entry['work_dir'], entry['output_dir'], apply=apply)
                entry['inspection'] = result
                if apply:
                    entry['state'] = result['cleanup']['state']
                    entry.pop('errors', None)
        except (OSError, ValueError, KeyError, TypeError) as error:
            entry.update(state='blocked', error=str(error))
        jobs.append(entry)
    counts = dict(Counter(entry['state'] for entry in jobs))
    return {'state': 'pending' if counts.get('pending') or counts.get('blocked') else 'complete',
            'keep_intermediates': keep, 'apply': apply, 'counts': counts, 'jobs': jobs}


def run_storage(path, kind, *, inspect=False, apply=False):
    if kind == 'build':
        from dataset import load
    else:
        from analysis import load
    return storage_report(path, load(path), kind, inspect=inspect, apply=apply)
