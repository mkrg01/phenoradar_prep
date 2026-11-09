"""Sample allocated bytes and quota-relevant inodes without following symlinks."""
import json
import fcntl
import gzip
import os
import signal
import stat
import subprocess
import time
import shutil
from pathlib import Path

from common import now, write_json, atomic_writer
from scratch_reservation import headroom


def measure(roots):
    """Deduplicate hard links and overlapping roots, including directory inodes."""
    result = dict(allocated_bytes=0, logical_bytes=0, inodes=0, files=0,
                  directories=0, symlinks=0, errors=[])
    seen = set()
    stack = [Path(root).absolute() for root in roots]
    while stack:
        path = stack.pop()
        try:
            info = path.lstat()
            identity = info.st_dev, info.st_ino
            if identity in seen:
                continue
            seen.add(identity)
            result['inodes'] += 1
            result['allocated_bytes'] += getattr(info, 'st_blocks', 0) * 512
            if stat.S_ISDIR(info.st_mode):
                result['directories'] += 1
                with os.scandir(path) as entries:
                    stack.extend(Path(entry.path) for entry in entries)
            elif stat.S_ISLNK(info.st_mode):
                result['symlinks'] += 1
            elif stat.S_ISREG(info.st_mode):
                result['files'] += 1
                result['logical_bytes'] += info.st_size
        except FileNotFoundError:
            pass  # A live job can remove a file between scandir and lstat.
        except OSError as error:
            result['errors'].append({'path': str(path), 'error': str(error)})
    return result


def maximum(old, new):
    return {key: max(old.get(key, 0), new[key]) for key in
            ('allocated_bytes', 'logical_bytes', 'inodes', 'files', 'directories', 'symlinks')}


def project_usage(root, interval):
    """Share one recent inventory across workers instead of scanning per sample."""
    root = Path(root).resolve()
    cache = root / 'work/storage_meter'
    cache.mkdir(parents=True, exist_ok=True)
    with (cache / '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        latest = cache / 'latest.json'
        try:
            saved = json.loads(latest.read_text())
            age = time.time() - saved['measured_at_seconds']
            if saved['root'] == str(root) and 0 <= age < interval:
                return saved['usage']
        except (OSError, ValueError, KeyError, TypeError):
            pass
        usage = measure([root])
        write_json(latest, {'root': str(root), 'measured_at_seconds': time.time(), 'usage': usage})
        return usage


def check_budget(usage, settings, fraction=1):
    if usage['errors']:
        raise ValueError(f"Storage inventory incomplete: {usage['errors'][:3]}")
    if settings.get('max_bytes') is not None and usage['allocated_bytes'] >= settings['max_bytes'] * fraction:
        raise ValueError('Project allocated-byte budget reached; stop/reduce work before resuming')
    if settings.get('max_inodes') is not None and usage['inodes'] >= settings['max_inodes'] * fraction:
        raise ValueError('Project inode budget reached; stop/reduce work before resuming')


def stop(process):
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        process.wait()
        return
    try:
        process.wait(timeout=15)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()
    # A shell can exit before background descendants. They must not keep
    # writing into scratch after its owner releases the reservation.
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass


def append_sample(receipt, sample, limit):
    series = receipt.with_suffix('.jsonl')
    line = (json.dumps(sample, separators=(',', ':')) + '\n').encode()
    if len(line) > limit:
        line = (json.dumps({'at': sample['at'], 'omitted': 'sample exceeds series_max_bytes'}) + '\n').encode()
    if series.exists() and series.stat().st_size + len(line) > limit:
        with series.open('rb') as source, atomic_writer(series.with_suffix('.jsonl.1.gz'), 'wb') as raw:
            with gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0) as zipped:
                shutil.copyfileobj(source, zipped)
        series.unlink()
    with series.open('ab') as output:
        output.write(line)


def compact_series(receipt):
    series = receipt.with_suffix('.jsonl')
    if series.is_file():
        with series.open('rb') as source, atomic_writer(series.with_suffix('.jsonl.gz'), 'wb') as raw:
            with gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0) as zipped:
                shutil.copyfileobj(source, zipped)
        series.unlink()


def monitored(command, *, cwd, env, shared_root, scratch, sample_roots, receipt, settings):
    """A sampling watchdog, not a filesystem quota or an exact peak profiler."""
    shared_root, receipt = Path(shared_root), Path(receipt)
    scratch = Path(scratch)
    record = {'schema_version': 1, 'state': 'running', 'started_at': now(),
              'interval_seconds': settings['interval_seconds'], 'command': command,
              'shared_peak': {}, 'scratch_peak': {}, 'sample_peak': {}, 'observations': 0,
              'job_id': os.environ.get('SLURM_JOB_ID') or os.environ.get('JOB_ID'), 'host': os.uname().nodename,
              'note': 'Sampled maxima are lower bounds; in-flight jobs can exceed a budget between scans.'}
    record['storage_limits'] = settings
    receipt.parent.mkdir(parents=True, exist_ok=True)
    process = None
    interrupted = False
    previous = {}

    def terminate(signum, frame):
        nonlocal interrupted
        interrupted = True
        raise InterruptedError(f'Job received signal {signum}')

    def observe():
        shared = project_usage(shared_root, settings['interval_seconds'])
        local = measure([scratch])
        sample = measure(sample_roots)
        for label, value in (('shared', shared), ('scratch', local), ('sample', sample)):
            if value['errors']:
                raise ValueError(f'{label} storage inventory incomplete: {value["errors"][:3]}')
            record[label + '_peak'] = maximum(record[label + '_peak'], value)
        record['observations'] += 1
        record['updated_at'] = now()
        write_json(receipt, record)
        append_sample(receipt, {'at': record['updated_at'], 'shared': shared,
                               'scratch': local, 'sample': sample}, settings.get('series_max_bytes', 1048576))
        check_budget(shared, settings, settings.get('stop_fraction', 1))
        for key, unit in (('allocated_bytes', 'bytes'), ('inodes', 'inodes')):
            ceiling = settings.get('scratch_max_' + unit)
            if ceiling is not None and local[key] >= ceiling:
                raise ValueError(f'Task scratch {unit} budget reached')
        headroom(scratch, settings)

    try:
        check_budget(project_usage(shared_root, settings['interval_seconds']), settings, settings['start_fraction'])
        for signum in (signal.SIGTERM, signal.SIGINT):
            previous[signum] = signal.signal(signum, terminate)
        observe()
        process = subprocess.Popen(command, cwd=cwd, env=env, start_new_session=True)
        record['child_started'] = True
        next_inventory = time.monotonic() + settings['interval_seconds']
        while True:
            try:
                code = process.wait(timeout=min(1, settings['interval_seconds']))
            except subprocess.TimeoutExpired:
                # Free space is cheap to poll between the larger inventories.
                headroom(scratch, settings)
                if time.monotonic() >= next_inventory:
                    observe()
                    next_inventory = time.monotonic() + settings['interval_seconds']
                continue
            observe()
            if code:
                raise subprocess.CalledProcessError(code, command)
            record['state'] = 'complete'
            return record
    except BaseException as error:
        record.update(state='interrupted' if interrupted else 'blocked' if process is None else 'failed', error=str(error))
        record['reason'] = ('signal' if interrupted else 'command_failure' if isinstance(error, subprocess.CalledProcessError)
                            else 'storage_admission' if process is None else 'storage_or_inventory_limit')
        raise
    finally:
        if process is not None:
            stop(process)
        record['finished_at'] = now()
        record['child_started'] = process is not None
        if process is not None:
            record['child_exit_code'] = process.returncode
        write_json(receipt, record)
        if settings.get('compact_logs', True):
            compact_series(receipt)
        for signum, handler in previous.items():
            signal.signal(signum, handler)
