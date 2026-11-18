"""Retain AGE worker diagnostics in one incrementally compressed tar per controller."""
import argparse
from contextlib import contextmanager
import fcntl
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import signal
import stat
import subprocess
import sys
import tarfile
import uuid
import zlib

from common import write_json

ARCHIVE_ERRORS = (OSError, ValueError, EOFError, tarfile.TarError, zlib.error)


@contextmanager
def _locked(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_symlink():
        raise ValueError(f'linked AGE log lock: {path}')
    with path.open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def archive_path(usage):
    parent = os.environ.get('JOB_ID') or uuid.uuid4().hex
    return Path(usage).absolute().parent / 'logs' / f'age-workers-{parent}.tar.gz'


def spool_paths(archive, signature):
    directory = Path(archive).with_suffix('').with_suffix('.pending')
    if directory.is_symlink():
        raise ValueError(f'linked AGE log spool: {directory}')
    directory.mkdir(parents=True, exist_ok=True)
    return [directory / f'{signature}.%j.{suffix}' for suffix in ('out', 'err')]


def _save(path, value):
    write_json(path, value, compact=True)
    with Path(path).open('rb') as handle:
        os.fsync(handle.fileno())


def _fingerprint(path):
    info = path.lstat()
    if not stat.S_ISREG(info.st_mode):
        raise ValueError(f'AGE log is not a regular file: {path}')
    return [info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns]


def _retire(pending, archive):
    # A committed append can be recovered even if the controller died before
    # publishing its worker receipt or retiring the individual files.
    status = Path(pending['status'])
    saved = json.loads(status.read_text())
    for entry in pending['files']:
        path = Path(entry['path'])
        if path.exists():
            if _fingerprint(path) != entry['fingerprint']:
                raise ValueError(f'AGE log changed after archival; retained: {path}')
            path.unlink()
            try:
                path.parent.rmdir()
            except OSError:
                pass
    saved['logs_archived'] = {'archive': str(archive), 'prefix': pending['prefix']}
    _save(status, saved)


def _load(archive, journal):
    if archive.is_symlink() or journal.is_symlink():
        raise ValueError('linked AGE log archive or journal')
    if journal.exists():
        state = json.loads(journal.read_text())
    else:
        if archive.exists() and archive.stat().st_size:
            raise ValueError(f'AGE log archive has no recovery journal: {archive}')
        state = {'schema_version': 1, 'size': 0, 'finalized': False, 'pending': None}
    actual = archive.stat().st_size if archive.exists() else 0
    if actual < state['size']:
        raise ValueError(f'AGE log archive was truncated: {archive}')
    if actual > state['size']:
        # No raw log is retired until its compressed append has been verified
        # and its committed size atomically published in this journal.
        with archive.open('r+b') as handle:
            handle.truncate(state['size'])
    if state['pending']:
        _retire(state['pending'], archive)
        state['pending'] = None
        _save(journal, state)
    return state


def _member(packed, name, source, size):
    info = tarfile.TarInfo(name)
    info.size = size
    info.mode = 0o600
    packed.write(info.tobuf())
    digest = hashlib.sha256()
    count = 0
    for block in iter(lambda: source.read(1024 * 1024), b''):
        count += len(block)
        packed.write(block)
        digest.update(block)
    if count != size:
        raise ValueError(f'AGE log size changed while archiving: {name}')
    packed.write(b'\0' * (-size % tarfile.BLOCKSIZE))
    return digest.hexdigest()


def _verify(archive, offset, expected):
    # Each append is a complete gzip member containing tar entries, without a
    # tar end marker between appends. Standard tar/gzip readers see one archive.
    with archive.open('rb') as raw:
        raw.seek(offset)
        with tarfile.open(fileobj=raw, mode='r|gz') as packed:
            found = {}
            for member in packed:
                if not member.isfile() or member.name in found:
                    raise ValueError('unexpected AGE log archive member')
                digest = hashlib.sha256()
                with packed.extractfile(member) as source:
                    for block in iter(lambda: source.read(1024 * 1024), b''):
                        digest.update(block)
                found[member.name] = digest.hexdigest()
            if found != expected:
                raise ValueError('AGE log archive content verification failed')


def archive_worker(status, outcome):
    status = Path(status).absolute()
    saved = json.loads(status.read_text())
    if not saved.get('log_archive'):
        return  # Historical receipts retain their existing behavior.
    if outcome.get('source') != 'qacct' or outcome['state'] not in {'COMPLETED', 'FAILED'}:
        return  # Eqw and disappearance do not establish that writers have stopped.
    archive = Path(saved['log_archive'])
    journal = archive.with_suffix('.json')
    with _locked(archive.with_suffix('.lock')):
        state = _load(archive, journal)
        if not journal.exists():
            _save(journal, state)
        saved = json.loads(status.read_text())
        if saved.get('logs_archived'):
            return
        prefix = saved['job_id']
        entries, missing = [], []
        for stream in ('stdout', 'stderr'):
            path = Path(saved[stream])
            if path.exists():
                entries.append({'path': str(path), 'fingerprint': _fingerprint(path),
                                'name': f'{prefix}/{stream}.txt'})
            else:
                missing.append(stream)
        if state['finalized']:
            state.update(size=state['data_size'], finalized=False)
            _save(journal, state)
        offset = state['size']
        metadata = json.dumps({**{k: saved.get(k) for k in ('job_id', 'controller_id', 'rule', 'stdout', 'stderr')},
                               'outcome': outcome, 'missing_streams': missing}, sort_keys=True).encode()
        expected = {}
        with archive.open('r+b' if archive.exists() else 'w+b') as raw:
            raw.truncate(offset)
            raw.seek(offset)
            with gzip.GzipFile(filename='', fileobj=raw, mode='wb', compresslevel=6, mtime=0) as packed:
                for entry in entries:
                    with open(entry['path'], 'rb') as source:
                        expected[entry['name']] = _member(packed, entry['name'], source, entry['fingerprint'][2])
                    if _fingerprint(Path(entry['path'])) != entry['fingerprint']:
                        raise ValueError('AGE log changed while archiving')
                name = f'{prefix}/job.json'
                expected[name] = _member(packed, name, io.BytesIO(metadata), len(metadata))
            raw.flush()
            os.fsync(raw.fileno())
            end = raw.tell()
        _verify(archive, offset, expected)
        state.update(size=end, pending={'status': str(status), 'prefix': prefix, 'files': entries})
        _save(journal, state)
        _retire(state['pending'], archive)
        state['pending'] = None
        _save(journal, state)


def finalize(archive):
    archive = Path(archive)
    journal = archive.with_suffix('.json')
    if not journal.exists():
        return
    with _locked(archive.with_suffix('.lock')):
        state = _load(archive, journal)
        if state['finalized']:
            return
        with archive.open('ab') as raw:
            with gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0) as packed:
                packed.write(b'\0' * (2 * tarfile.BLOCKSIZE))
            raw.flush()
            os.fsync(raw.fileno())
            end = raw.tell()
        state.update(data_size=state['size'], size=end, finalized=True)
        _save(journal, state)


def recover(directory, state_root):
    """Collect terminal workers, including those missed by a stopped controller."""
    from scheduler import age_states
    pending = {}
    errors = []
    for path in Path(state_root).glob('*.status.json'):
        saved = json.loads(path.read_text())
        archive = Path(saved.get('log_archive', '/'))
        if archive.parent == Path(directory) and not saved.get('logs_archived'):
            pending[saved['job_id']] = path
    if pending:
        for identifier, outcome in age_states(list(pending)).items():
            try:
                archive_worker(pending[identifier], outcome)
            except ARCHIVE_ERRORS as error:
                errors.append(str(error))
    for journal in Path(directory).glob('age-workers-*.tar.json'):
        try:
            finalize(journal.with_suffix('.gz'))
        except ARCHIVE_ERRORS as error:
            errors.append(str(error))
    if errors:
        raise ValueError('; '.join(errors[:5]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('run', 'recover'))
    parser.add_argument('--usage', required=True)
    parser.add_argument('--state', required=True)
    argv = sys.argv[1:]
    separator = argv.index('--') if '--' in argv else len(argv)
    args = parser.parse_args(argv[:separator])
    archive = archive_path(args.usage)
    def collect():
        try:
            recover(archive.parent, args.state)
            return True
        except (*ARCHIVE_ERRORS, subprocess.SubprocessError) as error:
            print(f'AGE log archival incomplete; raw logs and recovery journal retained: {error}', file=sys.stderr)
            return False
    collected = collect()
    if args.action == 'recover':
        raise SystemExit(0 if collected else 1)
    command = argv[separator + 1:]
    if not command:
        parser.error('a controller command is required after --')
    env = dict(os.environ, PHENORADAR_AGE_LOG_ARCHIVE=str(archive))
    process = subprocess.Popen(command, env=env)
    for signum in (signal.SIGINT, signal.SIGTERM):
        signal.signal(signum, lambda signum, frame: process.send_signal(signum))
    try:
        code = process.wait()
    finally:
        collect()
    raise SystemExit(code if code >= 0 else 128 - code)


if __name__ == '__main__':
    main()
