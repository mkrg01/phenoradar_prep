#!/usr/bin/env python3
"""Inspect or remove scratch from a successfully published ODB/Kofam job."""
import argparse
import fcntl
import json
import os
import shutil
import stat
import warnings
from contextlib import ExitStack
from pathlib import Path

from common import now, write_json


def owned_path(root, relative):
    """Reject traversal and symlink parents; a leaf symlink can be unlinked."""
    root, relative = Path(root).absolute(), Path(relative)
    if root.is_symlink() or relative.is_absolute() or not relative.parts or '..' in relative.parts:
        raise ValueError('cleanup requires a path strictly inside its managed directory')
    path = root / relative
    if path == root or any(parent.is_symlink() for parent in path.parents if parent != root and parent.is_relative_to(root)):
        raise ValueError(f'unsafe cleanup path: {path}')
    return path


def cleanup(root, relatives, receipt, *, errors=()):
    """Called under the job lock, only after durable results have been verified.

    Cleanup errors leave the computation successful and can be retried. Never
    follow symlinks into reference databases or user-supplied reads.
    """
    report = {'state': 'complete', 'at': now(), 'removed': [], 'errors': list(errors)}
    for relative in relatives:
        try:
            path = owned_path(root, relative)
            if path.is_symlink():
                path.unlink()
            elif path.is_dir():
                shutil.rmtree(path)
            elif path.exists():
                path.unlink()
            else:
                continue
            report['removed'].append(str(relative))
        except (OSError, ValueError) as error:
            report['errors'].append({'path': str(relative), 'error': str(error)})
    if report['errors']:
        report['state'] = 'pending'
        warnings.warn(f'Intermediate cleanup incomplete; see {receipt}', RuntimeWarning)
    try:
        write_json(receipt, report)
    except OSError as error:
        report["state"] = "pending"
        warnings.warn(f"Could not record cleanup: {receipt}: {error}", RuntimeWarning)
    return report


def scratch_usage(path):
    """Count each inode once; estimate bytes freed only for its last links."""
    inodes = {}
    path = Path(path)
    if path.is_symlink():
        return {'files': 0, 'allocated_bytes': 0, 'reclaimable_bytes': 0}
    walk = ([(path.parent, [], [path.name])] if path.is_file() and not path.is_symlink()
            else os.walk(path, followlinks=False))
    for directory, _, names in walk:
        for name in names:
            item = Path(directory) / name
            info = item.lstat()
            if not stat.S_ISREG(info.st_mode):
                continue
            entry = inodes.setdefault((info.st_dev, info.st_ino), [info.st_blocks * 512, info.st_nlink, 0])
            entry[2] += 1
    return {'files': sum(entry[2] for entry in inodes.values()),
            'allocated_bytes': sum(entry[0] for entry in inodes.values()),
            'reclaimable_bytes': sum(entry[0] for entry in inodes.values() if entry[1] == entry[2])}


def inspect_job(kind, work_dir, output_dir, apply=False):
    """Explicit legacy cleanup; default is a read-only, locked preview."""
    base, output = Path(work_dir).absolute(), Path(output_dir).absolute()
    if base.is_symlink() or output.is_symlink():
        raise ValueError('cleanup roots must not be symlinks')
    base, output = base.resolve(), output.resolve()
    if base == output or base.is_relative_to(output) or output.is_relative_to(base):
        raise ValueError('cleanup work and output directories must not overlap')
    with ExitStack() as stack:
        for path in (base / '.lock', output.parent / f'.{output.name}.lock'):
            # Old ODB jobs did not use an output lock. Their base lock suffices.
            if not path.exists() and path != base / '.lock':
                continue
            lock = stack.enter_context(path.open('r'))
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        status = json.loads((base / 'status.json').read_text())
        if status.get('state') != 'success':
            raise ValueError('cleanup requires a successful, inactive job')
        fingerprint = status['fingerprint']
        if kind == 'odb':
            from run_odb_chunk import completed_output
            record = completed_output(output, fingerprint, base.name)
            work = base / fingerprint
        else:
            from run_kofam import _completed, completed_work
            record = _completed(output, fingerprint)
            work = completed_work(record, base) if record is not None else None
        if record is None:
            raise ValueError('published results are missing or changed; scratch retained')
        if work.parent != base or work.is_symlink():
            raise ValueError('scratch must be a direct managed child')
        if work.exists():
            identity = work / 'identity.json'
            if identity.exists():
                if json.loads(identity.read_text()) != record['identity']:
                    raise ValueError('scratch identity differs from published results')
            else:
                # rmtree can remove identity.json before encountering a locked
                # or inaccessible file. A recorded partial cleanup is retryable.
                receipt = base / 'cleanup.json'
                previous = json.loads(receipt.read_text()) if receipt.is_file() else {}
                if previous.get('state') != 'pending' or not any(
                        error.get('path') == work.name for error in previous.get('errors', [])):
                    raise ValueError('scratch identity missing without a pending cleanup')
            def protect(value):
                if isinstance(value, dict):
                    if 'path' in value and Path(value['path']).resolve().is_relative_to(work):
                        raise ValueError('scratch contains a recorded input; refusing cleanup')
                    for child in value.values(): protect(child)
                elif isinstance(value, list):
                    for child in value: protect(child)
            protect(record['identity'])
        report = {'kind': kind, 'work': str(work), 'output': str(output),
                  'apply': apply, **scratch_usage(work)}
        if apply:
            report['cleanup'] = cleanup(base, [work.name], base / 'cleanup.json')
        return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kind', choices=['odb', 'kofam'])
    parser.add_argument('--work-dir', required=True, help='Job directory containing status.json and .lock')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--apply', action='store_true', help='Delete verified scratch; otherwise only report sizes')
    print(json.dumps(inspect_job(**vars(parser.parse_args())), indent=2))
