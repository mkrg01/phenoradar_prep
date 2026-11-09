"""Serialize cooperating tasks' scratch reservations on each host/filesystem."""
from contextlib import contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import stat
import uuid

from common import now


class ScratchAdmissionError(ValueError):
    """No child was started because storage headroom could not be reserved."""


def free_space(path):
    info = os.statvfs(path)
    return {'bytes': info.f_bavail * info.f_frsize,
            'inodes': info.f_favail if info.f_files else None}


def headroom(path, config):
    free = free_space(path)
    for unit in ('bytes', 'inodes'):
        minimum = config.get('scratch_min_free_' + unit, 0)
        if minimum and (free[unit] is None or free[unit] < minimum):
            raise ScratchAdmissionError(f'Insufficient scratch free {unit}: {free[unit]}; need {minimum}')
    return free


def coordination_directory(device):
    # TMPDIR can differ for every allocation on the same node. A stable host
    # directory coordinates /data1, /tmp, and job TMPDIRs on the same device.
    host = hashlib.sha256(os.uname().nodename.encode()).hexdigest()[:12]
    directory = Path('/tmp') / f'phenoradar-scratch-{os.getuid()}-{host}-{device:x}'
    directory.mkdir(mode=0o700, exist_ok=True)
    info = directory.lstat()
    if not stat.S_ISDIR(info.st_mode) or info.st_uid != os.getuid() or info.st_mode & 0o077:
        raise ScratchAdmissionError('Unsafe scratch reservation directory')
    return directory


@contextmanager
def reserve(local, config):
    local = Path(local).resolve()
    requested = {unit: config.get('scratch_reserve_' + unit, 0) for unit in ('bytes', 'inodes')}
    if not any(requested.values()):
        headroom(local, config)
        yield None
        return
    device = local.stat().st_dev
    directory = coordination_directory(device)
    guard = directory / '.lock'
    if guard.is_symlink():
        raise ScratchAdmissionError('Linked scratch reservation lock')
    handle, lease = None, None
    with guard.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        remaining = {'bytes': 0, 'inodes': 0}
        for path in sorted(directory.glob('lease-*.json')):
            if path.is_symlink():
                raise ScratchAdmissionError('Linked scratch reservation')
            with path.open('r+') as active:
                try:
                    fcntl.flock(active, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    record = json.load(active)
                    if record.get('device') != device:
                        raise ScratchAdmissionError('Scratch reservation device differs')
                    # Count each reservation in full against current free space.
                    # This deliberately avoids racy credits for files being removed
                    # or hardlinked from outside a task's scratch directory.
                    for unit in remaining:
                        value = record['reserved'][unit]
                        if type(value) is not int or value < 0:
                            raise ScratchAdmissionError('Invalid scratch reservation size')
                        remaining[unit] += value
                else:
                    # A dead process releases flock. Retire its lease, preserving
                    # any orphan scratch and charging it through actual free space.
                    path.unlink()
        free = headroom(local, config)
        for unit in remaining:
            need = requested[unit] + remaining[unit] + config.get('scratch_min_free_' + unit, 0)
            if need and (free[unit] is None or need > free[unit]):
                raise ScratchAdmissionError(f'Scratch reservation unavailable: need {need} {unit}, free {free[unit]}; no child started')
        lease = directory / ('lease-' + uuid.uuid4().hex + '.json')
        handle = lease.open('x+')
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            json.dump({'schema_version': 1, 'scratch': str(local), 'device': device,
                       'reserved': requested, 'pid': os.getpid(), 'created_at': now(),
                       'job_id': os.environ.get('SLURM_JOB_ID') or os.environ.get('JOB_ID')}, handle)
            handle.flush()
        except BaseException:
            handle.close()
            lease.unlink(missing_ok=True)
            raise
    try:
        yield lease
    finally:
        with guard.open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            lease.unlink(missing_ok=True)
            handle.close()
