"""Retire a completed build's regenerable ODB derivations, retaining publications."""
from contextlib import ExitStack
import fcntl
import json
from pathlib import Path
import re

from cleanup_work import cleanup, owned_path


def snapshot_paths(value, base):
    if isinstance(value, dict):
        if 'path' in value:
            yield (base / value['path']).resolve()
        for child in value.values():
            yield from snapshot_paths(child, base)
    elif isinstance(value, list):
        for child in value:
            yield from snapshot_paths(child, base)


def retire(build, protected, *, apply=False):
    build = Path(build).absolute()
    root = owned_path(build, 'work/cache/odb')
    if root.is_symlink():
        raise ValueError('linked private ODB cache')
    protected = [Path(p).absolute() for p in protected]
    # Published cache snapshots can still refer to a derivation. Such entries
    # remain live even if a database currently has a separate durable copy.
    for snapshot in root.glob('*/*/snapshot.json'):
        if snapshot.is_symlink() or any(p.is_symlink() for p in snapshot.parents if p.is_relative_to(root)):
            raise ValueError('linked mapping snapshot')
        protected.extend(snapshot_paths(json.loads(snapshot.read_text()), snapshot.parent))
    targets, retained = [], []
    with ExitStack() as locks:
        for name in ('.checks', '.partitions', '.tables'):
            directory = owned_path(build, root.relative_to(build) / name)
            if directory.is_symlink():
                raise ValueError('linked derived cache')
            if not directory.is_dir():
                continue
            for path in sorted(directory.iterdir()):
                if not re.fullmatch(r'[0-9a-f]{64}' + (r'\.json' if name == '.checks' else ''), path.name):
                    continue  # Never unlink lock names, including waiting writers' locks.
                if path.is_symlink() or any(p.is_symlink() for p in path.rglob('*')):
                    raise ValueError('linked derived cache member')
                if any(p == path or p.is_relative_to(path) for p in protected):
                    retained.append(str(path.relative_to(build)))
                    continue
                lock = path.with_suffix('.lock')
                if name != '.checks' and lock.exists():
                    if lock.is_symlink():
                        raise ValueError('linked derivation lock')
                    handle = locks.enter_context(lock.open('r'))
                    fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                targets.append(path.relative_to(build))
        report = {'targets': [str(t) for t in targets], 'retained': retained, 'apply': apply}
        if apply:
            report['cleanup'] = cleanup(build, targets, build / 'jobs/cleanup/derived-cache.json')
        return report
