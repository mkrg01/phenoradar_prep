"""Archive inactive successful native diagnostics, retaining products and observations."""
import hashlib
import json
import os
from pathlib import Path
import stat
import tarfile

from cleanup_work import owned_path
from common import atomic_writer, file_record, write_json


def files(work, protected, roots=('output', 'downloads')):
    protected = [Path(p).absolute() for p in protected]
    result = []
    for name in roots:
        root = owned_path(work, name)
        if root.is_symlink():
            raise ValueError(f'linked native diagnostic root: {root}')
        for directory, dirs, names in os.walk(root, followlinks=False):
            directory = Path(directory)
            dirs[:] = sorted(d for d in dirs if not (directory / d).is_symlink()
                             and (directory / d) != work / 'output/observations'
                             and not any(p == directory / d or p.is_relative_to(directory / d)
                                         for p in protected))
            for name in sorted(names):
                path = owned_path(work, (directory / name).relative_to(work))
                if (path in protected or path.is_symlink() or name.endswith('.lock')
                        or name.endswith(('.part', '.partial', '.tmp', '.fastq', '.fastq.gz', '.fq', '.fq.gz', '.sra'))):
                    continue
                if not stat.S_ISREG(path.lstat().st_mode):
                    raise ValueError(f'unsupported native diagnostic: {path}')
                result.append(path)
    return sorted(result)


def verify_archive(archive, entries):
    expected = {entry['path']: entry for entry in entries}
    with tarfile.open(archive, 'r:gz') as packed:
        seen = set()
        for member in packed:
            if not member.isfile() or member.name in seen or member.name not in expected:
                raise ValueError('unexpected or duplicate native archive member')
            digest = hashlib.sha256()
            with packed.extractfile(member) as source:
                for block in iter(lambda: source.read(1024 * 1024), b''):
                    digest.update(block)
            entry = expected[member.name]
            if member.size != entry['bytes'] or digest.hexdigest() != entry['sha256']:
                raise ValueError('native archive content differs from its inventory')
            seen.add(member.name)
        if seen != set(expected):
            raise ValueError('native archive has missing diagnostic files')


def compact(work, protected, *, apply=False):
    """Called under sample/native locks after all registered products verify.

    No scientific product or observation receipt moves. An interruption after
    archive commit can resume, verifying surviving source files before unlinking.
    """
    work = Path(work).absolute()
    archive = owned_path(work, 'diagnostics.tar.gz')
    receipt = owned_path(work, 'diagnostics.json')
    if archive.is_symlink() or receipt.is_symlink():
        raise ValueError('linked native archive or inventory')
    saved = json.loads(receipt.read_text()) if receipt.is_file() else None
    # Earlier archives remain readable/retryable. New archives also retire the
    # native content-digest index, which is a regenerable cache, not a product.
    roots = tuple(saved.get('roots', ('output', 'downloads'))) if saved else ('output', 'downloads', '.gg_cache')
    if any(name not in {'output', 'downloads', '.gg_cache'} for name in roots):
        raise ValueError('unsupported native archive roots')
    candidates = files(work, protected, roots)
    entries = [dict(file_record(p), path=str(p.relative_to(work))) for p in candidates]
    if saved:
        if file_record(archive) != saved['archive']:
            raise ValueError('native diagnostic archive changed')
        verify_archive(archive, saved['files'])
        expected = {entry['path']: entry for entry in saved['files']}
        if any(expected.get(entry['path']) != entry for entry in entries):
            raise ValueError('new or changed diagnostics after native archive commit')
    report = {'files': len(entries), 'logical_bytes': sum(e['bytes'] for e in entries),
              'archive': str(archive), 'apply': apply, 'state': 'available' if candidates else 'complete'}
    if not apply:
        return report
    if candidates and saved is None:
        # Existing uncommitted archives are replaced only after all sources verify.
        with atomic_writer(archive, 'wb') as raw:
            with tarfile.open(fileobj=raw, mode='w:gz', dereference=True) as packed:
                for path, entry in zip(candidates, entries):
                    packed.add(path, arcname=entry['path'], recursive=False)
        verify_archive(archive, entries)
        saved = {'schema_version': 1, 'roots': list(roots), 'archive': file_record(archive), 'files': entries}
        write_json(receipt, saved, compact=True)
    for path, entry in zip(candidates, entries):
        actual = dict(file_record(path), path=str(path.relative_to(work)))
        if actual != entry:
            raise ValueError('diagnostic changed before retirement')
        path.unlink()
    # Empty directories consume quota too; leave observation and product parents.
    for name in roots:
        for directory, dirs, names in os.walk(work / name, topdown=False, followlinks=False):
            path = Path(directory)
            if path.is_symlink() or path.is_relative_to(work / 'output/observations'):
                continue
            try:
                path.rmdir()
            except OSError:
                pass
    report['state'] = 'complete'
    return report
