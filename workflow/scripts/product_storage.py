"""Lossless gzip products with separate stored-byte and decoded-content identities."""
import gzip
import hashlib
from pathlib import Path

from common import atomic_writer, content_sha256


def decoded_sha256(path):
    digest = hashlib.sha256()
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def product_record(path):
    from dataset_assets import record, stat_identity
    entry = record(path)
    if str(path).endswith('.gz'):
        before = stat_identity(path)
        entry['content_sha256'] = decoded_sha256(path)
        if entry['stat'] != before or before != stat_identity(path):
            raise ValueError(f'file changed while registering content: {path}')
    return entry


def stage_product(entry, target, final=None):
    """Share matching encodings; convert legacy products without changing source files."""
    from dataset_assets import record, stage_file, stat_identity, verify
    source = verify(entry)
    target = Path(target)
    source_gzip, target_gzip = source.suffix == '.gz', target.suffix == '.gz'
    if source_gzip == target_gzip:
        saved = stage_file(entry, target, final)
        if source_gzip and 'content_sha256' not in saved:
            saved['content_sha256'] = decoded_sha256(target)
        return saved
    before = stat_identity(source)
    expected = content_sha256(entry)
    if source_gzip and 'content_sha256' not in entry:
        expected = decoded_sha256(source)
    if target.exists():
        saved = product_record(target)
        if content_sha256(saved) != expected:
            raise ValueError(f'refusing to replace different staged input: {target}')
    else:
        opener = gzip.open if source_gzip else open
        observed = hashlib.sha256()
        with opener(source, 'rb') as inp, atomic_writer(target, 'wb') as raw:
            zipped = gzip.GzipFile(filename='', fileobj=raw, mode='wb', compresslevel=6, mtime=0) if target_gzip else None
            output = zipped if zipped is not None else raw
            try:
                for chunk in iter(lambda: inp.read(1024 * 1024), b''):
                    observed.update(chunk)
                    output.write(chunk)
                if observed.hexdigest() != expected:
                    raise ValueError(f'product content changed while staging: {source}')
            finally:
                if zipped is not None:
                    zipped.close()
        saved = dict(record(target), content_sha256=expected) if target_gzip else record(target)
    if before != stat_identity(source):
        raise ValueError(f'registered file changed while staging: {source}')
    return dict(saved, path=str(final or target))
