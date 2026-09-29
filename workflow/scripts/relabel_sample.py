"""Relabel sample-bound sequence IDs without changing sequences or numeric text."""
import gzip
import hashlib
import os
import re
import tempfile
from pathlib import Path

from common import file_record
from sample_identity import SAFE


def text_open(path):
    return gzip.open(path, 'rt', newline='') if str(path).endswith('.gz') else open(path, newline='')


def renamed_gene(value, old, new):
    if not re.fullmatch(re.escape(old) + r'_g[0-9]+(?::[0-9]+-[0-9]+)?', value):
        raise ValueError(f'gene does not belong to sample {old}: {value}')
    return new + value[len(old):]


def lines(source, old, new, kind):
    """Change only the identifier field; keep all other bytes of decoded text."""
    header = False
    with text_open(source) as handle:
        for line in handle:
            if kind == 'fasta':
                if line.startswith('>'):
                    match = re.match(r'>(\S+)(.*)', line, re.DOTALL)
                    yield '>' + renamed_gene(match[1], old, new) + match[2]
                else:
                    yield line
                continue
            if line.startswith('#') or not line.strip():
                yield line
                continue
            values = line.split('\t')
            if kind in ('quant', 'mapping') and not header:
                expected = 'target_id' if kind == 'quant' else 'gene_id'
                if values[0] != expected:
                    raise ValueError(f'{kind} table must start with {expected}: {source}')
                header = True
                yield line
                continue
            index = 2 if kind == 'busco' else 0
            if kind == 'busco' and len(values) > 1 and values[1].strip() == 'Missing':
                yield line
                continue
            if len(values) <= index:
                raise ValueError(f'malformed {kind} row: {source}')
            values[index] = renamed_gene(values[index], old, new)
            yield '\t'.join(values)


def relabel(source, target, old, new, kind='fasta'):
    """Publish a new file atomically, retaining a matching existing native file.

    The decoded output digest is compared with an independent read-back. Existing
    targets are never overwritten unless their decoded contents already match.
    """
    source, target = Path(source), Path(target)
    if source.resolve() == target.resolve():
        raise ValueError('relabel requires a separate destination')
    if not SAFE.fullmatch(old) or not SAFE.fullmatch(new):
        raise ValueError('unsafe sample identifier')
    if kind not in ('fasta', 'busco', 'quant', 'mapping'):
        raise ValueError('unknown relabel format')
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix='.relabel-', dir=target.parent)
    expected = hashlib.sha256()
    count = 0
    try:
        with os.fdopen(fd, 'wb') as raw:
            output = gzip.GzipFile(filename='', fileobj=raw, mode='wb', compresslevel=1, mtime=0) if target.suffix == '.gz' else raw
            try:
                for line in lines(source, old, new, kind):
                    data = line.encode('utf-8')
                    expected.update(data)
                    output.write(data)
                    count += 1
            finally:
                if output is not raw: output.close()
        # Compare decoded bytes, since the native tool may use different gzip headers.
        candidate = target if target.exists() else Path(name)
        actual = hashlib.sha256()
        opener = gzip.open if target.suffix == '.gz' else open
        with opener(candidate, 'rb') as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b''):
                actual.update(chunk)
        if actual.digest() != expected.digest():
            raise ValueError(f'changed or conflicting relabelled output: {target}')
        if not target.exists(): os.replace(name, target)
    finally:
        if Path(name).exists(): Path(name).unlink()
    return {'source': file_record(source), 'output': file_record(target),
            'old_id': old, 'new_id': new, 'kind': kind, 'lines': count,
            'decoded_sha256': expected.hexdigest()}
