"""Reuse translations bound to CDS content and genetic code."""
import fcntl
import json
import os
import tempfile
from pathlib import Path
from contextlib import contextmanager

from common import file_record, now, write_json
from dataset_assets import link_file, record, verify
from translate_cds import translate


@contextmanager
def cache_lock(path):
    # Concurrent builds may request the same translation; wait for its publisher.
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def cache_path(cache_dir, cds, table):
    if type(table) is not int or table < 1:
        raise ValueError('translation table must be a positive integer')
    return Path(cache_dir) / cds['sha256'] / f'table_{table}'


def load_cached(path, cds, table):
    data = json.loads((path / 'receipt.json').read_text())
    if (data.get('schema_version') != 1 or data['cds_sha256'] != cds['sha256']
            or data['translation_table'] != table or data['protein']['path'] != 'protein.fa'):
        raise ValueError(f'translation cache identity differs: {path}')
    verify(dict(data['protein'], path=str(path / 'protein.fa')))
    return data


def publish(path, cds, protein, provenance):
    info = json.loads(Path(provenance).read_text())
    source = record(protein)
    if info['cds']['sha256'] != cds['sha256'] or info['protein']['sha256'] != source['sha256']:
        raise ValueError('translation provenance differs from CDS/protein')
    # Matching input/output hashes retain the translation's recorded validation.
    verify(cds)
    with tempfile.TemporaryDirectory(prefix='.protein-', dir=path.parent) as tmp:
        staged = Path(tmp) / 'entry'
        staged.mkdir()
        link_file(verify(source), staged / 'protein.fa')
        data = {'schema_version': 1, 'cds_sha256': cds['sha256'],
                'translation_table': info['translation_table'], 'seqkit': info['seqkit'],
                'sequences': info['sequences'],
                'protein': dict(record(staged / 'protein.fa'), path='protein.fa')}
        write_json(staged / 'receipt.json', data)
        os.rename(staged, path)
    return data


def register_translation(cds, protein, provenance, cache_dir, table=1):
    """Import a translation only when its recorded input and output hashes match."""
    cds = record(cds)
    info = json.loads(Path(provenance).read_text())
    if info['translation_table'] != table:
        raise ValueError('translation table differs from imported protein')
    if info['cds']['sha256'] != cds['sha256'] or info['protein']['sha256'] != file_record(protein)['sha256']:
        raise ValueError('translation provenance differs from CDS/protein')
    path = cache_path(cache_dir, cds, table)
    with cache_lock(path.with_suffix('.lock')):
        if path.exists():
            data = load_cached(path, cds, table)
            if data['protein']['sha256'] != info['protein']['sha256']:
                raise ValueError('conflicting translation for registered CDS/genetic code')
        else:
            publish(path, cds, protein, provenance)
    return path


def cached_translate(cds, output, provenance, cache_dir, seqkit='seqkit', table=1, threads=1):
    cds = record(cds)
    path = cache_path(cache_dir, cds, table)
    with cache_lock(path.with_suffix('.lock')):
        reused = path.exists()
        if reused:
            data = load_cached(path, cds, table)
        else:
            with tempfile.TemporaryDirectory(prefix='.translate-', dir=path.parent) as tmp:
                protein, receipt = Path(tmp) / 'protein.fa', Path(tmp) / 'translation.json'
                translate(cds['path'], protein, receipt, seqkit=seqkit, table=table, threads=threads)
                data = publish(path, cds, protein, receipt)
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.protein-', dir=output.parent) as tmp:
            staged = Path(tmp) / 'protein.fa'
            link_file(path / 'protein.fa', staged)
            os.replace(staged, output)
        verify(cds)
        write_json(provenance, {'created_at': now(), 'cds': file_record(cds['path']),
                    'protein': record(output), 'seqkit': data['seqkit'],
                    'translation_table': table, 'sequences': data['sequences'], 'reused': reused})
    return reused
