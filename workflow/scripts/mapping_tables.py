"""Immutable, species-local gene/OG tables; no combined mapping database."""
import argparse
from collections import OrderedDict, defaultdict
from contextlib import contextmanager
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
import csv
import fcntl
import gzip
import json
import os
from pathlib import Path
import shutil
import tempfile

from common import file_record, read_tsv, species_from_gene_id, write_json
from dataset_assets import SAFE, digest, link_file, record, stat_identity, verify
from translate_cds import fasta_ids


@contextmanager
def lock(path):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def relative_file(root, entry):
    relative = Path(entry['path'])
    path = (Path(root) / relative).resolve()
    if relative.is_absolute() or '..' in relative.parts or not path.is_relative_to(Path(root).resolve()):
        raise ValueError('mapping table path escapes snapshot')
    return dict(entry, path=str(path))


def checked(entry, cache_dir):
    """Hash unchanged imported files once, retaining their filesystem identity."""
    path = Path(entry['path']).resolve()
    receipt = Path(cache_dir) / '.checks' / (digest([str(path), entry['sha256']]) + '.json')
    if receipt.is_file():
        saved = json.loads(receipt.read_text())
        if saved['sha256'] != entry['sha256'] or saved['path'] != str(path):
            raise ValueError('mapping verification receipt differs')
        try: verify(saved)
        except ValueError as error: raise ValueError(f'ODB result changed after completion: {path}') from error
        if saved.get('stat') == stat_identity(path): return saved
    saved = record(path)
    if saved['sha256'] != entry['sha256']:
        raise ValueError(f'ODB result changed after completion: {path}')
    write_json(receipt, saved)
    return saved


def protein_record(path):
    path = Path(path).resolve()
    provenance = path.with_suffix('.json')
    if provenance.is_file():
        entry = dict(json.loads(provenance.read_text())['protein'], path=str(path))
        verify(entry)
        return dict(entry, stat=stat_identity(path))
    return record(path)


@lru_cache(maxsize=8)
def _read_tables(path, identity):
    path = Path(path)
    data = json.loads(path.read_text())
    if data.get('schema_version') != 2 or data.get('kind') != 'odb_tables':
        raise ValueError(f'expected species mapping tables: {path}')
    proteins = {p['species']: p for p in data['proteins']}
    if len(proteins) != len(data['proteins']) or set(proteins) != set(data['tables']):
        raise ValueError('mapping species differ from protein inventory')
    for name, entry in data['tables'].items():
        if (not SAFE.fullmatch(name) or entry['odb_species'] != name.replace('-', '_')
                or entry['protein_sha256'] != proteins[name]['sha256']):
            raise ValueError('mapping table identity differs from protein')
        relative_file(path.parent, entry['table'])
    return data


def load_tables(path, verify_files=True):
    path = Path(path).resolve()
    data = _read_tables(str(path), tuple(stat_identity(path)))
    if verify_files:
        for entry in data['tables'].values(): verify(relative_file(path.parent, entry['table']))
    return data


def read_species(snapshot, species):
    """Load at most one species, preserving unmapped genes and multiple OGs."""
    data = load_tables(snapshot, verify_files=False)
    if species not in data['tables']: raise ValueError(f'species absent from mappings: {species}')
    entry = data['tables'][species]
    path = verify(relative_file(Path(snapshot).parent, entry['table']))
    genes = defaultdict(list)
    with gzip.open(path, 'rt', newline='') as handle:
        rows = csv.DictReader(handle, delimiter='\t')
        if rows.fieldnames != ['gene_id', 'orthogroup']:
            raise ValueError(f'invalid species mapping columns: {path}')
        for row in rows:
            gene, og = row['gene_id'], row['orthogroup']
            if species_from_gene_id(gene) != species or og is None:
                raise ValueError(f'mapping gene species differs: {gene}')
            genes[gene]
            if og: genes[gene].append(og)
    return dict(genes), entry


def annotation_pairs(path):
    indices = None; saw_empty = False
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip(): continue
            fields = line.rstrip('\r\n').split('\t')
            if line.startswith('#'):
                if '#query' in fields and 'ODB_OG' in fields:
                    indices = fields.index('#query'), fields.index('ODB_OG')
                elif '#cluster_id' in fields and 'gene_id' in fields:
                    indices = fields.index('gene_id'), fields.index('#cluster_id')
                elif line.startswith('# No data found'): saw_empty = True
                continue
            if indices is None or len(fields) <= max(indices):
                raise ValueError(f'unrecognized ODB table at {path}:{number}')
            query, og = (fields[i] for i in indices)
            if not query or not og or og in {'-', 'NA', 'null'}:
                raise ValueError(f'empty query/OG at {path}:{number}')
            yield query, og
    if indices is None and not saw_empty:
        raise ValueError(f'ODB table has no recognized header: {path}')


def partition(annotation, members, cache_dir):
    """Split each legacy/chunk annotation once, with a bounded file-handle pool."""
    identity = digest([annotation['sha256'], sorted(members)])
    target = Path(cache_dir) / '.partitions' / identity
    with lock(target.with_suffix('.lock')):
        if (target / 'receipt.json').is_file(): return target
        print(f'Partitioning ODB annotations: {annotation["path"]}', flush=True)
        with tempfile.TemporaryDirectory(prefix='.partition-', dir=target.parent) as tmp:
            staging = Path(tmp); handles = OrderedDict(); counts = defaultdict(int)
            try:
                for gene, og in annotation_pairs(annotation['path']):
                    try: name = species_from_gene_id(gene)
                    except ValueError as error: raise ValueError(f'query does not belong to snapshot species: {gene}') from error
                    if name not in members: raise ValueError(f'query does not belong to snapshot species: {gene}')
                    handle = handles.pop(name, None)
                    if handle is None:
                        if len(handles) >= 64: handles.popitem(last=False)[1].close()
                        handle = open(staging / f'{name}.tsv', 'a')
                    handles[name] = handle
                    handle.write(f'{gene}\t{og}\n'); counts[name] += 1
            finally:
                for handle in handles.values(): handle.close()
            verify(annotation)
            shards = {s:dict(record(staging / f'{s}.tsv'), path=f'{s}.tsv') for s in counts}
            write_json(staging / 'receipt.json', {'annotation':annotation, 'counts':dict(counts), 'shards':shards})
            os.rename(staging, target)
    return target


def species_table(name, odb_name, protein, annotation, parts, cache_dir, version, node, shard=None):
    identity = digest([name, protein['sha256'], annotation['sha256'], version, node])
    target = Path(cache_dir) / '.tables' / identity
    with lock(target.with_suffix('.lock')):
        receipt = target / 'receipt.json'
        if receipt.is_file():
            entry = json.loads(receipt.read_text())
            if entry['protein_sha256'] != protein['sha256'] or entry['annotation_sha256'] != annotation['sha256'] or entry['odb_species'] != odb_name:
                raise ValueError('cached species mapping identity differs')
            verify(relative_file(target, entry['table']))
            return dict(entry, table=relative_file(target, entry['table']))
        genes = fasta_ids(protein['path'])
        if any(species_from_gene_id(g) != name for g in genes):
            raise ValueError(f'FASTA gene species differs: {name}')
        allowed = set(genes); pairs = defaultdict(set); rows = 0
        raw = parts / f'{name}.tsv'
        if raw.is_file():
            verify(relative_file(parts, shard))
            with raw.open() as handle:
                for line in handle:
                    gene, og = line.rstrip('\n').split('\t')
                    if gene not in allowed:
                        raise ValueError(f'query does not belong to this species FASTA inputs: {gene}')
                    pairs[gene].add(og); rows += 1
        unique = sum(map(len, pairs.values()))
        with tempfile.TemporaryDirectory(prefix='.table-', dir=target.parent) as tmp:
            staging = Path(tmp); output = staging / 'table.tsv.gz'
            # Deterministic bytes allow safe reuse and transport between projects.
            with output.open('wb') as raw_out, gzip.GzipFile(fileobj=raw_out, mode='wb', filename='', mtime=0, compresslevel=1) as handle:
                handle.write(b'gene_id\torthogroup\n')
                for gene in sorted(genes):
                    for og in sorted(pairs.get(gene, ())) or ['']:
                        handle.write(f'{gene}\t{og}\n'.encode())
            verify(protein)
            entry = {'odb_species':odb_name, 'protein_sha256':protein['sha256'],
                     'table':dict(record(output), path='table.tsv.gz'),
                     'qc':{'protein_genes':len(genes), 'input_annotation_rows':rows,
                           'unique_gene_og_pairs':unique, 'duplicate_pairs_removed':rows-unique,
                           'ambiguous_genes':sum(len(v)>1 for v in pairs.values())},
                     'annotation_sha256':annotation['sha256']}
            write_json(staging / 'receipt.json', entry)
            os.rename(staging, target)
        return dict(entry, table=relative_file(target, entry['table']))


def write_tables(outdir, entries, version='v12', node=3193, qc=None, references=()):
    """Publish only selected species, by linking validated immutable tables."""
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.tables-', dir=out) as tmp:
        staged = Path(tmp) / 'species'; staged.mkdir()
        saved = {}
        for name, entry in sorted(entries.items()):
            filename = entry['odb_species'] + '.tsv.gz'
            source = verify(entry['table']); destination = staged / filename
            link_file(source, destination)
            saved[name] = dict(entry, table=dict(entry['table'], path=f'species/{filename}', stat=stat_identity(destination)))
        if (out / 'species').exists(): shutil.rmtree(out / 'species')
        os.rename(staged, out / 'species')
        data = {'schema_version':2, 'kind':'odb_tables', 'version':version, 'node':node,
                'proteins':[{'species':s,'odb_species':e['odb_species'],'sha256':e['protein_sha256']} for s,e in saved.items()],
                'tables':saved, 'reference_sha256s':sorted(set(references)), 'qc':qc or {}}
        write_json(out / 'snapshot.json', data)
    return out / 'snapshot.json'


def subset(snapshot, species, outdir):
    data = load_tables(snapshot, verify_files=False)
    names = set(species)
    if not names or names - data['tables'].keys(): raise ValueError('selected species absent from mapping tables')
    entries = {s:dict(data['tables'][s], table=relative_file(Path(snapshot).parent, data['tables'][s]['table'])) for s in names}
    qc = {'mode':'existing', 'reused_species':len(names), 'mapped_species':0}
    return write_tables(outdir, entries, data['version'], data['node'], qc, data.get('reference_sha256s', []))


def collect(samples, chunks, chunk_dir, protein_dir, outdir, cache_dir, existing=None,
            version='v12', node=3193, source_plan=None, threads=1):
    from incremental_odb import load_snapshot
    if type(threads) is not int or threads < 1: raise ValueError('threads must be positive')
    species = {r['species']:r['odb_species'] for r in read_tsv(samples)}
    proteins = {s:protein_record(Path(protein_dir)/f'{odb}_protein.fa') for s,odb in species.items()}
    if source_plan:
        plan = json.loads(Path(source_plan).read_text())
        if plan['version'] != version or plan['node'] != node or set(plan['proteins']) != set(species):
            raise ValueError('ODB plan differs from selected inputs')
        if any(proteins[s]['sha256'] != plan['proteins'][s]['sha256'] for s in species):
            raise ValueError('protein changed after ODB planning')
        members = plan['sources']
    elif existing:
        members = [{'kind':'existing','root':str(existing),'species':sorted(species)}]
    else:
        members = [{'kind':'mapped', **c} for c in json.loads(Path(chunks).read_text())]
    if sorted(s for c in members for s in c['species']) != sorted(species):
        raise ValueError('chunk membership differs from selected species')
    entries = {}; references = set(); excluded = 0
    for chunk in members:
        if chunk['kind'] == 'existing':
            root = Path(chunk['root']); snapshot, original = load_snapshot(root, version, node)
            if chunk.get('snapshot_sha256') and file_record(root/'snapshot.json')['sha256'] != chunk['snapshot_sha256']:
                raise ValueError('ODB snapshot changed after planning')
            if set(chunk['species']) - original.keys(): raise ValueError('existing ODB snapshot has missing selected species')
            for name in chunk['species']:
                if original[name]['odb_species'] != species[name] or original[name]['sha256'] != proteins[name]['sha256']:
                    raise ValueError(f'protein differs from existing ODB input: {name}; remapping is required')
            references.update(snapshot.get('reference_sha256s', []))
            if snapshot.get('reference_sha256'): references.add(snapshot['reference_sha256'])
            if snapshot['schema_version'] == 2:
                for name in chunk['species']:
                    entry = snapshot['tables'][name]
                    entries[name] = dict(entry, table=relative_file(root, entry['table']))
                continue
            annotation = checked(dict(snapshot['annotations'], path=str(root/'annotations.tsv')), cache_dir)
        else:
            root = Path(chunk_dir)/chunk['chunk']; provenance = json.loads((root/'provenance.json').read_text())
            for result in provenance['results']: checked(dict(result,path=str(root/Path(result['path']).name)), cache_dir)
            annotation = next(r for r in provenance['results'] if r['path'].endswith('.og.annotations'))
            annotation = checked(dict(annotation,path=str(root/Path(annotation['path']).name)), cache_dir)
            original = {s:None for s in chunk['species']}
            references.add(provenance['identity']['reference']['sha256'])
        parts = partition(annotation, original, cache_dir)
        partitions = json.loads((parts/'receipt.json').read_text())
        counts = partitions['counts']
        excluded += sum(n for s,n in counts.items() if s not in chunk['species'])
        with ProcessPoolExecutor(max_workers=threads) as pool:
            pending = {pool.submit(species_table, name, species[name], proteins[name], annotation,
                                   parts, cache_dir, version, node, partitions['shards'].get(name)):name for name in chunk['species']}
            for index, future in enumerate(as_completed(pending), 1):
                entries[pending[future]] = future.result()
                if index % 100 == 0 or index == len(pending):
                    print(f'Mapping tables: {index}/{len(pending)} species ready', flush=True)
    kinds = {c['kind'] for c in members}
    qc = {'mode':'mixed' if len(kinds)>1 else 'existing' if kinds=={'existing'} else 'mapped',
          'reused_species':sum(len(c['species']) for c in members if c['kind']=='existing'),
          'mapped_species':sum(len(c['species']) for c in members if c['kind']=='mapped'),
          'excluded_annotation_rows':excluded}
    for key in ('input_annotation_rows','unique_gene_og_pairs','duplicate_pairs_removed','ambiguous_genes'):
        qc[key] = sum(e['qc'][key] for e in entries.values())
    return write_tables(outdir, entries, version, node, qc, references)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('samples','chunks','chunk-dir','protein-dir','outdir','cache-dir'): parser.add_argument('--'+key, required=True)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--existing'); parser.add_argument('--source-plan')
    parser.add_argument('--version',default='v12'); parser.add_argument('--node',type=int,default=3193)
    collect(**vars(parser.parse_args()))
