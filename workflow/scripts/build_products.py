#!/usr/bin/env python3
"""Publish verified build products and expose them read-only to analyses."""
import argparse
import json
import os
import sqlite3
import tempfile
from pathlib import Path

from common import now, read_tsv, write_json, write_tsv
from dataset_assets import digest, link_file, locked, record, verify
from translate_cds import fasta_ids


def load_complete(path, verify_files=True):
    from portable_build import completion_path, load_products
    path = completion_path(path)
    if not path.is_file(): raise ValueError(f'build is incomplete: missing {path}; finish build through mapping first')
    data = json.loads(path.read_text())
    if data.get('kind') != 'completed_build' or data.get('schema_version') not in (1, 2):
        raise ValueError('unsupported build completion record')
    payload = {k:v for k,v in data.items() if k != 'sha256'}
    if digest(payload) != data.get('sha256'): raise ValueError('build completion record changed')
    if data['schema_version'] == 2:
        return load_products(path, data, verify_files=verify_files)
    if verify_files:
        for entry in data['files']: verify(entry)
    return data


def complete(path):
    from dataset import load, materialize, item_products
    path = Path(path).resolve()
    with locked(path / '.complete.lock'):
        target = path / 'completed.json'
        if target.exists():
            load_complete(target)
            return target
        manifest = load(path, check_code=True)
        inputs = materialize(path)
        root = Path(manifest['root'])
        results = root / 'results' / manifest['analysis']['run_name']
        mapping = results / 'orthogroups/mapping'
        required = [mapping / name for name in ('mappings.sqlite','gene_orthogroups.tsv','merge_qc.json')]
        required.extend([results / 'metadata/samples.tsv', mapping / 'incremental_plan/plan.json'])
        missing = [str(p) for p in required if not p.is_file()]
        if missing: raise ValueError('build mapping incomplete: ' + ', '.join(missing))
        samples = read_tsv(results / 'metadata/samples.tsv')
        names = {i['species'] for i in manifest['items']}
        if len(samples) != len(names) or {r['species'] for r in samples} != names:
            raise ValueError('build mapping does not cover the complete metadata species set')
        sample_map = {r['species']:r for r in samples}
        plan = json.loads((mapping / 'incremental_plan/plan.json').read_text())
        if plan['node'] != manifest['analysis']['odb']['node'] or plan['version'] != 'v12' or set(plan['proteins']) != names:
            raise ValueError('mapping plan differs from build reference or membership')
        qc = json.loads((mapping / 'merge_qc.json').read_text())
        if qc['reused_species'] + qc['mapped_species'] != len(names): raise ValueError('mapping QC does not cover the build')
        reference_hashes = set()
        for source in qc['sources']:
            origin = json.loads(verify(source).read_text())
            reference_hashes.update(origin.get('reference_sha256s', []))
            known = origin.get('reference_sha256') or origin.get('identity', {}).get('reference', {}).get('sha256')
            if known: reference_hashes.add(known)
        products = {}
        with sqlite3.connect(required[0].as_uri() + '?mode=ro', uri=True) as db:
            if db.execute('PRAGMA quick_check').fetchone()[0] != 'ok': raise ValueError('invalid build mapping database')
            if {r[0] for r in db.execute('SELECT DISTINCT species FROM genes')} != names:
                raise ValueError('mapping database species differ from build')
            if db.execute('SELECT 1 FROM mappings m LEFT JOIN genes g ON m.query=g.query WHERE g.query IS NULL LIMIT 1').fetchone():
                raise ValueError('mapping database contains unknown genes')
            for item in manifest['items']:
                species = item['species']; run = item['row']['run']
                product = item_products(manifest, item)
                odb_species = species.replace('-', '_')
                protein = results / 'proteins' / f'{odb_species}_protein.fa'
                provenance = protein.with_suffix('.json')
                translation = json.loads(provenance.read_text())
                if (translation['cds']['sha256'] != product['reference']['cds']['sha256'] or
                    translation['protein']['sha256'] != record(protein)['sha256'] or
                    translation['translation_table'] != manifest['analysis']['translation']['table'] or
                    translation['protein']['sha256'] != plan['proteins'][species]['sha256']):
                    raise ValueError(f'build translation differs from reference: {species}')
                if {r[0] for r in db.execute('SELECT query FROM genes WHERE species=?', (species,))} != set(fasta_ids(protein)):
                    raise ValueError(f'mapping gene set differs from protein: {species}')
                if sample_map[species]['run'] != run: raise ValueError(f'mapping run differs: {species}')
                products[species] = {'row':item['row'], 'reference_id':product['reference']['reference_id'],
                    'odb_species':odb_species, 'counts':product['busco']['counts'],
                    'conditions':{stage:product[key].get('provenance', {}).get('condition')
                                  for stage,key in [('assembly','reference'),('busco','busco'),('quant','quant')]},
                    'cds':record(inputs / 'cds' / f'{species}_longestCDS.fa.gz'),
                    'busco':record(inputs / 'busco/full' / f'{species}.busco.full.tsv'),
                    'abundance':record(inputs / 'quant' / species / run / f'{run}_abundance.tsv'),
                    'protein':record(protein), 'translation':record(provenance)}
        files = [record(p) for p in required]
        files.extend(json.loads((path / 'input_receipt.json').read_text())['files'])
        files.extend(manifest['auxiliary'].values())
        files.extend(record(path / n) for n in ('build.json','pipeline.yaml','checksums.json','metadata.tsv'))
        for p in products.values(): files.extend([p['protein'], p['translation']])
        data = {'schema_version':1, 'kind':'completed_build', 'build_id':manifest['name'], 'created_at':now(),
                'root':str(root), 'input':str(inputs), 'fields':manifest['fields'],
                'translation':manifest['analysis']['translation'], 'lineage':manifest['analysis']['phylogeny']['lineage'],
                'odb':{'version':'v12','node':manifest['analysis']['odb']['node']},
                'mapping':record(required[0]), 'products':products, 'files':files,
                'odb_reference_sha256s':sorted(reference_hashes),
                'excluded_runs':manifest.get('excluded', [])}
        from portable_build import publish_products, publish_pointer
        return publish_pointer(path, publish_products(path, data))


def import_protein(completion, species, protein, provenance):
    # The driver validates all products once; each species job checks its own files.
    data = load_complete(completion, verify_files=False)
    selected = [p for p in data['products'].values() if p['odb_species'] == species]
    if len(selected) != 1: raise ValueError(f'protein absent from completed build: {species}')
    link_file(verify(selected[0]['protein']), protein)
    link_file(verify(selected[0]['translation']), provenance)


def subset_mapping(completion, samples, database, mappings, qc):
    data = load_complete(completion)
    rows = read_tsv(samples)
    names = {r['species'] for r in rows}
    if not names or names - data['products'].keys(): raise ValueError('analysis species absent from completed build')
    for row in rows:
        product = data['products'][row['species']]
        if row['run'] != product['row']['run'] or record(row['cds'])['sha256'] != product['cds']['sha256']:
            raise ValueError('analysis input differs from completed build')
    target = Path(database)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.build-mapping-', dir=target.parent)
    os.close(fd)
    try:
        with sqlite3.connect(temporary, uri=True) as db:
            db.execute('ATTACH DATABASE ? AS source', (Path(data['mapping']['path']).as_uri() + '?mode=ro',))
            db.executescript('CREATE TABLE genes (query TEXT PRIMARY KEY, species TEXT NOT NULL); '
                             'CREATE INDEX genes_species ON genes(species); '
                             'CREATE TABLE mappings (query TEXT NOT NULL REFERENCES genes(query), og TEXT NOT NULL, PRIMARY KEY(query,og)) WITHOUT ROWID; '
                             'CREATE TEMP TABLE selected (species TEXT PRIMARY KEY);')
            db.executemany('INSERT INTO selected VALUES (?)', ((s,) for s in sorted(names)))
            db.execute('INSERT INTO genes SELECT g.query,g.species FROM source.genes g JOIN selected s ON g.species=s.species')
            db.execute('INSERT INTO mappings SELECT m.query,m.og FROM source.mappings m JOIN genes g ON m.query=g.query')
            count = db.execute('SELECT count(*) FROM mappings').fetchone()[0]
            ambiguous = db.execute('SELECT count(*) FROM (SELECT query FROM mappings GROUP BY query HAVING count(*)>1)').fetchone()[0]
            write_tsv(mappings, ['#query','ODB_OG'], ({'#query':q,'ODB_OG':og} for q,og in db.execute('SELECT query,og FROM mappings ORDER BY query,og')))
        os.replace(temporary, target)
        write_json(qc, {'created_at':now(), 'mode':'completed_build', 'build':record(completion),
                       'reused_species':len(names), 'mapped_species':0, 'unique_gene_og_pairs':count,
                       'ambiguous_genes':ambiguous, 'species':sorted(names)})
    finally:
        if Path(temporary).exists(): Path(temporary).unlink()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['protein','mapping'])
    parser.add_argument('--completion', required=True)
    for name in ('species','protein','provenance','samples','database','mappings','qc'): parser.add_argument('--' + name)
    args = parser.parse_args()
    if args.action == 'protein': import_protein(args.completion,args.species,args.protein,args.provenance)
    else: subset_mapping(args.completion,args.samples,args.database,args.mappings,args.qc)
