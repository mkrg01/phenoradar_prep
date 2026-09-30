#!/usr/bin/env python3
"""Publish verified species products and expose selected tables to analyses."""
import argparse
import json
from pathlib import Path

from common import atomic_writer, now, read_tsv, write_json
from dataset_assets import digest, link_file, locked, record, verify
from mapping_tables import load_tables, protein_record, relative_file, subset
from layout import run_layout


def load_complete(path, verify_files=True):
    from portable_build import completion_path, load_products
    path = completion_path(path)
    data = json.loads(path.read_text())
    if data.get('kind') != 'completed_build' or data.get('schema_version') not in (3, 4, 5):
        raise ValueError('unsupported completed products: use a schema-3/4 products bundle or schema-5 database, or prepare a new build from metadata')
    if digest({k:v for k,v in data.items() if k != 'sha256'}) != data.get('sha256'):
        raise ValueError('build completion record changed')
    return load_products(path, data, verify_files=verify_files)


def complete(path):
    from dataset import load, materialize, item_products
    from portable_build import publish_products, publish_pointer
    path = Path(path).resolve()
    with locked(path / '.complete.lock'):
        target = path / 'completed.json'
        if target.exists():
            load_complete(target); return target
        manifest = load(path, check_code=True)
        inputs = materialize(path)
        results = Path(manifest['root']) / run_layout(manifest['analysis'])[0]
        with_expression = bool(manifest['analysis'].get('output_root'))
        mapping = results / 'orthogroups/mapping/snapshot.json'
        if not mapping.is_file(): raise ValueError(f'build mapping incomplete: missing {mapping}')
        tables = load_tables(mapping)
        samples = read_tsv(results / 'metadata/samples.tsv')
        names = {i['species'] for i in manifest['items']}
        if len(samples) != len(names) or {r['species'] for r in samples} != names or set(tables['tables']) != names:
            raise ValueError('build mapping does not cover the complete metadata species set')
        if tables['node'] != manifest['analysis']['odb']['ncbi_tax_id'] or tables['version'] != 'v12':
            raise ValueError('mapping reference differs from build')
        sample_map = {r['species']:r for r in samples}
        input_files = json.loads((path / 'input_receipt.json').read_text())['files']
        inventory = {entry['path']:entry for entry in input_files}
        products = {}
        for item in manifest['items']:
            species = item['species']; run = item['row']['run']
            product = item_products(manifest, item)
            odb_species = species.replace('-', '_')
            protein = results / 'proteins' / f'{odb_species}_protein.fa'
            provenance = protein.with_suffix('.json')
            translation = json.loads(provenance.read_text())
            protein_entry = protein_record(protein)
            if (translation['cds']['sha256'] != product['reference']['cds']['sha256'] or
                translation['protein']['sha256'] != protein_entry['sha256'] or
                translation['translation_table'] != manifest['analysis']['translation']['table'] or
                protein_entry['sha256'] != tables['tables'][species]['protein_sha256']):
                raise ValueError(f'build translation differs from reference: {species}')
            if sample_map[species]['run'] != run: raise ValueError(f'mapping run differs: {species}')
            products[species] = {'row':item['row'], 'reference_id':product['reference']['reference_id'],
                'odb_species':odb_species, 'counts':product['busco']['counts'],
                'conditions':{stage:product[key].get('provenance', {}).get('condition')
                              for stage,key in [('assembly','reference'),('busco','busco'),('quant','quant')]},
                'raw_inputs':product['quant'].get('provenance', {}).get('raw_inputs', {}),
                'cds':inventory[str(inputs / 'cds' / f'{species}_longestCDS.fa.gz')],
                'busco':inventory[str(inputs / 'busco/full' / f'{species}.busco.full.tsv')],
                'abundance':inventory[str(inputs / 'quant' / species / run / f'{run}_abundance.tsv')],
                'protein':protein_entry, 'translation':record(provenance)}
            if with_expression:
                expression = results / 'orthogroups/expression/runs' / f'{run}.tsv'
                qc_path = expression.with_suffix('.qc.json')
                if not expression.is_file() or not qc_path.is_file():
                    raise ValueError(f'build expression incomplete: {run}; run the database target before complete')
                expression_entry = record(expression)
                qc = json.loads(qc_path.read_text())
                expected_table = tables['tables'][species]['table']
                if (qc.get('species') != species or qc.get('run') != run or
                    qc.get('multimap') != 'error' or
                    qc.get('abundance', {}).get('sha256') != products[species]['abundance']['sha256'] or
                    qc.get('mapping_table', {}).get('sha256') != expected_table['sha256'] or
                    qc.get('expression', {}).get('sha256') != expression_entry['sha256']):
                    raise ValueError(f'build expression differs from inputs: {species}')
                products[species].update(expression=expression_entry, expression_qc=record(qc_path))
        files = [*input_files, record(mapping), record(results / 'metadata/samples.tsv')]
        files.extend(relative_file(mapping.parent, e['table']) for e in tables['tables'].values())
        files.extend(manifest['auxiliary'].values())
        files.extend(record(path / n) for n in ('build.json','pipeline.yaml','checksums.json','metadata.tsv'))
        for p in products.values(): files.extend([p['protein'], p['translation']])
        data = {'schema_version':5 if with_expression else 4, 'kind':'completed_build', 'build_id':manifest['name'], 'created_at':now(),
                'input':str(inputs), 'fields':manifest['fields'],
                'translation':manifest['analysis']['translation'], 'lineage':manifest['analysis']['phylogeny']['lineage'],
                'odb':{'version':'v12','node':manifest['analysis']['odb']['ncbi_tax_id']},
                'mapping':record(mapping), 'products':products, 'files':files,
                'excluded_runs':manifest.get('excluded', [])}
        if with_expression:
            # Preserve the schema-5 provenance field for existing bundles.
            data['tpm'] = {'multimap': 'error'}
        return publish_pointer(path, publish_products(path, data))


def import_protein(completion, species, protein, provenance):
    data = load_complete(completion, verify_files=False)
    selected = [p for p in data['products'].values() if p['odb_species'] == species]
    if len(selected) != 1: raise ValueError(f'protein absent from completed build: {species}')
    link_file(verify(selected[0]['protein']), protein)
    with atomic_writer(provenance, 'wb') as handle:
        handle.write(verify(selected[0]['translation']).read_bytes())


def subset_mapping(completion, samples, outdir):
    data = load_complete(completion, verify_files=False)
    rows = read_tsv(samples); names = {r['species'] for r in rows}
    if not names or names - data['products'].keys(): raise ValueError('analysis species absent from completed build')
    for row in rows:
        product = data['products'][row['species']]
        if row['run'] != product['row']['run']:
            raise ValueError('analysis input differs from completed build')
        verify(dict(product['cds'], path=row['cds']))
    return subset(verify(data['mapping']), names, outdir)


def import_expression(completion, samples, run, output, qc):
    data = load_complete(completion, verify_files=False)
    if data.get('tpm', {}).get('multimap') != 'error':
        raise ValueError('database TPM policy must be error')
    rows = [row for row in read_tsv(samples) if row['run'] == run]
    if len(rows) != 1: raise ValueError(f'expected exactly one selected sample for {run}')
    row = rows[0]
    product = data['products'].get(row['species'])
    if product is None or product['row']['run'] != run:
        raise ValueError('expression sample is absent from the database')
    verify(dict(product['abundance'], path=row['abundance']))
    link_file(verify(product['expression']), output)
    # A separate small receipt prevents output touches from changing another input's mtime.
    with atomic_writer(qc, 'wb') as handle:
        handle.write(verify(product['expression_qc']).read_bytes())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['protein','mapping','expression'])
    parser.add_argument('--completion', required=True)
    for name in ('species','protein','provenance','samples','outdir','run','output','qc'): parser.add_argument('--' + name)
    args = parser.parse_args()
    if args.action == 'protein': import_protein(args.completion,args.species,args.protein,args.provenance)
    elif args.action == 'mapping': subset_mapping(args.completion,args.samples,args.outdir)
    else: import_expression(args.completion,args.samples,args.run,args.output,args.qc)
