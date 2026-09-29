"""Self-contained completed products with paths relative to their manifest."""
import copy
import json
import os
import shutil
import tempfile
from pathlib import Path

from common import write_json, write_tsv
from dataset_assets import digest, link_file, record, stat_identity, verify

PRODUCT_FILES = ('cds', 'busco', 'abundance', 'protein', 'translation')
EXPRESSION_FILES = ('expression', 'expression_qc')


def inside_bundle(root, relative):
    relative = Path(relative)
    if relative.is_absolute() or '..' in relative.parts or relative == Path('.'):
        raise ValueError(f'invalid relative product path: {relative}')
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f'product path escapes bundle: {relative}')
    return path


def completion_path(path):
    path = Path(path).resolve()
    if path.is_dir():
        if (path / 'manifest.json').is_file():
            path /= 'manifest.json'
        elif (path / 'completed.json').is_file():
            path /= 'completed.json'
        elif (path / 'products/manifest.json').is_file():
            path /= 'products/manifest.json'
        else:
            path /= 'database/manifest.json'
    if not path.is_file():
        raise ValueError(f'build is incomplete: missing {path}; finish build through database first')
    data = json.loads(path.read_text())
    if data.get('kind') == 'completed_build_pointer':
        if data.get('schema_version') != 1 or digest({k:v for k,v in data.items() if k != 'sha256'}) != data.get('sha256'):
            raise ValueError('build completion record changed')
        entry = dict(data['manifest'], path=str(inside_bundle(path.parent, data['manifest']['path'])))
        path = verify(entry)
    return path


def load_products(path, data, verify_files=True):
    """Bind portable paths at read time; the on-disk manifest stays unchanged."""
    root = Path(path).resolve().parent
    # Metadata rows can be large; copy only records whose paths are rebound.
    bound = dict(data)
    bound['files'] = [dict(entry) for entry in data['files']]
    bound['products'] = {species:dict(product) for species,product in data['products'].items()}
    files = {}
    for entry in bound['files']:
        relative = entry['path']
        if relative in files:
            raise ValueError('duplicate product manifest path')
        entry['path'] = str(inside_bundle(root, relative))
        files[relative] = entry
        if verify_files: verify(entry)

    def bind(entry):
        stored = files.get(entry['path'])
        if stored is None or any(entry.get(k) != stored.get(k) for k in ('sha256', 'bytes')):
            raise ValueError('product is absent from manifest inventory')
        return copy.deepcopy(stored)

    bound['mapping'] = bind(bound['mapping'])
    for species, product in bound['products'].items():
        if data['schema_version'] >= 4:
            from sample_identity import annotate
            if annotate(product['row']) != product['row'] or product['row']['analysis_sample_id'] != species:
                raise ValueError('product sample identity differs from metadata')
        if product['row'].get('analysis_sample_id', product['row']['scientific_name'].replace(' ', '_')) != species:
            raise ValueError('product species differs from metadata')
        keys = PRODUCT_FILES + (EXPRESSION_FILES if data['schema_version'] >= 5 else ())
        for key in keys:
            if key not in product: raise ValueError(f'incomplete database product: {species}: {key}')
            product[key] = bind(product[key])
    for required in ('metadata.tsv', 'busco/summary.tsv', 'excluded_accessions.tsv', 'excluded_runs.tsv'):
        if required not in files: raise ValueError(f'incomplete product bundle: {required}')
    if data['schema_version'] >= 5 and data.get('tpm', {}).get('multimap') not in {'error', 'drop', 'split'}:
        raise ValueError('database requires a valid TPM ambiguity policy')
    bound['input'] = str(root)
    bound['bundle_root'] = str(root)
    return bound


def input_entries(completed):
    root = Path(completed['input'])
    for entry in completed['files']:
        path = Path(entry['path'])
        if not path.is_relative_to(root): continue
        relative = path.relative_to(root)
        if completed['schema_version'] == 1 or relative == Path('metadata.tsv') or relative.parts[0] in {'cds', 'busco', 'quant'}:
            yield entry, relative


def publish_products(build, data):
    """Publish a relocatable bundle after the ordinary completion checks pass."""
    from build_products import load_complete
    build = Path(build).resolve()
    target = build / ('database' if data['schema_version'] >= 5 else 'products')
    if target.exists():
        existing = load_complete(target)
        if existing['build_id'] != data['build_id']:
            raise ValueError('existing products belong to a different build')
        if any(existing[key] != data[key] for key in ('lineage', 'translation', 'odb', 'excluded_runs')):
            raise ValueError('existing products have different build settings')
        if existing.get('tpm') != data.get('tpm'):
            raise ValueError('existing database has a different TPM policy')
        if set(existing['products']) != set(data['products']):
            raise ValueError('existing products differ from completed build')
        for species, product in data['products'].items():
            if any(product.get(key) != existing['products'][species].get(key) for key in ('row', 'conditions')):
                raise ValueError('existing products have different metadata or conditions')
            for key in ('cds', 'busco', 'abundance', 'protein') + (('expression',) if data['schema_version'] >= 5 else ()):
                if product[key]['sha256'] != existing['products'][species][key]['sha256']:
                    raise ValueError('existing products differ from completed build')
            if product['translation']['sha256'] != existing['products'][species]['source_translation_sha256']:
                raise ValueError('existing translation differs from completed build')
            if data['schema_version'] >= 5 and product['expression_qc']['sha256'] != existing['products'][species]['source_expression_qc_sha256']:
                raise ValueError('existing expression QC differs from completed build')
        from mapping_tables import load_tables
        old_tables = load_tables(existing['mapping']['path'])['tables']
        new_tables = load_tables(data['mapping']['path'])['tables']
        if {s:e['table']['sha256'] for s,e in old_tables.items()} != {s:e['table']['sha256'] for s,e in new_tables.items()}:
            raise ValueError('existing mapping differs from completed build')
        return target / 'manifest.json'
    staging = Path(tempfile.mkdtemp(prefix='.' + target.name + '-', dir=build))
    files = {}

    def add(entry, relative):
        source = verify(entry)
        destination = staging / relative
        link_file(source, destination)
        saved = dict(entry, path=str(relative), stat=stat_identity(destination))
        files[str(relative)] = saved
        return copy.deepcopy(saved)

    def created(relative):
        saved = dict(record(staging / relative), path=str(relative))
        files[str(relative)] = saved
        return copy.deepcopy(saved)

    try:
        for entry, relative in input_entries(data): add(entry, relative)
        for name in ('source_metadata.tsv', 'excluded_accessions.tsv', 'excluded_runs.tsv'):
            source = build / name
            if source.exists(): add(record(source), name)
        if 'excluded_accessions.tsv' not in files:
            write_tsv(staging / 'excluded_accessions.tsv', ['accession', 'reason'], [])
            created('excluded_accessions.tsv')
        if 'excluded_runs.tsv' not in files:
            write_tsv(staging / 'excluded_runs.tsv', ['species', 'run', 'reason'], data.get('excluded_runs', []))
            created('excluded_runs.tsv')
        # These are historical records; their original paths are never dereferenced.
        for name in ('build.json', 'pipeline.yaml', 'checksums.json'):
            source = build / name
            if source.exists(): add(record(source), Path('provenance') / name)
        products = copy.deepcopy(data['products'])
        for species, product in products.items():
            for key in ('cds', 'busco', 'abundance'):
                relative = Path(product[key]['path']).relative_to(data['input'])
                product[key] = copy.deepcopy(files[str(relative)])
            protein_relative = Path('proteins') / Path(product['protein']['path']).name
            product['protein'] = add(product['protein'], protein_relative)
            product['source_translation_sha256'] = product['translation']['sha256']
            translation = json.loads(verify(product['translation']).read_text())
            translation['cds'] = {k:v for k,v in product['cds'].items() if k != 'stat'}
            translation['protein'] = {k:v for k,v in product['protein'].items() if k != 'stat'}
            translation_relative = protein_relative.with_suffix('.json')
            write_json(staging / translation_relative, translation)
            product['translation'] = created(translation_relative)
            if data['schema_version'] >= 5:
                run = product['row']['run']
                relative = Path('expression/runs') / f'{run}.tsv'
                product['expression'] = add(product['expression'], relative)
                product['source_expression_qc_sha256'] = product['expression_qc']['sha256']
                qc = json.loads(verify(product['expression_qc']).read_text())
                qc['abundance'] = {k:v for k,v in product['abundance'].items() if k != 'stat'}
                qc['expression'] = {k:v for k,v in product['expression'].items() if k != 'stat'}
                qc['mapping_table']['path'] = f"odb/species/{product['odb_species']}.tsv.gz"
                qc['mapping_table'].pop('stat', None)
                relative = relative.with_suffix('.qc.json')
                write_json(staging / relative, qc)
                product['expression_qc'] = created(relative)
        from mapping_tables import load_tables, relative_file, write_tables
        source_mapping = verify(data['mapping'])
        snapshot = load_tables(source_mapping)
        entries = {name:dict(entry, table=relative_file(source_mapping.parent, entry['table']))
                   for name,entry in snapshot['tables'].items()}
        mapping_path = write_tables(staging / 'odb', entries, snapshot['version'], snapshot['node'],
                                    snapshot.get('qc'), snapshot.get('reference_sha256s', []))
        published = load_tables(mapping_path, verify_files=False)
        for entry in published['tables'].values():
            relative = 'odb/' + entry['table']['path']
            files[relative] = dict(entry['table'], path=relative)
        mapping = created('odb/snapshot.json')
        portable = {k:copy.deepcopy(data[k]) for k in
                    ('kind','build_id','created_at','fields','translation','lineage','odb','excluded_runs')}
        if data['schema_version'] >= 5: portable['tpm'] = copy.deepcopy(data['tpm'])
        portable.update(schema_version=data["schema_version"], products=products, mapping=mapping, files=list(files.values()))
        portable['sha256'] = digest(portable)
        write_json(staging / 'manifest.json', portable)
        load_complete(staging)
        os.rename(staging, target)
    finally:
        if staging.exists(): shutil.rmtree(staging)
    return target / 'manifest.json'


def publish_pointer(build, manifest):
    build = Path(build).resolve()
    data = {'schema_version':1, 'kind':'completed_build_pointer',
            'manifest':dict(record(manifest), path=str(Path(manifest).relative_to(build)))}
    data['sha256'] = digest(data)
    write_json(build / 'completed.json', data)
    return build / 'completed.json'
