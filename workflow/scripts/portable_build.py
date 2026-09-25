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
        else:
            path /= 'products/manifest.json'
    if not path.is_file():
        raise ValueError(f'build is incomplete: missing {path}; finish build through mapping first')
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
    bound['odb_snapshot'] = bind(bound['odb_snapshot'])
    for species, product in bound['products'].items():
        if product['row']['scientific_name'].replace(' ', '_') != species:
            raise ValueError('product species differs from metadata')
        for key in PRODUCT_FILES: product[key] = bind(product[key])
    for required in ('metadata.tsv', 'busco/summary.tsv', 'excluded_accessions.tsv', 'excluded_runs.tsv'):
        if required not in files: raise ValueError(f'incomplete product bundle: {required}')
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
    target = build / 'products'
    if target.exists():
        existing = load_complete(target)
        if existing['build_id'] != data['build_id']:
            raise ValueError('existing products belong to a different build')
        if any(existing[key] != data[key] for key in ('lineage', 'translation', 'odb', 'excluded_runs')):
            raise ValueError('existing products have different build settings')
        if set(existing['products']) != set(data['products']):
            raise ValueError('existing products differ from completed build')
        for species, product in data['products'].items():
            if any(product.get(key) != existing['products'][species].get(key) for key in ('row', 'conditions')):
                raise ValueError('existing products have different metadata or conditions')
            for key in PRODUCT_FILES[:-1]:
                if product[key]['sha256'] != existing['products'][species][key]['sha256']:
                    raise ValueError('existing products differ from completed build')
            if product['translation']['sha256'] != existing['products'][species]['source_translation_sha256']:
                raise ValueError('existing translation differs from completed build')
        if existing['mapping']['sha256'] != data['mapping']['sha256']:
            raise ValueError('existing mapping differs from completed build')
        return target / 'manifest.json'
    staging = Path(tempfile.mkdtemp(prefix='.products-', dir=build))
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
        mapping = add(data['mapping'], 'odb/mappings.sqlite')
        mapping_root = Path(data['mapping']['path']).parent
        annotations = add(record(mapping_root / 'gene_orthogroups.tsv'), 'odb/annotations.tsv')
        if (mapping_root / 'merge_qc.json').exists():
            qc_path = mapping_root / 'merge_qc.json'
            add(record(qc_path), 'provenance/mapping_qc.json')
            for index, entry in enumerate(json.loads(qc_path.read_text()).get('sources', [])):
                add(entry, f'provenance/mapping_sources/{index:04d}.json')
        snapshot = {'schema_version':1, **data['odb'],
                    'proteins':[{'species':s, 'odb_species':p['odb_species'],
                                 **{k:v for k,v in p['protein'].items() if k != 'stat'}}
                                for s,p in sorted(products.items())],
                    'annotations':{k:('annotations.tsv' if k == 'path' else v) for k,v in annotations.items() if k != 'stat'},
                    'source_build':data['build_id']}
        if data.get('odb_reference_sha256s'):
            snapshot['reference_sha256s'] = data['odb_reference_sha256s']
        write_json(staging / 'odb/snapshot.json', snapshot)
        odb_snapshot = created('odb/snapshot.json')
        portable = {k:copy.deepcopy(data[k]) for k in
                    ('kind','build_id','created_at','fields','translation','lineage','odb','excluded_runs')}
        portable.update(schema_version=2, products=products, mapping=mapping,
                        odb_snapshot=odb_snapshot, files=list(files.values()))
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


def register_products(source, store, cache_dir, lineage, translation, node, excluded_runs=()):
    """Use a copied completed build as the seed for another project's builds."""
    from build_products import load_complete
    from dataset_assets import identities, register_reference, register_busco, register_quant
    from incremental_odb import import_snapshot
    source = completion_path(source)
    completed = load_complete(source)
    if completed['schema_version'] != 2:
        raise ValueError('register --products requires a portable products directory')
    if completed['lineage'] != lineage or completed['translation'] != translation or completed['odb'] != {'version':'v12','node':node}:
        raise ValueError('portable build lineage/translation/ODB settings differ from build.yaml')
    _, items = identities(Path(completed['input']) / 'metadata.tsv')
    if {i['species'] for i in items} != set(completed['products']):
        raise ValueError('portable metadata differs from product membership')
    results = []
    source_hash = completed['sha256']
    for item in items:
        if item['row']['run'] in excluded_runs:
            results.append({'species':item['species'], 'status':'excluded'}); continue
        product = completed['products'][item['species']]
        if item['row'] != product['row']:
            raise ValueError('portable run metadata differs from product')
        def provenance(stage):
            p = {'source':'portable_build', 'build_sha256':source_hash}
            condition = product.get('conditions', {}).get(stage)
            if condition: p['condition'] = condition
            return p
        ref = register_reference(store, item, product['cds']['path'], provenance('assembly'))
        register_busco(store, ref, full=product['busco']['path'], lineage=lineage, provenance=provenance('busco'))
        register_quant(store, ref, item, product['abundance']['path'], provenance('quant'))
        results.append({'species':item['species'], 'status':'registered'})
    odb = import_snapshot(Path(completed['odb_snapshot']['path']).parent, cache_dir, node=node)
    return results, odb
