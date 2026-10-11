"""Select matching completed database products and freeze them in a new build."""
import copy
import json
from pathlib import Path

from common import write_json
from dataset_assets import (ProductIntegrityError, counts, product_file,
                            product_json, record, stage_file, stat_identity, verify, verify_many)
from mapping_tables import load_tables, relative_file, write_tables
from portable_build import completion_path
from verification_cache import active_cache

PRODUCTS = ('cds', 'busco', 'abundance', 'protein', 'translation', 'expression', 'expression_qc')
SAMPLE_FIELDS = ('scientific_name', 'taxid', 'species_id', 'analysis_sample_id', 'run',
                 'private_file', 'lib_layout', 'read1_path', 'read2_path')


def software_provenance(entry):
    """Index the software pins without repeatedly parsing historical metadata."""
    path = verify(entry)
    before = stat_identity(path)
    cache = active_cache()
    kind = 'build-software-v1'
    indexed = cache.document(path, kind, before, entry['sha256']) if cache is not None else None
    if indexed is None:
        build = json.loads(path.read_text())
        gg = build.get('config', {}).get('genegalleon', {})
        indexed = {'pins': {k: gg.get(k) for k in ('image_uri', 'version', 'revision')},
                   'software_lock': build.get('software_lock')}
        if before != stat_identity(path):
            raise ValueError(f'file changed while reading build provenance: {path}')
        if cache is not None:
            cache.remember_document(path, kind, before, entry['sha256'], indexed)
    return indexed


def normalize(root, value):
    if value is None:
        return []
    if value == 'auto':
        # Resolve discovery now; preparation records these sources and snapshots
        # products locally, so an existing build never changes its reuse inputs.
        value = [str(p) for p in sorted(Path(root, 'results').glob('*/build.json'))]
    paths = [value] if isinstance(value, str) else value
    if not isinstance(paths, list) or any(not isinstance(p, str) or not p.strip() for p in paths):
        raise ValueError('reuse_from must be null, auto, a build/database path, or a list of paths')
    result = []
    for value in paths:
        path = (Path(root) / value).resolve()
        if not path.is_relative_to(Path(root).resolve()):
            raise ValueError(f'reuse_from must be inside the project for container mounts: {path}')
        if path.name == 'build.json' and path.is_file():
            pass
        elif path.is_dir() and (path / 'build.json').is_file() and not (path / 'completed.json').exists():
            path /= 'build.json'
        else:
            path = completion_path(path)
        if not path.is_relative_to(Path(root).resolve()):
            raise ValueError(f'reuse_from manifest must be inside the project: {path}')
        if str(path) not in result:
            result.append(str(path))
    return result


def validate_conditions(item, product, database, cfg):
    name = item['species']
    if database['translation'] != cfg['translation'] or database['lineage'] != cfg['busco']['lineage']:
        raise ValueError(f'reuse_from genetic code or BUSCO lineage differs: {name}')
    if database['odb'] != {'version': 'v12', 'node': cfg['odb']['ncbi_tax_id']}:
        raise ValueError(f'reuse_from ODB version/node differs: {name}')
    for key in SAMPLE_FIELDS:
        old, new = product['row'].get(key, ''), item['row'].get(key, '')
        if old and new and old != new:
            raise ValueError(f'reuse_from sample metadata differs: {name}: {key}')
    reference = product['cds']['sha256']
    if product.get('reference_id', reference) != reference:
        raise ValueError(f'reuse_from CDS reference identity differs: {name}')
    requested = item['row'].get('reference_id')
    if requested and requested != reference:
        raise ValueError(f'reuse_from reference_id differs: {name}')
    for stage, condition in product.get('conditions', {}).items():
        if condition and condition != cfg['conditions'].get(stage):
            raise ValueError(f'{stage} settings differ from reuse_from product: {name}; use reuse_from: null for a fresh build')


def validate_product(item, product, database, tables, cfg, verify_files=True, mapping_root=None):
    validate_conditions(item, product, database, cfg)
    name, reference = item['species'], product['cds']['sha256']
    for key in PRODUCTS:
        product_file(product[key], verify_files)
    counts(product['counts'])
    translation = product_json(product['translation'], verify_files=False)
    mapping = tables['tables'][name]
    if (translation['cds']['sha256'] != reference or
            translation['protein']['sha256'] != product['protein']['sha256'] or
            translation['translation_table'] != cfg['translation']['table'] or
            mapping['protein_sha256'] != product['protein']['sha256'] or
            mapping['odb_species'] != item['odb_species'] or
            product['odb_species'] != item['odb_species']):
        raise ValueError(f'reuse_from translation/mapping identity differs: {name}')
    root = mapping_root if mapping_root is not None else Path(database['mapping']['path']).parent.resolve()
    table = relative_file(root, mapping['table'], resolved_root=root)
    product_file(table, verify_files)
    qc = product_json(product['expression_qc'], verify_files=False)
    if (qc.get('species') != name or qc.get('run') != item['row']['run'] or
            qc.get('multimap') != 'error' or
            qc.get('abundance', {}).get('sha256') != product['abundance']['sha256'] or
            qc.get('mapping_table', {}).get('sha256') != table['sha256'] or
            qc.get('expression', {}).get('sha256') != product['expression']['sha256']):
        raise ValueError(f'reuse_from expression identity differs: {name}')
    return dict(mapping, table=table), translation


def select(paths, items, cfg, verify_files=True, rebuilds=None):
    """Read only selected samples; duplicate sources must agree on scientific products."""
    from build_products import load_complete
    if not items:
        return {}, {}, []
    wanted = {item['species']: item for item in items}
    selected, errors, sources = {}, {}, []
    stage_paths = [path for path in paths if Path(path).name == 'build.json']
    for path in paths:
        if path in stage_paths:
            continue
        source = record(path)
        data = load_complete(path, verify_files=False)
        verify(source)
        sources.append(source)
        names = sorted(wanted.keys() & data['products'].keys())
        if not names:
            continue
        software_lock = None
        provenance = Path(data['bundle_root']) / 'provenance/build.json'
        inventory = {entry['path']: entry for entry in data['files']}
        if str(provenance) in inventory:
            build = software_provenance(inventory[str(provenance)])
            if all(value == cfg['genegalleon'].get(k) for k, value in build['pins'].items()):
                software_lock = build['software_lock']
        tables = load_tables(verify(data['mapping']), verify_files=False)
        mapping_root = Path(data['mapping']['path']).parent.resolve()
        if (tables['version'], tables['node']) != ('v12', cfg['odb']['ncbi_tax_id']):
            raise ValueError(f'reuse_from mapping version/node differs: {path}')
        candidates = []
        for name in names:
            if name in errors:
                continue
            try:
                product = data['products'][name]
                validate_conditions(wanted[name], product, data, cfg)
                mapping = tables['tables'][name]
                mapping = dict(mapping, table=relative_file(mapping_root, mapping['table'], resolved_root=mapping_root))
                candidates.append((name, product, mapping))
            except ProductIntegrityError as error:
                if rebuilds is None or wanted[name]['row'].get('reference_id'):
                    errors[name] = str(error)
                else:
                    rebuilds[name] = str(error)
                selected.pop(name, None)
            except (ValueError, OSError, KeyError) as error:
                errors[name] = str(error)
                selected.pop(name, None)
        failures = {}
        if verify_files:
            entries = [(name, entry) for name, product, mapping in candidates
                       for entry in [*(product[key] for key in PRODUCTS), mapping['table']]]
            for (name, _), error in zip(entries, verify_many(entry for _, entry in entries)):
                if error is not None:
                    failures.setdefault(name, error)
        for name, product, mapping in candidates:
            try:
                if name in failures:
                    raise failures[name]
                mapping, translation = validate_product(wanted[name], product, data, tables, cfg, False,
                                                        mapping_root=mapping_root)
                identity = {'files': {k: product[k]['sha256'] for k in ('cds', 'busco', 'abundance', 'protein', 'expression')},
                            'mapping': mapping['table']['sha256'], 'counts': product['counts'],
                            'conditions': {stage: product.get('conditions', {}).get(stage) for stage in ('assembly', 'busco', 'quant')},
                            'reference': sorted(tables.get('reference_sha256s', []))}
                previous = selected.get(name)
                if previous and previous['identity'] != identity:
                    raise ValueError(f'conflicting reuse_from databases for {name}: {previous["source"]["path"]} and {path}')
                if not previous:
                    selected[name] = {'product': product, 'mapping': mapping, 'translation': translation,
                                      'identity': identity, 'source': source, 'software_lock': software_lock}
            except ProductIntegrityError as error:
                if rebuilds is None or wanted[name]['row'].get('reference_id'):
                    errors[name] = str(error)
                else:
                    rebuilds[name] = str(error)
                selected.pop(name, None)
            except (ValueError, OSError, KeyError) as error:
                errors[name] = str(error)
                selected.pop(name, None)
    references = {r for entry in selected.values() for r in entry['identity']['reference']}
    if len(references) > 1:
        raise ValueError('reuse_from databases use different ODB reference snapshots')
    from build_stage_reuse import select as select_stages
    select_stages(stage_paths, items, cfg, selected, errors, sources,
                  verify_files=verify_files, rebuilds=rebuilds)
    if rebuilds is not None:
        for name in selected:
            rebuilds.pop(name, None)
    return selected, errors, sources


def stage(staging, target, selected, cfg):
    """Snapshot validated products locally so retries do not depend on source databases."""
    from aggregate_tpm import cache_identity
    base = Path(staging) / 'work/cache'
    final = Path(target) / 'work/cache'
    tables = {}
    references = set()
    files = {}

    def saved(entry, relative):
        path = base / relative
        saved = stage_file(entry, path, final / relative)
        files[name].append(saved)
        return saved

    for name, chosen in selected.items():
        files[name] = []
        if chosen.get('kind') == 'stages':
            from build_stage_reuse import stage as stage_products
            files[name] = stage_products(staging, target, name, chosen)
            continue
        product = chosen['product']
        row = product['row']
        refid = product['cds']['sha256']
        prefix = Path('products') / name / refid
        cds = saved(product['cds'], prefix / 'cds.fa.gz')
        full = saved(product['busco'], prefix / ('busco.tsv.gz' if product['busco']['path'].endswith('.gz') else 'busco.tsv'))
        abundance = saved(product['abundance'], prefix / 'abundance.tsv')
        def provenance(stage):
            return {'source': 'database', 'database': chosen['source'],
                    'condition': product.get('conditions', {}).get(stage),
                    'raw_inputs': product.get('raw_inputs', {})}
        write_json(base / prefix / 'reference.json', {'schema_version': 2, 'species': name,
                   'taxid': row['taxid'], 'species_id': row['species_id'], 'run': row['run'],
                   'reference_id': refid, 'cds': cds, 'provenance': provenance('assembly')})
        write_json(base / prefix / 'busco.json', {'schema_version': 1, 'reference_id': refid,
                   'lineage': cfg['busco']['lineage'], 'counts': product['counts'], 'full': full,
                   'short': None, 'provenance': provenance('busco')})
        write_json(base / prefix / 'quant' / (row['run'] + '.json'), {'schema_version': 1,
                   'reference_id': refid, 'run': row['run'], 'abundance': abundance,
                   'sample': {k: row.get(k, '') for k in SAMPLE_FIELDS[5:]}, 'provenance': provenance('quant')})
        protein_dir = Path('proteins') / refid / f'table_{cfg["translation"]["table"]}'
        protein = saved(product['protein'], protein_dir / 'protein.fa')
        info = chosen['translation']
        write_json(base / protein_dir / 'receipt.json', {'schema_version': 1, 'cds_sha256': refid,
                   'translation_table': cfg['translation']['table'], 'seqkit': info['seqkit'],
                   'sequences': info['sequences'], 'protein': dict(protein, path='protein.fa')})
        tables[name] = chosen['mapping']
        references.update(chosen['identity']['reference'])
        identity = cache_identity({'species': name, 'run': row['run']}, abundance, chosen['mapping']['table'])
        expression_dir = Path('expression') / identity
        expression = saved(product['expression'], expression_dir / 'tpm.tsv')
        qc = copy.deepcopy(json.loads(verify(product['expression_qc']).read_text()))
        qc.update(abundance=abundance, expression=expression)
        qc['mapping_table'] = dict(chosen['mapping']['table'],
            path=str(final / 'odb' / f'v12_{cfg["odb"]["ncbi_tax_id"]}' / 'reused/species' / (product['odb_species'] + '.tsv.gz')))
        write_json(base / expression_dir / 'qc.json', qc)
        qc_entry = record(base / expression_dir / 'qc.json')
        write_json(base / expression_dir / 'receipt.json', {'identity': identity,
            'files': {'tpm.tsv': dict(expression, path='tpm.tsv'), 'qc.json': dict(qc_entry, path='qc.json')},
            'reused_from': chosen['source']})
        for relative in (prefix / 'reference.json', prefix / 'busco.json', prefix / 'quant' / (row['run'] + '.json'),
                         protein_dir / 'receipt.json', expression_dir / 'qc.json', expression_dir / 'receipt.json'):
            files[name].append(dict(record(base / relative), path=str(final / relative)))
    if tables:
        snapshot = write_tables(base / 'odb' / f'v12_{cfg["odb"]["ncbi_tax_id"]}' / 'reused', tables,
                                node=cfg['odb']['ncbi_tax_id'], references=references)
        saved_tables = load_tables(snapshot, verify_files=False)
        destination = final / 'odb' / f'v12_{cfg["odb"]["ncbi_tax_id"]}' / 'reused'
        for name, table in saved_tables['tables'].items():
            files[name].append(relative_file(destination, table['table']))
    return files
