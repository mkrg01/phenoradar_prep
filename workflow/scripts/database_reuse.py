"""Select matching completed database products and freeze them in a new build."""
import copy
import json
from pathlib import Path

from common import write_json
from dataset_assets import counts, link_file, record, stat_identity, verify
from mapping_tables import load_tables, relative_file, write_tables
from portable_build import completion_path

PRODUCTS = ('cds', 'busco', 'abundance', 'protein', 'translation')
SAMPLE_FIELDS = ('scientific_name', 'taxid', 'species_id', 'analysis_sample_id', 'run',
                 'private_file', 'lib_layout', 'read1_path', 'read2_path')


def normalize(root, value):
    if value is None:
        return []
    paths = [value] if isinstance(value, str) else value
    if not isinstance(paths, list) or any(not isinstance(p, str) or not p.strip() for p in paths):
        raise ValueError('reuse_from must be null, a database path, or a list of database paths')
    result = []
    for value in paths:
        path = (Path(root) / value).resolve()
        if not path.is_relative_to(Path(root).resolve()):
            raise ValueError(f'reuse_from must be inside the project for container mounts: {path}')
        path = completion_path(path)
        if not path.is_relative_to(Path(root).resolve()):
            raise ValueError(f'reuse_from manifest must be inside the project: {path}')
        if str(path) not in result:
            result.append(str(path))
    return result


def validate_product(item, product, database, tables, cfg):
    name = item['species']
    if database['schema_version'] < 4:
        raise ValueError('reuse_from requires sample identities (schema-4 or newer database)')
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
    for entry in product.get('raw_inputs', {}).values():
        if Path(entry['path']).exists(): verify(entry)
    for key in PRODUCTS:
        verify(product[key])
    counts(product['counts'])
    translation = json.loads(Path(product['translation']['path']).read_text())
    mapping = tables['tables'][name]
    if (translation['cds']['sha256'] != reference or
            translation['protein']['sha256'] != product['protein']['sha256'] or
            translation['translation_table'] != cfg['translation']['table'] or
            mapping['protein_sha256'] != product['protein']['sha256'] or
            mapping['odb_species'] != item['odb_species'] or
            product['odb_species'] != item['odb_species']):
        raise ValueError(f'reuse_from translation/mapping identity differs: {name}')
    table = relative_file(Path(database['mapping']['path']).parent, mapping['table'])
    verify(table)
    if 'expression' in product:
        verify(product['expression']); verify(product['expression_qc'])
        qc = json.loads(Path(product['expression_qc']['path']).read_text())
        if (qc.get('species') != name or qc.get('run') != item['row']['run'] or
                qc.get('multimap') != 'error' or
                qc.get('abundance', {}).get('sha256') != product['abundance']['sha256'] or
                qc.get('mapping_table', {}).get('sha256') != table['sha256'] or
                qc.get('expression', {}).get('sha256') != product['expression']['sha256']):
            raise ValueError(f'reuse_from expression identity differs: {name}')
    return dict(mapping, table=table), translation


def select(paths, items, cfg):
    """Read only selected samples; duplicate sources must agree on scientific products."""
    from build_products import load_complete
    wanted = {item['species']: item for item in items}
    selected, errors, sources = {}, {}, []
    for path in paths:
        source = record(path)
        data = load_complete(path, verify_files=False)
        verify(source)
        sources.append(source)
        names = sorted(wanted.keys() & data['products'].keys())
        if not names:
            continue
        tables = load_tables(verify(data['mapping']), verify_files=False)
        if (tables['version'], tables['node']) != ('v12', cfg['odb']['ncbi_tax_id']):
            raise ValueError(f'reuse_from mapping version/node differs: {path}')
        for name in names:
            if name in errors:
                continue
            try:
                product = data['products'][name]
                mapping, translation = validate_product(wanted[name], product, data, tables, cfg)
                identity = {'files': {k: product[k]['sha256'] for k in ('cds', 'busco', 'abundance', 'protein')},
                            'mapping': mapping['table']['sha256'], 'counts': product['counts'],
                            'conditions': {stage: product.get('conditions', {}).get(stage) for stage in ('assembly', 'busco', 'quant')},
                            'reference': sorted(tables.get('reference_sha256s', []))}
                previous = selected.get(name)
                expression = product.get('expression', {}).get('sha256')
                old_expression = previous['product'].get('expression', {}).get('sha256') if previous else None
                if previous and (previous['identity'] != identity or
                                 (expression and old_expression and expression != old_expression)):
                    raise ValueError(f'conflicting reuse_from databases for {name}: {previous["source"]["path"]} and {path}')
                if not previous or (expression and not old_expression):
                    selected[name] = {'product': product, 'mapping': mapping, 'translation': translation,
                                      'identity': identity, 'source': source}
            except (ValueError, OSError, KeyError) as error:
                errors[name] = str(error)
                selected.pop(name, None)
    references = {r for entry in selected.values() for r in entry['identity']['reference']}
    if len(references) > 1:
        raise ValueError('reuse_from databases use different ODB reference snapshots')
    return selected, errors, sources


def stage(staging, target, selected, cfg):
    """Snapshot validated products locally so retries do not depend on source databases."""
    from aggregate_tpm import cache_identity
    base = Path(staging) / 'work/cache'
    final = Path(target) / 'work/cache'
    tables = {}
    references = set()

    def saved(entry, relative):
        path = base / relative
        link_file(verify(entry), path)
        return dict(entry, path=str(final / relative), stat=stat_identity(path))

    for name, chosen in selected.items():
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
        if 'expression' in product:
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
    if tables:
        write_tables(base / 'odb' / f'v12_{cfg["odb"]["ncbi_tax_id"]}' / 'reused', tables,
                     node=cfg['odb']['ncbi_tax_id'], references=references)
