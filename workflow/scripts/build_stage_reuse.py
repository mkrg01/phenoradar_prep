"""Reuse verified sample stages without requiring a completed database."""
import copy
from pathlib import Path

from common import write_json
from dataset_assets import link_file, record, resolve, stat_identity, verify

STAGES = {'assembly': 'reference', 'busco': 'busco', 'quant': 'quant'}
IDENTITY_FIELDS = ('scientific_name', 'taxid', 'species_id', 'analysis_sample_id', 'run',
                   'private_file', 'lib_layout', 'read1_path', 'read2_path')


def fingerprints(entry):
    if entry.get('kind') != 'stages':
        product = entry['product']
        return {key: product[key]['sha256'] for key in ('cds', 'busco', 'abundance')}
    products = entry['products']
    result = {'cds': products['reference']['cds']['sha256']}
    if products['busco']:
        result['busco'] = products['busco']['full']['sha256']
    if products['quant']:
        result['abundance'] = products['quant']['abundance']['sha256']
    return result


def merge(selected, name, entry):
    previous = selected.get(name)
    if previous is None:
        selected[name] = entry
        return
    old, new = fingerprints(previous), fingerprints(entry)
    if any(old[key] != new[key] for key in old.keys() & new.keys()):
        raise ValueError(f'conflicting reuse_from stages for {name}; choose a specific source or reuse_from: null')
    if previous.get('kind') == 'stages':
        if entry.get('kind') != 'stages':
            selected[name] = entry
        else:
            for key in ('busco', 'quant'):
                if previous['products'][key] is None:
                    previous['products'][key] = entry['products'][key]


def select(paths, items, cfg, selected, errors, sources):
    from dataset import load
    wanted = {item['species']: item for item in items}
    for path in paths:
        source = record(path)
        build = load(Path(path).parent)
        old_items = {item['species']: item for item in build['items']}
        sources.append(source)
        for name in sorted(wanted.keys() & old_items.keys()):
            if name in errors:
                continue
            item, old = wanted[name], old_items[name]
            # Different inputs/settings are ordinary cache misses. Damaged or
            # conflicting matching products are reported rather than reused.
            if any(old['row'].get(k, '') != item['row'].get(k, '') for k in IDENTITY_FIELDS):
                continue
            if build['config'].get('conditions', {}).get('assembly') != cfg['conditions']['assembly']:
                continue
            try:
                products = resolve(build['config']['store'], old,
                                   build['config']['busco']['lineage'], need_full=True)
                if products['reference'] is None:
                    continue
                requested = item['row'].get('reference_id')
                if requested and requested != products['reference']['reference_id']:
                    raise ValueError(f'reuse_from reference_id differs: {name}')
                for stage, key in STAGES.items():
                    product = products[key]
                    if product and product.get('provenance', {}).get('condition') != cfg['conditions'][stage]:
                        products[key] = None
                if products['reference'] is None:
                    continue
                for product in (products['reference'], products['quant']):
                    for raw in (product or {}).get('provenance', {}).get('raw_inputs', {}).values():
                        if Path(raw['path']).exists():
                            verify(raw)
                for product in products.values():
                    if product:
                        product.setdefault('provenance', {})['reused_from'] = source
                merge(selected, name, {'kind': 'stages', 'products': products,
                                       'row': old['row'], 'source': source,
                                       'software_lock': build.get('software_lock')})
            except (ValueError, OSError, KeyError) as error:
                errors[name] = str(error)
                selected.pop(name, None)
        verify(source)


def stage(staging, target, name, chosen):
    """Copy receipts and hardlink verified outputs into this build's own cache."""
    products = chosen['products']
    reference = products['reference']
    prefix = Path('work/cache/products') / name / reference['reference_id']

    def saved(entry, filename):
        path = Path(staging) / prefix / filename
        link_file(verify(entry), path)
        return dict(entry, path=str(Path(target) / prefix / filename), stat=stat_identity(path))

    for key, filename in (('reference', 'reference.json'), ('busco', 'busco.json'),
                          ('quant', 'quant/' + chosen['row']['run'] + '.json')):
        product = products[key]
        if product is None:
            continue
        receipt = copy.deepcopy(product)
        receipt.setdefault('provenance', {}).setdefault('reused_from', chosen['source'])
        if key == 'reference':
            receipt['cds'] = saved(receipt['cds'], 'cds.fa.gz')
        elif key == 'busco':
            for field in ('full', 'short'):
                if receipt.get(field):
                    suffix = '.gz' if receipt[field]['path'].endswith('.gz') else ''
                    receipt[field] = saved(receipt[field], 'busco.' + field + suffix)
        else:
            receipt['abundance'] = saved(receipt['abundance'], 'abundance.tsv')
        write_json(Path(staging) / prefix / filename, receipt)
