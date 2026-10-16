"""Reuse completed samples across software updates; resume partial work only unchanged."""
import copy
from pathlib import Path

from common import write_json, content_sha256
from product_storage import stage_product
from dataset_assets import ProductIntegrityError, record, resolve, stage_file, verify

STAGES = {'assembly': 'reference', 'busco': 'busco', 'quant': 'quant'}
IDENTITY_FIELDS = ('scientific_name', 'taxid', 'species_id', 'analysis_sample_id', 'run',
                   'private_file', 'lib_layout', 'read1_path', 'read2_path')


def fingerprints(entry):
    if entry.get('kind') != 'stages':
        product = entry['product']
        return {key: content_sha256(product[key]) for key in ('cds', 'busco', 'abundance')}
    products = entry['products']
    result = {'cds': products['reference']['cds']['sha256']}
    if products['busco']:
        result['busco'] = products['busco']['full']['sha256']
    if products['quant']:
        result['abundance'] = content_sha256(products['quant']['abundance'])
    return result


def merge(selected, name, entry):
    previous = selected.get(name)
    if previous is None:
        selected[name] = entry
        return
    old_complete = previous.get('kind') != 'stages' or all(previous['products'][key] for key in STAGES.values())
    new_complete = entry.get('kind') != 'stages' or all(entry['products'][key] for key in STAGES.values())
    old_processing = (previous.get('processing') if previous.get('kind') == 'stages'
                      else previous['product'].get('processing'))
    new_processing = (entry.get('processing') if entry.get('kind') == 'stages'
                      else entry['product'].get('processing'))
    same_producer = old_processing is not None and old_processing == new_processing
    if old_complete != new_complete and not same_producer:
        if new_complete:
            selected[name] = entry
        return
    old, new = fingerprints(previous), fingerprints(entry)
    if any(old[key] != new[key] for key in old.keys() & new.keys()):
        raise ValueError(f'conflicting reuse_from stages for {name}; choose a specific source or reuse_from: null')
    if previous.get('kind') == 'stages':
        if entry.get('kind') != 'stages':
            selected[name] = entry
        else:
            # An older complete sample must replace a newer partial sample as
            # a unit; never manufacture one sample from different releases.
            if old_complete:
                return
            for key in ('busco', 'quant'):
                if previous['products'][key] is None:
                    previous['products'][key] = entry['products'][key]


def select(paths, items, cfg, selected, errors, sources, verify_files=True, rebuilds=None):
    from dataset import load, check_conditions
    from sample_processing import scientific_conditions, native_processing
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
            if scientific_conditions(build['config']) != scientific_conditions(cfg):
                continue
            try:
                products = resolve(build['config']['store'], old,
                                   build['config']['busco']['lineage'], need_full=True,
                                   verify_files=verify_files)
                if products['reference'] is None:
                    continue
                check_conditions(products, build['config'].get('conditions'))
                complete = all(products[key] for key in STAGES.values())
                if not complete and build['config'].get('conditions') != cfg['conditions']:
                    continue
                requested = item['row'].get('reference_id')
                if requested and requested != products['reference']['reference_id']:
                    raise ValueError(f'reuse_from reference_id differs: {name}')
                processing = native_processing(products, build)
                if processing and processing.get('compatibility') != scientific_conditions(cfg):
                    continue
                for product in products.values():
                    if product:
                        product.setdefault('provenance', {})['reused_from'] = source
                merge(selected, name, {'kind': 'stages', 'products': products,
                                       'row': old['row'], 'source': source,
                                       'processing': processing,
                                       'software_lock': build.get('software_lock')})
            except ProductIntegrityError as error:
                if rebuilds is None or item['row'].get('reference_id'):
                    errors[name] = str(error)
                else:
                    rebuilds[name] = str(error)
                selected.pop(name, None)
            except (ValueError, OSError, KeyError) as error:
                errors[name] = str(error)
                selected.pop(name, None)
        verify(source)


def stage(staging, target, name, chosen, cfg):
    """Copy receipts and hardlink verified outputs into this build's own cache."""
    products = chosen['products']
    reference = products['reference']
    prefix = Path('work/cache/products') / name / reference['reference_id']
    files = []

    def saved(entry, filename):
        path = Path(staging) / prefix / filename
        save = stage_product if filename == 'abundance.tsv.gz' else stage_file
        saved = save(entry, path, Path(target) / prefix / filename)
        files.append(saved)
        return saved

    for key, filename in (('reference', 'reference.json'), ('busco', 'busco.json'),
                          ('quant', 'quant/' + chosen['row']['run'] + '.json')):
        product = products[key]
        if product is None:
            continue
        receipt = copy.deepcopy(product)
        receipt.setdefault('provenance', {}).setdefault('reused_from', chosen['source'])
        if chosen.get('processing'):
            receipt['provenance']['processing'] = copy.deepcopy(chosen['processing'])
        receipt['provenance']['accepted_condition'] = cfg['conditions'][
            next(stage for stage, field in STAGES.items() if field == key)]
        if key == 'reference':
            receipt['cds'] = saved(receipt['cds'], 'cds.fa.gz')
        elif key == 'busco':
            for field in ('full', 'short'):
                if receipt.get(field):
                    suffix = '.gz' if receipt[field]['path'].endswith('.gz') else ''
                    receipt[field] = saved(receipt[field], 'busco.' + field + suffix)
        else:
            receipt['abundance'] = saved(receipt['abundance'], 'abundance.tsv.gz')
        write_json(Path(staging) / prefix / filename, receipt)
        files.append(dict(record(Path(staging) / prefix / filename),
                          path=str(Path(target) / prefix / filename)))
    return files
