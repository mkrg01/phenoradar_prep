"""Compact successful native diagnostics after the completed database verifies."""
import json
from pathlib import Path

from build_products import load_complete
from dataset import item_products, load, workspace
from native_archive import compact
from sample_review import jobs, sample_guard
from storage_meter import measure


def run(path, *, apply=False):
    path = Path(path).resolve()
    manifest = load(path, check_code=True)
    database = load_complete(path)
    before = measure([path])
    if before['errors']:
        raise ValueError('cannot compact an incompletely inventoried build')
    activity = jobs(path, manifest)
    if any(entry['activity'] for entry in activity.values()):
        raise ValueError('sample jobs are active or unresolved; compaction deferred')
    controller = path / 'jobs/database.json'
    if controller.is_file() and json.loads(controller.read_text()).get('state') != 'complete':
        raise ValueError('database controller has not recorded completion')
    report = {'apply': apply, 'before': before, 'samples': []}
    for item in manifest['items']:
        work = workspace(path, manifest, item)
        if not work.exists():
            continue
        row = {'sample': item['species'], 'state': 'retained'}
        with sample_guard(path, manifest, item, native=True):
            receipt = path / 'jobs/status' / f"{item['species']}.sample.json"
            job = json.loads(receipt.read_text()) if receipt.is_file() else {}
            # Observations omit wrapper failures such as the mandatory versions
            # command. Legacy receipts without a wrapper exit must be retained.
            if (manifest['config'].get('storage', {}).get('keep_intermediates', False)
                    or job.get('state') not in {'complete', 'reused'} or job.get('native_exit_code') != 0):
                row['reason'] = 'debug retention, incomplete native evidence, or failed native attempt'
                report['samples'].append(row)
                continue
            products = item_products(manifest, item)
            if not all(products.get(key) for key in ('reference', 'busco', 'quant')):
                raise ValueError('registered products missing; diagnostics retained')
            protected = [products['reference']['cds']['path'], products['quant']['abundance']['path']]
            protected.extend(products['busco'][key]['path'] for key in ('full', 'short')
                             if products['busco'].get(key))
            protected.extend(entry['path'] for entry in manifest.get('raw_inputs', {}).get(item['species'], {}).values())
            row.update(compact(work, protected, apply=apply))
        report['samples'].append(row)
    from derived_cache_cleanup import retire
    if manifest['config'].get('storage', {}).get('keep_intermediates', False):
        report['derived_cache'] = {'state': 'retained', 'reason': 'keep_intermediates', 'apply': apply}
    else:
        report['derived_cache'] = retire(path, [entry['path'] for entry in database.get('files', [])], apply=apply)
    load_complete(path)
    report['after'] = measure([path])
    report['database_verified'] = True
    return report
