#!/usr/bin/env python3
"""Import a schema-4 database and proven legacy OG expression without analysis jobs.

Source bundles and downstream results stay intact. Large files are hard-linked;
only expression sample labels and portable provenance records are rewritten.
"""
import argparse
import copy
import json
import math
import os
import shutil
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from common import now, write_json
from dataset_assets import digest, record, stat_identity, verify
from mapping_tables import load_tables
from portable_build import completion_path, load_products, publish_pointer


def expression_evidence(species, product, table, run_dir, receipts, sources, multimap):
    run = product['row']['run']
    path = run_dir / 'orthogroups/expression/runs' / f'{run}.tsv'
    qc_path = path.with_suffix('.qc.json')
    qc = json.loads(qc_path.read_text())
    receipt_path = receipts / f'{species}.json'
    receipt = json.loads(receipt_path.read_text())
    old = receipt['old_id']
    if (receipt['sample_id'] != species or receipt['run'] != run or
            qc['species'] != old or qc['run'] != run or qc['multimap'] != multimap):
        raise ValueError(f'expression identity or policy differs: {species}')
    conversions = receipt['conversions']
    for key, target in [('quant', product['abundance']), ('protein', product['protein']),
                        ('mapping', table['table'])]:
        change = conversions[key]
        if (change['old_id'] != old or change['new_id'] != species or
                change['output']['sha256'] != target['sha256']):
            raise ValueError(f'{key} conversion differs from database: {species}')
    if conversions['quant']['source']['sha256'] != qc['abundance']['sha256']:
        raise ValueError(f'expression abundance differs from converted input: {species}')
    annotation = table.get('annotation_sha256')
    source = sources.get(annotation)
    if source is None or source['proteins'].get(old) != conversions['protein']['source']['sha256']:
        raise ValueError(f'expression mapping source differs from database: {species}')
    if receipt['mapping']['annotation_sha256'] != annotation:
        raise ValueError(f'conversion mapping annotation differs: {species}')
    if qc.get('mapping_table', {}).get('sha256', conversions['mapping']['source']['sha256']) != conversions['mapping']['source']['sha256']:
        raise ValueError(f'expression mapping table differs: {species}')
    if qc['protein_genes'] != table['qc']['protein_genes']:
        raise ValueError(f'expression protein count differs: {species}')
    if 'expression' in qc:
        verify(dict(qc['expression'], path=str(path)))
    # Record legacy inputs now, then verify them again during conversion.
    return {'species': species, 'run': run, 'old': old, 'expression': record(path),
            'qc': record(qc_path), 'receipt': record(receipt_path)}


def rewrite_expression(task):
    evidence, database, abundance, table = task
    database = Path(database)
    source = verify(evidence['expression'])
    qc = json.loads(verify(evidence['qc']).read_text())
    species, run, old = (evidence[k] for k in ('species', 'run', 'old'))
    relative = Path('expression/runs') / f'{run}.tsv'
    output = database / relative
    output.parent.mkdir(parents=True, exist_ok=True)
    groups, original, normalized = set(), [], []
    with source.open() as src, output.open('x') as dest:
        if src.readline().rstrip('\r\n') != 'species\trun\torthogroup\ttpm_sum\ttpm':
            raise ValueError(f'invalid expression columns: {source}')
        dest.write('species\trun\torthogroup\ttpm_sum\ttpm\n')
        for line in src:
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 5 or fields[:2] != [old, run] or not fields[2] or fields[2] in groups:
                raise ValueError(f'invalid expression identity or duplicate OG: {source}')
            values = [float(v) for v in fields[3:]]
            if any(not math.isfinite(v) or v < 0 for v in values):
                raise ValueError(f'invalid expression values: {source}')
            if not math.isclose(values[1], values[0] / qc['retained_tpm'] * 1e6, rel_tol=1e-10, abs_tol=1e-10):
                raise ValueError(f'expression normalization differs from QC: {source}')
            groups.add(fields[2]); original.append(values[0]); normalized.append(values[1])
            # Keep numeric text exactly as originally written.
            dest.write(species + '\t' + '\t'.join(fields[1:]) + '\n')
    if (len(groups) != qc['orthogroups'] or not groups or
            not math.isclose(math.fsum(original), qc['retained_tpm'], rel_tol=1e-9) or
            not math.isclose(math.fsum(normalized), 1e6, rel_tol=1e-9)):
        raise ValueError(f'expression totals differ from QC: {source}')
    verify(evidence['expression']); verify(evidence['qc']); verify(evidence['receipt'])
    expression = dict(record(output), path=str(relative))
    qc.update(species=species, abundance={k: v for k, v in abundance.items() if k != 'stat'},
              mapping_table={k: v for k, v in table.items() if k != 'stat'},
              expression={k: v for k, v in expression.items() if k != 'stat'})
    qc['migration'] = {'kind': 'sample_id_only', 'source_expression': evidence['expression'],
                       'source_qc': evidence['qc'], 'conversion_receipt': evidence['receipt']}
    qc_relative = relative.with_suffix('.qc.json')
    write_json(database / qc_relative, qc)
    return species, expression, dict(record(database / qc_relative), path=str(qc_relative))


def migrate(source, destination, legacy_run, receipts, workers=1):
    destination, legacy_run, receipts = map(lambda p: Path(p).resolve(), (destination, legacy_run, receipts))
    if destination.exists():
        raise ValueError(f'destination already exists: {destination}; choose a new path')
    if type(workers) is not int or workers < 1:
        raise ValueError('workers must be a positive integer')
    source = completion_path(source)
    source_record = record(source)
    data = json.loads(source.read_text())
    if data.get('schema_version') != 4 or data.get('kind') != 'completed_build':
        raise ValueError('migration requires a schema-4 bundle with sample IDs')
    if digest({k: v for k, v in data.items() if k != 'sha256'}) != data.get('sha256'):
        raise ValueError('source manifest changed')
    print('Validating source database', flush=True)
    bound = load_products(source, data)
    mapping = load_tables(bound['mapping']['path'])
    if set(mapping['tables']) != set(data['products']):
        raise ValueError('mapping does not cover database samples exactly')
    run_provenance = record(legacy_run / 'run.json')
    run_config = json.loads(verify(run_provenance).read_text())['config']
    if (run_config['translation'] != data['translation'] or
            run_config['odb']['node'] != data['odb']['node'] or
            run_config['phylogeny']['lineage'] != data['lineage']):
        raise ValueError('legacy expression build settings differ')
    multimap = run_config['tpm']['multimap']
    if multimap not in {'error', 'drop', 'split'}:
        raise ValueError('invalid legacy TPM policy')
    merge_record = record(legacy_run / 'orthogroups/mapping/merge_qc.json')
    merge = json.loads(verify(merge_record).read_text())
    sources = {}
    for entry in merge['sources']:
        snapshot = json.loads(verify(entry).read_text())
        if snapshot['node'] != data['odb']['node'] or snapshot['version'] != data['odb']['version']:
            raise ValueError('legacy mapping reference differs')
        sources[snapshot['annotations']['sha256']] = {
            'record': entry, 'proteins': {p['species']: p['sha256'] for p in snapshot['proteins']}}
    print(f'Checking expression provenance for {len(data["products"])} samples', flush=True)
    evidence = []
    for i, (species, product) in enumerate(data['products'].items(), 1):
        evidence.append(expression_evidence(species, product, mapping['tables'][species],
                                            legacy_run, receipts, sources, multimap))
        if i % 500 == 0: print(f'Expression provenance: {i}/{len(data["products"])}', flush=True)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.stat().st_dev != destination.parent.stat().st_dev:
        raise ValueError('migration requires the same filesystem for hard links')
    staging = Path(tempfile.mkdtemp(prefix=f'.{destination.name}.migrate-', dir=destination.parent))
    database = staging / 'database'; database.mkdir()
    try:
        for i, entry in enumerate(bound['files'], 1):
            relative = Path(entry['path']).relative_to(source.parent)
            target = database / relative; target.parent.mkdir(parents=True, exist_ok=True)
            os.link(verify(entry), target)
            if i % 5000 == 0: print(f'Linked files: {i}/{len(bound["files"])}', flush=True)
        provenance = database / 'provenance/layout_migration'; provenance.mkdir(parents=True)
        for name, entry in [('source_manifest.json', source_record), ('legacy_run.json', run_provenance),
                            ('legacy_mapping_qc.json', merge_record)]:
            os.link(verify(entry), provenance / name)
        tasks = [(e, str(database), data['products'][e['species']]['abundance'],
                  dict(mapping['tables'][e['species']]['table'],
                       path='odb/' + mapping['tables'][e['species']]['table']['path'])) for e in evidence]
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for i, (species, expression, qc) in enumerate(pool.map(rewrite_expression, tasks, chunksize=1), 1):
                data['products'][species].update(expression=expression, expression_qc=qc)
                if i % 100 == 0 or i == len(tasks): print(f'Expression labels: {i}/{len(tasks)}', flush=True)
        # New links change inode ctime: capture final identities after all links exist.
        inventory = {}
        for entry in data['files']:
            relative = entry['path']
            inventory[relative] = dict(entry, stat=stat_identity(database / relative))
        mapping = copy.deepcopy(mapping)
        for entry in mapping['tables'].values():
            entry['table'] = dict(inventory['odb/' + entry['table']['path']], path=entry['table']['path'])
        write_json(database / data['mapping']['path'], mapping)
        relative = data['mapping']['path']
        inventory[relative] = dict(record(database / relative), path=relative)
        for product in data['products'].values():
            for key in ('expression', 'expression_qc'):
                inventory[product[key]['path']] = product[key]
            for key in ('cds', 'busco', 'abundance', 'protein', 'translation', 'expression', 'expression_qc'):
                product[key] = inventory[product[key]['path']]
        report = {'created_at': now(), 'source': source_record, 'legacy_run': run_provenance,
                  'samples': len(tasks), 'source_schema': 4, 'destination_schema': 5,
                  'scientific_recomputation': False, 'expression_numeric_values_preserved': True,
                  'source_retained': True, 'destination': str(destination),
                  'implementation': record(__file__)}
        write_json(provenance / 'migration.json', report)
        for path in provenance.iterdir():
            relative = str(path.relative_to(database))
            inventory[relative] = dict(record(path), path=relative)
        data.update(schema_version=5, tpm={'multimap': multimap}, files=list(inventory.values()),
                    mapping=inventory[data['mapping']['path']])
        data.pop('sha256'); data['sha256'] = digest(data)
        write_json(database / 'manifest.json', data)
        print('Validating migrated database', flush=True)
        load_products(database / 'manifest.json', data)
        verify(source_record); verify(run_provenance); verify(merge_record)
        for directory in ('work', 'logs', 'downstream'):
            (staging / directory).mkdir()
        write_json(staging / 'logs/migration.json', report)
        os.rename(staging, destination)
        publish_pointer(destination, destination / 'database/manifest.json')
    finally:
        if staging.exists(): shutil.rmtree(staging)
    print(f'Migrated {len(tasks)} samples: {destination / "database"}', flush=True)
    return destination


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'destination', 'legacy-run', 'receipts'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--workers', type=int, default=1)
    migrate(**vars(parser.parse_args()))
