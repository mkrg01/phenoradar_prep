"""Tiny portable databases representing completed native pipeline products."""
import json
from pathlib import Path

from common import write_json
from aggregate_tpm import aggregate
from dataset_assets import record, resolve
from mapping_tables import load_tables, relative_file, write_tables
from mapping_fixtures import make_mapping
from portable_build import publish_products
from translate_cds import fasta_ids


def database_from_stages(root, store, items, name='source'):
    build = Path(root) / 'resources' / name
    build.mkdir(parents=True, exist_ok=True)
    products = {}
    genes, pairs = [], []
    for item in items:
        species, run = item['species'], item['row']['run']
        stages = resolve(store, item, 'embryophyta_odb12', need_full=True)
        ref, busco, quant = [stages[k] for k in ('reference', 'busco', 'quant')]
        protein = build / 'proteins' / (item['odb_species'] + '_protein.fa')
        protein.parent.mkdir(exist_ok=True)
        ids = fasta_ids(ref['cds']['path'], compressed=True)
        protein.write_text(''.join(f'>{gene}\nMK*\n' for gene in ids))
        translation = protein.with_suffix('.json')
        write_json(translation, {'cds': ref['cds'], 'protein': record(protein), 'seqkit': 'test',
                                'translation_table': 1, 'sequences': len(ids)})
        genes.extend((gene, species) for gene in ids)
        pairs.extend((gene, f'OG{i}') for i, gene in enumerate(ids[:2], 1))
        products[species] = {'row': item['row'], 'odb_species': item['odb_species'],
            'reference_id': ref['reference_id'], 'counts': busco['counts'],
            'conditions': {s: p.get('provenance', {}).get('condition') for s,p in
                           [('assembly',ref),('busco',busco),('quant',quant)]},
            'raw_inputs': quant.get('provenance', {}).get('raw_inputs', {}),
            'cds': ref['cds'], 'busco': busco['full'], 'abundance': quant['abundance'],
            'protein': record(protein), 'translation': record(translation)}
    mapping = make_mapping(build / 'mapping/snapshot.json', genes, pairs)
    snapshot = load_tables(mapping)
    entries = {s: dict(e, protein_sha256=products[s]['protein']['sha256'],
                        table=relative_file(mapping.parent, e['table'])) for s,e in snapshot['tables'].items()}
    mapping = write_tables(build / 'bound_mapping', entries)
    # Materialize the stage products using the production publication layout.
    input_root = build / 'input'
    from dataset_assets import link_file
    from common import write_tsv
    files = []
    for species, p in products.items():
        for key, relative in [('cds', f'cds/{species}_longestCDS.fa.gz'),
                              ('busco', f'busco/full/{species}.busco.full.tsv'),
                              ('abundance', f'quant/{species}/{p["row"]["run"]}/{p["row"]["run"]}_abundance.tsv')]:
            target = input_root / relative
            link_file(p[key]['path'], target)
            p[key] = record(target); files.append(p[key])
    rows = [i['row'] for i in items]
    write_tsv(input_root/'metadata.tsv', list(rows[0]), rows)
    from dataset_assets import COUNTS
    write_tsv(input_root/'busco/summary.tsv', ['Species', *COUNTS],
              [dict(Species=s, **p['counts']) for s,p in products.items()])
    files += [record(input_root/'metadata.tsv'), record(input_root/'busco/summary.tsv')]
    samples = build / 'samples.tsv'
    write_tsv(samples, ['species', 'run', 'abundance'],
              [dict(species=s, run=p['row']['run'], abundance=p['abundance']['path'])
               for s,p in products.items()])
    for product in products.values():
        run = product['row']['run']
        expression = build / 'expression' / f'{run}.tsv'
        qc = expression.with_suffix('.qc.json')
        aggregate(samples, run, mapping, expression, qc)
        product.update(expression=record(expression), expression_qc=record(qc))
    data = {'schema_version': 5, 'kind': 'completed_build', 'build_id': name, 'created_at': 'fixture',
            'input': str(input_root), 'fields': list(rows[0]), 'translation': {'table': 1},
            'lineage': 'embryophyta_odb12', 'odb': {'version': 'v12', 'node': 3193},
            'mapping': record(mapping), 'products': products, 'files': files, 'excluded_runs': [],
            'tpm': {'multimap': 'error'}}
    return publish_products(build, data).parent
