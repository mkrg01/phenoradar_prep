#!/usr/bin/env python3
"""One-time, resumable species-ID to species_accession migration.

Writes a separate dataset, registry and ODB cache. Never edits source artifacts.
The destination can be removed before cutover; its per-sample receipts support
resuming after interruption. This tool is independent of normal build execution.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path

from common import now, read_tsv, write_json, write_tsv
from dataset_assets import (COUNTS, digest, identities, locked, record, verify,
                            register_reference, register_busco, register_quant)
from relabel_sample import relabel
from protein_cache import register_translation
from mapping_tables import (load_tables, relative_file, read_species, write_tables,
                            species_table)


def plan(metadata, store, products, legacy_results, old_cache, output):
    fields, items = identities(metadata)
    output = Path(output).resolve()
    store = Path(store).resolve()
    bundle = Path(products).resolve()
    # The old bundle format is intentionally consumed only by this migration.
    data = json.loads((bundle / 'manifest.json').read_text())
    if data.get('schema_version') != 3 or data.get('kind') != 'completed_build':
        raise ValueError('migration requires a schema-3 source products bundle')
    if digest({k:v for k,v in data.items() if k != 'sha256'}) != data.get('sha256'):
        raise ValueError('source manifest checksum differs')
    snapshot = load_tables(bundle / data['mapping']['path'])
    parts = list((Path(old_cache) / '.partitions').glob('*/receipt.json'))
    jobs, absent, names = [], [], set()
    for item in items:
        old = item['row']['species_id']
        if old in names:
            raise ValueError('legacy migration requires an unambiguous single run per old species')
        names.add(old)
        refs = list((store / old).glob('*/reference.json'))
        if not refs:
            absent.append(item)
            continue
        requested = item['row'].get('reference_id')
        if requested: refs = [p for p in refs if p.parent.name == requested]
        if len(refs) != 1:
            raise ValueError(f'ambiguous legacy reference: {old}')
        ref = json.loads(refs[0].read_text())
        bus = json.loads((refs[0].parent / 'busco.json').read_text())
        quant = refs[0].parent / 'quant' / (item['row']['run'] + '.json')
        if not quant.exists():
            raise ValueError(f'legacy run/reference association unavailable: {old}')
        quant = json.loads(quant.read_text())
        if ref['taxid'] != item['row']['taxid'] or quant['reference_id'] != ref['reference_id']:
            raise ValueError(f'legacy identity differs: {old}')
        product = data['products'].get(old)
        if product:
            if product['row']['run'] != item['row']['run'] or product['reference_id'] != ref['reference_id']:
                raise ValueError(f'bundle and registry differ: {old}')
            protein = dict(product['protein'], path=str(bundle / product['protein']['path']))
            translation = json.loads((bundle / product['translation']['path']).read_text())
        else:
            path = Path(legacy_results) / 'proteins' / (old.replace('-', '_') + '_protein.fa')
            translation = json.loads(path.with_suffix('.json').read_text())
            protein = dict(translation['protein'], path=str(path.resolve()))
        if translation['cds']['sha256'] != ref['reference_id']:
            raise ValueError(f'protein translation belongs to different CDS: {old}')
        entry = snapshot['tables'].get(old)
        if entry:
            if entry['protein_sha256'] != protein['sha256']:
                raise ValueError(f'mapping protein differs: {old}')
            entry = dict(entry, table=relative_file(bundle / 'odb', entry['table']))
        else:
            # Previously excluded samples still have verified legacy annotation shards.
            candidates = []
            for part in parts:
                receipt = json.loads(part.read_text())
                if old in receipt.get('shards', {}): candidates.append((part, receipt))
            if len(candidates) != 1:
                raise ValueError(f'need one verified legacy mapping shard: {old}')
            part, receipt = candidates[0]
            verify(relative_file(part.parent, receipt['shards'][old]))
            entry = species_table(old, old.replace('-', '_'), protein, receipt['annotation'],
                                  part.parent, output / 'legacy_tables', data['odb']['version'],
                                  data['odb']['node'], receipt['shards'][old])
        jobs.append({'item':item, 'old_id':old, 'reference':ref, 'busco':bus, 'quant':quant,
                     'protein':protein, 'translation':translation, 'mapping':entry,
                     'output':str(output), 'lineage':data['lineage'], 'table':data['translation']['table']})
    registered = {p.parent.parent.name for p in store.glob('*/*/reference.json')}
    if registered - names:
        raise ValueError(f'registered species absent from migration metadata: {sorted(registered - names)}')
    return fields, items, jobs, data, snapshot


def migrate_one(job):
    out = Path(job['output'])
    item, old = job['item'], job['old_id']
    sample, run = item['species'], item['row']['run']
    receipt = out / 'receipts' / (sample + '.json')
    fingerprint = digest(job)
    with locked(receipt.with_suffix('.lock')):
        if receipt.exists():
            saved = json.loads(receipt.read_text())
            if saved['source_identity'] != fingerprint:
                raise ValueError(f'migration source changed: {sample}')
            for entry in saved['files']: verify(entry)
            return saved
        root = out / 'input'
        paths = {'cds':root/'cds'/f'{sample}_longestCDS.fa.gz',
                 'busco':root/'busco/full'/f'{sample}.busco.full.tsv',
                 'quant':root/'quant'/sample/run/f'{run}_abundance.tsv',
                 'protein':root/'proteins'/f'{item["odb_species"]}_protein.fa',
                 'mapping':out/'tables'/f'{item["odb_species"]}.tsv.gz'}
        sources = {'cds':job['reference']['cds'], 'busco':job['busco']['full'],
                   'quant':job['quant']['abundance'], 'protein':job['protein'],
                   'mapping':job['mapping']['table']}
        formats = {'cds':'fasta','protein':'fasta','busco':'busco','quant':'quant','mapping':'mapping'}
        conversions = {}
        for kind, source in sources.items():
            conversions[kind] = relabel(verify(source), paths[kind], old, sample, formats[kind])
            verify(source)
        provenance = {'source':'sample_id_migration', 'old_sample_id':old,
                      'run':run, 'source_reference_id':job['reference']['reference_id'],
                      'source_identity':fingerprint}
        store = out/'dataset_assets'
        ref = register_reference(store,item,paths['cds'],dict(job['reference'].get('provenance',{}),**provenance))
        register_busco(store,ref,summary=job['busco']['counts'],full=paths['busco'],lineage=job['lineage'],
                       provenance=dict(job['busco'].get('provenance',{}),**provenance))
        register_quant(store,ref,item,paths['quant'],dict(job['quant'].get('provenance',{}),**provenance))
        translation = {k:job['translation'][k] for k in ('seqkit','translation_table','sequences')}
        translation.update(created_at=now(),cds=record(paths['cds']),protein=record(paths['protein']),migration=provenance)
        translation_path = paths['protein'].with_suffix('.json')
        write_json(translation_path,translation)
        register_translation(paths['cds'],paths['protein'],translation_path,store/'.proteins',job['table'])
        entry = dict(job['mapping'],odb_species=item['odb_species'],
                     protein_sha256=translation['protein']['sha256'],table=record(paths['mapping']))
        # Validate every new gene identifier and verify the stored mapping QC counts.
        from collections import defaultdict
        from relabel_sample import text_open
        genes = defaultdict(set)
        with text_open(paths['mapping']) as handle:
            next(handle)
            for line in handle:
                gene, og = line.rstrip('\n').split('\t')
                genes[gene]
                if og: genes[gene].add(og)
        from translate_cds import fasta_ids
        if set(genes) != set(fasta_ids(paths['protein'])):
            raise ValueError(f'mapping/protein IDs differ after migration: {sample}')
        if sum(map(len,genes.values())) != entry['qc']['unique_gene_og_pairs']:
            raise ValueError(f'mapping assignments changed: {sample}')
        saved = {'sample_id':sample,'old_id':old,'run':run,'reference_id':ref['reference_id'],
                 'source_identity':fingerprint,'files':[record(p) for p in [*paths.values(),translation_path]],
                 'mapping':entry,'counts':job['busco']['counts'],'conversions':conversions}
        write_json(receipt,saved)
        return saved


def migrate(metadata, store, products, legacy_results, old_cache, output, workers=4, limit=None):
    output = Path(output).resolve()
    output.mkdir(parents=True,exist_ok=True)
    with locked(output/'.migration.lock'):
        fields,items,jobs,data,snapshot = plan(metadata,store,products,legacy_results,old_cache,output)
        if limit is not None: jobs=jobs[:limit]
        write_tsv(output/'sample_id_map.tsv',['old_id','analysis_sample_id','run','taxid'],
                  [{'old_id':i['row']['species_id'],'analysis_sample_id':i['species'],'run':i['row']['run'],'taxid':i['row']['taxid']} for i in items])
        write_json(output/'plan.json',{'created_at':now(),'samples':len(items),'registered':len(jobs),
                   'source_metadata':record(metadata),'source_products':record(Path(products)/'manifest.json'),
                   'limited':limit is not None})
        print(f'Migrating {len(jobs)} registered samples with {workers} workers',flush=True)
        results = {}
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(migrate_one,j) for j in jobs]
            for future in as_completed(futures):
                saved=future.result();results[saved['sample_id']]=saved
                n=len(results)
                if n==1 or n%25==0 or n==len(jobs): print(f'Completed {n}/{len(jobs)}: {saved["sample_id"]}',flush=True)
        if limit is not None: return results
        if 'reference_id' not in fields: fields.append('reference_id')
        for item in items:
            item['row']['reference_id']=results[item['species']]['reference_id'] if item['species'] in results else ''
        write_tsv(output/'input/metadata.tsv',fields,[i['row'] for i in items])
        write_tsv(output/'input/busco/summary.tsv',['Species',*COUNTS],
                  [{'Species':s,**r['counts']} for s,r in sorted(results.items())])
        mapping=write_tables(output/'mapping',{s:r['mapping'] for s,r in results.items()},
                             data['odb']['version'],data['odb']['node'],
                             {'mode':'sample_id_migration','reused_species':len(results),'mapped_species':0},
                             snapshot.get('reference_sha256s',[]))
        from incremental_odb import import_snapshot
        cache=import_snapshot(mapping.parent,output/'odb_cache',data['odb']['version'],data['odb']['node'])
        write_json(output/'completed.json',{'created_at':now(),'samples':len(items),'registered':len(results),
                   'metadata':record(output/'input/metadata.tsv'),'mapping':record(mapping),'cache':str(cache),
                   'id_map':record(output/'sample_id_map.tsv'),'source_metadata':record(metadata),
                   'source_products':record(Path(products)/'manifest.json'),
                   'source_unchanged':True,'scientific_recomputation':False})
        return results


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('metadata','store','products','legacy-results','old-cache','output'):
        parser.add_argument('--'+name,required=True)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--limit',type=int,help='Pilot only; no completed migration is published')
    args=vars(parser.parse_args())
    if args['workers']<1 or (args['limit'] is not None and args['limit']<1): parser.error('workers/limit must be positive')
    migrate(**{k.replace('-','_'):v for k,v in args.items()})
