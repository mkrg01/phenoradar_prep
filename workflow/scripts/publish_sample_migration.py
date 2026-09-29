#!/usr/bin/env python3
"""Publish a completed build from a validated one-time ID migration.

No scientific executables are invoked. Default configuration is left unchanged
until the caller has reviewed the completed bundle and its validation report.
"""
import argparse
import json
from pathlib import Path

import yaml

from common import now, read_tsv, write_json
from dataset_assets import link_file, record, verify
from dataset import plan, prepare, load, materialize
from build_products import complete, load_complete
from mapping_tables import subset
from prepare_metadata import prepare as prepare_metadata


def publish(root, migration, name, config='config/build.yaml'):
    root=Path(root).resolve();migration=Path(migration).resolve()
    if not migration.is_relative_to(root):
        raise ValueError('migration must be inside the project')
    done=json.loads((migration/'completed.json').read_text())
    for key in ('metadata','mapping','id_map','source_metadata','source_products'):verify(done[key])
    cfg=yaml.safe_load((root/config).read_text())
    # Keep the curated metadata as the ongoing user input. Build derives sample IDs.
    if (root/cfg['metadata']).resolve()!=Path(done['source_metadata']['path']).resolve():
        raise ValueError('build metadata differs from migration source')
    cfg['name']=name
    cfg['store']=str((migration/'dataset_assets').relative_to(root))
    cfg['odb']['cache_dir']=str((migration/'odb_cache').relative_to(root))
    local=migration/'build.yaml'
    text=yaml.safe_dump(cfg,sort_keys=False)
    if local.exists() and local.read_text()!=text:
        raise ValueError('migration build settings changed')
    local.write_text(text)
    *_,report=plan(root,local)
    active=[r for r in report if r['assembly']!='excluded']
    pending=[r for r in active if any(r[s]!='reuse' for s in ('assembly','busco','quant'))]
    if pending:
        raise ValueError('migration does not cover all selected build samples: '+json.dumps(pending[:10]))
    print(f'Validated reuse of assembly/BUSCO/quant for {len(active)} samples',flush=True)
    build=root/'builds'/name
    if not build.exists():prepare(root,name,local)
    frozen=load(build,check_code=True)
    if frozen['config']['store']!=str((migration/'dataset_assets').resolve()):
        raise ValueError('existing build uses a different registry')
    if (build/'completed.json').exists():
        data=load_complete(build)
        return build
    inputs=materialize(build)
    results=root/'results'/('build_'+name)
    prepare_metadata(inputs/'metadata.tsv',inputs/'busco/summary.tsv',inputs/'cds',inputs/'quant',
                     root/'resources/taxonomy/taxa.sqlite',results/'metadata',threshold=0,missing_taxonomy='allow')
    samples=read_tsv(results/'metadata/samples.tsv')
    for i,row in enumerate(samples,1):
        source=migration/'input/proteins'/f'{row["odb_species"]}_protein.fa'
        target=results/'proteins'/source.name
        link_file(source,target)
        link_file(source.with_suffix('.json'),target.with_suffix('.json'))
        if i%500==0:print(f'Staged proteins {i}/{len(samples)}',flush=True)
    subset(verify(done['mapping']),{r['species'] for r in samples},results/'orthogroups/mapping')
    print('Publishing and validating completed products',flush=True)
    complete(build)
    data=load_complete(build)
    if len(data['products'])!=len(active) or data['schema_version']!=4:
        raise ValueError('completed sample membership or schema differs')
    write_json(migration/'build_validation.json',{
        'created_at':now(),'build':str(build),'completed':record(build/'completed.json'),
        'sample_count':len(data['products']),'biological_species_count':len({p['row']['species_id'] for p in data['products'].values()}),
        'schema_version':data['schema_version'],'scientific_recomputation':False,
        'assembly_busco_quant_reused':len(active),'mapping_reused':len(active),
        'prepared_config':record(local)})
    print(f'Completed sample build: {build}',flush=True)
    return build


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',default='.')
    parser.add_argument('--migration',required=True)
    parser.add_argument('--name',required=True)
    parser.add_argument('--config',default='config/build.yaml')
    publish(**vars(parser.parse_args()))
