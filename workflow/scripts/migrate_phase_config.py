#!/usr/bin/env python3
"""Convert old dataset/config settings without changing artifacts or caches."""
import argparse
import copy
import json
from pathlib import Path
import yaml
from phase_config import deep_merge, read_yaml


def migrate(root, legacy_dataset, legacy_config=None, build_output='config/build.local.yaml', analysis_output='config/analysis.local.yaml'):
    root=Path(root).resolve()
    def path(value): return (root/Path(value)).resolve()
    source=path(legacy_dataset)
    if source.is_dir(): source /= 'dataset.json'
    if source.suffix=='.json':
        manifest=json.loads(source.read_text())
        cfg=manifest['config']
        pipeline=manifest.get('selection_analysis',manifest['analysis'])
        metadata=source.parent/'metadata.tsv'
    else:
        cfg=read_yaml(source)
        pipeline=read_yaml(path(legacy_config)) if legacy_config else read_yaml(root/'workflow/pipeline_defaults.yaml')
        if cfg.get('analysis_config'): pipeline=deep_merge(pipeline,read_yaml(path(cfg['analysis_config'])))
        metadata=cfg['metadata']
    build=read_yaml(root/'config/build.yaml')
    build['excluded_accessions']=cfg.get('excluded_accessions')
    build.update(metadata=str(metadata),store=cfg['store'],genegalleon=copy.deepcopy(cfg['genegalleon']),translation=pipeline['translation'],busco={'lineage':pipeline['phylogeny']['lineage']})
    build['odb']={k:v for k,v in pipeline['odb'].items() if k!='incremental'}
    build['odb']['chunk_size']=cfg.get('odb_chunk_size',20)
    old=cfg['slurm']; slurm=build['slurm']
    for key in ('partition','account','concurrency','array_size'):
        if key in old: slurm[key]=old[key]
    slurm['jobs']=old.get('downstream_jobs',10)
    slurm['stages']=dict(copy.deepcopy(old['stages']))
    slurm['stages']['controller']=slurm['stages'].pop('downstream')
    if old.get('downstream_profile'):
        profile_path=path(old['downstream_profile'])/'config.yaml'
        if profile_path.exists():
            profile=read_yaml(profile_path)
            defaults=profile.get('default-resources',{})
            slurm['default_resources']={k:v for k,v in defaults.items() if k in {'mem_mb','runtime'}}
            for key in ('partition','account'):
                if not slurm.get(key) and defaults.get('slurm_'+key): slurm[key]=defaults['slurm_'+key]
            slurm['rules']=copy.deepcopy(profile.get('set-resources',{}))
            for rule,cpus in profile.get('set-threads',{}).items(): slurm['rules'].setdefault(rule,{})['cpus']=cpus
    analysis=read_yaml(root/'config/analysis.yaml')
    for key in ('seed','trait','selection','tpm','alignment','kegg','phylogeny','exclude_species'):
        analysis[key]=copy.deepcopy(pipeline[key])
    analysis['phylogeny'].pop('lineage',None)
    inputs=path(pipeline.get('input_root','input'))
    analysis['inputs']={key:str(inputs/filename) if (inputs/filename).is_file() else None for key,filename in
                        [('species_trait','species_trait.tsv'),('species_list','species_list.txt'),('calibrations','calibrations.tsv')]}
    analysis['slurm']={k:copy.deepcopy(v) for k,v in slurm.items() if k not in {'concurrency','array_size','stages'}}
    analysis['slurm']['stages']={'controller':copy.deepcopy(slurm['stages']['controller'])}
    for rule in ('odb_map','prepare_odb_reference'): analysis['slurm']['rules'].pop(rule,None)
    outputs=[path(build_output),path(analysis_output)]
    if len(set(outputs))!=2 or any(p.exists() for p in outputs): raise ValueError('migration outputs must be distinct new files; existing files are not overwritten')
    for output,value in zip(outputs,(build,analysis)):
        output.parent.mkdir(parents=True,exist_ok=True)
        with output.open('x') as handle: yaml.safe_dump(value,handle,sort_keys=False)
    return outputs


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',default='.')
    parser.add_argument('--legacy-dataset',required=True,help='old dataset YAML or datasets/<name>/dataset.json')
    parser.add_argument('--legacy-config',help='old resolved/flat Snakemake configuration')
    parser.add_argument('--build-output',default='config/build.local.yaml')
    parser.add_argument('--analysis-output',default='config/analysis.local.yaml')
    for path in migrate(**vars(parser.parse_args())): print(path)
