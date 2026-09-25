#!/usr/bin/env python3
"""Freeze and run independent analyses of a completed build."""
import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml
from build_products import load_complete
from portable_build import completion_path, input_entries
from common import now, read_tsv, write_json
from configuration import validate_analysis, validate_keys
from dataset import absolute, implementation, inside, load_execution
from dataset_assets import SAFE, digest, link_file, locked, record, verify
from phase_config import deep_merge, execution_settings, read_yaml, validate_slurm, write_profile

TARGETS = ('all','alignments','kegg','phylogeny','phylogeny_prepare','taxonomy_check',
           'contrast_pairs','phylogeny_calibrations','timetree','phenoradar_inputs')
ANALYSIS_KEYS = {'build','inputs','seed','trait','selection','tpm','alignment','kegg','phylogeny','exclude_species','slurm'}


def settings(root, config, build=None):
    root = Path(root).resolve()
    override = read_yaml(config)
    if set(override) - ANALYSIS_KEYS: raise ValueError('unknown analysis settings: ' + ', '.join(sorted(set(override) - ANALYSIS_KEYS)))
    cfg = deep_merge(read_yaml(root / 'config/analysis.yaml'), override)
    if 'lineage' in cfg['phylogeny']: raise ValueError('BUSCO lineage belongs to build.yaml')
    if set(cfg['inputs']) - {'species_trait','species_list','calibrations'}: raise ValueError('unknown analysis input')
    validate_slurm(cfg['slurm'], build=False)
    source = inside(root, absolute(root, build or cfg['build']))
    source = inside(root, completion_path(source))
    completed = load_complete(source)
    if completed['schema_version'] == 1 and Path(completed['root']).resolve() != root:
        raise ValueError('completed build must belong to this project for container mounts')
    base = read_yaml(root / 'workflow/pipeline_defaults.yaml')
    resolved = deep_merge(base, {k:v for k,v in cfg.items() if k not in {'build','inputs','slurm'}})
    resolved['translation'] = completed['translation']
    resolved['phylogeny']['lineage'] = completed['lineage']
    resolved['odb'].update(node=completed['odb']['node'], incremental=False, existing_results=None)
    resolved['build_manifest'] = str(source)
    validate_keys(resolved); validate_analysis(resolved)
    threshold = resolved['selection']['busco_threshold']
    if type(threshold) not in (int,float) or not 0 <= threshold <= 1: raise ValueError('selection.busco_threshold must be between 0 and 1')
    names = set(completed['products'])
    requested = names
    inputs = cfg['inputs']
    if resolved['selection']['species_list']:
        if not inputs.get('species_list'): raise ValueError('selection.species_list requires inputs.species_list')
        values = absolute(root, inputs['species_list']).read_text().splitlines()
        if not values or len(values) != len(set(values)) or set(values) - names:
            raise ValueError('species list must contain unique species from the completed build')
        requested = set(values)
    exclusions = resolved['exclude_species']
    if not isinstance(exclusions,list) or any(not isinstance(s,str) for s in exclusions) or len(exclusions) != len(set(exclusions)) or set(exclusions) - names:
        raise ValueError('exclude_species must contain unique species from the completed build')
    requested = requested - set(exclusions)
    report = []
    for name, product in sorted(completed['products'].items()):
        counts = product['counts']
        fraction = (counts['busco_cds_single'] + counts['busco_cds_duplicated']) / counts['busco_cds_total']
        reason = 'exclude_species' if name in exclusions else 'outside_species_list' if name not in requested else 'below_busco_threshold' if fraction < threshold else ''
        report.append({'species':name, 'selected':not bool(reason), 'reason':reason, 'busco_complete_fraction':fraction})
    if not any(r['selected'] for r in report): raise ValueError('no species passed analysis selection')
    for key,value in inputs.items():
        if value and not absolute(root,value).is_file(): raise ValueError(f'analysis input missing: {key}: {value}')
    return cfg, resolved, source, completed, sorted(requested), report


def prepare(root, name, config, build=None):
    root = Path(root).resolve()
    if not SAFE.fullmatch(name): raise ValueError('analysis name must be a simple directory name')
    if name.startswith('build_'): raise ValueError('analysis names beginning with build_ are reserved for build outputs')
    target = root / 'analyses' / name
    if target.exists() or (root / 'results' / name).exists(): raise ValueError('analysis/run name already exists; resume it or choose a new name')
    cfg, resolved, source, completed, requested, report = settings(root, config, build)
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.prepare-', dir=target.parent))
    try:
        inputs = staging / 'input'
        for entry, relative in input_entries(completed):
            link_file(verify(entry), inputs / relative)
        for key, filename in [('species_trait','species_trait.tsv'), ('calibrations','calibrations.tsv')]:
            if cfg['inputs'].get(key): shutil.copy2(absolute(root,cfg['inputs'][key]), inputs / filename)
        # One frozen selection combines the user list and explicit exclusions.
        (inputs / 'species_list.txt').write_text(''.join(s + '\n' for s in requested))
        resolved['selection']['species_list'] = True
        resolved['exclude_species'] = []
        resolved['run_name'] = name
        resolved['input_root'] = str(target / 'input')
        (staging / 'pipeline.yaml').write_text(yaml.safe_dump(resolved, sort_keys=False))
        (staging / 'analysis.yaml').write_text(yaml.safe_dump(cfg, sort_keys=False))
        records = [dict(record(p), path=str(target / p.relative_to(staging))) for p in sorted(staging.rglob('*')) if p.is_file()]
        manifest = {'schema_version':1, 'name':name, 'root':str(root), 'created_at':now(), 'config':cfg,
                    'pipeline':resolved, 'build':record(source), 'files':records, 'implementation':implementation(root), 'selection':report}
        write_json(staging / 'analysis.json', manifest)
        write_json(staging / 'checksums.json', {'analysis.json':record(staging / 'analysis.json')['sha256']})
        os.rename(staging,target)
    finally:
        if staging.exists(): shutil.rmtree(staging)
    return target


def load(path, check_code=False):
    path = Path(path).resolve()
    hashes = json.loads((path / 'checksums.json').read_text())
    if record(path / 'analysis.json')['sha256'] != hashes['analysis.json']: raise ValueError('frozen analysis changed')
    manifest = json.loads((path / 'analysis.json').read_text())
    verify(manifest['build']); load_complete(manifest['build']['path'])
    for entry in manifest['files']: verify(entry)
    if check_code:
        for entry in manifest['implementation']: verify(entry)
    return manifest


def status(path):
    manifest = load(path)
    done = Path(path) / 'completed.json'
    if done.exists():
        for entry in json.loads(done.read_text())['files']: verify(entry)
    return {'name':manifest['name'], 'build':manifest['build']['path'], 'state':'complete' if done.exists() else 'prepared',
            'selection':manifest['selection'],
            'submissions':[json.loads(p.read_text()) for p in sorted((Path(path)/'jobs').glob('submission_*.json')) if not p.name.endswith('.resources.json')]}


def run(path, target='all', execution=None, local=False, cores=1, mem_mb=8000):
    path = Path(path).resolve()
    manifest = load(path, check_code=True)
    root = Path(manifest['root'])
    validate_analysis(manifest['pipeline'], [target])
    slurm = load_execution(execution) if execution else manifest['config']['slurm']
    profile = write_profile(path / 'jobs' / (Path(execution).stem if execution else 'local') / 'profile', slurm)
    command = [str(root/'run_pipeline.sh')]
    if local: command += ['--cores',str(cores),'--resources',f'mem_mb={mem_mb}']
    else: command += ['--slurm','--profile',str(profile),'--jobs',str(slurm['jobs'])]
    command += ['--configfile',str(path/'pipeline.yaml')]
    with locked(path/'.run.lock'):
        subprocess.run([*command,'--',target], cwd=root, check=True)
        if target == 'all':
            subprocess.run([str(root/'run_pipeline.sh'),'--cores','1','--resources','mem_mb=4000',
                            '--configfile',str(path/'pipeline.yaml'),'--','phenoradar_inputs'], cwd=root, check=True,
                           env={k:v for k,v in os.environ.items() if not k.startswith('SLURM_')})
            output = root/'results'/manifest['name']
            write_json(path/'completed.json', {'created_at':now(), 'build':manifest['build'],
                       'files':[record(p) for p in sorted(output.rglob('*')) if p.is_file()]})


def submit(path, target='all', dry_run=False, resources=None):
    path = Path(path).resolve()
    manifest = load(path, check_code=True)
    validate_analysis(manifest['pipeline'], [target])
    slurm = execution_settings(manifest['config']['slurm'], resources, build=False)
    jobs = path/'jobs'; jobs.mkdir(exist_ok=True); (jobs/'logs').mkdir(exist_ok=True)
    with locked(jobs/'.submit.lock'):
        previous = [p for p in jobs.glob('submission_*.json') if not p.name.endswith('.resources.json')]
        if not dry_run:
            ids = []
            for p in previous:
                old = json.loads(p.read_text())
                if old['state'] in {'submitting','unknown'}: raise ValueError(f'unresolved submission in {p}; inspect Slurm before retrying')
                if old.get('job_id'): ids.append(old['job_id'])
            if ids:
                import pwd
                queued = subprocess.check_output(['squeue','--noheader','--user',pwd.getpwuid(os.getuid()).pw_name,'--format=%i'],text=True).splitlines()
                if any(j.strip().split('_',1)[0] in ids for j in queued): raise ValueError('analysis still has queued/running jobs')
        number = 1 + max([int(p.name.split('_')[1].split('.')[0]) for p in jobs.glob('submission_*.resources.json')] + [0])
        execution = jobs/f'submission_{number:04d}.resources.json'
        write_json(execution, {'slurm':slurm,'sha256':digest(slurm)})
        script = jobs/f'{number:04d}_analysis.sh'
        invocation = [sys.executable,str(Path(manifest['root'])/'workflow/scripts/analysis.py'),'run','--analysis',str(path),
                      '--target',target,'--execution',str(execution)]
        script.write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec ' + shlex.join(invocation) + '\n')
        script.chmod(0o755)
        job = slurm['stages']['controller']
        cmd = ['sbatch','--parsable','--nodes=1','--ntasks=1',f"--cpus-per-task={job['cpus']}",f"--mem={job['mem_mb']}M",f"--time={job['time']}",
               '--chdir='+manifest['root'],'--job-name='+manifest['name'],f'--output={jobs}/logs/{number:04d}_%j.out',f'--error={jobs}/logs/{number:04d}_%j.err']
        for key in ('partition','account'):
            if slurm.get(key): cmd.append(f'--{key}={slurm[key]}')
        cmd.append(str(script))
        if dry_run:
            print(shlex.join(cmd)); return cmd
        receipt = jobs/f'submission_{number:04d}.json'
        result = {'created_at':now(),'target':target,'command':cmd,'resources':record(execution),'state':'submitting','job_id':None}
        write_json(receipt,result)
        try: identifier = subprocess.check_output(cmd,text=True).strip().split(';')[0]
        except (OSError,subprocess.CalledProcessError):
            result['state']='rejected'; write_json(receipt,result); raise
        if not identifier.isdigit():
            result['state']='unknown'; write_json(receipt,result); raise ValueError('unrecognized sbatch result; inspect queue before retrying')
        result.update(state='submitted',job_id=identifier); write_json(receipt,result)
        return cmd


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    for name in ('plan','prepare'):
        p=sub.add_parser(name); p.add_argument('--root',default='.'); p.add_argument('--config',default='config/analysis.yaml'); p.add_argument('--build')
        if name=='prepare': p.add_argument('--name',required=True)
    for name in ('status','submit','run'):
        p=sub.add_parser(name); p.add_argument('--analysis',required=True)
        if name!='status': p.add_argument('--target',choices=TARGETS,default='all')
        if name=='submit': p.add_argument('--dry-run',action='store_true'); p.add_argument('--resources')
        if name=='run':
            p.add_argument('--execution'); p.add_argument('--local',action='store_true'); p.add_argument('--cores',type=int,default=1); p.add_argument('--mem-mb',type=int,default=8000)
    args=parser.parse_args()
    if args.command in {'plan','prepare'}:
        root=Path(args.root).resolve(); config=absolute(root,args.config)
        if args.command=='plan': print(json.dumps(settings(root,config,args.build)[-1],indent=2))
        else: print(prepare(root,args.name,config,args.build))
    elif args.command=='status': print(json.dumps(status(args.analysis),indent=2))
    elif args.command=='submit': submit(args.analysis,args.target,args.dry_run,args.resources)
    else: run(args.analysis,args.target,args.execution,args.local,args.cores,args.mem_mb)


if __name__=='__main__':
    try: main()
    except (ValueError,OSError,subprocess.CalledProcessError) as error:
        print(f'analysis: {error}',file=sys.stderr); sys.exit(1)
