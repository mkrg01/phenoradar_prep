#!/usr/bin/env python3
"""Compare one production BUSCO chain with and without Slurm grouping."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

def identity(path):
    s = path.stat()
    return s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('project', 'workdir', 'image', 'fasta', 'manifest'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--require-identical-trees', action='store_true',
                        help='also require byte-identical VeryFastTree output')
    args = parser.parse_args()
    root, work = args.project.resolve(), args.workdir.resolve()
    sources = [args.image.resolve(), args.fasta.resolve(), args.manifest.resolve()]
    protected = {p: identity(p) for p in sources}
    work.mkdir(parents=True, exist_ok=False)
    rules = (root / 'workflow/rules/phylogeny.smk').read_text()
    rules = 'def busco_marker_mem_mb' + rules.split('def busco_marker_mem_mb', 1)[1]
    rules = rules.split('rule merge_busco_gene_trees:', 1)[0]
    rules = rules.replace('"../envs/phylogeny.yaml"', repr(str(root / 'workflow/envs/phylogeny.yaml')))
    real_sbatch = os.environ.get('PHENORADAR_REAL_SBATCH') or shutil.which('sbatch')
    if not real_sbatch or Path(real_sbatch).resolve() == root / 'workflow/slurm/sbatch':
        raise ValueError('a real sbatch executable is required')
    products = ('alignments/raw/marker.faa', 'alignments/raw/marker.json',
                'alignments/marker.faa', 'alignments/marker.json', 'alignments/marker.columns.tsv',
                'gene_trees/marker.nwk', 'gene_trees/marker.json')
    results, commands = {}, {}
    for mode in ('grouped', 'separate'):
        case = work / mode
        out = case / 'results/test/phylogeny/representatives'
        (out / 'markers').mkdir(parents=True)
        (out / 'plan').mkdir()
        shutil.copy2(sources[1], out / 'markers/marker.faa')
        shutil.copy2(sources[2], out / 'plan/species.tsv')
        (out / 'markers/.snakemake_timestamp').touch()
        prefix = f'''import json
import sys
sys.path.insert(0, {str(root / "workflow/scripts")!r})
from layout import REPRESENTATIVES as REPRESENTATIVES_REL
from slurm_submit import configure_scheduler
onstart:
    configure_scheduler(workflow)
containerized: {str(sources[0])!r}
OUT = "results/test"
LOG = "logs/test"
SCRIPTS = {str(root / "workflow/scripts")!r}
PHYLO_RUN = OUT + "/{{phylo_branch}}"
PHY = {{"min_taxa": 4, "min_protein_length": 100, "trimal_mode": "gappyout"}}
PYTHON = "python"
RANDOM_SEED = 12345
rule all:
    input: OUT + "/phylogeny/representatives/gene_trees/marker.nwk"
rule collect_busco_markers:
    output: fasta=directory(PHYLO_RUN + "/markers")
'''
        (case / 'Snakefile').write_text(prefix + (rules if mode == 'grouped' else
                                                 rules.replace('    group: "busco_marker"\n', '')))
        env = dict(os.environ, PHENORADAR_REAL_SBATCH=real_sbatch,
                   PHENORADAR_SUBMIT_PYTHON=sys.executable,
                   PHENORADAR_SUBMIT_STATE=str(case / 'submissions'),
                   PATH=str(root / 'workflow/slurm') + os.pathsep + os.environ['PATH'],
                   XDG_CACHE_HOME=str(work / 'cache'))
        cmd = [sys.executable, '-m', 'snakemake', '--snakefile', str(case / 'Snakefile'),
               '--directory', str(case), '--executor', 'slurm', '--jobs', '1',
               '--deployment-method', 'apptainer', 'conda', '--conda-base-path', '/opt/conda',
               '--apptainer-args', '--cleanenv --bind ' + str(root.parent),
               '--slurm-no-account', '--slurm-status-command', 'sacct',
               '--slurm-init-seconds-before-status-checks', '5', '--seconds-between-status-checks', '5',
               '--default-resources', 'mem_mb=8000', 'runtime=10',
               '--latency-wait', '30', '--scheduler', 'greedy', '--printshellcmds', 'all']
        print(json.dumps({'phase': 'running', 'mode': mode}), flush=True)
        with (work / (mode + '.log')).open('w') as log:
            subprocess.run(cmd, cwd=case, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        worker_ids = re.findall(r'SLURM jobid (\d+)', (work / (mode + '.log')).read_text())
        if len(worker_ids) != (1 if mode == 'grouped' else 3):
            raise ValueError(f'unexpected {mode} worker count: {worker_ids}')
        for name in products:
            if not (out / name).is_file() or not (out / name).stat().st_size:
                raise ValueError(f'missing/empty product: {out / name}')
        qc = json.loads((out / 'gene_trees/marker.json').read_text())
        if qc['status'] != 'retained':
            raise ValueError(f'marker did not exercise the real tools: {qc}')
        results[mode] = {'workers': worker_ids, 'taxa': qc['taxa'], 'sites': qc['sites']}
        commands[mode] = cmd, env
    for name in ('alignments/raw/marker.faa', 'alignments/marker.faa',
                 'alignments/marker.columns.tsv'):
        rel = 'results/test/phylogeny/representatives/' + name
        if (work / 'grouped' / rel).read_bytes() != (work / 'separate' / rel).read_bytes():
            raise ValueError(f'grouped output differs: {name}')
    tree_dir = 'results/test/phylogeny/representatives/gene_trees'
    tree_reports = [json.loads((work / mode / tree_dir / 'marker.json').read_text())
                    for mode in ('grouped', 'separate')]
    left, right = tree_reports
    if (left['command'][:-1] != right['command'][:-1] or
            left['alignment']['sha256'] != right['alignment']['sha256'] or
            left['executable']['sha256'] != right['executable']['sha256'] or
            any(left[key] != right[key] for key in ('status', 'taxa', 'sites'))):
        raise ValueError('tree inputs, tools, arguments or retained taxa/sites differ')
    tree_identical = ((work / 'grouped' / tree_dir / 'marker.nwk').read_bytes() ==
                      (work / 'separate' / tree_dir / 'marker.nwk').read_bytes())
    if args.require_identical_trees and not tree_identical:
        raise ValueError('grouped tree differs; check tool reproducibility independently')
    out = work / 'grouped/results/test/phylogeny/representatives'
    completed = {out / name: identity(out / name) for name in products}
    cmd, env = commands['grouped']
    with (work / 'resume.log').open('w') as log:
        subprocess.run(cmd, cwd=work / 'grouped', env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    if 'SLURM jobid' in (work / 'resume.log').read_text():
        raise ValueError('completed group was submitted again')
    if any(identity(p) != saved for p, saved in {**protected, **completed}.items()):
        raise ValueError('source inputs/image or completed products changed')
    result = {'phase': 'verified_complete', 'cases': results, 'identical_alignment_outputs': 3,
              'tree_inputs_tools_arguments_identical': True, 'tree_bytes_identical': tree_identical,
              'protected_unchanged': True, 'resume_submissions': 0}
    (work / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
