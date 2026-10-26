"""Audit and upgrade an unfinished v0.4.14 analysis to the batched workflow."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import shlex
import sys

from analysis import load
from analysis_build_index import create_index
from build_products import load_complete
from common import now, read_tsv, write_json
from dataset import implementation
from dataset_assets import locked, record, verify
from mapping_tables import load_tables, relative_file
from verification_cache import verification_cache

PREDECESSOR = 'dddd331b5a2770bbe7709370cf3695834cc261eb'
CHANGED = {'VERSION', 'run_pipeline.sh', 'workflow/rules/alignments.smk',
           'workflow/rules/build_inputs.smk', 'workflow/rules/tpm.smk',
           'workflow/scripts/analysis.py', 'workflow/scripts/build_products.py',
           'workflow/scripts/dataset.py', 'workflow/scripts/phase_config.py'}
ADDED = {'workflow/scripts/alignment_batches.py', 'workflow/scripts/analysis_build_index.py',
         'workflow/scripts/slurm_submit.py', 'workflow/scripts/upgrade_analysis.py', 'workflow/slurm/sbatch'}


def verified_imports(path, manifest, completed):
    """Validate source-equivalent imports before resetting only their scheduler metadata."""
    root = Path(manifest['root'])
    out = root / manifest['pipeline']['output_root']
    reusable = []
    for species, product in completed['products'].items():
        for kind, first, second, report_key, source_key in [
            ('import_build_protein', out / 'proteins' / f"{product['odb_species']}_protein.fa.gz",
             out / 'proteins' / f"{product['odb_species']}_protein.json", 'protein', 'translation'),
            ('import_build_expression', out / 'orthogroups/expression/runs' / f"{product['row']['run']}.tsv.gz",
             out / 'orthogroups/expression/runs' / f"{product['row']['run']}.qc.json", 'expression', 'expression_qc'),
        ]:
            if not first.exists() or not second.exists():
                continue
            verify(dict(product[report_key], path=str(first)))
            original = json.loads(verify(product[source_key]).read_text())
            saved = json.loads(second.read_text())
            if ({k: v for k, v in original.items() if k != report_key} !=
                    {k: v for k, v in saved.items() if k != report_key} or
                    any(saved[report_key][k] != product[report_key][k] for k in ('bytes', 'sha256'))):
                raise ValueError(f'cannot reuse changed imported provenance: {second}')
            for output in (first, second):
                # Protein outputs were absolute; expression outputs were relative.
                reusable.append((kind, str(output) if kind.endswith('protein') else str(output.relative_to(root))))
    mapping = out / 'orthogroups/mapping/snapshot.json'
    if mapping.exists():
        current = load_tables(mapping, verify_files=False)
        original = load_tables(completed['mapping']['path'], verify_files=False)
        selected = {r['species'] for r in read_tsv(out / 'metadata/samples.tsv')}
        if set(current['tables']) != selected:
            raise ValueError('cannot reuse mapping with a changed selection')
        for species, entry in current['tables'].items():
            source = original['tables'][species]
            if (any(entry[k] != source[k] for k in ('odb_species', 'protein_sha256')) or
                    any(entry['table'][k] != source['table'][k] for k in ('bytes', 'sha256'))):
                raise ValueError(f'cannot reuse changed mapping: {species}')
            verify(relative_file(mapping.parent, entry['table']))
        reusable.extend(('import_build_mapping', str(p.relative_to(root)))
                        for p in (mapping, mapping.parent / 'species'))
    return reusable


def upgrade(path, apply=False):
    path = Path(path).resolve()
    manifest = load(path)
    root = Path(manifest['root'])
    if (path / 'completed.json').exists():
        raise ValueError('completed analyses are immutable')
    if (root / 'VERSION').read_text().strip() != '0.4.15':
        raise ValueError('this audited upgrade requires v0.4.15')
    subprocess.run(['git', 'diff', '--quiet', 'HEAD'], cwd=root, check=True)
    new = implementation(root)
    if {x['path']: x['sha256'] for x in new} == {x['path']: x['sha256'] for x in manifest['implementation']}:
        return {'state': 'already_upgraded'}
    previous = {}
    for entry in manifest['implementation']:
        relative = str(Path(entry['path']).relative_to(root))
        old = subprocess.check_output(['git', 'show', f'{PREDECESSOR}:{relative}'], cwd=root)
        if hashlib.sha256(old).hexdigest() != entry['sha256']:
            raise ValueError(f'frozen code does not match audited v0.4.14: {relative}')
        previous[relative] = entry['sha256']
    for entry in new:
        relative = str(Path(entry['path']).relative_to(root))
        if relative not in previous:
            if relative not in ADDED:
                raise ValueError(f'unaudited new source: {relative}')
        elif previous[relative] != entry['sha256'] and relative not in CHANGED:
            raise ValueError(f'unaudited source change: {relative}')
    if set(previous) - {str(Path(x['path']).relative_to(root)) for x in new}:
        raise ValueError('audited upgrade cannot remove source files')
    completed = load_complete(manifest['build']['path'])
    reuse = verified_imports(path, manifest, completed)
    result = {'state': 'ready', 'verified_import_outputs': len(reuse), 'from_commit': PREDECESSOR}
    if not apply:
        return result
    from snakemake.persistence import Persistence
    before = record(path / 'analysis.json')['sha256']
    with locked(path / '.run.lock'), locked(path / 'jobs/.submit.lock'):
        if record(path / 'analysis.json')['sha256'] != before:
            raise ValueError('analysis changed during upgrade verification')
        queue = subprocess.check_output(['squeue', '-h', '-u', str(__import__('os').getuid()),
                                         '-o', '%i|%Z'], text=True)
        if any(line.strip().endswith('|' + str(root)) for line in queue.splitlines()):
            raise ValueError('repository still has queued/running jobs; inspect them before upgrading')
        backup = path / 'jobs/code_migration_0003'
        backup.mkdir(exist_ok=False)
        for name in ('analysis.json', 'checksums.json', 'run.json'):
            if (path / name).exists():
                shutil.copy2(path / name, backup / name)
        persistence = Persistence(nolock=True, path=root / '.snakemake')
        metadata = []
        for rule, output in reuse:
            saved = persistence.metadata(output)
            if saved and saved.get('rule') == rule:
                metadata.append({'output': output, 'metadata': saved})
        write_json(backup / 'import_metadata.json', metadata)
        create_index(manifest['build']['path'], path / 'build_inputs', completed,
                     {r['species'] for r in manifest['selection'] if r['selected']})
        binary = Path(sys.executable).with_name('snakemake')
        if not binary.is_file():
            raise ValueError('run the upgrade with the Snakemake host environment Python')
        launcher = backup / 'snakemake'
        launcher.write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec ' + shlex.quote(str(binary)) +
                            ' --rerun-triggers mtime input params code "$@"\n')
        launcher.chmod(0o755)
        for item in metadata:
            persistence.cleanup_metadata(item['output'])
            persistence.cleanup_metadata(item['output'])
        manifest['implementation'] = new
        manifest['compatible_launcher'] = record(launcher)
        manifest['files'].extend(record(p) for p in sorted((path / 'build_inputs').rglob('*.json')))
        event = {'created_at': now(), 'from_commit': PREDECESSOR,
                 'to_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
                 'backup': str(backup), 'verified_import_outputs': len(reuse),
                 'reason': 'Batch scheduling and indexed imports; scientific algorithms and environments unchanged.'}
        manifest.setdefault('code_migrations', []).append(event)
        write_json(path / 'analysis.json', manifest)
        write_json(path / 'checksums.json', {'analysis.json': record(path / 'analysis.json')['sha256']})
        result.update(state='upgraded', backup=str(backup), reset_import_metadata=len(metadata))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis', required=True)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    with verification_cache(root / '.cache/verification.sqlite'):
        print(json.dumps(upgrade(args.analysis, args.apply), indent=2))
