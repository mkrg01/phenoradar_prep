"""Audit reusable outputs and explicitly migrate unfinished analysis conditions."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import shlex
import sys
import re
from types import SimpleNamespace

import yaml

from analysis import load, submission_guard
from analysis_build_index import create_index
from build_products import load_complete
from common import content_sha256, now, read_tsv, write_json
from dataset import implementation
from dataset_assets import locked, record, stat_identity, verify
from mapping_tables import load_tables, relative_file
from product_storage import decoded_sha256
from verification_cache import active_cache, verification_cache

PREDECESSOR = 'dddd331b5a2770bbe7709370cf3695834cc261eb'
BATCHED = 'e7a0e1343f909666df140ea26c8d9f50fcdcaa0e'
OCCUPANCY = '2729d40481974a1ffa1db73e9d6982d3f279d342'
CONTAINER = '35c93ecbe3e43e48c4cac142fff39174847056fd'
COMPRESSED_IMPORTS = '641fcfd0cc1bb7a44ff2fa25a087ea3f714e43ee'
CHANGED = {'VERSION', 'run_pipeline.sh', 'workflow/rules/alignments.smk',
           'workflow/rules/build_inputs.smk', 'workflow/rules/tpm.smk',
           'workflow/scripts/analysis.py', 'workflow/scripts/build_products.py',
           'workflow/scripts/dataset.py', 'workflow/scripts/phase_config.py',
           'workflow/scripts/configuration.py', 'workflow/pipeline_defaults.yaml',
           'workflow/scripts/filter_species.py', 'workflow/scripts/upgrade_analysis.py',
           'workflow/container.smk', 'workflow/scripts/dataset_software.py'}
ADDED = {'workflow/scripts/alignment_batches.py', 'workflow/scripts/analysis_build_index.py',
         'workflow/scripts/slurm_submit.py', 'workflow/scripts/upgrade_analysis.py', 'workflow/slurm/sbatch',
         'workflow/scripts/alignment_selection.py', 'workflow/scripts/container_pull.py'}
BATCHED_CHANGES = {'VERSION', 'workflow/rules/alignments.smk', 'workflow/scripts/configuration.py',
                   'workflow/pipeline_defaults.yaml', 'workflow/scripts/filter_species.py',
                   'workflow/scripts/upgrade_analysis.py', 'workflow/container.smk',
                   'workflow/scripts/dataset_software.py'}
CONTAINER_CHANGES = {'VERSION', 'workflow/container.smk', 'workflow/scripts/dataset_software.py',
                     'workflow/scripts/upgrade_analysis.py'}
SCHEDULER_PREDECESSOR = '092e4a409fea7dc008f1629a8e7870a8ad43c873'
SCHEDULER_REPAIR = 'd2fa8dfc30901007ce5fb5460e9152c14b92c4e0'
SCHEDULER_CHANGES = {'VERSION', 'workflow/Snakefile', 'workflow/scripts/phase_config.py',
                     'workflow/scripts/slurm_submit.py', 'workflow/scripts/versioning.py',
                     'workflow/scripts/analysis.py', 'workflow/scripts/analysis_inspection.py',
                     'workflow/scripts/upgrade_analysis.py'}


def upgrade_scheduler(path, manifest, apply=False, min_sample_fraction=None):
    """Migrate audited controller fixes while preserving the saved tool image."""
    root = Path(manifest['root'])
    before, after = occupancy_condition(manifest, min_sample_fraction)
    if before != after:
        raise ValueError('scheduler-only migration cannot change analysis conditions')
    subprocess.run(['git', 'diff', '--quiet', 'HEAD'], cwd=root, check=True)
    current = implementation(root)
    previous = {str(Path(x['path']).relative_to(root)): x['sha256'] for x in manifest['implementation']}
    observed = {str(Path(x['path']).relative_to(root)): x['sha256'] for x in current}
    if observed == previous:
        return {'state': 'already_upgraded'}
    if set(previous) != set(observed):
        raise ValueError('scheduler-only migration cannot add or remove workflow sources')
    candidates = [revision for revision in (SCHEDULER_PREDECESSOR, SCHEDULER_REPAIR)
                  if hashlib.sha256(subprocess.check_output(['git', 'show', f'{revision}:VERSION'],
                                                            cwd=root)).hexdigest() == previous.get('VERSION')]
    if len(candidates) != 1:
        raise ValueError('frozen code must match audited v0.5.3 or v0.5.4')
    revision = candidates[0]
    source_version = subprocess.check_output(['git', 'show', f'{revision}:VERSION'], cwd=root, text=True).strip()
    for relative, checksum in previous.items():
        saved = subprocess.check_output(['git', 'show', f'{revision}:{relative}'], cwd=root)
        if hashlib.sha256(saved).hexdigest() != checksum:
            raise ValueError(f'frozen code differs from audited v{source_version}: {relative}')
        if observed[relative] != checksum and relative not in SCHEDULER_CHANGES:
            raise ValueError(f'scheduler migration cannot change scientific code or environments: {relative}')
    submission_guard(path, manifest)
    result = {'state': 'ready', 'from_commit': revision,
              'changed_sources': sorted(p for p in observed if observed[p] != previous[p]),
              'reset_metadata': 0}
    if not apply:
        return result
    before_manifest = record(path / 'analysis.json')['sha256']
    with locked(path / '.run.lock'), locked(path / 'jobs/.submit.lock'):
        if record(path / 'analysis.json')['sha256'] != before_manifest:
            raise ValueError('analysis changed during scheduler migration')
        submission_guard(path, manifest)
        numbers = [int(p.name.rsplit('_', 1)[1]) for p in (path / 'jobs').glob('code_migration_*')
                   if re.fullmatch(r'code_migration_[0-9]+', p.name)]
        backup = path / f'jobs/code_migration_{max(numbers, default=0) + 1:04d}'
        backup.mkdir(exist_ok=False)
        for name in ('analysis.json', 'checksums.json', 'run.json', 'pipeline.yaml', 'analysis.yaml'):
            if (path / name).exists():
                shutil.copy2(path / name, backup / name)
        from versioning import IMAGE_REPOSITORY
        image = manifest.get('compatible_container_image', f'docker://{IMAGE_REPOSITORY}:v{source_version}')
        manifest['implementation'] = current
        manifest['compatible_container_image'] = image
        event = {'created_at': now(), 'from_commit': revision,
                 'to_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
                 'backup': str(backup), 'container_image': image,
                 'changed_sources': result['changed_sources'],
                 'reason': 'Scheduler-only repair; scientific conditions, tool image, '
                           'outputs and Snakemake output metadata are retained.'}
        manifest.setdefault('code_migrations', []).append(event)
        write_json(path / 'analysis.json', manifest)
        write_json(path / 'checksums.json', {'analysis.json': record(path / 'analysis.json')['sha256']})
        return dict(result, state='upgraded', backup=str(backup), container_image=image)


def audit_source(root, saved, current):
    previous = {str(Path(x['path']).relative_to(root)): x['sha256'] for x in saved}
    candidates = [revision for revision in (PREDECESSOR, BATCHED, OCCUPANCY, CONTAINER, COMPRESSED_IMPORTS)
                  if hashlib.sha256(subprocess.check_output(
                      ['git', 'show', f'{revision}:VERSION'], cwd=root)).hexdigest() == previous.get('VERSION')]
    if len(candidates) != 1:
        raise ValueError('frozen code must match audited v0.4.14 through v0.4.18')
    revision = candidates[0]
    allowed = {PREDECESSOR: CHANGED, BATCHED: BATCHED_CHANGES, OCCUPANCY: CONTAINER_CHANGES,
               CONTAINER: {'VERSION', 'workflow/scripts/upgrade_analysis.py'},
               COMPRESSED_IMPORTS: {'VERSION', 'workflow/scripts/upgrade_analysis.py'}}[revision]
    for relative, checksum in previous.items():
        old = subprocess.check_output(['git', 'show', f'{revision}:{relative}'], cwd=root)
        if hashlib.sha256(old).hexdigest() != checksum:
            raise ValueError(f'frozen code differs from audited predecessor: {relative}')
    for entry in current:
        relative = str(Path(entry['path']).relative_to(root))
        if relative not in previous:
            if relative not in ADDED:
                raise ValueError(f'unaudited new source: {relative}')
        elif previous[relative] != entry['sha256'] and relative not in allowed:
            raise ValueError(f'unaudited source change: {relative}')
    if set(previous) - {str(Path(x['path']).relative_to(root)) for x in current}:
        raise ValueError('audited upgrade cannot remove source files')
    return revision


def occupancy_condition(manifest, requested):
    """Omitting the new option preserves the legacy all-OG condition explicitly."""
    from alignment_selection import validate_fraction
    before = manifest['pipeline']['alignment'].get('min_sample_fraction', 0)
    after = before if requested is None else validate_fraction(requested)
    return before, after


def verified_import_product(source, target, saved):
    """Check exact copies by stored bytes and converted imports by decoded bytes."""
    if (Path(source['path']).suffix == '.gz') == (target.suffix == '.gz'):
        verify(dict(source, path=str(target)))
        return source
    # Importers deliberately gzip historical plain products. Bind the saved
    # receipt to the actual file before trusting its decoded-content identity.
    before = stat_identity(target)
    verify(dict(saved, path=str(target)))
    cache = active_cache()
    indexed = cache.document(target, 'import-decoded-v1', before, saved['sha256']) if cache else None
    observed = indexed['sha256'] if indexed else decoded_sha256(target)
    if before != stat_identity(target):
        raise ValueError(f'import changed while verifying content: {target}')
    expected = content_sha256(source)
    if str(source['path']).endswith('.gz') and 'content_sha256' not in source:
        expected = decoded_sha256(verify(source))
    if observed != expected or ('content_sha256' in saved and saved['content_sha256'] != observed):
        raise ValueError(f'cannot reuse changed imported content: {target}')
    if cache is not None and indexed is None:
        cache.remember_document(target, 'import-decoded-v1', before, saved['sha256'], {'sha256': observed})
    return dict(saved, bytes=before[2])


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
            original = json.loads(verify(product[source_key]).read_text())
            saved = json.loads(second.read_text())
            imported = verified_import_product(product[report_key], first, saved[report_key])
            if ({k: v for k, v in original.items() if k != report_key} !=
                    {k: v for k, v in saved.items() if k != report_key} or
                    any(saved[report_key][k] != imported[k] for k in ('bytes', 'sha256'))):
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


def upgrade(path, apply=False, min_sample_fraction=None):
    path = Path(path).resolve()
    manifest = load(path)
    root = Path(manifest['root'])
    if (path / 'completed.json').exists():
        raise ValueError('completed analyses are immutable')
    if (root / 'VERSION').read_text().strip() in {'0.5.4', '0.5.5'}:
        return upgrade_scheduler(path, manifest, apply, min_sample_fraction)
    if (root / 'VERSION').read_text().strip() != '0.4.19':
        raise ValueError('this audited upgrade requires v0.4.19')
    subprocess.run(['git', 'diff', '--quiet', 'HEAD'], cwd=root, check=True)
    new = implementation(root)
    before_fraction, after_fraction = occupancy_condition(manifest, min_sample_fraction)
    same_code = {x['path']: x['sha256'] for x in new} == {
        x['path']: x['sha256'] for x in manifest['implementation']}
    if same_code and before_fraction == after_fraction:
        return {'state': 'already_upgraded'}
    revision = (subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip()
                if same_code else audit_source(root, manifest['implementation'], new))
    completed = load_complete(manifest['build']['path'])
    reuse = verified_imports(path, manifest, completed)
    result = {'state': 'ready', 'verified_import_outputs': len(reuse), 'from_commit': revision,
              'alignment.min_sample_fraction': {'before': before_fraction, 'after': after_fraction}}
    if not apply:
        return result
    from snakemake.persistence import Persistence
    from snakemake.io import IOFile
    binary = Path(sys.executable).with_name('snakemake')
    if not binary.is_file():
        raise ValueError('run the upgrade with the Snakemake host environment Python')
    before = record(path / 'analysis.json')['sha256']
    with locked(path / '.run.lock'), locked(path / 'jobs/.submit.lock'):
        if record(path / 'analysis.json')['sha256'] != before:
            raise ValueError('analysis changed during upgrade verification')
        queue = subprocess.check_output(['squeue', '-h', '-u', str(__import__('os').getuid()),
                                         '-o', '%i|%Z'], text=True)
        current_job = __import__('os').environ.get('SLURM_JOB_ID')
        if any(line.strip().endswith('|' + str(root)) and line.split('|')[0] != current_job
               for line in queue.splitlines()):
            raise ValueError('repository still has queued/running jobs; inspect them before upgrading')
        numbers = [int(p.name.rsplit('_', 1)[1]) for p in (path / 'jobs').glob('code_migration_*')
                   if re.fullmatch(r'code_migration_[0-9]+', p.name)]
        # Snakemake 9.8 requires the checksum setting even for metadata-only
        # operations. No workflow/DAG evaluation or input hashing is needed here.
        metadata_dag = SimpleNamespace(workflow=SimpleNamespace(
            dag_settings=SimpleNamespace(max_checksum_file_size=0)))
        persistence = Persistence(nolock=True, path=root / '.snakemake', dag=metadata_dag)
        metadata = []
        for rule, output in reuse:
            saved = persistence.metadata(IOFile(output))
            if saved and saved.get('rule') == rule:
                metadata.append({'output': output, 'metadata': saved})
        backup = path / f'jobs/code_migration_{max(numbers, default=2) + 1:04d}'
        backup.mkdir(exist_ok=False)
        for name in ('analysis.json', 'checksums.json', 'run.json', 'pipeline.yaml', 'analysis.yaml'):
            if (path / name).exists():
                shutil.copy2(path / name, backup / name)
        write_json(backup / 'import_metadata.json', metadata)
        if not (path / 'build_inputs/index.json').exists():
            create_index(manifest['build']['path'], path / 'build_inputs', completed,
                         {r['species'] for r in manifest['selection'] if r['selected']})
        launcher = backup / 'snakemake'
        launcher.write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec ' + shlex.quote(str(binary)) +
                            ' --rerun-triggers mtime input params code "$@"\n')
        launcher.chmod(0o755)
        for item in metadata:
            # cleanup_metadata short-circuits after removing an incomplete
            # record; the second call also clears the corresponding metadata.
            output = IOFile(item['output'])
            persistence.cleanup_metadata(output)
            persistence.cleanup_metadata(output)
        manifest['implementation'] = new
        manifest['compatible_launcher'] = record(launcher)
        for section in ('config', 'pipeline'):
            manifest[section]['alignment']['min_sample_fraction'] = after_fraction
        (path / 'pipeline.yaml').write_text(yaml.safe_dump(manifest['pipeline'], sort_keys=False))
        (path / 'analysis.yaml').write_text(yaml.safe_dump(manifest['config'], sort_keys=False))
        files = {entry['path']: entry for entry in manifest['files']}
        for p in [path / 'pipeline.yaml', path / 'analysis.yaml', *sorted((path / 'build_inputs').rglob('*.json'))]:
            entry = record(p)
            files[entry['path']] = entry
        manifest['files'] = list(files.values())
        event = {'created_at': now(), 'from_commit': revision,
                 'to_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
                 'backup': str(backup), 'verified_import_outputs': len(reuse),
                 'analysis_conditions': result['alignment.min_sample_fraction'],
                 'reason': 'Explicit OG sample-occupancy selection; retained OG alignment algorithms, '
                           'tool environments and other analysis conditions unchanged.'}
        manifest.setdefault('code_migrations', []).append(event)
        write_json(path / 'analysis.json', manifest)
        write_json(path / 'checksums.json', {'analysis.json': record(path / 'analysis.json')['sha256']})
        result.update(state='upgraded', backup=str(backup), reset_import_metadata=len(metadata))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis', required=True)
    parser.add_argument('--apply', action='store_true')
    parser.add_argument('--min-sample-fraction', type=float,
                        help='explicitly change OG sample occupancy; omitted preserves the saved condition')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    with verification_cache(root / '.cache/verification.sqlite'):
        print(json.dumps(upgrade(args.analysis, args.apply, args.min_sample_fraction), indent=2))
