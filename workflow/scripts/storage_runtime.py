"""Common node scratch, log compaction and storage measurements for every backend."""
import argparse
from contextlib import contextmanager
import csv
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import shlex
import sys
import tempfile
import uuid

import yaml

from common import mount_csv, now, write_json
from storage_meter import measure, monitored
from scratch_reservation import ScratchAdmissionError, headroom, reserve


DEFAULTS = {'keep_intermediates': False, 'scratch_root': 'auto', 'compact_logs': True,
            'max_bytes': None, 'max_inodes': None, 'interval_seconds': 60,
            'scratch_min_free_bytes': 0, 'scratch_min_free_inodes': 0,
            'scratch_max_bytes': None, 'scratch_max_inodes': None,
            'scratch_reserve_bytes': 0, 'scratch_reserve_inodes': 0,
            'allow_shared_scratch': True, 'start_fraction': .85, 'stop_fraction': .95,
            'series_max_bytes': 1048576, 'task_limits': {}, 'projection_factor': 1.3}
TASK_KEYS = {'scratch_min_free_bytes', 'scratch_min_free_inodes', 'scratch_max_bytes',
             'scratch_max_inodes', 'scratch_reserve_bytes', 'scratch_reserve_inodes'}


def settings(storage):
    unknown = set(storage) - DEFAULTS.keys()
    if unknown:
        raise ValueError(f'unknown storage settings: {sorted(unknown)}')
    result = dict(DEFAULTS, **storage)
    for key in ('keep_intermediates', 'compact_logs', 'allow_shared_scratch'):
        if type(result[key]) is not bool:
            raise ValueError(f'storage.{key} must be true or false')
    for key in ('max_bytes', 'max_inodes', 'scratch_max_bytes', 'scratch_max_inodes'):
        if result[key] is not None and (type(result[key]) is not int or result[key] < 1):
            raise ValueError(f'storage.{key} must be null or a positive integer')
    for key in ('interval_seconds', 'scratch_min_free_bytes', 'scratch_min_free_inodes',
                'scratch_reserve_bytes', 'scratch_reserve_inodes', 'series_max_bytes'):
        if type(result[key]) is not int or result[key] < (1 if key == 'interval_seconds' else 0):
            raise ValueError(f'invalid storage.{key}')
    if result['series_max_bytes'] < 1024:
        raise ValueError('storage.series_max_bytes must be at least 1024')
    if not all(type(result[k]) in (int, float) for k in ('start_fraction', 'stop_fraction')) or not (
            0 < result['start_fraction'] < result['stop_fraction'] <= 1):
        raise ValueError('storage requires 0 < start_fraction < stop_fraction <= 1')
    for unit in ('bytes', 'inodes'):
        reservation, ceiling = result['scratch_reserve_' + unit], result['scratch_max_' + unit]
        if reservation:
            if ceiling is None:
                result['scratch_max_' + unit] = reservation
            elif ceiling > reservation:
                raise ValueError('scratch reservation must cover its scratch ceiling')
    profiles = result['task_limits']
    if not isinstance(profiles, dict):
        raise ValueError('storage.task_limits must be a mapping')
    import re
    for name, values in profiles.items():
        if not isinstance(name, str) or not re.fullmatch(r'[A-Za-z0-9_]+', name) or not isinstance(values, dict) or set(values) - TASK_KEYS:
            raise ValueError('invalid storage.task_limits entry')
        settings(dict(result, task_limits={}, **values))
    if type(result['projection_factor']) not in (int, float) or not 1 <= result['projection_factor'] <= 10:
        raise ValueError('storage.projection_factor must be between 1 and 10')
    root = result['scratch_root']
    if not isinstance(root, str) or (root not in {'auto', 'env', 'workspace'} and not Path(root).is_absolute()) or any(c in root for c in ':,\n\r'):
        raise ValueError('storage.scratch_root must be auto, env, workspace or an absolute path without bind delimiters')
    return result


def task_settings(storage, label, command):
    """AGE rule names and Slurm's pinned Snakemake --target-jobs select profiles."""
    config = settings(storage)
    names = {label}
    if label == 'workflow':
        words = list(command)
        if len(words) >= 3 and words[1] == '-c':
            words = shlex.split(words[2])
        for index, word in enumerate(words):
            if word == '--target-jobs':
                for target in words[index + 1:]:
                    if target.startswith('-'):
                        break
                    names.add(target.split(':', 1)[0])
            elif word.startswith('--target-jobs='):
                names.add(word.split('=', 1)[1].split(':', 1)[0])
    matches = [config['task_limits'][name] for name in names if name in config['task_limits']]
    for key in TASK_KEYS:
        values = [profile[key] for profile in matches if key in profile]
        if values:
            config[key] = None if any(v is None for v in values) else max(values)
    return settings(dict(config, task_limits={}))


def storage_policy(storage, root):
    """Operational retry policy changes storage only, never frozen science inputs."""
    filename = os.environ.get('PHENORADAR_STORAGE_POLICY')
    if not filename:
        return storage, None
    path = (Path(root) / filename).resolve()
    raw = path.read_bytes()
    document = yaml.safe_load(raw)
    if not isinstance(document, dict) or set(document) != {'storage'} or not isinstance(document['storage'], dict):
        raise ValueError('Storage retry policy must contain only a storage mapping')
    result = copy.deepcopy(storage)
    for key, value in document['storage'].items():
        if key == 'task_limits' and isinstance(value, dict):
            profiles = result.setdefault('task_limits', {})
            for label, values in value.items():
                if not isinstance(values, dict):
                    raise ValueError('storage.task_limits entries must be mappings')
                profiles[label] = dict(profiles.get(label, {}), **values)
        else:
            result[key] = value
    settings(result)
    return result, {'path': str(path), 'sha256': hashlib.sha256(raw).hexdigest()}


@contextmanager
def task_scratch(root, config):
    root = Path(root).resolve()
    requested = config['scratch_root']
    if config['keep_intermediates']:
        requested = 'workspace'
    inherited = os.environ.get('PHENORADAR_TASK_SCRATCH')
    candidates = ([Path('/data1'), Path(os.environ.get('TMPDIR', '/tmp')), Path('/tmp')]
                  if requested == 'auto' else [Path(os.environ.get('TMPDIR', ''))] if requested == 'env'
                  else [] if requested == 'workspace' else [Path(requested)])
    failures = []
    for candidate in candidates:
        if not candidate.is_absolute() or not candidate.is_dir():
            continue
        path = candidate.resolve()
        if inherited and path.is_relative_to(Path(inherited)):
            continue  # Submission-host scratch is never an execution-node path.
        if path.is_relative_to(root) or (requested == 'auto' and path.stat().st_dev == root.stat().st_dev):
            continue
        if os.access(path, os.W_OK | os.X_OK):
            failures.append(path)
    if requested == 'workspace' or requested == 'auto' and config['allow_shared_scratch']:
        parent = root / 'work/scratch'
        parent.mkdir(parents=True, exist_ok=True)
        failures.append(parent)
    errors = []
    for parent in dict.fromkeys(failures):
        try:
            headroom(parent, config)
            local = Path(tempfile.mkdtemp(prefix='phenoradar-task-', dir=parent))
        except (OSError, ScratchAdmissionError) as error:
            errors.append(f'{parent}: {error}')
            continue
        entered = False
        try:
            with reserve(local, config):
                entered = True
                write_json(local / 'owner.json', {'root': str(root), 'job_id': os.environ.get('SLURM_JOB_ID') or os.environ.get('JOB_ID'), 'at': now()})
                yield local
            return
        except ScratchAdmissionError as error:
            if entered:
                raise
            errors.append(str(error))
        finally:
            if not entered or not config['keep_intermediates']:
                shutil.rmtree(local)
    raise ScratchAdmissionError('No scratch location satisfies storage policy: ' + (
        '; '.join(errors) or 'no eligible node-local path; shared fallback is disabled or not selected'))


def environment(local, config):
    env = dict(os.environ)
    old = env.get('PHENORADAR_TASK_SCRATCH')
    if os.environ.get('PHENORADAR_SCHEDULER') == 'slurm' and os.environ.get('SLURM_JOB_ID'):
        env.pop('PHENORADAR_TASK_CPUS', None)
        env.pop('PHENORADAR_TASK_MEM_MB', None)
    env.update(PHENORADAR_TASK_SCRATCH=str(local), TMPDIR=str(local / 'tmp'),
               PHENORADAR_STORAGE_ACTIVE='1',
               PYTHONPYCACHEPREFIX=str(local / 'pycache'),
               PHENORADAR_COMPACT_ODB_LOGS='1' if config['compact_logs'] else '0')
    for name in ('tmp', 'pycache'):
        (local / name).mkdir()
    slots = env.get('SLURM_CPUS_PER_TASK') or env.get('PHENORADAR_TASK_CPUS') or env.get('NSLOTS')
    if slots is not None:
        cpus = int(slots)
        if cpus < 1:
            raise ValueError('Invalid task CPU allocation')
        for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
            env[key] = str(cpus)
            for prefix in ('APPTAINERENV_', 'SINGULARITYENV_'):
                env[prefix + key] = str(cpus)
    # --cleanenv still honors APPTAINERENV_*. Mounts are evaluated by the host
    # runtime, so a different node can choose its own path after scheduling.
    rows = []
    for variable in ('APPTAINER_MOUNT', 'SINGULARITY_MOUNT'):
        for line in env.get(variable, '').splitlines():
            if not line.strip():
                continue
            fields = next(csv.reader([line]))
            options = dict(field.split('=', 1) for field in fields if '=' in field)
            if old and options.get('src') == old:
                continue
            if fields not in rows:
                rows.append(fields)
    rows.append(['type=bind', 'src=' + str(local), 'dst=' + str(local)])
    env['APPTAINER_MOUNT'] = mount_csv(rows)
    env['SINGULARITY_MOUNT'] = env['APPTAINER_MOUNT']
    for prefix in ('APPTAINER', 'SINGULARITY'):
        env[prefix + '_TMPDIR'] = str(local / 'tmp')
        for key in ('PHENORADAR_TASK_SCRATCH', 'PHENORADAR_COMPACT_ODB_LOGS'):
            env[prefix + 'ENV_' + key] = env[key]
    return env


def report(build, expected_samples=None):
    """Measure all retained work, with a verified-completion gate on projection."""
    from dataset import load
    build = Path(build).resolve()
    manifest = load(build)
    root = Path(manifest['root'])
    config = settings(manifest['config'].get('storage', {}))
    count = len(manifest['items'])
    complete = (build / 'completed.json').is_file()
    if complete:
        from build_products import load_complete
        load_complete(build)
    retained, project, database = [measure([p]) for p in (build, root, build / 'database')]
    metrics = [json.loads(p.read_text()) for p in (build / 'jobs/storage').glob('*.json')]
    peaks = {label: {key: max([0, *[m.get(label + '_peak', {}).get(key, 0) for m in metrics]])
                     for key in ('allocated_bytes', 'inodes')} for label in ('shared', 'scratch', 'sample')}
    errors = retained['errors'] + project['errors'] + database['errors']
    if expected_samples is not None and (type(expected_samples) is not int or expected_samples < count):
        raise ValueError('expected_samples must be at least the measured build sample count')
    projection = None
    if expected_samples and complete and count and not errors:
        projection = {}
        variable_factor = expected_samples / count * config['projection_factor']
        # Retained work includes logs, raw ODB hits/annotations and all caches.
        for key, limit in (('allocated_bytes', 'max_bytes'), ('inodes', 'max_inodes')):
            estimate = max(0, project[key] - retained[key]) + math.ceil(retained[key] * variable_factor)
            # Include sampled transient shared-space growth as well as final size.
            transient = max(0, peaks['shared'][key] - project[key])
            estimate += math.ceil(transient * config['projection_factor'])
            projection[key] = {'estimate': estimate, 'budget': config[limit],
                               'within_budget': estimate < config[limit] if config[limit] else None}
    data = {'at': now(), 'build': str(build), 'samples': count, 'complete': complete,
            'retained_build': retained, 'database': database, 'project': project,
            'sampled_peaks': peaks, 'expected_samples': expected_samples, 'projection': projection,
            'attempts': len(metrics), 'errors': errors,
            'limitations': ['Periodic inventories can miss short peaks.',
                           'Project scope excludes other users and external inputs/caches.',
                           'A small pilot and safety factor do not guarantee full-build usage.',
                           'Scratch peaks depend on chunk size and worker concurrency.']}
    write_json(build / 'reports/storage.json', data)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['run'])
    parser.add_argument('--pipeline', required=True)
    parser.add_argument('--root', required=True)
    parser.add_argument('--usage', required=True)
    parser.add_argument('--label', required=True)
    parser.add_argument('--sample-map')
    parser.add_argument('--task-variable')
    argv = sys.argv[1:]
    separator = argv.index('--') if '--' in argv else len(argv)
    args = parser.parse_args(argv[:separator])
    if os.environ.get('PE_HOSTFILE'):
        hosts = {line.split()[0] for line in Path(os.environ['PE_HOSTFILE']).read_text().splitlines() if line.strip()}
        if len(hosts) > 1:
            raise ValueError('This workflow requires all requested CPU slots on one node')
    command = argv[separator + 1:]
    if not command:
        parser.error('a command is required after --')
    storage, policy = storage_policy(yaml.safe_load(Path(args.pipeline).read_text()).get('storage', {}), args.root)
    config = task_settings(storage, args.label, command)
    config['policy_source'] = policy
    sample_roots = []
    if args.sample_map:
        samples = json.loads(Path(args.sample_map).read_text())
        sample_roots = [samples[os.environ[args.task_variable]]]
    receipt = Path(args.usage) / (args.label + '_' + uuid.uuid4().hex + '.json')
    try:
        with task_scratch(args.root, config) as local:
            monitored(command, cwd=args.root, env=environment(local, config), shared_root=args.root,
                      scratch=local, sample_roots=sample_roots, receipt=receipt, settings=config)
    except ScratchAdmissionError as error:
        if not receipt.exists():
            write_json(receipt, {'schema_version': 1, 'state': 'blocked', 'reason': 'scratch_admission',
                                 'child_started': False, 'error': str(error), 'finished_at': now(),
                                 'storage_limits': config, 'label': args.label})
        raise


if __name__ == '__main__':
    main()
