"""Prepare an immutable BUSCO download once and mount it into native sample jobs."""
import csv
import fcntl
from functools import lru_cache
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from common import now, write_json
from dataset_assets import ProductIntegrityError, UncachedVerification, digest, record, verify


DESTINATION = '/usr/local/db/busco_downloads'


def settings(root, config):
    if set(config) - {'lineage', 'cache_dir'}:
        raise ValueError('unknown build.busco settings')
    if not isinstance(config.get('lineage'), str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', config['lineage']):
        raise ValueError('busco.lineage must be a dataset name')
    cache = config.get('cache_dir', 'resources/busco')
    if cache is not None:
        if not isinstance(cache, str) or not cache.strip() or '\n' in cache or '\r' in cache:
            raise ValueError('busco.cache_dir must be null or a directory path without newlines')
        cache = Path(root) / cache
        config['cache_dir'] = str(cache.resolve())
    else:
        config['cache_dir'] = None


@lru_cache(maxsize=1)
def runtime():
    # Match the normal GeneGalleon preference. A singularity compatibility binary
    # supplied by Apptainer is supported too; SingularityCE and Docker are not.
    if (os.environ.get('GG_CONTAINER_RUNTIME', 'auto') != 'auto'
            or os.environ.get('GG_CONTAINER_DOCKER_IMAGE') or os.environ.get('GG_WRAPPER_IMAGE')):
        raise ValueError('shared references require the Apptainer runtime; set busco.cache_dir and taxonomy.cache_dir to null to disable')
    command = shutil.which('singularity') or shutil.which('apptainer')
    if not command:
        raise ValueError('shared references require Apptainer')
    version = subprocess.check_output([command, '--version'], text=True).strip()
    help_text = subprocess.check_output([command, 'exec', '--help'], text=True)
    if 'apptainer' not in version.lower() or '--mount' not in help_text:
        raise ValueError('shared references require Apptainer with APPTAINER_MOUNT support; '
                         'set busco.cache_dir and taxonomy.cache_dir to null to disable')
    return command


def tree_files(directory, label='BUSCO'):
    if Path(directory).is_symlink():
        raise ProductIntegrityError(f'{label} reference must not be a symlink: {directory}')
    files = []
    for path in sorted(Path(directory).rglob('*')):
        if path.is_symlink():
            raise ProductIntegrityError(f'{label} reference must not contain symlinks: {path}')
        if path.is_file():
            files.append(path)
    return files


def read_reference(root, lineage, expected=None):
    root = Path(root)
    marker = root / 'reference.json'
    if expected is not None:
        if Path(expected['receipt']['path']) != marker:
            raise ProductIntegrityError('BUSCO reference path differs from the frozen snapshot')
        verify(expected['receipt'])
    if not marker.is_file():
        raise ProductIntegrityError(f'incomplete BUSCO reference: {root}; use a new cache directory')
    data = json.loads(marker.read_text())
    if data.get('schema_version') != 1 or data.get('lineage') != lineage:
        raise ProductIntegrityError('BUSCO reference identity differs from the requested lineage')
    inventory = data['files']
    if not inventory or digest(inventory) != data.get('sha256'):
        raise ProductIntegrityError('BUSCO reference inventory changed')
    if expected is not None and data['sha256'] != expected['sha256']:
        raise ProductIntegrityError('BUSCO reference differs from the frozen snapshot')
    directory = root / 'dataset'
    paths = tree_files(directory)
    if [str(p.relative_to(directory)) for p in paths] != [entry['path'] for entry in inventory]:
        raise ProductIntegrityError('BUSCO reference files are missing or unexpected')
    for entry in inventory:
        verify(dict(entry, path=str(directory / entry['path'])))
    return {'lineage': lineage, 'sha256': data['sha256'], 'dataset': str(directory),
            'dataset_version': data['dataset_version'], 'receipt': record(marker)}


def prepare(config, software, *, expected=None, cached_only=False):
    """Pin the first complete download. Never update or repair a published cache."""
    if not config.get('cache_dir'):
        return None
    lineage = config['lineage']
    root = Path(config['cache_dir']) / lineage
    if expected is not None or root.exists():
        return read_reference(root, lineage, expected)
    if cached_only:
        raise UncachedVerification('uncached BUSCO reference preparation')
    command = runtime()
    root.parent.mkdir(parents=True, exist_ok=True)
    with (root.parent / (lineage + '.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if root.exists():
            return read_reference(root, lineage, expected)
        staging = Path(tempfile.mkdtemp(prefix='.' + lineage + '-', dir=root.parent))
        try:
            env = dict(os.environ)
            for prefix in ('APPTAINER', 'SINGULARITY'):
                for key in ('BIND', 'BINDPATH', 'MOUNT'):
                    env.pop(prefix + '_' + key, None)
                env[prefix + '_TMPDIR'] = str(staging / 'tmp')
            (staging / 'tmp').mkdir()
            # Activate the same base environment as GeneGalleon's transcriptome core.
            script = ('if [ -f /opt/conda/etc/profile.d/conda.sh ]; then '
                      'source /opt/conda/etc/profile.d/conda.sh; conda activate base; fi; '
                      'exec busco --download_path /busco/dataset --download "$1"')
            print(f'BUSCO reference: preparing shared {lineage} under {root}', flush=True)
            mount = io.StringIO()
            csv.writer(mount, lineterminator='').writerow(['type=bind', f'src={staging}', 'dst=/busco'])
            subprocess.run([command, 'exec', '--contain', '--mount', mount.getvalue(),
                            software['image'], 'bash', '-euc', script, 'busco-download', lineage],
                           env=env, check=True)
            dataset = staging / 'dataset'
            lineage_dir = dataset / 'lineages' / lineage
            cfg = lineage_dir / 'dataset.cfg'
            if not cfg.is_file():
                cfg = lineage_dir / 'info/dataset.cfg'
            if not cfg.is_file() or not cfg.stat().st_size or not any((lineage_dir / 'hmms').glob('*.hmm')):
                raise ValueError(f'incomplete BUSCO lineage download: {lineage}')
            inventory = [dict(path=str(p.relative_to(dataset)), bytes=p.stat().st_size, sha256=record(p)['sha256'])
                         for p in tree_files(dataset)]
            version = next((line.split('=', 1)[1].strip() for line in cfg.read_text().splitlines()
                            if '=' in line and line.split('=', 1)[0].strip() == 'creation_date'), 'unknown')
            write_json(staging / 'reference.json', {'schema_version': 1, 'lineage': lineage,
                       'created_at': now(), 'dataset_version': version,
                       'sha256': digest(inventory), 'files': inventory}, compact=True)
            shutil.rmtree(staging / 'tmp')
            staging.rename(root)
            return read_reference(root, lineage)
        finally:
            if staging.exists():
                shutil.rmtree(staging)


def mount_environment(env, reference, destination=DESTINATION, label='BUSCO'):
    """Append a CSV mount without losing unrelated mounts or allowing a shadow."""
    runtime()
    rows = []
    for variable in ('APPTAINER_MOUNT', 'SINGULARITY_MOUNT'):
        for line in env.get(variable, '').splitlines():
            if not line.strip():
                continue
            fields = next(csv.reader([line], strict=True))
            options = dict(field.split('=', 1) for field in fields if '=' in field)
            mounted = next((options[key] for key in ('dst', 'destination', 'dest', 'target') if key in options), None)
            mounted = os.path.normpath(mounted) if mounted else None
            if mounted and (mounted == destination or destination.startswith(mounted.rstrip('/') + '/')
                            or mounted.startswith(destination + '/')):
                raise ValueError(f'existing {variable} mount overlaps the shared {label} reference')
            if fields not in rows:
                rows.append(fields)
    rows.append(['type=bind', 'src=' + reference['dataset'], 'dst=' + destination, 'ro'])
    text = io.StringIO()
    csv.writer(text, lineterminator='\n').writerows(rows)
    env['APPTAINER_MOUNT'] = text.getvalue().rstrip('\n')
    env.pop('SINGULARITY_MOUNT', None)


def verify_mount(reference, software, env):
    """Fail before GeneGalleon can silently fall back to a per-sample download."""
    data = json.loads(verify(reference['receipt']).read_text())
    lineage = reference['lineage']
    candidates = {f'lineages/{lineage}/dataset.cfg', f'lineages/{lineage}/info/dataset.cfg'}
    cfg = next(entry for entry in data['files'] if entry['path'] in candidates)
    script = ('import hashlib, pathlib, sys; p=pathlib.Path(sys.argv[1]); '
              'sys.exit(0 if hashlib.sha256(p.read_bytes()).hexdigest()==sys.argv[2] '
              'else "shared BUSCO mount differs")')
    # GeneGalleon clears these variables before starting its container. Perform
    # the probe with the same binds, rather than ones the native job will discard.
    env = dict(env)
    for prefix in ('APPTAINER', 'SINGULARITY'):
        for key in ('BIND', 'BINDPATH'):
            env.pop(prefix + '_' + key, None)
    subprocess.run([runtime(), 'exec', software['image'], 'python', '-c', script,
                    DESTINATION + '/' + cfg['path'], cfg['sha256']], env=env, check=True)
