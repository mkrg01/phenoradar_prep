"""Pull SIF images with bounded compression and publish only complete images."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

from dataset_assets import locked


def budgets(cpus, memory_mb, env):
    if any(type(n) is not int or n < 1 for n in (cpus, memory_mb)):
        raise ValueError('container CPU and memory budgets must be positive integers')
    memory_bytes = memory_mb * 1000000
    if env.get('SLURM_JOB_ID'):
        allocated_cpus = int(env.get('SLURM_CPUS_PER_TASK', '1'))
        allocated_memory = int(env.get('SLURM_MEM_PER_NODE', '0'))
        if not allocated_memory:
            allocated_memory = int(env.get('SLURM_MEM_PER_CPU', '0')) * allocated_cpus
        if allocated_cpus < 1 or allocated_memory < 1:
            raise ValueError('container preparation requires a finite Slurm CPU/memory allocation')
        cpus = min(cpus, allocated_cpus)
        memory_bytes = min(memory_bytes, allocated_memory * 1048576)
    if memory_bytes < 4000000000:
        raise ValueError('allocate at least 4 GB for container preparation; 16 GB is recommended')
    # Compression buffers use at most one quarter of the job budget, capped at
    # 4 GiB. Leave room for unpacking, the runtime, Python and filesystem cache.
    return cpus, min(4096, memory_bytes // 1048576 // 4)


def memory_mib(value):
    match = re.fullmatch(r'([0-9]+)([KMG]?)', value.strip(), re.I)
    if not match:
        raise ValueError('cannot interpret the runtime mksquashfs mem limit: ' + value)
    number, unit = int(match[1]), match[2].upper()
    return number * {'': 1, 'K': 1 / 1024, 'M': 1, 'G': 1024}[unit]


def compression_config(runtime, env, destination, cpus, memory):
    settings = dict(line.split('=', 1) for line in subprocess.check_output(
        [runtime, 'buildcfg'], text=True, env=env).splitlines() if '=' in line)
    prefix = 'APPTAINER' if 'APPTAINER_CONF_FILE' in settings else 'SINGULARITY'
    source = (env.get(prefix + '_CONFIG_FILE') or env.get('SINGULARITY_CONFIG_FILE')
              or settings[prefix + '_CONF_FILE'])
    text = Path(source).read_text()
    active = dict(re.findall(r'^\s*(mksquashfs (?:procs|mem))\s*=\s*([^#\n]+)', text, re.M))
    old_cpus = int(active.get('mksquashfs procs', '0').strip())
    old_memory = active.get('mksquashfs mem', '').strip()
    # Keep any stricter administrator/user limits.
    cpus = min(cpus, old_cpus) if old_cpus else cpus
    memory = min(memory, int(memory_mib(old_memory))) if old_memory else memory
    if cpus < 1 or memory < 1:
        raise ValueError('invalid runtime compression limits')
    if settings.get(prefix + '_SUID_INSTALL') == '1' and os.geteuid() != 0:
        if not old_cpus or not old_memory or old_cpus > cpus or memory_mib(old_memory) > memory:
            raise ValueError('this setuid runtime requires administrator mksquashfs procs/mem limits; '
                             'use an unprivileged installation for job-specific limits')
        return source, cpus, memory
    text = re.sub(r'^\s*mksquashfs (?:procs|mem)\s*=.*$', '', text, flags=re.M)
    destination.write_text(text.rstrip() + f'\nmksquashfs procs = {cpus}\nmksquashfs mem = {memory}M\n')
    return str(destination), cpus, memory


def inspect(runtime, image, env):
    if not image.is_file() or not image.stat().st_size:
        raise ValueError('container image is missing or empty: ' + str(image))
    return json.loads(subprocess.check_output(
        [runtime, 'inspect', '--json', str(image)], text=True, env=env))


def pull_image(runtime, uri, destination, *, cpus=1, memory_mb=16000, architecture=None, env=None):
    """Use an invocation-local config; never change the installed runtime config."""
    env = dict(os.environ if env is None else env)
    cpus, memory = budgets(cpus, memory_mb, env)
    destination = Path(destination).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with locked(destination.with_name(destination.name + '.lock')):
        if destination.exists():
            inspect(runtime, destination, env)
            return destination
        with tempfile.TemporaryDirectory(prefix='.container-', dir=destination.parent) as temporary:
            staging = Path(temporary)
            config, cpus, memory = compression_config(runtime, env, staging / 'runtime.conf', cpus, memory)
            scratch = staging / 'tmp'
            scratch.mkdir()
            for prefix in ('APPTAINER', 'SINGULARITY'):
                env[prefix + '_CONFIG_FILE'] = config
                env[prefix + '_TMPDIR'] = str(scratch)
            env['TMPDIR'] = str(scratch)
            image = staging / 'image.sif'
            command = [runtime, 'pull']
            if architecture:
                command += ['--arch', architecture]
            command += [str(image), uri]
            print(f'container compression: {cpus} CPU(s), {memory} MiB buffer limit',
                  file=sys.stderr, flush=True)
            subprocess.run(command, env=env, check=True)
            inspect(runtime, image, env)
            image.replace(destination)
    return destination


def prepare(root, cpus=1, memory_mb=16000, cache_dir=None):
    from versioning import resolve_container_image
    root = Path(root).resolve()
    uri = resolve_container_image(root)
    runtime = shutil.which('apptainer') or shutil.which('singularity')
    if not runtime:
        raise ValueError('install Apptainer or Singularity to prepare the workflow container')
    # Match Snakemake 9.8's URL-keyed cache, including tag-based references.
    cache = Path(cache_dir).resolve() if cache_dir else root / '.snakemake/singularity'
    image = cache / (hashlib.md5(uri.encode()).hexdigest() + '.simg')
    pull_image(runtime, uri, image, cpus=cpus, memory_mb=memory_mb)
    subprocess.run([runtime, 'exec', '--cleanenv', str(image), 'conda', 'info', '--json'],
                   stdout=subprocess.DEVNULL, check=True)
    return image


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', default=Path(__file__).resolve().parents[2])
    parser.add_argument('--cpus', type=int, default=1)
    parser.add_argument('--memory-mb', type=int, default=16000)
    parser.add_argument('--cache-dir')
    args = parser.parse_args()
    prepare(args.root, args.cpus, args.memory_mb, args.cache_dir)
