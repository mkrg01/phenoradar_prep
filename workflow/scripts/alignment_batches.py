"""Weighted OG batches with bounded parallelism and per-OG recovery."""
import argparse
from contextlib import contextmanager
from contextvars import ContextVar
from concurrent.futures import ThreadPoolExecutor, as_completed
import errno
import heapq
import json
import os
from pathlib import Path
import shutil
import tempfile

from align_orthogroups import align, SAFE_OG
from common import atomic_writer, file_record, write_json
from dataset_assets import digest, record, stat_identity

_INSPECTION = ContextVar('alignment_batch_inspection', default=None)


@contextmanager
def inspect_batches():
    """Collect invalid receipts for a read-only planner instead of deleting them."""
    invalid = set()
    token = _INSPECTION.set(invalid)
    try:
        yield invalid
    finally:
        _INSPECTION.reset(token)


def plan(inputs, outdir, batch_size=500, heavy_bytes=16000000):
    if type(batch_size) is not int or batch_size < 1:
        raise ValueError('batch size must be positive')
    inputs, outdir = Path(inputs), Path(outdir)
    collection = file_record(inputs / 'provenance.json')
    groups = json.loads((inputs / 'provenance.json').read_text())['orthogroups']
    validate_groups(groups)
    weights = [(og, (inputs / f'{og}.faa').stat().st_size) for og in groups]
    heavy = [[item] for item in weights if item[1] >= heavy_bytes]
    small = sorted((item for item in weights if item[1] < heavy_bytes), key=lambda x: (-x[1], x[0]))
    count = (len(small) + batch_size - 1) // batch_size
    bins = [[] for _ in range(count)]
    heap = [(0, i) for i in range(count)]
    for item in small:
        total, i = heapq.heappop(heap)
        bins[i].append(item)
        if len(bins[i]) < batch_size:
            heapq.heappush(heap, (total + item[1], i))
    outdir.mkdir(parents=True, exist_ok=True)
    batches = []
    for i, members in enumerate(heavy + bins):
        name = f'batch_{i:06d}'
        payload = {'batch': name, 'collection_sha256': collection['sha256'],
                   'orthogroups': sorted(og for og, _ in members),
                   'input_bytes': sum(size for _, size in members)}
        payload['sha256'] = digest(payload)
        path = outdir / f'{name}.json'
        if not path.exists() or json.loads(path.read_text()) != payload:
            write_json(path, payload)
        batches.append(name)
    write_json(outdir / 'index.json', {'collection': collection, 'batches': batches,
                                      'orthogroups': len(groups), 'batch_size': batch_size})


def batch_targets(plan_dir, outdir, reports):
    """Invalidate derived receipts for missing/changed members in aggregate targets.

    This audit runs on aggregate DAG expansion, not on each batch worker. OG
    files are never deleted: a repaired batch reuses its remaining valid OGs.
    """
    plan_dir = Path(plan_dir)
    index = json.loads((plan_dir / 'index.json').read_text())
    targets = []
    for name in index['batches']:
        path = Path(outdir) / 'batches' / f'{name}.json'
        if path.exists():
            try:
                receipt = json.loads(path.read_text())
                batch = json.loads((plan_dir / f'{name}.json').read_text())
                if (receipt['plan_sha256'] != batch['sha256'] or
                        set(receipt['members']) != set(batch['orthogroups'])):
                    raise ValueError('batch membership changed')
                for og, files in receipt['members'].items():
                    expected = [Path(index['collection']['path']).parent / f'{og}.faa',
                                Path(outdir) / f'{og}.faa', Path(reports) / f'{og}.json']
                    if len(files) != len(expected) or any(
                            Path(entry['path']).resolve() != filename.resolve() or
                            entry['stat'] != stat_identity(filename)
                            for entry, filename in zip(files, expected)):
                        raise ValueError('batch member changed')
            except (OSError, ValueError, KeyError, TypeError):
                if _INSPECTION.get() is None:
                    path.unlink()
                else:
                    _INSPECTION.get().add(str(path))
        targets.append(str(path))
    return targets


def reusable(fasta, output, report, threads, executable):
    try:
        saved = json.loads(Path(report).read_text())
        if saved['method'] != 'singleton':
            if saved['threads'] != threads or saved['executable']['sha256'] != executable['sha256']:
                return False
        for path, key in [(fasta, 'input'), (output, 'alignment')]:
            current = file_record(path)
            if any(current[k] != saved[key][k] for k in ('bytes', 'sha256')):
                return False
        return True
    except (OSError, ValueError, KeyError, TypeError):
        return False


def move_alignment(source, destination):
    """Move without duplicating storage; publish cross-filesystem copies atomically."""
    source, destination = Path(source), Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        source.replace(destination)
    except OSError as error:
        if error.errno != errno.EXDEV:
            raise
        fd, name = tempfile.mkstemp(prefix=f'.{destination.name}.', dir=destination.parent)
        os.close(fd)
        try:
            shutil.copy2(source, name)
            os.replace(name, destination)
            source.unlink()
        finally:
            Path(name).unlink(missing_ok=True)


def copy_report(source, destination):
    with Path(source).open('rb') as original, atomic_writer(destination, 'wb') as target:
        shutil.copyfileobj(original, target)


def cache_directory(cache, inputs, outdir, reports):
    path = Path(cache)
    if any(path.resolve().is_relative_to(Path(folder).resolve()) for folder in (inputs, outdir, reports)):
        raise ValueError('alignment cache must be outside inputs, final alignments and reports')
    return path


def validate_groups(groups):
    if (not isinstance(groups, list) or any(not isinstance(og, str) or not SAFE_OG.fullmatch(og) for og in groups)
            or len(set(groups)) != len(groups)):
        raise ValueError('batch OG IDs must be unique, safe filename components')


def run_batch(batch, inputs, outdir, reports, receipt, cpus=12, threads_per_og=4, command='famsa', cache=None):
    if any(type(value) is not int or value < 1 for value in (cpus, threads_per_og)):
        raise ValueError('CPU counts must be positive')
    payload = json.loads(Path(batch).read_text())
    if digest({k: v for k, v in payload.items() if k != 'sha256'}) != payload.get('sha256'):
        raise ValueError('alignment batch plan changed')
    validate_groups(payload['orthogroups'])
    inputs, outdir, reports = Path(inputs), Path(outdir), Path(reports)
    cache = cache_directory(cache, inputs, outdir, reports) if cache is not None else None
    if file_record(inputs / 'provenance.json')['sha256'] != payload['collection_sha256']:
        raise ValueError('alignment collection changed')
    threads = min(cpus, threads_per_og)
    parallel = max(1, cpus // threads)
    binary = shutil.which(command)
    executable = file_record(binary) if binary else None

    def process(og):
        fasta, output, provenance = inputs / f'{og}.faa', outdir / f'{og}.faa', reports / f'{og}.json'
        reused = reusable(fasta, output, provenance, threads, executable)
        if not reused and cache is not None and reusable(
                fasta, cache / f'{og}.faa', cache / f'{og}.json', threads, executable):
            copy_report(cache / f'{og}.json', provenance)
            move_alignment(cache / f'{og}.faa', output)
            reused = True
        if not reused:
            align(fasta, output, provenance, threads=threads, command=command)
        if cache is not None:
            # A selected OG lives in the final directory, not in both locations.
            (cache / f'{og}.faa').unlink(missing_ok=True)
            (cache / f'{og}.json').unlink(missing_ok=True)
        return og, [record(p) for p in (fasta, output, provenance)], reused

    members, reused_count = {}, 0
    with ThreadPoolExecutor(max_workers=parallel) as pool:
        futures = [pool.submit(process, og) for og in payload['orthogroups']]
        try:
            for future in as_completed(futures):
                og, files, reused = future.result()
                members[og] = files
                reused_count += int(reused)
        except BaseException:
            for future in futures:
                future.cancel()
            raise
    write_json(receipt, {'plan_sha256': payload['sha256'], 'cpus': cpus,
                         'parallel': parallel, 'threads_per_og': threads,
                         'reused': reused_count, 'computed': len(members) - reused_count,
                         'members': members})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    p = sub.add_parser('plan')
    p.add_argument('--inputs', required=True); p.add_argument('--outdir', required=True)
    p.add_argument('--batch-size', type=int, default=500)
    p = sub.add_parser('batch')
    for flag in ('batch', 'inputs', 'outdir', 'reports', 'receipt'):
        p.add_argument('--' + flag, required=True)
    p.add_argument('--cpus', type=int, default=12)
    p.add_argument('--threads-per-og', type=int, default=4)
    p.add_argument('--cache')
    args = vars(parser.parse_args())
    action = args.pop('action')
    (plan if action == 'plan' else run_batch)(**args)
