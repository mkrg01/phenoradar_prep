"""Checkpoint taxonomy screen batches and run them locally or as scheduler arrays."""
import argparse
from collections import Counter, defaultdict
from contextlib import closing
import csv
from decimal import Decimal
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time

import angiosperms353 as qc
from angiosperms353_reference import BASE_URL, RELEASE, RANKS, KewReference, ReferenceTaxonomy, accessions, kew_families
from common import atomic_writer, now, text_open, write_json
from dataset_assets import SAFE, digest, locked, record, verify
from phase_config import DEFAULT_TAXONOMY_BATCH_SIZE, array_concurrency, execution_config, normalize_slurm, time_minutes
import scheduler
from storage_runtime import task_scratch, task_settings


def batch_path(output, index):
    return output / 'batches' / f'{index:04}'


def batches(plan):
    size = plan['settings']['batch_size']
    samples = plan['samples']
    for offset in range(0, len(samples), size):
        yield offset // size + 1, samples[offset:offset + size]


def tool_versions():
    return {name: subprocess.check_output(command, text=True).strip() for name, command in
            [('captus', ['captus', '--version']), ('mmseqs', ['mmseqs', 'version'])]}


def recover(samples, path, plan, work, log_dir, threads, ram, blat):
    runs = [row['run'] for row in samples]
    if path.exists():
        with text_open(path) as handle:
            saved = json.load(handle)
        if saved['plan'] != plan['id'] or saved['runs'] != runs or len(saved['markers']) != len(runs):
            raise ValueError(f'recovered markers belong to a different batch: {path}')
        return saved['markers']
    only_loci = set(plan['only_loci']) if plan['only_loci'] is not None else None
    files = qc.extract(samples, work, threads=threads, ram=ram, blat=blat, log_dir=log_dir)
    recovered = [{locus: sequence for locus, sequence in qc.markers(file).items()
                  if only_loci is None or locus in only_loci} for file, _ in zip(files, samples, strict=True)]
    with atomic_writer(path) as handle:
        json.dump(dict(plan=plan['id'], runs=runs, markers=recovered), handle, separators=(',', ':'))
    return recovered


def batch_result(output, index, samples, plan):
    path = batch_path(output, index) / 'completed.json'
    if not path.is_file():
        return None
    result = json.loads(path.read_text())
    if (result['plan'] != plan['id']
            or [row['run'] for row in result['summaries']] != [row['run'] for row in samples]):
        raise ValueError(f'incompatible batch result: {path}')
    verify(result['loci'])
    return result


def run_batch(output, plan, index, samples, threads, ram, storage):
    work = batch_path(output, index)
    with locked(work / '.lock'):
        if batch_result(output, index, samples, plan) is not None:
            return
        samples = [dict(row, source_runs=set(row['source_runs']), exclude_accessions=set(row['exclude_accessions']))
                   for row in samples]
        for entry in plan['implementation']:
            verify(entry)
        for row in samples:
            verify(row['cds'])
        blat = str(verify(plan['tools']['blat']))
        if tool_versions() != {name: version for name, version in plan['tools'].items() if name != 'blat'}:
            raise ValueError('screen tool versions differ from the prepared plan')
        families = kew_families(verify(plan['reference']['manifests'][-1]))
        taxonomy = ReferenceTaxonomy(verify(plan['taxonomy']))
        settings = plan['settings']
        with closing(taxonomy.db), task_scratch(Path(plan['root']), task_settings(storage, 'taxonomy_screen', [])) as scratch:
            marker_file = work / 'markers.json.gz'
            reused_markers = marker_file.is_file()
            started = time.monotonic()
            try:
                recovered = recover(samples, marker_file, plan, scratch, work, threads, ram, blat)
            except ValueError as error:
                with (work / 'captus.log').open('a') as log:
                    print(error, file=log)
                raise
            recovery_seconds = time.monotonic() - started
            by_locus = defaultdict(dict)
            for sample_index, markers in enumerate(recovered):
                for locus, sequence in markers.items():
                    by_locus[locus][sample_index] = sequence
            paths = {locus: verify(plan['reference']['loci'][locus]) for locus in sorted(by_locus)}
            ancestry = [{rank: value[0] for rank, value in taxonomy.lineage(int(row['taxid'])).items() if rank in RANKS}
                        for row in samples]
            started = time.monotonic()
            with atomic_writer(work / 'loci.tsv.gz') as handle:
                evidence = csv.DictWriter(handle, qc.LOCUS_FIELDS, delimiter='\t', lineterminator='\n')
                evidence.writeheader()
                summaries = qc.score_batch(samples, ancestry, by_locus, paths, taxonomy, families, scratch, work,
                                           threads, ram, settings['omit_same_species'], settings['min_loci'],
                                           Decimal(settings['max_registered_support']), evidence)
            for row in summaries:
                for field in ('registered', 'alternative'):
                    row[f'{field}_name'] = taxonomy.node(row[field])[1] if row[field] is not None else ''
            search_and_vote_seconds = time.monotonic() - started
        write_json(work / 'completed.json', dict(plan=plan['id'], summaries=summaries, reference=list(paths),
                   loci=record(work / 'loci.tsv.gz'), threads=threads, ram_gb=ram, reused_markers=reused_markers,
                   recovery_seconds=round(recovery_seconds, 3),
                   search_and_vote_seconds=round(search_and_vote_seconds, 3)), compact=True)


def finalize(output, plan):
    with locked(output / '.finalize.lock'):
        results, missing = [], []
        for index, samples in batches(plan):
            result = batch_result(output, index, samples, plan)
            if result is None:
                missing.append(index)
            else:
                results.append(result)
        if missing:
            raise ValueError(f'incomplete taxonomy screen batches: {missing}; rerun screen-taxonomy with the same output')
        loci = sorted({locus for result in results for locus in result['reference']})
        summaries = [row for result in results for row in result['summaries']]
        manifest = output / 'screen.json'
        if manifest.is_file():
            saved = json.loads(manifest.read_text())
            if saved.get('plan_id') != plan['id']:
                raise ValueError(f'incompatible taxonomy screen: {manifest}')
            for entry in [saved['loci'], saved['summary'], saved['exclusion_candidates'], *saved['plots'].values()]:
                verify(entry)
        else:
            reference = [*plan['reference']['manifests'], *(plan['reference']['loci'][locus] for locus in loci)]
            for entry in [plan['taxonomy'], *reference, *[row['cds'] for row in plan['samples']],
                          plan['build_manifest'], plan['database_manifest']]:
                if entry:
                    verify(entry)
            with atomic_writer(output / 'loci.tsv') as handle:
                csv.writer(handle, delimiter='\t', lineterminator='\n').writerow(qc.LOCUS_FIELDS)
                for result in results:
                    with text_open(result['loci']['path']) as source:
                        next(source)
                        shutil.copyfileobj(source, handle)
            files = qc.write_report(output, summaries, Decimal(plan['settings']['max_registered_support']))
            settings = dict(plan['settings'], batch_resources=[
                dict(batch=index, threads=result['threads'], ram_gb=result['ram_gb'])
                for index, result in enumerate(results, start=1)])
            write_json(manifest, dict(schema_version=3, kind='angiosperms353_qc', created_at=now(), plan_id=plan['id'],
                       build_manifest=plan['build_manifest'], database_manifest=plan['database_manifest'],
                       selection_config=plan['selection_config'], taxonomy=plan['taxonomy'], tools=plan['tools'],
                       reference=dict(release=RELEASE, url=BASE_URL, files=reference), settings=settings,
                       samples={row['run']: dict(taxid=str(row['taxid']), cds=row['cds'],
                                excluded_accessions=row['exclude_accessions'], source_runs=row['source_runs'])
                                for row in plan['samples']}, loci=record(output / 'loci.tsv'), **files))
        return dict(output=str(output), samples=len(summaries), loci=len(loci), batches=len(results),
                    statuses=dict(Counter(row['status'] for row in summaries)))


def prepare_plan(output, samples, cache, taxonomy, build, database, selection_config,
                 root, settings, only_loci, blat):
    sample_rows, seen = [], set()
    for row in samples:
        run = row['run']
        if not SAFE.fullmatch(run) or run in seen or not str(row['taxid']).isdigit() or int(row['taxid']) <= 0:
            raise ValueError(f'unsafe or repeated run, or invalid taxid: {run}')
        seen.add(run)
        species_name = row.get('scientific_name', '').replace('_', ' ').casefold()
        if settings['omit_same_species'] and not species_name:
            raise ValueError('--omit-same-species requires scientific_name for every sample')
        sample_rows.append(dict(run=run, taxid=str(row['taxid']), cds=row['cds'], species_name=species_name,
                                source_runs=sorted(accessions(run) | accessions(row.get('source_runs', ''))),
                                exclude_accessions=sorted({entry.upper() for entry in
                                                          re.split(r'[,;\s]+', row.get('exclude_accessions', '')) if entry})))
    executable = shutil.which(blat or 'blat')
    if executable is None:
        raise ValueError(f'BLAT executable not found: {blat or "blat"}')
    inputs = dict(root=str(root), cache=str(Path(cache).resolve()), samples=sample_rows,
                  build_manifest=build, database_manifest=database, taxonomy=record(taxonomy),
                  only_loci=sorted(only_loci) if only_loci is not None else None, settings=settings,
                  tools=dict(tool_versions(), blat=record(executable)),
                  implementation=[record(Path(__file__).with_name(name)) for name in
                                  ('angiosperms353.py', 'angiosperms353_reference.py', 'taxonomy_screen_jobs.py')])
    identity = digest(inputs)
    path = output / 'plan.json'
    if path.is_file():
        saved = json.loads(path.read_text())
        if saved['id'] != identity:
            raise ValueError('existing taxonomy screen has different inputs, software or scientific settings; use a new output')
        return saved
    if any(p.name != '.screen.lock' for p in output.iterdir()):
        raise ValueError(f'output contains files from another run: {output}')
    reference = KewReference(cache)
    manifests = reference.manifests()
    if only_loci is not None and not only_loci:
        raise ValueError('validation loci must be a nonempty subset of Kew Angiosperms353')
    paths = reference.loci(only_loci)
    plan = dict(inputs, id=identity, selection_config=selection_config,
                reference=dict(manifests=[record(path) for path in manifests],
                               loci={locus: record(path) for locus, path in paths.items()}))
    write_json(path, plan, compact=True)
    return plan


def job_resources(config, threads, ram, cluster):
    execution = normalize_slurm(execution_config(config))
    job = {'cpus': 4, 'mem_gb': 16, 'time': '08:00:00'} | execution.get('stages', {}).get('taxonomy_screen', {})
    if threads is not None:
        job['cpus'] = threads
    if ram is None:
        ram = min(8, job['mem_gb'] - 2) if cluster else 8
    if (type(job['cpus']) is not int or job['cpus'] < 1 or type(ram) is not int or ram < 2
            or type(job['mem_gb']) is not int or job['mem_gb'] < 1):
        raise ValueError('screen requires at least one CPU and at least 2 GB of tool memory')
    if cluster and ram > job['mem_gb'] - 2:
        raise ValueError('--ram must leave at least 2 GB within the taxonomy_screen scheduler memory request')
    time_minutes(job['time'])
    return execution, job, ram


def submit(output, plan, missing, execution, job, ram, storage):
    from build_dispatch import check_submissions

    root = Path(plan['root'])
    jobs = output / 'jobs'
    jobs.mkdir(exist_ok=True)
    settings = scheduler.resolve(execution)
    limits = dict(execution, concurrency=execution.get('concurrency', execution.get('jobs')),
                  stages={'taxonomy_screen': job})
    if limits['concurrency'] is not None and (type(limits['concurrency']) is not int or limits['concurrency'] < 1):
        raise ValueError('execution.total_limits.jobs must be a positive integer or null')
    for name, value in limits.get('total_limits', {}).items():
        if value is not None and (type(value) is not int or value < 1):
            raise ValueError(f'execution.total_limits.{name} must be a positive integer or null')
    cap = array_concurrency(limits, 'taxonomy_screen')
    limit = scheduler.array_size(settings) if missing else 1
    number = 1 + max((int(path.stem.split('_')[1]) for path in jobs.glob('submission_*.json')), default=0)
    receipt = jobs / f'submission_{number:04}.json'
    state = dict(created_at=now(), jobs=[])
    script_base = [sys.executable, str(Path(__file__).resolve()), str(output)]

    def enqueue(label, script, dependencies, indices=None, resources=job):
        log = jobs / f'{number:04}_{label}'
        template = str(log) + ('_%A_%a' if indices else '_%j')
        command = scheduler.batch_command(settings, resources, root, f'taxonomy_{label}',
                                          template + '.out', template + '.err', script=script,
                                          indices=indices,
                                          concurrency=min(cap, len(indices)) if cap is not None and indices else None,
                                          afterany=dependencies)
        entry = dict(stage=label, scheduler=settings['scheduler'], state='submitting', command=command)
        state['jobs'].append(entry)
        write_json(receipt, state)
        try:
            response = subprocess.check_output(command, cwd=root, text=True)
            job_id = scheduler.parse_id(response, settings['scheduler'])
        except (OSError, subprocess.SubprocessError, ValueError):
            entry['state'] = 'unknown'
            write_json(receipt, state)
            raise
        entry.update(state='submitted', job_id=job_id)
        write_json(receipt, state)
        return job_id

    check_submissions(jobs)
    arrays = []
    # Chain arrays to preserve the combined resource cap on clusters with small array limits.
    for offset in range(0, len(missing), limit):
        indices = missing[offset:offset + limit]
        label = f'batch_{offset // limit + 1:03}'
        script = jobs / f'{number:04}_{label}.sh'
        task_file = script.with_suffix('.json')
        write_json(task_file, dict(indices=indices, threads=job['cpus'], ram=ram, storage=storage))
        command = shlex.join([*script_base, 'batch', str(task_file)])
        script.write_text(f'#!/bin/bash\nset -eu\nexec {command}\n')
        script.chmod(0o755)
        arrays.append(enqueue(label, script, arrays[-1:], range(1, len(indices) + 1)))
    # AGE hold_jid does not require success. The finalizer checks all batch receipts.
    finalizer = jobs / f'{number:04}_finalize.sh'
    finalizer.write_text('#!/bin/bash\nset -eu\nexec ' + shlex.join([*script_base, 'finalize']) + '\n')
    finalizer.chmod(0o755)
    final_job = dict(cpus=1, mem_gb=min(job['mem_gb'], 8), time=job['time'])
    final_id = enqueue('finalize', finalizer, arrays[-1:], resources=final_job)
    return dict(status='submitted', output=str(output), scheduler=settings['scheduler'], arrays=arrays,
                finalize_job=final_id, pending_batches=len(missing), concurrency=cap)


def screen(samples, *, cache, taxonomy, output, build=None, database=None, selection_config=None, min_loci=50,
           registered_support=None, support=None, threads=None, ram=None, omit_same_species=False, only_loci=None,
           blat=None, batch_size=DEFAULT_TAXONOMY_BATCH_SIZE, root=None, config=None, cluster=False):
    try:
        import matplotlib  # noqa: F401
    except ImportError as error:
        raise ValueError('the Angiosperms353 environment must include matplotlib') from error
    cutoff = qc.registered_threshold(min_loci, registered_support, support)
    if type(batch_size) is not int or batch_size < 1 or not samples:
        raise ValueError('screen requires samples and a positive batch size')
    output = Path(output).resolve()
    root = Path(root or output.parent).resolve()
    config = config or {}
    execution, job, ram = job_resources(config, threads, ram, cluster)
    storage = config.get('storage', {})
    task_settings(storage, 'taxonomy_screen', [])
    settings = dict(min_loci=min_loci, max_registered_support=str(cutoff), batch_size=batch_size,
                    rank_selection='closest_ancestor_with_reference', vote_denominator='comparable_loci',
                    voting='common_identity_split_all_winners', min_bases=qc.MIN_BASES, coverage=qc.COVERAGE,
                    evalue=qc.EVALUE, captus_scipio_passes=1, mmseqs_kmer=10,
                    omit_same_species=omit_same_species, restricted_loci=only_loci is not None)
    with locked(output / '.screen.lock'):
        plan = prepare_plan(output, samples, cache, taxonomy, build, database, selection_config,
                            root, settings, only_loci, blat)
        pending = [(index, rows) for index, rows in batches(plan) if batch_result(output, index, rows, plan) is None]
        if (output / 'screen.json').is_file():
            return finalize(output, plan)
        if cluster:
            return submit(output, plan, [index for index, _ in pending], execution, job, ram, storage)
        if (output / 'jobs').exists():
            from build_dispatch import check_submissions
            check_submissions(output / 'jobs')
        for index, rows in pending:
            run_batch(output, plan, index, rows, job['cpus'], ram, storage)
        return finalize(output, plan)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    commands = parser.add_subparsers(dest='command', required=True)
    commands.add_parser('finalize')
    worker = commands.add_parser('batch')
    worker.add_argument('task_file', type=Path)
    args = parser.parse_args()
    plan = json.loads((args.output / 'plan.json').read_text())
    if args.command == 'batch':
        task = int(os.environ.get('SLURM_ARRAY_TASK_ID') or os.environ['SGE_TASK_ID'])
        spec = json.loads(args.task_file.read_text())
        if not 1 <= task <= len(spec['indices']):
            raise ValueError(f'invalid taxonomy screen array index: {task}')
        index = spec['indices'][task - 1]
        run_batch(args.output, plan, index, dict(batches(plan))[index], spec['threads'], spec['ram'], spec['storage'])
    else:
        print(json.dumps(finalize(args.output, plan), indent=2))


if __name__ == '__main__':
    main()
