"""Allocated-node comparison on copies of small/median/large real OG inputs."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import psutil

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
from alignment_batches import plan
from common import file_record, write_json


def benchmark(inputs, outdir, image, environment):
    inputs, outdir = Path(inputs), Path(outdir)
    names = json.loads((inputs / 'provenance.json').read_text())['orthogroups']
    weighted = sorted(((inputs / f'{og}.faa').stat().st_size, og) for og in names)
    selected = sorted({weighted[max(0, min(len(weighted)-1, int(q * (len(weighted)-1)) + offset))][1]
                       for q in (.05, .50, .90, .99) for offset in (-1, 0, 1)})
    copies = outdir / 'inputs'; copies.mkdir(parents=True, exist_ok=True)
    for og in selected:
        shutil.copy2(inputs / f'{og}.faa', copies / f'{og}.faa')
    write_json(copies / 'provenance.json', {'orthogroups': selected,
                                          'source': file_record(inputs / 'provenance.json')})
    plans = outdir / 'plan'; plan(copies, plans)
    index = json.loads((plans / 'index.json').read_text())
    runner = Path(__file__).resolve().parents[1] / 'workflow/scripts/alignment_batches.py'
    measurements = []
    for label, cpus, threads in [('serial_4', 4, 4), ('parallel_3x4', 12, 4), ('parallel_6x2', 12, 2)]:
        folder = outdir / label; folder.mkdir(exist_ok=True)
        began, peak = time.monotonic(), 0
        with (folder / 'run.log').open('w') as log:
            for batch in index['batches']:
                command = ['apptainer', 'exec', '--cleanenv', '--bind', str(Path.cwd().parent), str(image),
                           '/usr/bin/env', f'PATH={environment}/bin:/usr/bin:/bin', 'OPENBLAS_NUM_THREADS=1',
                           f'{environment}/bin/python', str(runner), 'batch', '--batch', str(plans / f'{batch}.json'),
                           '--inputs', str(copies), '--outdir', str(folder / 'alignments'),
                           '--reports', str(folder / 'reports'), '--receipt', str(folder / f'{batch}.json'),
                           '--cpus', str(cpus), '--threads-per-og', str(threads)]
                process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
                parent = psutil.Process(process.pid)
                while process.poll() is None:
                    try:
                        family = [parent, *parent.children(recursive=True)]
                        rss = sum(p.memory_info().rss for p in family if p.is_running())
                        peak = max(peak, rss)
                    except psutil.Error:
                        pass
                    time.sleep(.05)
                if process.returncode:
                    raise RuntimeError(f'{label} failed; see {folder / "run.log"}')
        measurements.append({'mode': label, 'seconds': time.monotonic() - began,
                             'peak_rss_mb': peak / 1000000})
    baseline = outdir / 'serial_4/alignments'
    for measurement in measurements:
        candidate = outdir / measurement['mode'] / 'alignments'
        measurement['identical_alignments'] = all((baseline / f'{og}.faa').read_bytes() ==
                                                 (candidate / f'{og}.faa').read_bytes() for og in selected)
    if not measurements[1]['identical_alignments']:
        raise ValueError('batching changed four-thread FAMSA outputs')
    result = {'orthogroups': selected, 'input_bytes': {og: (copies / f'{og}.faa').stat().st_size for og in selected},
              'measurements': measurements}
    write_json(outdir / 'result.json', result)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('inputs', 'outdir', 'image', 'environment'):
        parser.add_argument('--' + name, required=True)
    benchmark(**vars(parser.parse_args()))
