"""Select alignment OGs by the fraction of selected sample IDs represented."""
import argparse
from collections import Counter
from decimal import Decimal, ROUND_CEILING
import json
import math
from pathlib import Path
import shutil
import tempfile

from alignment_batches import cache_directory, copy_report, move_alignment, validate_groups
from common import atomic_writer, file_record, now, read_tsv, write_json, write_tsv
from mapping_tables import load_tables, read_species


def validate_fraction(value):
    if type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError('alignment.min_sample_fraction must be a number between 0 and 1')
    return value


def minimum_samples(fraction, total):
    validate_fraction(fraction)
    return int((Decimal(str(fraction)) * total).to_integral_value(rounding=ROUND_CEILING))


def check_record(path, saved):
    current = file_record(path)
    if any(current[key] != saved[key] for key in ('bytes', 'sha256')):
        raise ValueError(f'alignment selection input changed: {path}')
    return current


def measure(inputs, samples, mapping, outdir):
    """Count each OG once per sample, reading each selected mapping table once."""
    inputs, outdir = Path(inputs), Path(outdir)
    collection = json.loads((inputs / 'provenance.json').read_text())
    groups = collection['orthogroups']
    validate_groups(groups)
    check_record(samples, collection['samples'])
    check_record(mapping, collection['mapping'])
    # The historical 'species' column is the canonical computational sample ID,
    # including its run accession in current databases. Repeated rows/copies do
    # not add samples; different sample IDs of one biological species do.
    sample_ids = {row['species'] for row in read_tsv(samples)}
    if not sample_ids:
        raise ValueError('alignment selection requires selected samples')
    tables = load_tables(mapping, verify_files=False)
    if set(tables['tables']) != sample_ids:
        raise ValueError('alignment selection mapping differs from selected samples')
    counts = Counter()
    for sample in sorted(sample_ids):
        assignments, _ = read_species(mapping, sample)
        present = {og for names in assignments.values() for og in names}
        counts.update(present)
    if set(counts) != set(groups):
        raise ValueError('alignment collection OGs differ from selected mappings')
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / 'occupancy.tsv'
    write_tsv(path, ['orthogroup', 'samples_with_sequence', 'total_samples'],
              [dict(orthogroup=og, samples_with_sequence=counts[og], total_samples=len(sample_ids))
               for og in groups])
    write_json(outdir / 'provenance.json', {
        'created_at': now(), 'source_collection': file_record(inputs / 'provenance.json'),
        'samples': file_record(samples), 'mapping': file_record(mapping),
        'total_samples': len(sample_ids), 'candidate_orthogroups': len(groups),
        'counting_unit': 'distinct selected analysis_sample_id (samples.tsv species column)',
        'occupancy': file_record(path),
    })


def select(inputs, occupancy, outdir, min_sample_fraction=0.5):
    """Publish only selected input links; keep the original collected FASTAs."""
    validate_fraction(min_sample_fraction)
    inputs, occupancy, outdir = Path(inputs), Path(occupancy), Path(outdir)
    measured = json.loads((occupancy / 'provenance.json').read_text())
    check_record(inputs / 'provenance.json', measured['source_collection'])
    for key in ('samples', 'mapping'):
        check_record(measured[key]['path'], measured[key])
    check_record(occupancy / 'occupancy.tsv', measured['occupancy'])
    rows = read_tsv(occupancy / 'occupancy.tsv')
    groups = json.loads((inputs / 'provenance.json').read_text())['orthogroups']
    validate_groups(groups)
    if len(rows) != len(groups) or {row['orthogroup'] for row in rows} != set(groups):
        raise ValueError('alignment occupancy table differs from the collection')
    total = measured['total_samples']
    required = minimum_samples(min_sample_fraction, total)
    retained, report = [], []
    for row in rows:
        count = int(row['samples_with_sequence'])
        if int(row['total_samples']) != total or not 1 <= count <= total:
            raise ValueError('invalid alignment sample count')
        accepted = count >= required
        if accepted:
            retained.append(row['orthogroup'])
        report.append(dict(row, sample_fraction=count / total, selected=int(accepted)))
    settings = {'min_sample_fraction': min_sample_fraction, 'total_samples': total,
                'minimum_samples': required, 'candidate_orthogroups': len(groups),
                'selected_orthogroups': len(retained), 'excluded_orthogroups': len(groups) - len(retained),
                'counting_unit': measured['counting_unit']}
    outdir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.og-selection-', dir=outdir.parent) as temporary:
        staging = Path(temporary) / 'selected'
        staging.mkdir()
        for og in retained:
            original = (inputs / f'{og}.faa').resolve(strict=True)
            (staging / f'{og}.faa').symlink_to(original)
        write_tsv(staging / 'occupancy.tsv',
                  ['orthogroup', 'samples_with_sequence', 'total_samples', 'sample_fraction', 'selected'], report)
        selection_record = dict(file_record(staging / 'occupancy.tsv'),
                                path=str((outdir / 'occupancy.tsv').resolve()))
        write_json(staging / 'provenance.json', {
            'created_at': now(), 'orthogroups': retained, 'source_collection': measured['source_collection'],
            'samples': measured['samples'], 'mapping': measured['mapping'],
            'selection': settings, 'occupancy': selection_record,
        })
        if outdir.exists():
            shutil.rmtree(outdir)
        staging.replace(outdir)
    print(json.dumps(settings), flush=True)
    return settings


def finish(inputs, outdir, reports, cache=None):
    """Publish selected OGs only, moving other alignments to a reusable work cache."""
    inputs, outdir, reports = Path(inputs), Path(outdir), Path(reports)
    cache = cache_directory(cache if cache is not None else inputs.parent / 'cache', inputs, outdir, reports)
    selected = json.loads((inputs / 'provenance.json').read_text())
    for key in ('source_collection', 'samples', 'mapping'):
        check_record(selected[key]['path'], selected[key])
    groups = selected['orthogroups']
    validate_groups(groups)
    alignments = []
    for og in groups:
        report_path = reports / f'{og}.json'
        report = json.loads(report_path.read_text())
        expected = outdir / f'{og}.faa'
        if Path(report['alignment']['path']).name != expected.name:
            raise ValueError(f'misplaced OG alignment: {og}')
        check_record(inputs / f'{og}.faa', report['input'])
        current = check_record(expected, report['alignment'])
        alignments.append({'orthogroup': og, 'alignment': current, 'provenance': file_record(report_path)})
    occupancy = inputs / 'occupancy.tsv'
    check_record(occupancy, selected['occupancy'])
    rows = read_tsv(occupancy)
    if {r['orthogroup'] for r in rows if r['selected'] == '1'} != set(groups):
        raise ValueError('selected alignment inventory differs from occupancy')
    # Verify every selected alignment before moving any former output. Keep its
    # report with the cached file so later reuse checks input, output and aligner.
    expected_names = {f'{og}.faa' for og in groups}
    candidates = {r['orthogroup'] + '.faa' for r in rows}
    outdir.mkdir(parents=True, exist_ok=True)
    for path in outdir.glob('*.faa'):
        if path.name not in candidates:
            path.unlink()  # Preserve the existing policy for obsolete OGs.
        elif path.name not in expected_names:
            report = reports / f'{path.stem}.json'
            cached_report = cache / report.name
            if report.is_file():
                copy_report(report, cached_report)
            else:
                cached_report.unlink(missing_ok=True)
            move_alignment(path, cache / path.name)
    for path in cache.glob('*.faa'):
        if path.name not in candidates:
            path.unlink()
            path.with_suffix('.json').unlink(missing_ok=True)
    with occupancy.open() as source, atomic_writer(outdir / 'occupancy.tsv') as target:
        shutil.copyfileobj(source, target)
    write_json(outdir / 'provenance.json', {
        'created_at': now(), 'collection': file_record(inputs / 'provenance.json'),
        'selection': selected['selection'], 'occupancy': file_record(outdir / 'occupancy.tsv'),
        'alignments': alignments,
    })


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='action', required=True)
    for action, flags in {'measure': ('inputs', 'samples', 'mapping', 'outdir'),
                          'select': ('inputs', 'occupancy', 'outdir'),
                          'finish': ('inputs', 'outdir', 'reports')}.items():
        command = commands.add_parser(action)
        for flag in flags:
            command.add_argument('--' + flag, required=True)
        if action == 'select':
            command.add_argument('--min-sample-fraction', type=float, default=0.5)
        if action == 'finish':
            command.add_argument('--cache')
    args = vars(parser.parse_args())
    globals()[args.pop('action')](**args)
