"""Sample occupancy, copy handling, boundary counts and reusable OG outputs."""
import json
import sys

import pytest

from alignment_batches import plan, run_batch
from alignment_selection import measure, select, finish, minimum_samples
from align_orthogroups import collect, align
from common import read_tsv, write_tsv
from mapping_fixtures import make_mapping


@pytest.fixture
def occupancy_inputs(tmp_path):
    names = ['Plant_R1', 'Plant_R2', 'Other_R3', 'Empty_R4']
    rows = [dict(species=name, odb_species=name, run='R' + str(i + 1)) for i, name in enumerate(names)]
    # Repeating the same canonical sample row must not change its weight.
    write_tsv(tmp_path / 'samples.tsv', list(rows[0]), rows + [rows[0]])
    genes = [(name + '_g1', name) for name in names] + [('Plant_R1_g2', 'Plant_R1'),
                                                          ('Plant_R1_g3', 'Plant_R1')]
    pairs = [(name + '_g1', 'OGthree') for name in names[:3]]
    pairs += [(name + '_g1', 'OGhalf') for name in names[:2]]
    pairs += [('Plant_R1_g2', 'OGlow'), ('Plant_R1_g3', 'OGlow')]
    mapping = make_mapping(tmp_path / 'mapping/snapshot.json', genes, pairs)
    proteins = tmp_path / 'proteins'; proteins.mkdir()
    for name in names:
        (proteins / f'{name}_protein.fa').write_text(''.join(
            f'>{gene}\nMAA\n' for gene, sample in genes if sample == name))
    raw = tmp_path / 'inputs'
    collect(tmp_path / 'samples.tsv', mapping, proteins, raw)
    occupancy = tmp_path / 'occupancy'
    measure(raw, tmp_path / 'samples.tsv', mapping, occupancy)
    return raw, occupancy, tmp_path / 'selected'


def test_default_fraction_counts_samples_not_copies_or_biological_species(occupancy_inputs):
    raw, occupancy, selected = occupancy_inputs
    counts = {r['orthogroup']: r for r in read_tsv(occupancy / 'occupancy.tsv')}
    assert {name: int(row['samples_with_sequence']) for name, row in counts.items()} == {
        'OGthree': 3, 'OGhalf': 2, 'OGlow': 1}
    assert {r['total_samples'] for r in counts.values()} == {'4'}
    original = {p: p.stat().st_mtime_ns for p in raw.iterdir()}
    result = select(raw, occupancy, selected)
    assert result['minimum_samples'] == 2 and result['total_samples'] == 4
    assert result['selected_orthogroups'] == 2 and result['excluded_orthogroups'] == 1
    assert {p.stem for p in selected.glob('*.faa')} == {'OGthree', 'OGhalf'}
    assert all(p.is_symlink() for p in selected.glob('*.faa'))
    assert all(p.stat().st_mtime_ns == saved for p, saved in original.items())


@pytest.mark.parametrize('fraction,total,minimum', [(0, 5582, 0), (0.5, 5582, 2791),
    (0.5, 3, 2), (0.07, 100, 7), (0.14, 100, 14), (1, 3, 3)])
def test_inclusive_boundary_uses_exact_decimal_ceiling(fraction, total, minimum):
    assert minimum_samples(fraction, total) == minimum


def test_zero_and_one_fractions_and_empty_selection(occupancy_inputs):
    raw, occupancy, selected = occupancy_inputs
    assert select(raw, occupancy, selected, 0)['selected_orthogroups'] == 3
    assert select(raw, occupancy, selected, 1)['selected_orthogroups'] == 0
    plans = selected.parent / 'plans'
    plan(selected, plans)
    assert json.loads((plans / 'index.json').read_text())['batches'] == []
    out = selected.parent / 'alignments'
    finish(selected, out, selected.parent / 'reports')
    saved = json.loads((out / 'provenance.json').read_text())
    assert saved['alignments'] == [] and saved['selection']['min_sample_fraction'] == 1


def test_threshold_changes_reuse_alignments_and_keep_excluded_results(occupancy_inputs, tmp_path):
    raw, occupancy, selected = occupancy_inputs
    out, reports, plans = tmp_path / 'out', tmp_path / 'reports', tmp_path / 'plans'
    famsa = tmp_path / 'famsa'
    famsa.write_text(f'#!{sys.executable}\nfrom pathlib import Path\nimport sys\n'
                     'Path(sys.argv[-1]).write_bytes(Path(sys.argv[-2]).read_bytes())\n')
    famsa.chmod(0o755)
    # A completed result excluded by the new default must survive untouched.
    align(raw / 'OGlow.faa', out / 'OGlow.faa', reports / 'OGlow.json', threads=4, command=str(famsa))
    old = (out / 'OGlow.faa').stat().st_mtime_ns
    for fraction, expected, reused in [(0.5, {'OGhalf', 'OGthree'}, 0),
                                       (0.75, {'OGthree'}, 1),
                                       (0.25, {'OGlow', 'OGhalf', 'OGthree'}, 3)]:
        select(raw, occupancy, selected, fraction)
        plan(selected, plans)
        receipt = out / 'batches/batch_000000.json'
        run_batch(plans / 'batch_000000.json', selected, out, reports, receipt,
                  cpus=4, command=str(famsa))
        assert json.loads(receipt.read_text())['reused'] == reused
        finish(selected, out, reports)
        saved = json.loads((out / 'provenance.json').read_text())
        assert {r['orthogroup'] for r in saved['alignments']} == expected
        assert saved['selection']['min_sample_fraction'] == fraction
        assert (out / 'OGlow.faa').stat().st_mtime_ns == old


def test_selection_rejects_changed_sources(occupancy_inputs):
    raw, occupancy, selected = occupancy_inputs
    (raw / 'provenance.json').write_text('{}')
    with pytest.raises(ValueError, match='input changed'):
        select(raw, occupancy, selected)


def test_selection_exports_omit_preserved_below_threshold_results(occupancy_inputs, tmp_path):
    from analysis_validation import discover
    raw, occupancy, selected = occupancy_inputs
    select(raw, occupancy, selected, 1)
    source = tmp_path / 'source'
    out = source / 'orthogroups/alignments'
    out.mkdir(parents=True)
    (out / 'OGlow.faa').write_bytes((raw / 'OGlow.faa').read_bytes())
    finish(selected, out, tmp_path / 'reports')
    samples = [{**r, 'scientific_name': r['species'], 'taxid': str(i + 1)}
               for i, r in enumerate(read_tsv(tmp_path / 'samples.tsv')[:4])]
    write_tsv(source / 'metadata/samples.tsv', list(samples[0]), samples)
    write_tsv(source / 'metadata/species_metadata.tsv', ['species'],
              [dict(species=r['species']) for r in samples])
    observed = discover(source)
    assert observed['sections']['alignments'] == 'ready'
    assert all(not str(p).endswith('OGlow.faa') for p in observed['files'])
