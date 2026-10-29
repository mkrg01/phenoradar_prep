"""One ready-to-use metadata table, with explicit provenance for every pair source."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from common import read_tsv, write_tsv
from consumer_metadata import create, validate


@pytest.fixture
def cohort(tmp_path):
    base = tmp_path / 'metadata/species_metadata_base.tsv'
    taxonomy = tmp_path / 'metadata/metadata_high_busco.tsv'
    output = tmp_path / 'metadata/species_metadata.tsv'
    samples = ['Plant_A_run_1', 'Plant_A_run_2', 'Plant_B_run_3', 'Plant_C_run_4']
    write_tsv(base, ['species', 'habit', 'contrast_pair_id', 'family'],
              [dict(species=s, habit=t, contrast_pair_id='stale', family='Family A')
               for s, t in zip(samples, ['0', '0', '1', ''])])
    write_tsv(taxonomy, ['species', 'family', 'order', 'genus'],
              [dict(species=s, family='Family A', order='Order A', genus=s.split('_run')[0]) for s in samples])
    paths = {}
    for branch, pair_values in [('representatives', ['7', '7', '7', '']), ('all', ['2', '', '2', ''])]:
        path = tmp_path / f'phylogeny/{branch}/contrast/species_metadata.tsv'
        write_tsv(path, ['species', 'habit', 'contrast_pair_id'],
                  [dict(species=s, habit=t, contrast_pair_id=p)
                   for s, t, p in zip(samples, ['0', '0', '1', ''], pair_values)])
        paths[branch] = path
    return dict(base=base, taxonomy=taxonomy, output=output, trait='habit'), paths


@pytest.mark.parametrize('branches', [[], ['representatives'], ['all', 'representatives']])
def test_explicit_sources_keep_cohort_taxonomy_missing_values_and_branch_identity(cohort, branches):
    kwargs, paths = cohort
    rows = create(**kwargs, pair_tables={branch: paths[branch] for branch in branches})
    assert len(rows) == 4 and [r['species'] for r in rows] == [r['species'] for r in read_tsv(kwargs['base'])]
    assert all(r['family'] == 'Family A' and r['order'] == 'Order A' and r['class'] == '' for r in rows)
    assert rows[-1]['habit'] == ''
    if len(branches) == 1:
        assert [r['contrast_pair_id'] for r in rows] == ['7', '7', '7', '']
    else:
        assert all(r['contrast_pair_id'] == '' for r in rows)
    if len(branches) > 1:
        assert [r['contrast_pair_id_all'] for r in rows] == ['2', '', '2', '']
        assert [r['contrast_pair_id_representatives'] for r in rows] == ['7', '7', '7', '']
    assert validate(kwargs['output'], 'habit') == rows
    report = json.loads(kwargs['output'].with_suffix('.json').read_text())
    assert set(report['sources']['contrasts']) == set(branches)
    assert report['sources']['base']['path'] == 'metadata/species_metadata_base.tsv'


def test_partial_contrast_tables_leave_other_samples_unpaired(cohort):
    kwargs, paths = cohort
    path = paths['representatives']
    rows = read_tsv(path)[:1]
    write_tsv(path, list(rows[0]), rows)
    result = create(**kwargs, pair_tables={'representatives': path})
    assert [r['contrast_pair_id'] for r in result] == ['7', '', '', '']


@pytest.mark.parametrize('problem', ['unknown_branch', 'unknown_sample', 'trait_mismatch', 'duplicate',
                                     'taxonomy_mismatch', 'conflicting_order', 'missing_taxonomy'])
def test_invalid_joins_do_not_replace_previous_outputs(cohort, problem):
    kwargs, paths = cohort
    create(**kwargs)
    previous = {p: p.read_bytes() for p in [kwargs['output'], kwargs['output'].with_suffix('.json')]}
    pairs = {'representatives': paths['representatives']}
    pair_rows = read_tsv(paths['representatives'])
    if problem == 'unknown_branch': pairs = {'arbitrary': paths['representatives']}
    elif problem == 'unknown_sample': pair_rows[0]['species'] = 'Unknown'
    elif problem == 'trait_mismatch': pair_rows[0]['habit'] = '1'
    elif problem == 'duplicate': pair_rows.append(dict(pair_rows[0]))
    else:
        rows = read_tsv(kwargs['taxonomy'])
        if problem == 'taxonomy_mismatch': rows[0]['family'] = 'Other family'
        elif problem == 'missing_taxonomy': rows.pop()
        else: rows.append({**rows[0], 'order': 'Conflicting order'})
        write_tsv(kwargs['taxonomy'], list(rows[0]), rows)
    write_tsv(paths['representatives'], list(pair_rows[0]), pair_rows)
    with pytest.raises(ValueError): create(**kwargs, pair_tables=pairs)
    assert all(p.read_bytes() == data for p, data in previous.items())


@pytest.mark.parametrize('key', ['base', 'taxonomy', 'pairs', 'output'])
def test_verification_detects_changed_sources_and_output(cohort, key):
    kwargs, paths = cohort
    create(**kwargs, pair_tables={'representatives': paths['representatives']})
    path = paths['representatives'] if key == 'pairs' else kwargs[key]
    path.write_text(path.read_text().replace('Family A', 'Family B').replace('\t7', '\t8').replace('stale', 'changed'))
    with pytest.raises(ValueError, match='checksum mismatch'): validate(kwargs['output'], 'habit')


def test_production_rule_waits_for_selected_contrasts_and_ignores_stale_branches(cohort):
    kwargs, paths = cohort
    root = kwargs['output'].parent.parent
    snakemake = os.environ.get('SNAKEMAKE_BIN') or shutil.which('snakemake')
    if not snakemake:
        pytest.skip('Snakemake required')
    repo = Path(__file__).resolve().parents[1]
    code = (repo / 'workflow/rules/metadata.smk').read_text()
    rule = code[code.index('def consumer_pair_tables():'):code.index('rule record_run:')]
    for path in paths.values():
        path.rename(path.with_name('source.tsv'))
    workflow = root / 'workflow/Snakefile'
    workflow.parent.mkdir()
    (root / 'envs').mkdir()
    shutil.copyfile(repo / 'workflow/envs/analysis.yaml', root / 'envs/analysis.yaml')
    config = root / 'config.json'
    workflow.write_text(f'''
import json
from pathlib import Path
configfile: {str(config)!r}
OUT = {str(root)!r}
META = OUT + '/metadata'
LOG = OUT + '/logs'
SCRIPTS = {str(repo / 'workflow/scripts')!r}
PYTHON = {sys.executable!r}
CONSUMER_METADATA = 'metadata/species_metadata.tsv'
BASE_METADATA = 'metadata/species_metadata_base.tsv'
PHY = config['phylogeny']
TREE_BRANCHES = {{name: 'phylogeny/' + name for name in ['all', 'phenotyped', 'representatives']}}
rule all:
    input: META + '/species_metadata.tsv'
rule saved_contrasts:
    input: OUT + '/phylogeny/{{branch}}/contrast/source.tsv'
    output: OUT + '/phylogeny/{{branch}}/contrast/species_metadata.tsv'
    shell: 'cp {{input:q}} {{output:q}}'
''' + rule)
    for enabled, branches in [(True, ['representatives']), (True, ['all', 'representatives']),
                              (False, ['all', 'representatives'])]:
        config.write_text(json.dumps({'trait': 'habit', 'phylogeny':
                          {'trees': branches, 'contrast_pairs': {'enabled': enabled}}}))
        result = subprocess.run([snakemake, '--snakefile', str(workflow), '--cores', '1'],
                                cwd=root, text=True, capture_output=True, timeout=60)
        assert result.returncode == 0, result.stdout + result.stderr
        ready = validate(kwargs['output'], 'habit')
        if not enabled:
            assert all(r['contrast_pair_id'] == '' for r in ready)
            assert not any(k.startswith('contrast_pair_id_') for k in ready[0])
        elif len(branches) == 1:
            assert ready[0]['contrast_pair_id'] == '7'
        else:
            assert ready[0]['contrast_pair_id_all'] == '2'
            assert ready[0]['contrast_pair_id_representatives'] == '7'
