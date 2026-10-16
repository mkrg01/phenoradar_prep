"""Analysis plans stay readable while preserving complete selection reports."""
import json
from pathlib import Path

import pytest
import yaml

import analysis
from build_products import load_complete
from common import read_tsv
from test_datasets import dataset_project, imported
from test_submit_cli import invoke


@pytest.fixture
def plan_project(dataset_project):
    root = dataset_project
    build = imported(root)
    config = root / 'config/analysis.yaml'
    cfg = yaml.safe_load(config.read_text())
    cfg['build'] = str(build.relative_to(root))
    cfg['inputs'] = {'species_trait': None, 'species_list': None, 'calibrations': None}
    cfg['selection']['busco_threshold'] = 0.5
    cfg['exclude_species'] = []
    cfg['alignment']['enabled'] = False
    cfg['kegg']['enabled'] = False
    cfg['phylogeny']['trees'] = []
    for section in ('dating', 'taxonomy_check', 'contrast_pairs'):
        cfg['phylogeny'][section]['enabled'] = False
    config.write_text(yaml.safe_dump(cfg))
    return root, build, config, cfg


@pytest.mark.parametrize('selection', ['busco', 'species_list', 'explicit'])
def test_plan_summarizes_selection_and_preserves_every_sample(plan_project, monkeypatch, capsys, selection):
    root, build, config, cfg = plan_project
    reasons = {'Alpha_plant_A1': '', 'Beta_sp-X_B1': '', 'Gamma_plant_G1': 'below_busco_threshold'}
    if selection == 'species_list':
        (root / 'input/subset.txt').write_text('Beta_sp-X\nGamma_plant\n')
        cfg['inputs']['species_list'] = 'input/subset.txt'
        reasons['Alpha_plant_A1'] = 'outside_species_list'
    elif selection == 'explicit':
        cfg['exclude_species'] = ['Beta_sp-X']
        reasons['Beta_sp-X_B1'] = 'exclude_species'
    config.write_text(yaml.safe_dump(cfg))
    # The analysis report must not overwrite the build plan.
    build_report = root / 'reports/plan/plan.tsv'
    build_report.parent.mkdir(parents=True)
    build_report.write_text('saved build plan\n')

    invoke(monkeypatch, analysis, 'plan', '--root', root, '--build', build)
    output = capsys.readouterr().out
    directory = root / 'reports/analysis/plan'
    rows = read_tsv(directory / 'plan.tsv')
    excluded = read_tsv(directory / 'excluded.tsv')
    assert len(rows) == 3
    assert list(rows[0]) == ['species', 'selected', 'reason', 'busco_complete_fraction']
    assert {r['species']: r['reason'] for r in rows} == reasons
    assert {r['species']: r['selected'] for r in rows} == {
        name: str(not reason) for name, reason in reasons.items()}
    assert {r['species']: float(r['busco_complete_fraction']) for r in rows} == {
        'Alpha_plant_A1': 0.8, 'Beta_sp-X_B1': 0.6, 'Gamma_plant_G1': 0.4}
    assert excluded == [r for r in rows if r['selected'] == 'False']
    assert f'Samples: {3 - len(excluded)} selected; {len(excluded)} excluded' in output
    assert f'Explicit exclusions: {int(selection == "explicit")}' in output
    assert f'Outside species list: {int(selection == "species_list")}' in output
    assert 'Below BUSCO threshold: 1' in output
    assert 'BUSCO threshold: 50%' in output
    assert 'Analyses (all target): OG expression, PhenoRadar inputs' in output
    assert 'Phylogeny trees:' not in output
    assert f'Build: {build}' in output and f'Config: {config}' in output
    assert f'Details: {directory / "plan.tsv"}' in output
    assert f'Excluded: {directory / "excluded.tsv"}' in output
    assert not any(name in output for name in reasons)
    assert build_report.read_text() == 'saved build plan\n'
    assert not (root / 'results').exists()


@pytest.mark.parametrize('absolute_outdir', [False, True])
def test_plan_json_keeps_original_payload_and_supports_custom_paths(plan_project, monkeypatch, capsys, tmp_path, absolute_outdir):
    root, build, _, cfg = plan_project
    config = root / 'custom-analysis.yaml'
    cfg['selection']['busco_threshold'] = 0.7
    config.write_text(yaml.safe_dump(cfg))
    directory = tmp_path / 'external-reports' if absolute_outdir else root / 'preview'
    outdir = directory if absolute_outdir else 'preview'
    expected = analysis.settings(root, config, build)[-1]

    # Resolve build from the custom config, and relative report paths from --root.
    monkeypatch.chdir(tmp_path)
    invoke(monkeypatch, analysis, 'plan', '--root', root, '--config', config.name,
           '--outdir', outdir, '--json')
    rows = json.loads(capsys.readouterr().out)
    assert rows == expected
    assert len(read_tsv(directory / 'plan.tsv')) == len(rows)
    assert len(read_tsv(directory / 'excluded.tsv')) == 2
    assert not (root / 'reports/analysis/plan').exists()


def test_plan_overwrites_reports_and_keeps_header_when_none_excluded(plan_project, monkeypatch, capsys):
    root, _, config, cfg = plan_project
    invoke(monkeypatch, analysis, 'plan', '--root', root)
    capsys.readouterr()
    directory = root / 'reports/analysis/plan'
    assert len(read_tsv(directory / 'excluded.tsv')) == 1

    cfg['selection']['busco_threshold'] = 0
    config.write_text(yaml.safe_dump(cfg))
    invoke(monkeypatch, analysis, 'plan', '--root', root)
    output = capsys.readouterr().out
    assert 'Samples: 3 selected; 0 excluded' in output
    assert all(row['selected'] == 'True' for row in read_tsv(directory / 'plan.tsv'))
    assert (directory / 'excluded.tsv').read_text() == 'species\tselected\treason\tbusco_complete_fraction\n'


@pytest.mark.parametrize('trees, extra_branches', [(['representatives'], False), (['all', 'phenotyped'], True)])
def test_plan_names_enabled_analyses_without_predicting_tree_tip_counts(plan_project, monkeypatch, capsys, trees, extra_branches):
    root, _, config, cfg = plan_project
    cfg['alignment']['enabled'] = True
    cfg['kegg']['enabled'] = True
    cfg['phylogeny']['trees'] = trees
    cfg['phylogeny']['contrast_pairs']['enabled'] = True
    for section in ('dating', 'taxonomy_check'):
        cfg['phylogeny'][section]['enabled'] = extra_branches
    config.write_text(yaml.safe_dump(cfg))

    invoke(monkeypatch, analysis, 'plan', '--root', root)
    output = capsys.readouterr().out
    assert 'Samples: 2 selected; 1 excluded' in output
    assert 'OG expression, OG alignments, KEGG, phylogeny' in output
    assert 'contrast pairs, PhenoRadar inputs' in output
    assert ('dating, taxonomy check' in output) is extra_branches
    assert 'Phylogeny trees: ' + ', '.join(trees) in output


def test_plan_with_thousands_of_samples_keeps_stdout_short(plan_project, monkeypatch, capsys):
    root, build, config, cfg = plan_project
    rows = [{'species': f'Plant_{i}_RUN{i}', 'selected': True, 'reason': '',
             'busco_complete_fraction': 0.8} for i in range(5583)]
    # Supply many selection rows without constructing thousands of input artifacts.
    monkeypatch.setattr(analysis, 'settings', lambda *args: (cfg, cfg, build, {}, [], rows))
    invoke(monkeypatch, analysis, 'plan', '--root', root)
    output = capsys.readouterr().out
    assert 'Samples: 5583 selected; 0 excluded' in output
    assert len(output.splitlines()) <= 14
    assert 'Plant_' not in output
    assert len(read_tsv(root / 'reports/analysis/plan/plan.tsv')) == 5583


@pytest.mark.parametrize('options', [(), ('--json',)])
def test_plan_still_rejects_damaged_database_products(plan_project, monkeypatch, capsys, options):
    root, build, _, _ = plan_project
    product = load_complete(build, verify_files=False)['products']['Alpha_plant_A1']['protein']
    Path(product['path']).write_text('damaged protein\n')
    with pytest.raises(ValueError, match='registered file changed'):
        invoke(monkeypatch, analysis, 'plan', '--root', root, *options)
    assert capsys.readouterr().out == ''
    assert not (root / 'reports/analysis/plan').exists()
