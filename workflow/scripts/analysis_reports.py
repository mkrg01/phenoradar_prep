"""Compact analysis previews with detailed sample selection in TSV files."""
from collections import Counter
from pathlib import Path

from common import write_tsv


PLAN_FIELDS = ['species', 'selected', 'reason', 'busco_complete_fraction']


def write_plan(rows, directory):
    directory = Path(directory)
    write_tsv(directory / 'plan.tsv', PLAN_FIELDS, rows)
    write_tsv(directory / 'excluded.tsv', PLAN_FIELDS,
              [row for row in rows if not row['selected']])


def print_plan(rows, config, build, config_path, directory):
    selected = sum(row['selected'] for row in rows)
    reasons = Counter(row['reason'] for row in rows if not row['selected'])
    print(f'Build: {build}')
    print(f'Config: {config_path}')
    print(f'Samples: {selected} selected; {len(rows) - selected} excluded')
    print(f'BUSCO threshold: {config["selection"]["busco_threshold"] * 100:g}%')
    for label, reason in [('Explicit exclusions', 'exclude_species'),
                          ('Outside species list', 'outside_species_list'),
                          ('Below BUSCO threshold', 'below_busco_threshold')]:
        print(f'  {label}: {reasons[reason]}')
    analyses = ['OG expression']
    if config['alignment']['enabled']:
        analyses.append('OG alignments')
    if config['kegg']['enabled']:
        analyses.append('KEGG')
    phylogeny = config['phylogeny']
    if phylogeny['trees']:
        analyses.append('phylogeny')
    for section, label in [('dating', 'dating'), ('taxonomy_check', 'taxonomy check'),
                           ('contrast_pairs', 'contrast pairs')]:
        if phylogeny[section]['enabled']:
            analyses.append(label)
    analyses.append('PhenoRadar inputs')
    print('Analyses (all target): ' + ', '.join(analyses))
    if phylogeny['trees']:
        print('Phylogeny trees: ' + ', '.join(phylogeny['trees']))
    print(f'Details: {Path(directory) / "plan.tsv"}')
    print(f'Excluded: {Path(directory) / "excluded.tsv"}')
