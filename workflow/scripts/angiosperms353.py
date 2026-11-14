"""Screen database samples for discordance with Kew Angiosperms353 references."""
from collections import Counter, defaultdict
from contextlib import contextmanager
import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

from angiosperms353_reference import INTERMEDIATE_RANKS, LEGACY_RANKS, LOCUS, RANKS, accessions, fasta, reference_rows
from common import now, read_tsv, write_json, write_tsv
from dataset_assets import record, verify
from phase_config import DEFAULT_REGISTERED_SUPPORT

MIN_BASES = 150
COVERAGE = '0.70'
EVALUE = '1e-10'
CIGAR = re.compile(r'([0-9]+)([MID=X])')
COMPLEMENT = str.maketrans('ACGT', 'TGCA')
SUMMARY_FIELDS = ['run', 'taxid', 'resolved_taxid', 'resolved_species_taxid', 'resolved_species_name',
                  'rank', 'registered', 'status', 'extracted_loci',
                  'comparable_loci', 'registered_votes', 'tied_loci', 'registered_name', 'alternative', 'alternative_name',
                  'alternative_votes', 'support', 'alternative_counts', 'alternative_loci',
                  'voted_loci', 'alternative_votes_total', 'registered_support', 'discordance_support']
LOCUS_FIELDS = ['run', 'rank', 'locus', 'vote', 'winner', 'registered_reference',
                'alternative_reference', 'shared_bases', 'registered_matches', 'alternative_matches',
                'alternative_winners', 'alternative_candidates']
PLOT_NAMES = ('registered_support.svg', 'registered_support_by_rank.svg')


def registered_threshold(min_loci, registered_support=None, support=None):
    if type(min_loci) is not int or min_loci < 1:
        raise ValueError('minimum loci must be positive')
    if registered_support is not None and support is not None:
        raise ValueError('specify either registered support or alternative support')
    try:
        if support is not None:
            alternative = Decimal(str(support))
            cutoff = Decimal(1) - alternative if alternative.is_finite() else alternative
        else:
            cutoff = Decimal(str(DEFAULT_REGISTERED_SUPPORT if registered_support is None else registered_support))
    except InvalidOperation as error:
        raise ValueError('support must be numeric') from error
    if not cutoff.is_finite() or not 0 <= cutoff <= 1:
        raise ValueError('registered support must be between zero and one')
    return cutoff


def load_panel(path):
    path = Path(path).resolve()
    rows = read_tsv(path)
    if (not rows or not {'run', 'taxid', 'cds'} <= rows[0].keys()
            or any(None in row or any(value is None for value in row.values()) or not row['cds'] for row in rows)):
        raise ValueError('sample TSV requires run, taxid and cds columns')
    return [dict(row, cds=record(path.parent / row['cds'])) for row in rows]


def build_panel(path, runs=None):
    from dataset import item_products, load
    from sample_review import sample_guard

    path = Path(path).resolve()
    if path.name == 'build.json':
        path = path.parent
    build = load(path)
    items = {item['row']['run']: item for item in build['items']}
    requested = set(items if runs is None else runs)
    if not requested or requested - items.keys():
        raise ValueError(f'no samples selected or runs absent from build: {sorted(requested - items.keys())}')
    samples = []
    for run in sorted(requested):
        item = items[run]
        with sample_guard(path, build, item):
            reference = item_products(build, item).get('reference')
        if reference is None:
            raise ValueError(f'{run}: no completed CDS assembly; select runs with --runs')
        samples.append(dict(item['row'], cds=reference['cds']))
    return samples, record(path / 'build.json')


def database_panel(path, runs=None):
    from build_products import load_complete
    from portable_build import completion_path

    manifest = record(completion_path(path))
    database = load_complete(manifest['path'], verify_files=False)
    products = {product['row']['run']: product for product in database['products'].values()}
    if len(products) != len(database['products']):
        raise ValueError('database contains duplicate run accessions')
    requested = set(products if runs is None else runs)
    if not requested or requested - products.keys():
        raise ValueError(f'no samples selected or runs absent from database: {sorted(requested - products.keys())}')
    verify(manifest)
    return [dict(products[run]['row'], cds=products[run]['cds']) for run in sorted(requested)], manifest


def extract(samples, work, *, threads, ram, blat, log_dir):
    inputs = work / 'cds'
    inputs.mkdir()
    by_source, names = {}, []
    for row in samples:
        source = verify(row['cds'])
        if source not in by_source:
            name = by_source[source] = f'sample{len(by_source)}'
            suffix = '.fna.gz' if source.suffix == '.gz' else '.fna'
            (inputs / f'{name}{suffix}').symlink_to(source)
        names.append(by_source[source])
    output = work / 'captus'
    command = ['captus', 'extract', '-f', str(inputs), '-a', str(work / 'assemblies'),
               '-o', str(output), '-n', 'Angiosperms353', '--disable_stitching', '--ignore_depth',
               '--max_paralogs', '-1', '--max_loci_scipio_x2', '0', '--threads', str(threads), '--ram', str(ram), '--show_less']
    command.extend(['--blat_path', blat])
    with (log_dir / 'captus.log').open('wb') as log:
        try:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True,
                           env=dict(os.environ, TMPDIR=str(work)))
        except subprocess.CalledProcessError as error:
            raise ValueError(f'Captus failed (exit {error.returncode}); see {log.name}') from error
    paths = [output / f'{name}__captus-ext/01_coding_NUC/NUC_coding_NT.fna' for name in names]
    # Captus 1.5.1 can exit zero after a failed BLAT/Scipio subprocess.
    for path in output.glob('*__captus-ext/01_coding_NUC/**/*.log'):
        check_extraction_log(path)
    missing = [path for path in set(paths) if not path.parent.is_dir()]
    if missing:
        raise ValueError(f'Captus did not process {len(missing)} assemblies; inspect {log_dir / "captus.log"}')
    return paths


def check_extraction_log(path):
    with path.open(errors='replace') as handle:
        for line in handle:
            if re.search(r'assertion failed|segmentation fault|dyld\[|cannot execute|can.t locate .*\.pm|out of memory|killed by signal', line, re.I):
                raise ValueError(f'Captus dependency failed in {path}: {line.strip()}')


def markers(path):
    if not path.is_file():
        return {}
    candidates = defaultdict(set)
    for header, sequence in fasta(path):
        parts = header.split()[0].split('__')
        if len(parts) not in {2, 3} or not LOCUS.fullmatch(parts[1]):
            raise ValueError(f'unexpected Captus locus header: {header}')
        candidates[parts[1]].add(sequence)
    return {locus: sequence for locus, copies in candidates.items() if len(copies) == 1
            for sequence in copies if len(sequence) >= MIN_BASES}


@dataclass(frozen=True)
class Hit:
    target: int
    score: int
    start: int
    end: int
    query: str
    reference: str

    def bases(self):
        position = self.start - 1
        step = 1 if self.start <= self.end else -1
        bases = {}
        for query, reference in zip(self.query, self.reference, strict=True):
            if query != '-':
                if query in 'ACGT' and reference in 'ACGT':
                    bases[position] = query == reference
                position += step
        if position != self.end - 1 + step:
            raise ValueError('MMseqs2 query alignment disagrees with query coordinates')
        return bases


def source_label(ref):
    return f'{ref[0]}|{ref[1]}|{ref[2]}'


def vote(run, rank, locus, length, registered, hits, refs):
    original = hits.get(registered)
    row = dict(run=run, rank=rank, locus=locus, vote='no_registered_hit')
    if original is None:
        return row
    row['registered_reference'] = source_label(refs[original.target])
    other_hits = [(group, hit) for group, hit in hits.items() if group != registered]
    if not other_hits:
        return dict(row, vote='no_alternative_hit')
    best_score = max(hit.score for _, hit in other_hits)
    alternatives = [(group, hit) for group, hit in other_hits if hit.score == best_score]
    alternatives.sort()
    original_bases = original.bases()
    alternative_bases = {group: hit.bases() for group, hit in alternatives}
    shared = set(original_bases)
    for bases in alternative_bases.values():
        shared.intersection_update(bases)
    if len(shared) < max(MIN_BASES, (7 * length + 9) // 10):
        return dict(row, vote='insufficient_overlap')

    old = sum(original_bases[i] for i in shared)
    candidates = {group: dict(reference=source_label(refs[hit.target]),
                              matches=sum(alternative_bases[group][i] for i in shared))
                  for group, hit in alternatives}
    new = max(candidate['matches'] for candidate in candidates.values())
    best = [group for group, candidate in candidates.items() if candidate['matches'] == new]
    row.update(alternative_reference=candidates[best[0]]['reference'], shared_bases=len(shared),
               registered_matches=old, alternative_matches=new,
               alternative_candidates=json.dumps(candidates, sort_keys=True))
    if old > new:
        return dict(row, vote='registered', winner=registered)
    if old < new:
        return dict(row, vote='alternative', winner=best[0] if len(best) == 1 else None,
                    alternative_winners=best)
    return dict(row, vote='tie', alternative_winners=best)


@contextmanager
def search_locus(locus, queries, reference_path, assignments, work, threads, ram, log_dir):
    """Search the complete same-locus database without a fixed top-k cap."""
    with tempfile.TemporaryDirectory(prefix=f'{locus}-', dir=work) as temp:
        path = Path(temp)
        refs = []
        with (path / 'target.fna').open('w') as handle:
            for species, repository, accession, ranks, sequence in reference_rows(reference_path, locus, assignments):
                handle.write(f'>{len(refs)}\n{sequence}\n')
                refs.append((species, repository, accession, ranks, accessions(accession), sequence))
        if not refs:
            raise ValueError(f'Kew locus {locus} has no taxonomically resolved reference sequences')
        with (path / 'query.fna').open('w') as handle:
            for i, sequence in queries.items():
                handle.write(f'>{i}\n{sequence}\n')
        with (log_dir / 'mmseqs.log').open('ab') as log:
            log.write(f'\nAngiosperms353 locus {locus}: {len(queries)} queries, {len(refs)} references\n'.encode())
            log.flush()
            subprocess.run(['mmseqs', 'easy-search', str(path / 'query.fna'), str(path / 'target.fna'),
                            str(path / 'hits.tsv'), str(path / 'tmp'), '--search-type', '3',
                            '-k', '10', '--max-seqs', str(len(refs)), '--strand', '2',
                            '--format-output', 'query,target,raw,qstart,qend,tstart,tend,cigar',
                            '-e', EVALUE, '-c', COVERAGE, '--cov-mode', '2', '--threads', str(threads),
                            '--shuffle', '0', '--split-memory-limit', f'{ram}G', '-v', '2'],
                           stdout=log, stderr=subprocess.STDOUT, check=True, env=dict(os.environ, TMPDIR=str(path)))
        yield refs, path / 'hits.tsv'


def alignment_from_cigar(query, target, qstart, qend, tstart, tend, cigar):
    """Decode the nucleotide traceback against the original FASTA records."""
    qstep = 1 if qstart <= qend else -1
    tstep = 1 if tstart <= tend else -1
    query = query[min(qstart, qend) - 1:max(qstart, qend)]
    target = target[min(tstart, tend) - 1:max(tstart, tend)]
    if qstep < 0:
        query = query.translate(COMPLEMENT)[::-1]
    if tstep < 0:
        target = target.translate(COMPLEMENT)[::-1]
    operations = CIGAR.findall(cigar)
    if not operations or ''.join(count + op for count, op in operations) != cigar:
        raise ValueError(f'invalid MMseqs2 nucleotide CIGAR: {cigar!r}')
    qpos = tpos = 0
    qaln, taln = [], []
    for count, op in operations:
        size = int(count)
        if op == 'D':
            qaln.append('-' * size)
        else:
            qaln.append(query[qpos:qpos + size])
            qpos += size
        if op == 'I':
            taln.append('-' * size)
        else:
            taln.append(target[tpos:tpos + size])
            tpos += size
    if (qpos != abs(qend - qstart) + 1 or tpos != abs(tend - tstart) + 1
            or qpos != len(query) or tpos != len(target)):
        raise ValueError('MMseqs2 nucleotide CIGAR disagrees with alignment coordinates')
    return ''.join(qaln), ''.join(taln)


def reference_allowed(sample, reference, omit_same_species):
    species, _, accession, _, source_runs, _ = reference
    return not (source_runs & sample['source_runs'] or accession.upper() in sample['exclude_accessions']
                or (omit_same_species and species.replace('_', ' ').casefold() == sample['species_name']))


def closest_reference_rank(lineage, refs, sample, omit_same_species):
    """Find the nearest named ancestor with an unmasked reference sequence."""
    wanted = set(lineage.values())
    if not wanted:
        return None
    closest = next(iter(lineage.values()))
    present = set()
    for reference in refs:
        if reference_allowed(sample, reference, omit_same_species):
            present.update(wanted.intersection(reference[3]))
            if closest in present:
                return next(iter(lineage))
    return next((rank for rank, taxid in lineage.items() if taxid in present), None)


def read_hits(path, queries, refs, samples, *, omit_same_species, ranks_by_query=None):
    # Keep only the best alignment per taxon for each query. This bounds memory
    # independently of the number of specimens in well-sampled taxa.
    kept = defaultdict(lambda: {rank: {} for rank in RANKS})
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            q, t, score, start, end, tstart, tend, cigar = line.rstrip('\r\n').split('\t')
            index, target_id = int(q), int(t)
            ref = refs[target_id]
            if index not in queries:
                raise ValueError(f'unexpected MMseqs2 query index: {index}')
            if not reference_allowed(samples[index], ref, omit_same_species):
                continue
            score, start, end = int(score), int(start), int(end)
            hit = None
            considered = RANKS if ranks_by_query is None else ranks_by_query[index]
            for rank in considered:
                group = ref[3][RANKS.index(rank)]
                if group is None:
                    continue
                groups = kept[index][rank]
                previous = groups.get(group)
                if previous is None or (score, -target_id) > (previous.score, -previous.target):
                    if hit is None:
                        qaln, taln = alignment_from_cigar(queries[index], ref[5], start, end,
                                                           int(tstart), int(tend), cigar)
                        hit = Hit(target_id, score, start, end, qaln, taln)
                    groups[group] = hit
    return kept


def summarize(sample, rank, registered, rows, extracted, min_loci, support):
    alternatives = Counter()
    registered_votes = Fraction(0)
    alternative_loci = ties = voted = 0
    for row in rows:
        if row['vote'] == 'registered':
            registered_votes += 1
            voted += 1
        elif row['vote'] in {'alternative', 'tie'}:
            winners = row['alternative_winners']
            registered_tied = row['vote'] == 'tie'
            share = Fraction(1, len(winners) + registered_tied)
            for winner in winners:
                alternatives[winner] += share
            if registered_tied:
                registered_votes += share
                ties += 1
            alternative_loci += 1
            voted += 1
    best, votes = min(alternatives.items(), key=lambda item: (-item[1], item[0]), default=(None, Fraction(0)))
    status = decision_status(rank, registered_votes, voted, min_loci, support)
    # Decide using exact fractions before formatting report values.
    def number(weight):
        return int(weight) if weight.denominator == 1 else float(weight)

    alternative_votes_total = voted - registered_votes
    registered_support = Decimal(registered_votes.numerator) / registered_votes.denominator / voted if voted else None
    return dict(run=sample['run'], taxid=sample['taxid'], rank=rank, registered=registered,
                status=status, extracted_loci=extracted, comparable_loci=voted,
                registered_votes=number(registered_votes), tied_loci=ties,
                alternative=best, alternative_votes=number(votes),
                support=str(Decimal(votes.numerator) / votes.denominator / voted) if voted else None,
                alternative_counts=json.dumps({group: number(weight) for group, weight in alternatives.items()},
                                              sort_keys=True),
                alternative_loci=alternative_loci, voted_loci=voted,
                alternative_votes_total=number(alternative_votes_total),
                registered_support=str(registered_support) if voted else None,
                discordance_support=str(1 - registered_support) if voted else None)


def decision_status(rank, registered_votes, voted, min_loci, cutoff):
    if not rank:
        return 'no_reference'
    if voted < min_loci:
        return 'insufficient_comparison'
    if registered_votes < voted * Fraction(cutoff):
        return 'discordant'
    return 'no_discordance'


def exclusion_reason(row):
    rank = row['rank']
    return f"angiosperms353_{rank}_discordant_{row['registered']}_outside_registered_{rank}"


def plot_support(output, summaries, cutoff):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    def histogram(ax, rows, title):
        assessed = [row for row in rows if row['status'] in {'no_discordance', 'discordant'}]
        scores = [float(row['registered_support']) * 100 for row in assessed]
        ax.hist(scores, bins=range(0, 101, 5), color='grey')
        percent = cutoff * 100
        ax.axvline(float(percent), color='red', linestyle='--', label=f'Threshold: {percent.normalize():f}%')
        flagged = sum(row['status'] == 'discordant' for row in assessed)
        ax.set(xlim=(0, 100), xticks=range(0, 101, 20), ylabel='Number of samples',
               title=f'{title} · {len(assessed)} assessed, {flagged} flagged, {len(rows) - len(assessed)} unassessed')
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(loc='upper center')

    fig, ax = plt.subplots(figsize=(7, 4))
    histogram(ax, summaries, 'All ranks')
    ax.set_xlabel('Registered taxon support (%)')
    fig.tight_layout()
    fig.savefig(output / PLOT_NAMES[0])
    plt.close(fig)

    by_rank = defaultdict(list)
    for row in summaries:
        by_rank[row['rank'] or 'Unassigned'].append(row)
    fig, axes = plt.subplots(len(by_rank), 1, sharex=True, squeeze=False, figsize=(7, 3 * len(by_rank)))
    for ax, (rank, rows) in zip(axes[:, 0], sorted(by_rank.items()), strict=True):
        histogram(ax, rows, rank)
    axes[-1, 0].set_xlabel('Registered taxon support (%)')
    fig.tight_layout()
    fig.savefig(output / PLOT_NAMES[1])
    plt.close(fig)


def write_report(output, summaries, cutoff):
    plot_support(output, summaries, cutoff)
    write_tsv(output / 'summary.tsv', SUMMARY_FIELDS, summaries)
    write_tsv(output / 'exclusion_candidates.tsv', ['accession', 'taxid', 'reason'],
              [dict(accession=row['run'], taxid=row['taxid'], reason=exclusion_reason(row))
               for row in summaries if row['status'] == 'discordant'])
    return dict(summary=record(output / 'summary.tsv'), plots={name: record(output / name) for name in PLOT_NAMES},
                exclusion_candidates=record(output / 'exclusion_candidates.tsv'))


def report(screen_path, output, *, registered_support, selection_config=None):
    """Reapply a threshold to frozen votes without running sequence searches."""
    source = record(screen_path)
    data = json.loads(verify(source).read_text())
    if data.get('kind') != 'angiosperms353_qc' or data.get('schema_version') != 3:
        raise ValueError('replotting requires a version 3 Angiosperms353 screen')
    min_loci = data['settings']['min_loci']
    cutoff = registered_threshold(min_loci, registered_support)
    summaries = read_tsv(verify(data['summary']))
    if not summaries:
        raise ValueError('Angiosperms353 screen contains no sample summaries')
    votes = Counter()
    with verify(data['loci']).open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['vote'] == 'registered':
                votes[row['run']] += 1
            elif row['vote'] == 'tie':
                votes[row['run']] += Fraction(1, 1 + len(json.loads(row['alternative_winners'])))
    for row in summaries:
        row['status'] = decision_status(row['rank'], votes[row['run']], int(row['comparable_loci']), min_loci, cutoff)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    files = write_report(output, summaries, cutoff)
    verify(source)
    if selection_config:
        verify(selection_config)
    updated = dict(data, **files, created_at=now(), source_screen=source, selection_config=selection_config,
                   settings=dict(data['settings'], max_registered_support=str(cutoff)))
    write_json(output / 'screen.json', updated)
    return dict(output=str(output), samples=len(summaries), statuses=dict(Counter(row['status'] for row in summaries)))


def score_batch(batch, ancestry, by_locus, paths, assignments, work, log_dir, threads, ram,
                omit_same_species, min_loci, support, evidence):
    rows, extracted, selected = defaultdict(list), Counter(), {}
    ranks = [tuple(lineage) for lineage in ancestry]
    for locus, path in paths.items():
        queries = by_locus[locus]
        with search_locus(locus, queries, path, assignments, work, threads, ram, log_dir) as (refs, hits_path):
            ranks_by_query = {}
            for index in queries:
                current = ranks[index]
                if current and selected.get(index) != current[0]:
                    rank = closest_reference_rank(ancestry[index], refs, batch[index], omit_same_species)
                    if rank in current:
                        selected[index] = rank
                        ranks[index] = current[:current.index(rank) + 1]
                # Keep lower-rank evidence until the closest available rank is known.
                ranks_by_query[index] = ranks[index]
            hits = read_hits(hits_path, queries, refs, batch, omit_same_species=omit_same_species,
                             ranks_by_query=ranks_by_query)
            for index, sequence in queries.items():
                extracted[index] += 1
                for rank in ranks_by_query[index]:
                    result = vote(batch[index]['run'], rank, locus, len(sequence), ancestry[index][rank],
                                  hits.get(index, {}).get(rank, {}), refs)
                    rows[index, rank].append(result)
    summaries = []
    for index, sample in enumerate(batch):
        rank = selected.get(index)
        results = rows[index, rank]
        evidence.writerows(results)
        summaries.append(summarize(sample, rank, ancestry[index].get(rank), results,
                                   extracted[index], min_loci, support))
    return summaries


def screen(samples, **kwargs):
    from taxonomy_screen_jobs import screen as run_screen
    return run_screen(samples, **kwargs)


def load_screen(path, build_path):
    path = Path(path).resolve()
    data = json.loads(path.read_text())
    version = data.get('schema_version')
    if data.get('kind') != 'angiosperms353_qc' or version not in (1, 2, 3) or not data.get('build_manifest'):
        raise ValueError('review requires a completed Angiosperms353 screen of a frozen build')
    if Path(verify(data['build_manifest'])).resolve() != (Path(build_path) / 'build.json').resolve():
        raise ValueError('Angiosperms353 screen belongs to a different build')
    if data['settings']['restricted_loci'] or data['settings']['omit_same_species']:
        raise ValueError('validation-only Angiosperms353 screens cannot be used for sample review')
    expected_ranks = {1: set(LEGACY_RANKS), 2: set(INTERMEDIATE_RANKS), 3: {'', *RANKS}}[version]
    summaries = defaultdict(dict)
    for row in read_tsv(verify(data['summary'])):
        if row['rank'] not in expected_ranks or row['rank'] in summaries[row['run']]:
            raise ValueError('duplicate or invalid rank in Angiosperms353 screen')
        summaries[row['run']][row['rank']] = row
    if set(summaries) != set(data['samples']) or any(
            (len(ranks) != 1 if version == 3 else set(ranks) != expected_ranks) for ranks in summaries.values()):
        raise ValueError('incomplete Angiosperms353 screen summary')
    return record(path), data['samples'], summaries
