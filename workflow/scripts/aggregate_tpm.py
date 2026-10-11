#!/usr/bin/env python3
"""Sum and renormalize OG TPM per run, rejecting genes assigned to multiple OGs."""
import argparse
import csv
import math
import json
import os
import tempfile
from collections import defaultdict
from pathlib import Path

from common import file_record, now, species_from_gene_id, write_json, write_tsv, text_open, content_sha256
from mapping_tables import read_species, load_tables, relative_file, lock
from dataset_assets import SAFE, digest, verify
from product_storage import product_record, stage_product
from sample_table import sample_table


def selected_sample(samples, run):
    rows = sample_table(samples)[2].get(run, [])
    if len(rows) != 1:
        raise ValueError(f"expected exactly one manifest row for run {run}")
    return rows[0]


def _aggregate(samples, run, mapping, output, qc):
    sample = selected_sample(samples, run)
    values = {}
    with text_open(sample["abundance"]) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"target_id", "tpm"}.issubset(reader.fieldnames or []):
            raise ValueError("abundance table requires target_id and tpm")
        for row in reader:
            gene = row["target_id"]
            if not gene or gene in values:
                raise ValueError(f"empty or duplicate target_id in {run}: {gene}")
            value = float(row["tpm"])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"invalid TPM in {run}: {gene}")
            values[gene] = value
    total = math.fsum(values.values())
    if not values or total <= 0:
        raise ValueError(f"{run}: abundance has no positive TPM")
    genes, entry = read_species(mapping, sample['species'])
    foreign = next((gene for gene in values if species_from_gene_id(gene) != sample['species']), None)
    if foreign:
        raise ValueError(f"{run}: target ID belongs to a different species: {foreign}")
    table_record = file_record(relative_file(Path(mapping).parent, entry['table'])['path'])
    mapping = {gene:groups for gene,groups in genes.items() if groups and gene in values}
    protein_genes = len(genes)
    quantified_proteins = len(genes.keys() & values.keys())
    ambiguous = [gene for gene, groups in mapping.items() if len(groups) > 1]
    if ambiguous:
        raise ValueError(f"{run}: {len(ambiguous)} genes map to multiple OGs; examples {ambiguous[:5]}")
    summed, included = defaultdict(float), set()
    for gene, groups in mapping.items():
        group, = groups
        summed[group] += values[gene]
        included.add(gene)
    matched_tpm = math.fsum(values[gene] for gene in mapping)
    included_tpm = math.fsum(values[gene] for gene in included)
    if included_tpm <= 0:
        raise ValueError(f"{run}: no positive TPM maps to retained OGs; check FASTA/query IDs")
    result = [{"species": sample["species"], "run": run, "orthogroup": group,
               "tpm_sum": value, "tpm": value / included_tpm * 1e6}
              for group, value in sorted(summed.items())]
    write_tsv(output, ["species", "run", "orthogroup", "tpm_sum", "tpm"], result)
    write_json(qc, {"created_at": now(), "species": sample["species"], "run": run, "multimap": "error",
                    "targets": len(values), "protein_genes": protein_genes, "quantified_proteins": quantified_proteins,
                    "mapped_targets": len(mapping), "ambiguous_targets": len(ambiguous), "retained_targets": len(included),
                    "total_tpm": total, "mapped_tpm": matched_tpm, "mapped_tpm_fraction": matched_tpm / total,
                    "retained_tpm": included_tpm, "retained_tpm_fraction": included_tpm / total,
                    "orthogroups": len(result), "abundance": product_record(sample["abundance"]),
                    "mapping_table": table_record, "expression": product_record(output)})


def cache_identity(sample, abundance, table):
    return digest([sample['species'], sample['run'], content_sha256(abundance), table['sha256'], 'error',
                   [file_record(Path(__file__).with_name(name))['sha256']
                    for name in ('aggregate_tpm.py', 'mapping_tables.py', 'common.py')]])


def aggregate(samples, run, mapping, output, qc, *, cache_dir=None):
    if not cache_dir:
        return _aggregate(samples, run, mapping, output, qc)
    sample = selected_sample(samples, run)
    abundance = product_record(sample['abundance'])
    entry = load_tables(mapping, verify_files=False)['tables'][sample['species']]
    table = relative_file(Path(mapping).parent, entry['table'])
    verify(table)
    identity = cache_identity(sample, abundance, table)
    cache = Path(cache_dir) / identity
    with lock(cache.with_suffix('.lock')):
        if not (cache / 'receipt.json').is_file():
            with tempfile.TemporaryDirectory(prefix='.expression-', dir=cache.parent) as tmp:
                staging = Path(tmp)
                filename = 'tpm.tsv.gz' if str(output).endswith('.gz') else 'tpm.tsv'
                _aggregate(samples, run, mapping, staging / filename, staging / 'qc.json')
                verify(abundance); verify(table)
                files = {name:dict(product_record(staging / name), path=name) for name in (filename, 'qc.json')}
                write_json(staging / 'receipt.json', {'identity':identity, 'files':files})
                os.rename(staging, cache)
        receipt = json.loads((cache / 'receipt.json').read_text())
        if receipt['identity'] != identity: raise ValueError('expression cache identity differs')
        filename = 'tpm.tsv.gz' if 'tpm.tsv.gz' in receipt['files'] else 'tpm.tsv'
        published = {}
        published['tpm'] = stage_product(dict(receipt['files'][filename], path=str(cache / filename)), output)
        report = json.loads(verify(dict(receipt['files']['qc.json'], path=str(cache / 'qc.json'))).read_text())
        report.update(abundance=abundance, expression=published['tpm'])
        write_json(qc, report)
        from dataset_assets import record
        published['qc'] = record(qc)
        return published


def aggregate_many(samples, runs, mapping, outdir, *, cache_dir=None, compressed=False):
    """Share manifest/mapping indexes and keep each run's independent checkpoint."""
    runs = list(runs)
    if not runs or len(set(runs)) != len(runs) or any(not SAFE.fullmatch(run) for run in runs):
        raise ValueError('expression batch requires unique, safe, nonempty runs')
    for run in runs:
        selected_sample(samples, run)
    load_tables(mapping, verify_files=False)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    (out / 'receipt.json').unlink(missing_ok=True)
    # Checkpoints outlive a failed Snakemake directory output even without a configured cache.
    cache_dir = cache_dir or out.parent / '.expression-cache'
    files, errors = {}, []
    for run in runs:
        try:
            suffix = '.tsv.gz' if compressed else '.tsv'
            published = aggregate(samples, run, mapping, out / (run + suffix),
                                  out / f'{run}.qc.json', cache_dir=cache_dir)
            files[run] = {key: dict(entry, path=Path(entry['path']).name)
                          for key, entry in published.items()}
        except (OSError, ValueError) as error:
            errors.append(f'{run}: {error}')
    if errors:
        raise ValueError('expression batch failed:\n' + '\n'.join(errors))
    write_json(out / 'receipt.json', {'schema_version': 1, 'kind': 'expression_batch',
                                    'runs': runs, 'files': files})


def publish_run(batch_dir, run, output, qc):
    """Expose the existing per-run paths using a small verified batch receipt."""
    batch = Path(batch_dir)
    receipt = json.loads((batch / 'receipt.json').read_text())
    if (receipt.get('schema_version') != 1 or receipt.get('kind') != 'expression_batch'
            or run not in receipt['runs']):
        raise ValueError(f'run absent from completed expression batch: {run}')
    entries = {key: relative_file(batch, entry) for key, entry in receipt['files'][run].items()}
    for key in ('tpm', 'qc'):
        verify(entries[key])
    expression = stage_product(entries['tpm'], output)
    report = json.loads(Path(entries['qc']['path']).read_text())
    report['expression'] = expression
    write_json(qc, report)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument('--run')
    selection.add_argument('--runs', nargs='+')
    for flag in ["samples", "mapping", "output", "qc", "outdir", "batch-dir"]:
        parser.add_argument(f"--{flag}")
    parser.add_argument("--cache-dir")
    parser.add_argument('--gzip', action='store_true', help='Compress batched TPM tables')
    args = parser.parse_args()
    required = ('run', 'output', 'qc') if args.batch_dir else (
        ('samples', 'mapping', 'outdir') if args.runs else ('samples', 'mapping', 'output', 'qc'))
    if any(getattr(args, key) is None for key in required):
        parser.error('required arguments: ' + ', '.join('--' + key for key in required))
    if args.batch_dir:
        publish_run(args.batch_dir, args.run, args.output, args.qc)
    elif args.runs:
        aggregate_many(args.samples, args.runs, args.mapping, args.outdir,
                       cache_dir=args.cache_dir, compressed=args.gzip)
    else:
        aggregate(args.samples, args.run, args.mapping, args.output, args.qc, cache_dir=args.cache_dir)
