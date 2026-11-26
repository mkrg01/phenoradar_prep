#!/usr/bin/env python3
"""Map a reference protein FASTA to OrthoDB, preserving its gene identifiers."""

import argparse
import json
from pathlib import Path

from common import atomic_writer, file_record, write_json
from mapping_tables import annotation_pairs, lock, write_gene_mapping
from run_odb_chunk import run as run_odb
from translate_cds import fasta_ids


def map_reference_proteins(proteins, reference, output_dir, jobs=1, command='ODB-mapper', prefix=''):
    protein, reference, output = (Path(path).resolve() for path in (proteins, reference, output_dir))
    if jobs < 1:
        raise ValueError('jobs must be at least 1')
    label = 'chunk_001'
    with lock(output / '.mapping.lock'):
        inputs = {'protein': file_record(protein), 'reference': file_record(reference)}
        identity = output / 'inputs.json'
        if identity.is_file():
            if json.loads(identity.read_text()) != inputs:
                raise ValueError('mapping inputs changed; use a new output directory')
        genes = set(fasta_ids(protein))
        ref = json.loads(reference.read_text())
        write_json(identity, inputs)
        manifest = output / 'proteins.fs'
        with atomic_writer(manifest) as handle:
            handle.write(f'{protein}\n')

        odb_dir = output / 'odb'
        mapped = run_odb(manifest, reference, odb_dir, output / '.work', label,
                         command=command, prefix=prefix, version=ref['version'], node=ref['node'],
                         jobs=jobs, batch_size=jobs, gzip_results=True)
        annotation = odb_dir / f'{label}.og.annotations.gz'
        pairs = {}
        rows = 0
        for gene, og in annotation_pairs(annotation):
            if gene not in genes:
                raise ValueError(f'ODB annotation query is absent from the protein FASTA: {gene}')
            pairs.setdefault(gene, set()).add(og)
            rows += 1
        if inputs != {'protein': file_record(protein), 'reference': file_record(reference)}:
            raise ValueError('protein FASTA or reference changed during mapping')

        table = output / 'gene_to_orthogroup.tsv.gz'
        with atomic_writer(table, 'wb') as handle:
            write_gene_mapping(handle, genes, pairs)
        unique = sum(map(len, pairs.values()))
        result = {
            'schema_version': 1,
            **inputs,
            'orthodb': {'version': ref['version'], 'node': ref['node'],
                        'fingerprint': mapped['fingerprint'],
                        'provenance': file_record(odb_dir / 'provenance.json')},
            'table': file_record(table),
            'qc': {'protein_genes': len(genes), 'mapped_genes': len(pairs),
                   'unmapped_genes': len(genes) - len(pairs),
                   'ambiguous_genes': sum(len(groups) > 1 for groups in pairs.values()),
                   'unique_gene_og_pairs': unique, 'duplicate_pairs_removed': rows - unique},
        }
        write_json(output / 'mapping_provenance.json', result)
        return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proteins', required=True, help='Protein FASTA (.gz accepted) with unique gene IDs')
    parser.add_argument('--reference', required=True, help='Prepared OrthoDB reference.json')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--jobs', type=int, default=1)
    parser.add_argument('--command', default='ODB-mapper')
    parser.add_argument('--prefix', default='', help='Optional ODB-mapper environment prefix')
    map_reference_proteins(**vars(parser.parse_args()))
