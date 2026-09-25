# Analysis has no producers for translation, ODB references or ODB-mapper jobs.
from build_products import load_complete
COMPLETED_BUILD = load_complete(config['build_manifest'], verify_files=False)
BUILD_PROTEINS = {p['odb_species']:p for p in COMPLETED_BUILD['products'].values()}

rule import_build_protein:
    input:
        completion=config['build_manifest'],
        protein=lambda wc: BUILD_PROTEINS[wc.species]['protein']['path'],
        provenance=lambda wc: BUILD_PROTEINS[wc.species]['translation']['path'],
        code=f'{SCRIPTS}/build_products.py',
        helpers=[f'{SCRIPTS}/portable_build.py', f'{SCRIPTS}/dataset_assets.py']
    output:
        protein=f'{PROTEINS}/{{species}}_protein.fa',
        provenance=f'{PROTEINS}/{{species}}_protein.json'
    conda: '../envs/analysis.yaml'
    resources: mem_mb=1000
    shell:
        '{PYTHON:q} {input.code:q} protein --completion {input.completion:q} --species {wildcards.species:q} '
        '--protein {output.protein:q} --provenance {output.provenance:q}'

rule import_build_mapping:
    input:
        completion=config['build_manifest'],
        mapping=COMPLETED_BUILD['mapping']['path'],
        samples=f'{META}/samples.tsv',
        proteins=lambda wc: sorted({f'{PROTEINS}/{r["odb_species"]}_protein.fa' for r in sample_rows(wc)}),
        code=f'{SCRIPTS}/build_products.py',
        helpers=[f'{SCRIPTS}/portable_build.py', f'{SCRIPTS}/dataset_assets.py']
    output:
        snapshot=f'{MAPPING}/snapshot.json',
        tables=directory(f'{MAPPING}/species')
    params: out=MAPPING
    conda: '../envs/analysis.yaml'
    resources: mem_mb=8000
    shell:
        '{PYTHON:q} {input.code:q} mapping --completion {input.completion:q} --samples {input.samples:q} '
        '--outdir {params.out:q}'
