# Analysis has no producers for translation, ODB references or ODB-mapper jobs.
from analysis_build_index import BuildInputs
from common import sha256
BUILD_INDEX = f'{OUT}/build_inputs/index.json'
BUILD_INPUTS = BuildInputs(config['build_manifest'], BUILD_INDEX)
BUILD_IMPORT_SIGNATURE = [sha256(Path(SCRIPTS) / name) for name in ('build_products.py', 'analysis_build_index.py')]

rule import_build_protein:
    input:
        completion=config['build_manifest'],
        protein=lambda wc: BUILD_INPUTS.product(wc.species, 'odb_species')['protein']['path'],
        provenance=lambda wc: BUILD_INPUTS.product(wc.species, 'odb_species')['translation']['path'],
        code=ancient(f'{SCRIPTS}/build_products.py'),
        helpers=[f'{SCRIPTS}/portable_build.py', f'{SCRIPTS}/dataset_assets.py', f'{SCRIPTS}/product_storage.py']
    output:
        protein=f'{PROTEINS}/{{species}}_protein.fa.gz',
        provenance=f'{PROTEINS}/{{species}}_protein.json'
    conda: '../envs/analysis.yaml'
    params: index=BUILD_INDEX, implementation=BUILD_IMPORT_SIGNATURE
    resources: mem_mb=4000
    shell:
        '{PYTHON:q} {input.code:q} protein --completion {input.completion:q} --species {wildcards.species:q} '
        '--protein {output.protein:q} --provenance {output.provenance:q} --index {params.index:q}'

rule import_build_mapping:
    input:
        completion=config['build_manifest'],
        mapping=lambda wc: BUILD_INPUTS.data['mapping']['path'],
        samples=f'{META}/samples.tsv',
        proteins=lambda wc: sorted({f'{PROTEINS}/{r["odb_species"]}_protein.fa.gz' for r in sample_rows(wc)}),
        code=ancient(f'{SCRIPTS}/build_products.py'),
        helpers=[f'{SCRIPTS}/portable_build.py', f'{SCRIPTS}/dataset_assets.py', f'{SCRIPTS}/product_storage.py']
    output:
        snapshot=f'{MAPPING}/snapshot.json',
        tables=directory(f'{MAPPING}/species')
    params: out=MAPPING, index=BUILD_INDEX, implementation=BUILD_IMPORT_SIGNATURE
    conda: '../envs/analysis.yaml'
    resources: mem_mb=8000
    shell:
        '{PYTHON:q} {input.code:q} mapping --completion {input.completion:q} --samples {input.samples:q} '
        '--outdir {params.out:q} --index {params.index:q}'
