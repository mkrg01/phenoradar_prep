"""Sample-level generating conditions and a publication-friendly processing table."""
import copy
from pathlib import Path

from common import write_tsv


FIELDS = ['sample_id', 'scientific_name', 'run_accession', 'genegalleon_version',
          'assembly_method', 'busco_lineage', 'busco_complete_percent',
          'phenoradar_prep_version', 'orthodb_release', 'orthodb_clade']
PRODUCT_KEYS = ('reference', 'busco', 'quant')


def software_version(root):
    path = Path(root) / 'VERSION'
    return path.read_text().strip() if path.is_file() else 'unknown'


def scientific_conditions(cfg):
    from dataset import FIXED_GENEGALLEON_SETTINGS
    return {'settings': {**cfg['genegalleon'].get('settings', {}), **FIXED_GENEGALLEON_SETTINGS},
            'translation': cfg['translation'], 'busco_lineage': cfg['busco']['lineage']}


def generating_conditions(manifest):
    from dataset_assets import digest
    cfg = manifest['config']
    gg = cfg['genegalleon']
    lock = manifest.get('software_lock') or {}
    source = lock.get('source', {})
    identity = source.get('identity', {})
    downloaded = source.get('kind') == 'downloaded_source'
    version = source.get('version') or identity.get('version') or (
        gg.get('version') if downloaded else None) or 'unknown'
    revision = identity.get('revision') or (gg.get('revision') if downloaded else None)
    # A local checkout need not be the release named in the build configuration.
    # Keep its content identity in the same display column as the version.
    repository = Path(gg.get('repository') or '.')
    source_digest = digest([(str(Path(e['path']).relative_to(repository)), e['sha256'])
                            for e in source.get('files', [])]) if source.get('files') else None
    if source.get('kind') == 'local_source':
        version += '+local.' + source_digest[:12] if source_digest else '+local'
    return {'genegalleon_version': version, 'revision': revision,
            'image_uri': lock.get('container', {}).get('identity', {}).get('uri'),
            'source_sha256': source_digest,
            'container_sha256': [e['sha256'] for e in lock.get('container', {}).get('files', [])],
            'compatibility': scientific_conditions(cfg),
            'assembly_method': gg.get('settings', {}).get('assembly_method', 'unknown')}


def native_processing(products, manifest=None):
    """Never infer a generating release from the configuration of a reuse build."""
    present = [products[key] for key in PRODUCT_KEYS if products.get(key)]
    recorded = [p.get('provenance', {}).get('processing') for p in present]
    if recorded and all(recorded):
        if any(value != recorded[0] for value in recorded[1:]):
            raise ValueError('sample contains different GeneGalleon generating conditions; rebuild the sample')
        return copy.deepcopy(recorded[0])
    # Historical native receipts prove their producer only when every stage
    # agrees with the source build's actual software receipts and conditions.
    if manifest is not None and present and not any(recorded):
        conditions = manifest['config'].get('conditions', {})
        stages = [('assembly', 'reference'), ('busco', 'busco'), ('quant', 'quant')]
        if manifest.get('genegalleon') and all(
            p.get('provenance', {}).get('source') == 'genegalleon'
            and p['provenance'].get('software') == manifest['genegalleon']
            and p['provenance'].get('condition') == conditions.get(stage)
            for stage, key in stages if (p := products.get(key))):
            return generating_conditions(manifest)
    return None


def completed_processing(product, source_config=None):
    processing = copy.deepcopy(product.get('processing'))
    if processing is None and source_config is not None and any(product.get('conditions', {}).values()):
        # This proves the accepted scientific settings, not the generating release.
        processing = {'genegalleon_version': 'unknown', 'assembly_method': 'unknown',
                      'compatibility': scientific_conditions(source_config)}
    return processing


def write_summary(path, data):
    rows = []
    for name, product in sorted(data['products'].items()):
        processing = product.get('processing') or {}
        counts = product['counts']
        total = counts['busco_cds_total']
        complete = counts['busco_cds_single'] + counts['busco_cds_duplicated']
        node = data['odb']['node']
        rows.append({'sample_id': name, 'scientific_name': product['row']['scientific_name'],
                     'run_accession': product['row']['run'],
                     'genegalleon_version': processing.get('genegalleon_version', 'unknown'),
                     'assembly_method': processing.get('assembly_method', 'unknown'),
                     'busco_lineage': data['lineage'],
                     'busco_complete_percent': f'{100 * complete / total:.2f}' if total else 'unknown',
                     'phenoradar_prep_version': product.get('phenoradar_prep_version', 'unknown'),
                     'orthodb_release': data['odb']['version'],
                     'orthodb_clade': {3193: 'Embryophyta', 33090: 'Viridiplantae'}.get(
                         node, f'NCBI Taxonomy ID {node}')})
    write_tsv(path, FIELDS, rows)
