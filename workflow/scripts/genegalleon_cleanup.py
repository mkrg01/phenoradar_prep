"""Retire native work only after all three sample products are durable."""
import json
import os
from pathlib import Path

from cleanup_work import owned_path
from common import sha256
from dataset_assets import stat_identity, verify


POLICY_VERSION = 1
DIRECTORIES = (
    'downloads/busco_downloads', 'downloads/ete_taxonomy',
    'output/transcriptome_assembly/assembled_transcripts_with_isoforms',
    'output/transcriptome_assembly/longest_cds_transcript',
    'output/transcriptome_assembly/longest_cds',
    'output/transcriptome_assembly/amalgkit_quant',
    'output/transcriptome_assembly/amalgkit_merge',
)


def completed_targets(build, work, manifest, item, products, *, apply):
    """Return owned targets and retention reasons; no deletion happens here."""
    if not all(products.get(key) for key in ('reference', 'busco', 'quant')):
        return None
    receipt = Path(build) / 'jobs/status' / f"{item['species']}.sample.json"
    if not receipt.is_file() or json.loads(receipt.read_text()).get('state') not in {'complete', 'reused'}:
        return None
    # Recheck the actual registered bytes, including legacy native-path products.
    entries = [products['reference']['cds'], products['quant']['abundance']]
    entries.extend(products['busco'][key] for key in ('full', 'short') if products['busco'].get(key))
    protected = [verify(entry).resolve() for entry in entries]
    raw = manifest.get('raw_inputs', {}).get(item['species'], {})
    protected.extend(Path(entry['path']).resolve() for entry in raw.values())
    for key in ('busco_reference', 'taxonomy_reference'):
        reference = manifest.get(key)
        if reference:
            protected.extend((Path(reference['dataset']).resolve(), Path(reference['receipt']['path']).resolve()))
    targets, retained = [], []

    def add(relative):
        path = owned_path(build, work.relative_to(build) / relative)
        if not path.exists() and not path.is_symlink():
            return
        if any(p == path or p.is_relative_to(path) for p in protected):
            retained.append({'path': str(path.relative_to(build)), 'reason': 'registered product or original input'})
        else:
            targets.append(path.relative_to(build))

    for relative in DIRECTORIES:
        add(Path(relative))

    # A staged private read may be a hard link or a cross-filesystem copy.
    # Never retire the only available copy or a caller-owned input tree.
    private = []
    for key, entry in raw.items():
        original = Path(entry['path'])
        relative = Path('input/reads') / item['species'] / (key + ''.join(original.suffixes))
        staged = owned_path(build, work.relative_to(build) / relative)
        if not staged.exists() and not staged.is_symlink():
            continue
        reason = None
        if staged.is_symlink() or not staged.is_file():
            reason = 'staged read is not an owned regular file'
        elif not original.is_file():
            reason = 'original input is unavailable'
        elif original.resolve() == staged or original.resolve().is_relative_to(staged.parent):
            reason = 'original input is inside the staged read directory'
        elif any(p == staged for p in protected[:len(entries)]):
            reason = 'registered product'
        elif not os.path.samefile(original, staged):
            before = stat_identity(original), stat_identity(staged)
            if before[0][2] != before[1][2]:
                reason = 'original input differs from staged copy'
            elif apply:
                if sha256(original) != sha256(staged):
                    reason = 'original input differs from staged copy'
                if before != (stat_identity(original), stat_identity(staged)):
                    raise ValueError('private input changed during cleanup verification')
        if reason:
            retained.append({'path': str(staged.relative_to(build)), 'reason': reason})
        else:
            private.append(staged)
    if private:
        directory = private[0].parent
        if {p for p in directory.iterdir()} == set(private):
            targets.append(directory.relative_to(build))
        else:
            targets.extend(p.relative_to(build) for p in private)
    return targets, retained
