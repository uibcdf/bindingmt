"""Emit a frozen actual Interface 0.6 / Footprint 0.2 core-reader corpus."""

import hashlib
import json
import subprocess
from pathlib import Path


def build(destination):
    """Create new public records and a source/input/producer inventory once."""
    import topomt as tmt
    from topomt.provider_output import ProviderRun

    import bindingmt as bmt
    from bmt_test_support.pisa_case import ORIGINAL, source, verify

    destination = Path(destination)
    destination.mkdir()
    root = Path(__file__).resolve().parents[1]
    cases = {}
    for name in ['positive', 'empty', 'failed', 'unknown']:
        filename = (
            'pisa-empty.zip'
            if name == 'empty'
            else 'pisa-failed.zip'
            if name == 'failed'
            else 'pisa.zip'
        )
        path = ORIGINAL / filename
        molecular = (
            source(empty=name == 'empty') if name in ('positive', 'empty') else None
        )
        owner = (
            tmt.load_provider_output(path, molecular_system=molecular)
            if molecular is not None
            else ProviderRun.load(path)
        )
        cases[name] = bmt.BindingInterface.from_topomt(
            owner,
            question={'kind': 'association'},
            participants={
                'A': {
                    'source_id': 'frozen-original-pisa',
                    'source_atom_indices': list(range(20)),
                },
                'B': {
                    'source_id': 'frozen-original-pisa',
                    'source_atom_indices': list(range(20, 40)),
                },
            },
            grouping=[['A'], ['B']],
            artifact_ref={
                'uri': filename,
                'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            },
        )
        cases[name].save(destination / f'I-{name}.json')
    molecular = source()
    for name, participant, groups in [
        ('A', 'A', ['other-bonds']),
        ('B', 'B', ['other-bonds']),
        ('empty', 'A', ['h-bonds']),
        ('mixed', 'A', ['h-bonds', 'other-bonds']),
    ]:
        footprint = cases['positive'].project_reported_footprint(
            molecular,
            bundle_file=ORIGINAL / 'pisa.zip',
            finding_id='1',
            participant=participant,
            contact_groups=groups,
        )
        footprint.save(destination / f'F-{name}.json')
    manifest = {
        'scope': 'Actual emitted experimental records for independent core reading; no new original PISA execution.',
        'bindingmt_base': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=root, text=True
        ).strip(),
        'bindingmt_runtime_sha256': {
            str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted((root / 'bindingmt').rglob('*.py'))
        },
        'topomt_source': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=root.parent / 'topomt', text=True
        ).strip(),
        'molsysmt_source': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=root.parent / 'molsysmt', text=True
        ).strip(),
        'authoring_program_sha256': hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest(),
        'files': {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(destination.glob('*.json'))
        },
        'independent_expectations': verify()['expected'],
        'producer_state': 'Base commit plus exact actual runtime fingerprints; new source is not yet committed at production time.',
    }
    (destination / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest
