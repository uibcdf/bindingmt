"""Author the small walkthrough input through public MolSysMT tools once."""

import hashlib
import json
import subprocess
from pathlib import Path


def build(destination):
    """Create new files; never regenerate over an existing example dataset."""
    import molsysmt as msm
    import pyunitwizard as puw
    from molsysmt.native import Structures

    destination = Path(destination)
    destination.mkdir()
    root = Path(__file__).resolve().parents[1]
    parent = root / 'tests/data/interface_native_handoff/cb6af716/source.h5msm'
    parent_sha = hashlib.sha256(parent.read_bytes()).hexdigest()
    if parent_sha != '75ed110abb15515ed33806cdd03d9958eb231b3fbe2427dec0ffbc70aa287d1b':
        raise ValueError('The frozen parent source differs from the authored recipe.')
    system = msm.h5msm.read(str(parent))
    system.interactions = {}
    system.chemical_states.append_state()
    states = msm.convert(system.chemical_states, to_form='molsysmt.ChemicalStatesDict')
    states.data['states'][1]['connectivity_completeness'] = 'complete'
    system.chemical_states = msm.convert(states, to_form='molsysmt.ChemicalStates')
    msm.set(system, chemical_state=0, formal_charge=[-1, 1, 0, 0])
    msm.set(system, chemical_state=1, formal_charge=[0, 0, 1, -1])
    system.structures = Structures(
        coordinates=puw.quantity(
            [
                [[0.25, 0, 0], [0, 0, 0], [5, 0, 0], [9, 0, 0]],
                [[5, 0, 0], [9, 0, 0], [0, 0, 0], [0.25, 0, 0]],
                [[5, 0, 0], [9, 0, 0], [0, 0, 0], [1.5, 0, 0]],
            ],
            'nm',
        )
    )
    msm.set(system, structure_chemical_state_index=[0, 1, 1])
    path = destination / 'source.h5msm'
    msm.h5msm.write(system, str(path), format_version='1.0')
    manifest = {
        'case': 'bindingmt-authored-site-interface-v1',
        'scope': 'Four-ion software/input-consistency control; no biological claim.',
        'source': {
            'file': path.name,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'bytes': path.stat().st_size,
            'n_atoms': 4,
            'n_structures': 3,
            'structure_chemical_state_indices': [0, 1, 1],
            'parent_sha256': parent_sha,
        },
        'recipe': {
            'site_selection': 'formal_charge>0',
            'partner_selection': 'formal_charge<0',
            'contact_cutoff': '0.3 nm',
            'probe_radius': '0.14 nm',
            'n_sphere_points': 240,
        },
        'expected': [
            {
                'structure': 0,
                'state': 0,
                'region': [1],
                'partner': [0],
                'pairs': [[1, 0]],
                'footprint': [1],
                'burial': 'positive',
            },
            {
                'structure': 1,
                'state': 1,
                'region': [2],
                'partner': [3],
                'pairs': [[2, 3]],
                'footprint': [2],
                'burial': 'positive',
            },
            {
                'structure': 2,
                'state': 1,
                'region': [2],
                'partner': [3],
                'pairs': [],
                'footprint': [],
                'burial': 'zero',
            },
        ],
        'producer': {
            'bindingmt_base_commit': subprocess.check_output(
                ['git', 'rev-parse', 'HEAD'], cwd=root, text=True
            ).strip(),
            'molsysmt_commit': subprocess.check_output(
                ['git', 'rev-parse', 'HEAD'], cwd=root.parent / 'molsysmt', text=True
            ).strip(),
            'molsysmt_version': msm.__version__,
            'authoring_program_sha256': hashlib.sha256(
                Path(__file__).read_bytes()
            ).hexdigest(),
        },
    }
    with (destination / 'case.json').open('x') as stream:
        json.dump(manifest, stream, indent=2, allow_nan=False)
        stream.write('\n')
    return manifest
