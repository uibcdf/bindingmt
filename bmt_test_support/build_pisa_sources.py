"""Freeze native copies of the offered PISA sources through public MolSysMT."""

import hashlib
import json
import subprocess
from pathlib import Path


def build(destination):
    """Create a new destination; keep original captures and sources unchanged."""
    import molsysmt as msm
    import numpy as np
    import pyunitwizard as puw

    root = Path(__file__).resolve().parents[1]
    original = root / 'tests/data/interface_topomt_pisa/3a8ee979'
    protocol = json.loads((original / 'protocol.json').read_text())
    destination = Path(destination)
    destination.mkdir()
    rows = {}
    for stem, filename in [('positive', 'input.pdb'), ('empty', 'empty_input.pdb')]:
        source = original / filename
        payload = source.read_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        if digest != protocol['files'][filename]:
            raise ValueError('The offered original source differs from its protocol.')
        system = msm.convert(str(source), to_form='molsysmt.MolSys')
        output = destination / f'{stem}.h5msm'
        msm.h5msm.write(system, str(output), format_version='1.0')
        restored = msm.h5msm.read(str(output))
        np.testing.assert_array_equal(
            puw.get_value(msm.get(system, coordinates=True), to_unit='nm'),
            puw.get_value(msm.get(restored, coordinates=True), to_unit='nm'),
        )
        assert source.read_bytes() == payload
        rows[stem] = {
            'original_file': filename,
            'original_sha256': digest,
            'file': output.name,
            'sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
            'bytes': output.stat().st_size,
            'n_atoms': msm.get(restored, n_atoms=True),
            'n_structures': msm.get(restored, n_structures=True),
        }
    report = {
        'scope': 'New native copies of unchanged offered sources; no new PISA execution.',
        'bindingmt_base': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=root, text=True
        ).strip(),
        'molsysmt_source': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=root.parent / 'molsysmt', text=True
        ).strip(),
        'producer_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'sources': rows,
    }
    (destination / 'manifest.json').write_text(json.dumps(report, indent=2) + '\n')
    return report
