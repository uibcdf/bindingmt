# Inspect a preserved interface analysis

This development-source example uses the small, already saved PISA case supplied
by TopoMT. It recovers the original analysis, preserves its attribution and
projects the reported contact atoms of one participant. TopoMT and the original
native molecular source are required for recovery/projection. PISA execution,
network access and new docking are unnecessary.

Run from the checkout root in the development environment with its qualified
MolSysMT and optional TopoMT sources. Use new output paths: saving never overwrites.
These experimental capabilities postdate the `0.1.0` tag.

```python
from pathlib import Path

import bindingmt as bmt
import molsysmt as msm
import topomt as tmt

data = Path('tests/data/interface_topomt_pisa')
system = msm.h5msm.read(str(data / 'native_6deb0200/positive.h5msm'))
bundle = data / '3a8ee979/pisa.zip'
original = tmt.load_provider_output(bundle, molecular_system=system)

interface = bmt.BindingInterface.from_topomt(
    original,
    question={'kind': 'association'},
    participants={
        'A': {'source_id': 'original-pisa-control', 'source_atom_indices': list(range(20))},
        'B': {'source_id': 'original-pisa-control', 'source_atom_indices': list(range(20, 40))},
    },
    grouping=[['A'], ['B']],
)
print(interface.summary()['analyses'])
footprint = interface.project_reported_footprint(
    system,
    bundle_file=bundle,
    finding_id='1',
    participant='A',
    contact_groups=['other-bonds'],
)
print(footprint.source_atom_indices)
interface.save('preserved-interface.json')
footprint.save('reported-contact-atoms.json')
```

The expected A membership is
`[0, 1, 2, 3, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 18]`.
The footprint retains every contributing original contact/endpoint, its parent,
run, finding, criterion, quantities and scoped historical correspondence checks.
Select `participant='B'` for B's reported atoms. Selecting `['h-bonds']` on this
finding returns a completed empty footprint; a completed analysis with no
interface findings has no finding to project.

Saved records can be inspected later with only BindingMT's normal core
dependencies. Reading/copying/summary do not import TopoMT, open the source or
bundle, query MolSysMT, contact a server or endorse a different current structure.

```python
import bindingmt as bmt

saved = bmt.ContactFootprint.read('reported-contact-atoms.json', experimental=True)
print(saved.source_atom_indices)
print(saved.summary()['current_assessment'])  # 'not_performed'
```

For new use against the unchanged complete original source, explicitly call
`saved.check_context(system, bundle_file=bundle)`. It returns a new check receipt
and preserves the original record. A changed source, missing asset or unavailable
owner keeps its actual diagnosis. A correspondence match does not establish a
current binding site, biological binding, an energy prediction or a current PISA
finding. The first molecular profile uses the complete original axis and reference
frame/state; transferring to another source is outside it.

The [receiving contract](../../devguide/binding_interface_topomt_receiving.md)
records the exact sources, frozen producer identities, tests and limits. Native
MolSysMT contact/burial calculations remain available through the separate
[Site/interface walkthrough](../site_interface/README.md).
