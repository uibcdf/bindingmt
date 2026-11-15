# From a BindingSite to measured interface observations

This runnable source example connects the current public `BindingSite`, region,
`BindingInterface` contact/burial and `ContactFootprint` APIs. It needs an existing
compatible BindingMT/MolSysMT development environment; see the
[development setup](../../devguide/DEVELOPMENT.md). It downloads nothing and uses
no docking, visualization or external scientific engine. The records remain
explicitly experimental. The tag `0.1.0` predates this example; use the source
revision containing this directory, not that tag's API.

Copy this whole directory anywhere. From that directory, choose a new output
path whose parent already exists:

```bash
python walkthrough.py /tmp/bindingmt-example-1 --structure-index 1
python walkthrough.py /tmp/bindingmt-example-1 --read
```

The first command performs small native MolSysMT calculations and saves the
results. The second reads the saved history without reading molecular input or
recalculating. Each calculation needs a fresh output directory; an existing one
is refused. To compare all three structures, run the first command with indices
`0`, `1` and `2`, using a different output directory each time. On a shared
machine, run these sequentially with your environment's CPU/memory controls.
The calculation uses one native worker, Cartesian coordinates, no PBC, a 0.3 nm
contact cutoff, a 0.14 nm probe and 240 sphere points. Its explicit example work
bounds are 16 dense selection pairs and 1440 surface evaluations; those bounds
are not a RAM reservation.

## What the example shows

`source.h5msm` is an authored four-ion, three-structure, two-chemical-state
control. Its coordinates and charges deliberately distinguish source identity,
chemical state, a selected candidate region and actual observations. These are
software consistency controls, not a biological binding-site criterion or an
affinity prediction.

| Structure | Associated state | Candidate region atoms | Partner atoms | Contact pairs | Region footprint | Buried areas |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 0 | `[1]` | `[0]` | `[[1, 0]]` | `[1]` | Positive on both sides |
| 1 | 1 | `[2]` | `[3]` | `[[2, 3]]` | `[2]` | Positive on both sides |
| 2 | 1 | `[2]` | `[3]` | `[]` | `[]` | Zero on both sides |

All indices above are **source atom indices**, not row positions in an extracted
participant. The candidate rule is the authored `formal_charge>0`; the partner
rule is `formal_charge<0`, both evaluated on the chosen structure's associated
state. Between structures 0 and 1, that state changes and so do the selected
atoms. Between structures 1 and 2, the same atoms remain selected, but their
separation changes from 0.25 nm to 1.5 nm. A completed empty contact result is
therefore valid and differs from a failed calculation.

For structure 1, the displayed region is `[2]`, the observed footprint is `[2]`,
and the buried areas are approximately 0.401981 and 0.509150 nm². These are
side-specific sampled SASA losses under the stated model. The program displays
actual returned quantities after explicit conversion to nm²; saved records keep
full precision, negotiated units, numerical diagnostics and producer identities.
Neither contact counts nor buried area are binding free energies or validation
of an association claim.

## Follow the public calls

Read [`walkthrough.py`](walkthrough.py), the executable source of this guide:

1. `msm.h5msm.read` loads the native `molsysmt.MolSys`. The frozen manifest's
   SHA-256 identifies its bytes; `source.revision` carries that identity.
2. `BindingSite.from_selection(..., structure_indices=[frame],
   chemical_state='structure')` captures the selected members and their scope.
   `site.resolve_region(system, structure_index=frame)` resolves that fixed
   snapshot. It returns the public region mapping; the example passes its
   `projection.exported_target_atom_indices` to the interface operations.
   Molecular selection remains MolSysMT-owned.
3. `calculate_interface_contacts` and `calculate_interface_burial` receive the
   same explicitly chosen source/frame and disjoint participants. Contacts and
   SASA are computed by MolSysMT. They return **separate BindingInterface
   records**, each preserving its method, coverage, axes, units and provenance.
4. `contacts.project_footprint(participant='region')` derives the actually
   observed region atoms from the contact result. It does not redefine the
   candidate site. For structure 2 this footprint is empty, while the site
   still selects atom 2.
5. Public `save`/`read` methods round-trip the site, both interfaces and the
   footprint. The public `BindingSite.region_from_dict` validates the separately
   saved region mapping. The example checks exact record preservation before
   completing its output.

A captured site is not a live selector: reusing the structure-1 snapshot as if
it belonged to structure 2 is refused. Running this example on another structure
creates another snapshot from the explicit rule. This illustrates one supported
fixed-structure workflow, not a dynamic pocket finder, ensemble aggregator,
geometric search-domain export or cross-component adoption profile.

## Saved files and reading without molecular input

| File | Meaning |
| --- | --- |
| `site.json` | Captured BindingSite, with its source, selection and structure/state scope |
| `region.json` | Public resolved-region record, including the site reference and maps |
| `contacts.json` | Native contact BindingInterface and completed contact coverage |
| `burial.json` | Native SASA/burial BindingInterface, side areas and numerical qualification |
| `footprint.json` | Derived ContactFootprint with its parent contact reference |
| `source.h5msm` | Exact copy of the frozen molecular input for deliberate later live work |
| `run.json` | Local example receipt linking the five result files and record identities |

`run.json` is written last, after scientific execution and public round-trip
checks. A failed calculation or save propagates its error and leaves no
completed receipt. Retain partial files for diagnosis and retry in a new
directory. This is not a multi-file transaction, a new library interchange
format, a signature or independent scientific verification. UUIDs differ between
fresh runs; compare observations, not whole files from independent runs.

The `--read` route checks the receipt's five file hashes and record references,
then uses existing public readers. It does **not** open either `source.h5msm` or
`case.json`, and it never reselects or calculates. To demonstrate that boundary,
copy `walkthrough.py` alone to another directory, keep the five JSON results and
`run.json`, and invoke it with `--read`. The saved molecular copy is optional
for that operation. Reading preserves historical producer identities; changing
application unit policy or installed code does not relabel them. Hash checks
catch mismatched bytes; they do not authenticate who created a result.

## Input provenance and qualification

[`case.json`](case.json) records the independently authored expected observations,
recipe, exact source hash, parent input and producing source/helper identities.
The runner uses its recipe and input digest, **not** its expected answers, to
compute results. The H5MSM 1.0 file is frozen; do not regenerate it to follow a
newer provider. Its authoring program is
[`build_site_interface_input.py`](../../bmt_test_support/build_site_interface_input.py),
which uses public MolSysMT operations and requires a new destination.

The [delivery record](../../devguide/archive/site_interface_walkthrough.md)
names exact qualified sources, unit policies, cold-reader controls and limits.
The [tests](../../tests/test_site_interface_walkthrough.py) execute this actual
copyable program, including a fresh reader process with molecular work blocked.
The files also have exact declarations for MolSysSuite's unshipped installed-test
support staging; that validates support inputs and does not establish a new
public package or installed-artifact qualification.
