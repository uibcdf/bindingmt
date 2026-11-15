---
summary: Resolve native interface participants on the chosen structure and its associated chemical state.
issue: uibcdf/bindingmt#117
status: resolved
opened: 2026-11-15
closed: 2026-11-15
severity: high
verification: reproduced
area: [core, scientific]
guard: tests/test_interface_selection_context.py
normative:
blocked_by: []
supersedes: []
---

# BindingInterface participant selection used an unrelated default context

## What

At source `2988dd16bc63086bc3db009e92dd776566865f91`, contacts and burial share
`admit_axes`; ionic characterization uses `_scope`. Both call public
`msm.select` without forwarding the chosen frame/state. MolSysMT correctly uses
its documented all-structure/reference-state defaults. Later calculations and
recorded state can therefore describe a different context from the participants.

The small authored two-frame/two-state control has reference-state positive/
negative atoms `[1]`/`[0]`, but frame 1 is associated with state 1 and needs
`[2]`/`[3]`. Its relevant pair is 0.25 nm apart. With `structure_index=1` and
`formal_charge>0`/`formal_charge<0`, all three operations select `[1]`/`[0]`.
Contacts and scoped ionic observations report completed zero instead of the
one observation obtained from the independently specified `[2]`/`[3]` axes.
Burial also uses the wrong participants; its surface-result count remains one
and must not be confused with a positive interface finding.

## How

Both selection calls now forward `structure_indices=[frame]` and
`chemical_state="structure"`. MolSysMT performs selection, associated-state
resolution and all molecular/numerical operations. An unresolved association
refuses through its native diagnostic before required calculation; there is no
reference-state fallback. Known single-state numeric selections retain their
ordinary behavior. Public docstrings name this scope explicitly.

The correction changes no signature, schema, atom mapping or scientific
definition. Records retain their actual resolved axes and original provenance.
Historical records remain readable without selection or reassessment; the fix
does not silently repair earlier results. Correcting an affected historical
calculation requires a separate explicit calculation on its intended input.
Caller source/topology/chemistry consistency and absence of concurrent edits
retain their existing limits; this is no observer or source-authentication claim.

## Why and measured evidence

The complete 34-case new module has **28 failures / 6 passes** on original
runtime and **34 passes** after correction. Final source qualification covers
**138 cases per normal/pm-fs policy**: new scope 34, complete contacts 29,
complete burial 51 and seven native-only ionic selectors covering 24 cases.
Final populations have no failures/errors/skips and use disjoint nodes.

Actual chemical, spatial and combined queries resolve the chosen frame/state;
owner call arguments, positive/empty findings, numeric frame 0/frame 1 controls,
unknown-association refusal, original MemoryError/OS ENOMEM/KeyboardInterrupt/
ordinary-error identity, Pint-to-unyt historical read and actual ionic/contact/
surface annotation composition are exercised. They are small authored
consistency controls, not biological, binding-affinity or method-accuracy evidence.
No PharmacophoreMT construction/matching is imported or executed.

Six actual JSON files from the two final policies also pass 12 fresh-process
Pint/unyt reads under pm/fs, with molecular query/build/calculation routes blocked,
no original source supplied and unchanged bytes/producers. These reads are
separate from the 138-case pytest population. The initial reader file glob counted
pytest's `current` symlink aliases; that harness failure is preserved. The final
reader admits the six distinct physical output files without changing their bytes.

The authored derivative is created through public MolSysMT operations and typed
ChemicalStatesDict conversion from the unchanged four-ion source asset, SHA-256
`75ed110abb15515ed33806cdd03d9958eb231b3fbe2427dec0ffbc70aa287d1b`.
Its new monatomic state explicitly declares complete empty connectivity.
Original fixture files/producers and their resource declaration retain their
bytes; no regenerated source is substituted for historical evidence. An initial
probe omitting the new state's connectivity declaration is preserved as input
preparation failure. The corrected probe reproduces all three wrong selections.

Runtime/test/resource/provider/native snapshots agree before/after final cells.
The unchanged source closure is MolSysMT `6deb0200cd7e384128eded46caaf751a57e65bc6`,
PyUnitWizard `cc51c8e600f78f108bc563edd3221bd8788ead86`, ArgDigest
`7be790b7334e16715fceebbfa8a2b400edf768bb`, SMonitor
`4a4cd1e415ea57e50c22925c36c2e727772fb3ab` and DepDigest
`5e93eaff314fe7cddf4eee89842977023b3a8b00`. Actual development import origins and
mapped MolSysMT native SHA-256
`2dbcf13d148cb3ae13a5908707fcdd61c75a2de5027e50c566ce07de81bdb038`
are checked separately; whole native build/ABI authenticity is not established.
Existing generated editable-version metadata is recorded separately.
Fetched newer provider commits remain unmerged/unqualified for this correction.

## What was refuted

This is a BindingMT orchestration defect, not a MolSysMT default-selection bug.
Explicit owner queries and explicit numeric controls agree with independent
membership expectations. A consumer selection parser, chemical repair, alternate
contact/SASA/ionic kernel, registry change or historical codec rewrite is
unnecessary. Declaring complete monatomic connectivity fixes only the initial
probe setup; it does not fix the reproduced consumer selection error.

## Scope and exclusions

The guard applies to the admitted native Cartesian contact/surface/ionic profile.
It adds no PBC, state transfer, arbitrary chemistry, pharmacophore semantics,
T/H/P consumer adoption or unrelated frozen-results repair. Earlier 368/1191-case
Site receiving keeps its original sources. Full-suite/hosted/installed/public/
platform/native-engine/final-pair and shared gates retain their owners.

Serial CPU6/nice10/one-worker cells use monitored reserve13GiB/RSS512MiB/AS2GiB/
CPU45s/wall180s controls. No new resource interruption or old heavy-hold restart
occurs. Six applicable workspace and fresh committed-export checks pass; the
export retires. Exact probes/XML/source/resource/origin/publication receipts
retain #50 custody under `.molsyssuite/interface-selection-context-112-20261115/`,
review 2026-11-22. No provider edit/install/build, tag/package/release or manual
remote CI occurs; skipped CI is not a pass.

## Acceptance criteria

The named guard verifies actual frame/state-dependent membership and numerical
consequences, not only forwarding syntax. It also receives unknown-state failure
before required children, direct owner-error identity and immutable historical
readback/composition. #112 retains its broader implementation/consumer scope;
closing #117 qualifies this independently bounded correction only.
