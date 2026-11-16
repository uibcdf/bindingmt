---
summary: Saved native contact/footprint readers inspect malformed second-operation options before checking their type.
issue: uibcdf/bindingmt#120
status: resolved
opened: 2026-11-15
closed: 2026-11-15
severity: medium
verification: reproduced
area: [core, persistence]
guard: tests/test_contact_record_options.py
normative:
blocked_by: []
supersedes: []
---

# Reject malformed saved contact operation options

## What

At source `abc0af924161ad698d5416f7f029d508e3435436`, the public native
BindingInterface and ContactFootprint readers raise an internal `AttributeError`
when the second saved required operation has scalar/null `options`. Both contact
calls require structured option objects, but endpoint-order admission checks
only the first one's type. The shared validator later calls `.get` on either.
No malformed record is accepted; this is a diagnostic/admission defect.

## How

Check the mapping type for each call before inspecting its fields, using the
existing `ArgumentError`/`calls` route. Retain exact valid native records, their
identity and quantity semantics. Do not add a schema or broaden scientific input.

Use actual public positive/completed-zero contact records and their independent
footprints; exercise both operation positions, missing/scalar/list options,
public from-dict and file reading. Preserve input bytes and caller unit policy,
block molecular work during saved admission and propagate original resource/
process failures. The new guard is written before the source correction.

## Why

The #118 walkthrough makes independent history reading a concrete supported
workflow. Invalid saved input should produce the existing structured field
error rather than expose a Python implementation detail. Scope follows #112's
D2 malformed-record and no-recalculation controls.

## What was refuted

The valid source records read correctly. These are data-only corruption probes,
not new calculations, completed-zero outcomes, provider execution failures or
resource interruptions. The original cutoff and source bytes remain unchanged.
A separate malformed QuantityRecord error belongs to PyUnitWizard; it must not
be repaired with a consumer codec or a blanket exception wrapper here.

## Acceptance criteria and limits

The addressed public routes reject malformed `calls` as `ArgumentError`, leave
input unchanged and retain existing valid historical evidence without molecular
queries. Scientific values, wire formats, producer identities, native bytes,
shared prefix and the 0.1.0 tag remain unchanged. #111/#112/#30 broader and
owner-input gates retain their independent scope. Qualification, custody and
publication evidence are appended before closure.

## Delivered correction and exact evidence

The shared contact validator now checks `isinstance(options, Mapping)` before
its option predicates for every required call. It uses the existing `_require`
and SMonitor `ArgumentError` without wrapping provider/process exceptions. Wire
versions, scientific values, operation definitions, payload discovery and public
signatures remain unchanged. Empty first-operation mapping still receives the
existing `participants` diagnostic because its endpoint axes are absent; the
fix does not reclassify that established mapping/axis error.

On unchanged runtime `abc0af924161ad698d5416f7f029d508e3435436`, the finalized
172-case regression selection has **48 failures / 124 passes**. Every failure
is the reproduced internal `AttributeError` for a nonmapping second-call option
value. After the correction, **268 complete affected cases pass per normal and
pm/fs policy**: the new 172 cases plus existing contact, interface-data and the
actual #118 walkthrough controls. No old interrupted selection is restarted.

The new controls exercise actual positive and completed-zero public native
records and independent footprints. They cover both operation positions,
missing/null/scalar/list/empty-map options, from-dict and file reading with valid
outer seals, unchanged input/file bytes, and molecular queries blocked during
saved admission. Positive/empty valid readers retain exact producers/values
under Pint nm/ps, Pint pm/fs and unyt pm/fs. Injected original `MemoryError`,
`OSError(ENOMEM)`, `KeyboardInterrupt` and `SystemExit` from the quantity provider
retain their exact identity, rather than becoming invalid-data or zero results.

The same **172 cases pass on official SDK-staged unshipped inputs**, against the
existing editable runtime, with the copied test outside the source tree and its
source `pythonpath` cleared. The SDK checks/binds/retires the namespace and
exact files; support trees/imports are retired. Existing declared #118 source
inputs are reused unchanged. No helper/data/dependency/declaration edits or
private namespace stager are needed. This is not installed-package qualification.

Actual composition stays MolSysMT `6deb0200cd7e384128eded46caaf751a57e65bc6`,
PyUnitWizard `cc51c8e600f78f108bc563edd3221bd8788ead86`, ArgDigest
`7be790b7334e16715fceebbfa8a2b400edf768bb`, SMonitor
`4a4cd1e415ea57e50c22925c36c2e727772fb3ab` and DepDigest
`5e93eaff314fe7cddf4eee89842977023b3a8b00`. The actual-origin/mapped-native
check retains SHA-256
`2dbcf13d148cb3ae13a5908707fcdd61c75a2de5027e50c566ce07de81bdb038`;
whole native build/ABI authenticity remains unassessed. Generated BindingMT
version bytes retain SHA-256
`9b7ab0e153e277e83521b0b6314f77c6cd29f81b19e315e4921de3d1a070d493` and
report `0.1.0+19.g0ab5102`; they do not replace the actual source identity.
Runtime/test/declaration/example bytes are checked before/after the selections.
Newly fetched provider commits remain unmerged/unqualified. The final delivery
commit/public readback are linked in [#120](https://github.com/uibcdf/bindingmt/issues/120).

Six workspace and fresh committed-export gates cover Ruff/format, indexes,
document syntax/targets, exact shared resources and Conda contract. The reporting
guard and actual-origin check also pass; the fresh export is retired. No new
tag/package/install/build or manual remote CI; skipped CI is not a pass. #30
retains complete hosted/full-suite/installed/public/platform gates. The finite
T/H/P owner requests still have no new supplied case in this inspection.

## Independent provider feedback and preserved attempts

[PyUnitWizard #158](https://github.com/uibcdf/pyunitwizard/issues/158) receives
the separately isolated `QuantityRecord.from_dict` malformed block metadata
case. Its documented `RecordError` is instead a `TypeError` for null `blocks`
or `block_digests`. The tiny direct reproduction imports no BindingMT/molecular
source. Refreshed `c1cf832063f61efda0e9d68031b1ed2d149a6cb4` has the same
`record.py` bytes as the executed provider. This local fix neither catches nor
repairs that provider boundary; provider receiving remains pending.

Keep the initial bounded structural probes and all setup/test attempts. The
first regression launcher selected a relative path from its task directory and
ran no tests. The next attempt had the 48 actual failures plus eight test-author
field expectations, and compared lazy backend loading as if it were policy
mutation. The finalized baseline corrects those harness issues, primes the
chosen backend, preserves the original diagnostics and confirms 48/124. No
attempt is converted into a pass or described as a resource interruption.
The first status CLI used isolated mode without its sibling tool modules; the
second selected BindingMT against a historical registry. The documented CLI
then refreshes the five registered providers; BindingMT is inspected directly.
These bounded setup failures precede scientific work and change no sibling files.

All cells use CPU 6/nice 10/native one, monitored 13 GiB available reserve,
RSS 512 MiB/AS 2 GiB/CPU 45 s/wall 180 s. RAM/pressure/swap/load preflight is
recorded; there is no Slurm allocation. Controller availability remains unknown
after the earlier sandbox socket denial. No new resource interruption occurs.
Raw attempts, exact inventories, export/staging/origin/publication receipts stay
under `.molsyssuite/record-admission-112-20261115/`, #50 custody, review
**2026-11-22**. Shared/human work and migration data are preserved.

Next receive the first ready T/H/P case and the separate owner fix for #158.
The new central latest-policy handoff is independently owned by BindingMT #119
and MolSysSuite #201; neither its receipt nor a version edit supplies scientific
or installed qualification.
