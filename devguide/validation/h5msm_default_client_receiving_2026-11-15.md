# Original H5MSM packet: default write and migration receiving

Owner: [BindingMT #73](https://github.com/uibcdf/bindingmt/issues/73).
Provider: [MolSysMT #256](https://github.com/uibcdf/molsysmt/issues/256).
Status: **the finite default-IO source profile passes under both unit policies**.
Execution date: 2026-11-15, America/Mexico_City; stored UTC timestamps are 2026-11-16.

## Selected scope and prior acceptance

The maintainer forwards MolSysMT's outstanding request and selects this small
storage-only receiving within the current work. This resumes this finite handoff,
without restarting the separately held scientific/preparation programme or the
later final MolSysMT/MolSysViewer pair validation.

Earlier [experimental P5](h5msm_p5_client_receiving_2026-11-10.md) and
[explicit-version 1.0/P6 receiving](h5msm_v10_client_receiving_2026-11-10.md)
remain valid within their original populations. They are not substituted for
this different default-argument profile. No earlier preparation, minimization,
PDBQT export, docking/search or molecular campaign is repeated here.

## Exact original inputs and independent expected data

Reuse the three original complete/projected/merged **0.6** fixtures and original
typed expected JSON from the
[requested manifest at 08442cb14](https://github.com/uibcdf/molsysmt/blob/08442cb1496c8f6698ccc17ae1f7d21566d8713b/devtools/data/h5msm_p5_review_20261110/manifest.json).
All nine packet files, including the existing unused PDBQT references, match
that exact Git source byte-for-byte. No input or oracle is regenerated.

| Original input | SHA-256 |
| --- | --- |
| `complete.h5msm` | `188d077fcc036218e7ca1e873ff257d368dfa2b81b3965b70b714692de665ff5` |
| `projected.h5msm` | `7440d9d237d9fa57672ce6fd2afedba3f73096c31419dd99ee1d4d2132c86160` |
| `merged.h5msm` | `55d33285b688021eb72749b8fe60f76f6dbeb0f7c7ec6f922688e83e4563ca5f` |

The existing frozen
[`h5msm_client_oracle.py`](../../bmt_test_support/h5msm_client_oracle.py)
observes public native data and compares it with the original typed JSON. Its
two documented transport normalizations remain unchanged; it implements no
writer, restoration, chemical preparation or alternate molecular representation.
Whole snapshots retain topology, chemical states, structure/state associations,
coordinates/structure metadata, mechanics/parameter columns, charge and type
assignment reports, preparation history and named interactions/coverage. Original
index labels, dtypes, maps, quantities and scientific producer strings are checked
inside these full observations, rather than only comparing atom counts.

## Actual source, import and binary decisions

| Identity | Observed value |
| --- | --- |
| BindingMT producing base | `c09c6f90302e3b726371b1869afa33749235179c`; exact new worker/oracle/guard bytes are archived separately. No BindingMT runtime implementation changes in this receiving. |
| Actual editable MolSysMT | `6deb0200cd7e384128eded46caaf751a57e65bc6` |
| Requested owner reference | `08442cb1496c8f6698ccc17ae1f7d21566d8713b` |
| Actual mapped MolSysMT ABI3 | `2dbcf13d148cb3ae13a5908707fcdd61c75a2de5027e50c566ce07de81bdb038` |
| Python | 3.14.10, Linux, existing `molsyssuite@uibcdf_3.14` |
| PyUnitWizard | `cc51c8e600f78f108bc563edd3221bd8788ead86` |
| ArgDigest | `7be790b7334e16715fceebbfa8a2b400edf768bb` |
| SMonitor | `4a4cd1e415ea57e50c22925c36c2e727772fb3ab` |
| DepDigest | `5e93eaff314fe7cddf4eee89842977023b3a8b00` |

The requested owner reference was fetched, with **no merge or source/native/
installation change**. Comparing the full package and Rust/build inputs against
the actual editable source finds only two changed package files:
`form/mmcif_PdbxContainers_DataContainer/_bond_state.py` and
`form/mmcif_PdbxContainers_DataContainer/to_molsysmt_MolSys.py`.
The H5MSM/default-conversion/native domains used here and Rust inputs are unchanged.
This records a deliberate reuse decision for this finite H5 profile, not receipt
of the separate mmCIF changes or all of the owner reference.

The actual ABI3 path/digest is checked against its live mapped file in each
producer. After normal package import, every public callable on the loaded Rust
module is prohibited for the selected IO/observation operations; none of those
guards is reached. `kernel_calls=[]` records that empty guard-event list, not
complete bootstrap instrumentation or binary build authentication.
Original optional engines, subprocesses and network calls are prohibited during
the IO/observation cells. Native builder-source authenticity, general ABI
qualification and the owner's different historical reference binary are not
claimed.

A separately labelled post-run check verifies actual editable import origins
and the mapped native again; it records h5py **3.16.0**, HDF5 **2.2.0** and
SHA-256 values of its observed HDF5 extension files. That check does not replay
the writer or retroactively supply original writer build authentication.

The new roots literally report `creator_version=0.23.0+904.gff25f1202` from the
existing generated build stamp. This is **not** the live source Git identity.
Original roots retain `creator_version=0+unknown` within migration metadata;
original scientific producer data inside the recovered systems remains exact.
No version stamp, bibliography or earlier scientific producer is rewritten.

## Public commands and complete results

For each original, each destination is absent and new. The literal public calls
omit `format_version`; the archived worker's AST independently confirms both call
sites, rather than trusting a prose label.

```python
molsys = msm.h5msm.read(origin)
msm.h5msm.write(molsys, native_destination)
msm.convert(origin, to_form=migration_destination)
recovered = msm.h5msm.read(destination)
```

Policies use Pint quantities and the explicit Pint parser:

- Ordinary: `nm, ps, K, mole, dalton, e, kJ/mol, radians`.
- Alternate: `pm, fs, K, mole, dalton, coulomb, J/mol, radians`.

| Original | Native default write, ordinary | Whole-file default migration, ordinary | Native default write, pm/fs | Whole-file default migration, pm/fs |
| --- | --- | --- | --- | --- |
| complete | PASS | PASS | PASS | PASS |
| projected | PASS | PASS | PASS | PASS |
| merged | PASS | PASS | PASS | PASS |

All **12 requested operations pass**. A first two-operation complete/ordinary
admission also passes; it overlaps the final normal selection and is not pooled
as additional receiving or as an independent scientific case. The three storage
fixtures and repetitions do not constitute twelve independent molecular results.

Every actual destination declares H5MSM **1.0**, reopens publicly and matches the
complete original expected snapshot exactly. Native writes leave their borrowed
systems unchanged. Every whole-file migration additionally has
`json.loads(target_root['migration_source']) == original_root` for **all original
root attributes**. Original input hashes survive every operation; caller unit
policy and complete before/after runtime/input fingerprints agree. New target
creation timestamps and file hashes are recorded individually, not mistaken for
changes in original scientific provenance or required to match across executions.

No failure or resource interruption occurs in this finite population. The
admission lasts 4.29 s; complete ordinary and pm/fs cells last 5.22 s and 5.43 s.
Maximum sampled group RSS is **296,218,624 bytes**. Cells are serial, CPU 6/nice 10,
all native/build workers one, AS 2 GiB/group RSS 512 MiB/CPU 45 s/wall 180 s,
with monitored **13 GiB available reserve**. No Slurm allocation is held; the
previous sandbox-denied controller probe is not evidence that Slurm is absent.
Free resources are not reserved, other work is preserved, and old heavy holds
remain unretried.

## Published records and independent inspection

The [primary manifest](data/h5msm_default_client_2026-11-15/manifest.json)
authenticates **76 compressed members**, 1,345,998 stored / 8,856,909 decoded bytes;
SHA-256 `5addb125362b9b3462ee98ab5b034fec157f4cd2dc1aee8b0b3e0885e445a849`.
It preserves exact original inputs/oracles, new worker/guard/observation code,
all admission/whole attempts, complete snapshots, actual output files, roots,
source identities and resource receipts. Original packet authorship remains
MolSysMT's; these are new consumer observations, not regenerated original data.

The [reader/origin companion](data/h5msm_default_client_reader_2026-11-15/manifest.json)
preserves **16 members**, including the independent reader and its receipts.
The fresh reader authenticates every compressed/decoded byte, checks the actual
stored HDF5 roots with h5py, all migration-source metadata, original input hashes,
literal omitted-version call sites and all snapshots against the original typed
oracles. It verifies all twelve distinct policy/source/route observations, keeping
the admission separate. MolSysMT/BindingMT/quantity/optional-engine imports and
molecular replay are prohibited; no preparation or native calculation is repeated.

Raw BindingMT #50 custody is `.molsyssuite/h5msm-default-client-20261115/`,
review **2026-11-22**. Original sources and borrowed environments remain intact.
The source-only records add no installed-test inputs, dependency/pin/schema or
runtime changes. Local documentation/report/resource checks and the required
fresh committed-source export retain separate receipts in that custody.

## Closure boundary and next action

This completes BindingMT's execution of the specifically requested
**default write/default whole-file migration** profile for MolSysMT #256.
MolSysMT owns provider acceptance and any stable-format freeze. Share this report
and its immutable publishing commit in BindingMT #73; retain previous P5/P6
acceptances without rerunning them.

No new package, build, installation, tag/release, hosted matrix or final installed
MolSysMT/MolSysViewer pair is qualified or launched. Other Python/platforms,
arbitrary H5MSM contents, scientific accuracy and the wider #73 objective retain
their own gates. BindingSite #111 and BindingInterface #112 remain at their
[current implementation checkpoint](../README.md); this storage receipt supplies
none of their outstanding regional/consumer owner inputs.
