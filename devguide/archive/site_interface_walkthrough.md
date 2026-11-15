---
summary: Connect the existing Site, resolved region, native interface observations and historical readers in a copyable source example.
issue: uibcdf/bindingmt#118
status: resolved
opened: 2026-11-15
closed: 2026-11-15
severity: medium
verification: reproduced
area: [core, documentation, scientific]
guard: tests/test_site_interface_walkthrough.py
normative:
blocked_by: []
supersedes: []
---

# BindingSite and BindingInterface reproducible walkthrough

## What

The maintainer accepts a complete example of the already delivered public
BindingSite/region and BindingInterface operations. The
[copyable walkthrough](../../examples/site_interface/README.md) and its
[executable program](../../examples/site_interface/walkthrough.py) now load a
native MolSys, choose an explicit structure/state, capture a candidate site,
resolve its region and calculate actual contact and burial observations through
MolSysMT. They project the observed contact footprint and save/read all records.

The frozen authored control has four atoms, three structures and two states.
Structure 0/state 0 selects region `[1]` and partner `[0]`; structure 1/state 1
selects `[2]`/`[3]`. Both pairs are 0.25 nm apart. Structure 2/state 1 retains
`[2]`/`[3]` but separates them by 1.5 nm. With the explicit 0.3 nm cutoff, it
has completed contacts `[]`, footprint `[]` and zero side burial, while its
candidate region remains `[2]`. Positive side losses for structure 1 are about
0.401981 and 0.509150 nm² under the stated model. These are software and input
consistency observations, not biological cases or affinity evidence.

## How

Public `BindingSite.from_selection` captures `formal_charge>0` on the chosen
structure and its associated chemical state. `resolve_region` returns the public
region mapping; its exported source atom indices become one interface
participant. The counterpart is MolSysMT's state-aware `formal_charge<0`
selection on the same structure. The example creates a new captured site for
each run; an old snapshot is never silently retargeted.

Contacts and burial remain separate native `BindingInterface` records. Public
`ContactFootprint` projection records what was observed, preserving the parent
contact reference. All molecular selection, H5MSM, contacts and SASA remain
MolSysMT-owned. There is no new runtime facade, scientific schema, alternate
molecular implementation or optional engine.

The program exclusively creates a new output directory. It saves five public
result records and an exact copy of the source, verifies public round-trip
preservation, then writes `run.json` last. That receipt links exact file bytes,
record IDs and logical digests; it is example orchestration, not a new library
interchange authority, signature or multi-file transaction. Calculation/save
errors retain identity; partial files remain available without a completed
receipt. A subsequent run requires a new directory.

The `--read` route opens only the five result JSONs and their receipt. It
checks exact files and references, then uses existing public readers. It needs
neither molecular source nor `case.json`, and does not recalculate or infer live
assurance. It preserves original producers and quantity records. Display uses
explicit conversion to nm²; computation and storage retain negotiated quantities.

## Exact input and producing sources

[case.json](../../examples/site_interface/case.json) separates the independent
expected observations from the recipe used by the runner. It retains the parent
input and original producing identities. No frozen historical input was changed.

| Item | Exact identity |
| --- | --- |
| BindingMT runtime base/input producer | `47090f8af425896e9e243047bcac537ebd2c68b3` |
| New frozen H5MSM 1.0 input | 133184 bytes; SHA-256 `eb3174c05f4a9ab0a562eeb1fb151e3ca17c4e149af49ceb74de965613babe4d` |
| Original parent source | SHA-256 `75ed110abb15515ed33806cdd03d9958eb231b3fbe2427dec0ffbc70aa287d1b` |
| Unshipped authoring program | SHA-256 `0f6698d94caf087986d241a2999bec1cb0e99c7737d1f9b22862dc8a0cbfc39d` |
| MolSysMT source | `6deb0200cd7e384128eded46caaf751a57e65bc6` |
| PyUnitWizard source | `cc51c8e600f78f108bc563edd3221bd8788ead86` |
| ArgDigest source | `7be790b7334e16715fceebbfa8a2b400edf768bb` |
| SMonitor source | `4a4cd1e415ea57e50c22925c36c2e727772fb3ab` |
| DepDigest source | `5e93eaff314fe7cddf4eee89842977023b3a8b00` |
| Actually loaded MolSysMT native extension | SHA-256 `2dbcf13d148cb3ae13a5908707fcdd61c75a2de5027e50c566ce07de81bdb038` |
| Shared publication/support SDK | `95f693be10c9205ab9b912cb53453710a207bf45` |
| Synchronized common guide | SHA-256 `0188fee9330856bb1d40bdd5739c275a73471b602df20ba3a3f0ac8476677678` |

Execution uses Python 3.14.10 in `molsyssuite@uibcdf_3.14`. The actual development
origins and mapped native path/digest are checked, without claiming whole native
build/ABI authenticity. Generated BindingMT version bytes remain independently
identified: reported `0.1.0+19.g0ab5102`, SHA-256
`9b7ab0e153e277e83521b0b6314f77c6cd29f81b19e315e4921de3d1a070d493`.
The input's MolSysMT producer reports `0.23.0+904.gff25f1202`; its actual source
commit above is authoritative complementary evidence. Neither report is silently
rewritten to a later receiving commit.

The source test receipt fingerprints runtime/test/example/helper/declaration
bytes before and after execution. The final delivery commit and public readback
are linked in [#118](https://github.com/uibcdf/bindingmt/issues/118). Providers
remain at the qualified commits above: fetched newer remote commits are
unmerged/unqualified and are not imported as part of this delivery.

## Acceptance and evidence

The durable [17-case guard](../../tests/test_site_interface_walkthrough.py)
executes the canonical example, not a duplicate implementation:

- Nine complete runs cover three source structures under Pint nm/ps, Pint pm/fs
  and unyt pm/fs. Independent manifest expectations check region/partner axes,
  associated state, actual pairs, observed footprint and positive/zero side
  areas. Exact round trips, unchanged source bytes and caller policy are checked.
- Two copied-CLI cases cover positive and completed-zero output. Each runs from
  another directory with Python isolated mode, refuses an existing destination,
  then reads in a fresh process after both original input files and the saved
  source copy are removed. Molecular selection/get/set/extract/convert/copy,
  H5MSM reading, contacts and SASA are blocked. Vina and MolSysViewer stay absent.
  Saved result bytes, producers and observations remain unchanged.
- Five controls cover captured-structure refusal, missing completion receipt,
  changed result bytes, original required-calculation failure and partial-save
  failure without a completed receipt. The final control checks exact frozen
  input length/hash and the actual authoring program identity.

All **17 source cases pass**. The same **17 cases pass using the official SDK's
staged unshipped inputs**, with the test module copied outside the source tree,
its source `pythonpath` setting cleared, the canonical example staged below
`bmt_test_support`, and namespace binding/byte validation/retirement performed
by the SDK. This remains support qualification against the existing editable
runtime, not installed package, public artifact, platform or full-suite evidence.

The exact [resources declaration](../../devtools/conda-build/resources.toml)
adds only the authoring helper and the canonical example's code/source/manifest.
Code targets stay in the unshipped namespace; prior fixtures and producer
fingerprints remain intact. Runtime package-discovery rules are unchanged.
Workspace and fresh committed-export gates cover Ruff, formatting, report
indexes, document syntax/targets, shared exact resources and the Conda contract.
The reporting-protocol guard and actual-origin check also pass. Fresh export
and staged-support temporary trees/imports are retired after checks.

## Attempts, resources and custody

All cells use CPU 6, nice 10, one native/build worker, a monitored 13 GiB
available-memory reserve, group RSS 512 MiB, AS 2 GiB, CPU 45 s and wall 180 s.
Preflight records about 43 GiB available, no memory pressure or active swap
traffic. There is no Slurm allocation. The sandbox denied the controller socket;
controller availability is unknown, not inferred down from that failed probe.
No owned cell has a resource interruption; earlier heavy holds remain unretried.

The first 15-case source attempt preserves nine test-authoring failures
(`contexts` is a serialized field, not a public object attribute), then its
corrected 15-case selection passes. The final expanded 17-case selection also
covers cold completed-zero reading and partial-save failure. One initial support
attempt rejects a caller's unsupported context key before staging; the corrected
SDK invocation passes. These are preserved failed attempts, not provider bugs or
resource holds. No scientific fallback or provider patch was needed.

Raw attempts, source/resource/native inventories and publication/export receipts
remain under `.molsyssuite/site-interface-walkthrough-20261115/`, BindingMT #50
custody, review **2026-11-22**. Peak RSS and memory minima are in each actual
`resource-outcome.json`; no checkpoint broadens the limits to fit a calculation.
The migration directory and human/shared environments are preserved.

## Scope, review and next step

The walkthrough exposes the existing public API's current composition: a region
is a validated mapping, contact and burial are independent interface records,
and the example links saved records through a local receipt. This delivery does
not propose or claim a merged runtime result API. Domain scientific defaults,
new site discovery, dynamic/ensemble aggregation, transformations, pharmacophore
construction, ionic observations and search-domain application are outside this
finite example. Existing native membership/geometry application keeps its own
qualification.

T/H/P receiving in [the owner checkpoint](../binding_site_owner_receiving.md)
remains independent: TopoMT #199, PharmacophoreMT #226 and PeptideMT #60 retain
unreceived owner inputs. PyUnitWizard #156 and PharmacophoreMT #202 retain their
separate limits. #111/#112 and #30 remain open. No new provider defect, common
contract or dependency route was introduced. No package install/build, new tag,
release or manual remote CI was performed; skipped CI is not a pass. Historical
1191/368/138-case populations retain their exact original sources and scope.

Next review the copied walkthrough as a user, use it to discuss concrete API
friction, and receive whichever finite owner case becomes ready. Admit that
case's actual profile and qualify its public route rather than claiming generic
cross-component interoperability from this authored control.
