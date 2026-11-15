# BindingSite: proposed first saved TopoMT receiving case

Owner: [BindingMT #111](https://github.com/uibcdf/bindingmt/issues/111).
Profile review: 2026-11-14; public-reader readiness: 2026-11-15.
**Public saved recovery/resave and changed-context refusal pass under both unit
policies. Owner admission remains [TopoMT #199](https://github.com/uibcdf/topomt/issues/199);
the BindingMT adapter and complete T0/T1 workflow remain unimplemented.**

This makes the input request in the [implementation plan](binding_site_implementation_plan.md#10-t-one-complete-saved-topomt-receiving-workflow)
concrete. It follows [TopoMT's reviewed offer](https://github.com/uibcdf/bindingmt/issues/111#issuecomment-6814170065).
It does not freeze the proposed shared profile or qualify native DFND, every
original provider, a scientific binding finding or an installed/public pairing.

## Candidate and exact input identities

Propose the existing Pocketeer structural control used by TopoMT's
[`tests/test_saved_provider_workflows.py`](https://github.com/uibcdf/topomt/blob/c972b13d91d19ac763f1f0ce0ae75ebf3712ae0b/tests/test_saved_provider_workflows.py).
Its `_archived_pocketeer` helper already selects the source and original bundle
from the public evidence archive. This is a candidate for reuse, not a request
to repeat preparation or pocket detection or to install the original engine.

Inspected TopoMT source: `c972b13d91d19ac763f1f0ce0ae75ebf3712ae0b`.
Archive: [`devguide/pocketeer/evidence/pocketeer_saved_output_public_receptor_20261106.zip`](https://github.com/uibcdf/topomt/blob/c972b13d91d19ac763f1f0ce0ae75ebf3712ae0b/devguide/pocketeer/evidence/pocketeer_saved_output_public_receptor_20261106.zip).
The outer archive is 182226 bytes, SHA-256
`a795eb8778f861ddff39e1cfc86284b87c6c819b1b53b788ffb33865a413f1cb`.

| Exact member | Bytes | SHA-256 |
| --- | --- | --- |
| `structural-control/source.h5msm` | 365600 | `6312e2b0b3c25f08df22302f380c4cd52dfdf98a8becb3d2796515b91e2d0ca8` |
| `structural-control/original.zip` | 41202 | `3b926b3f780ced981b1e5b17d060aa65b9b292b2ac9063837dbc9092123df042` |

The original manifest declares provider `pocketeer`, backend `library`, run ID
`22ee8bc47f7e6a64a6193955dd92c93b981fe3aecc34714a9d4f8cc589d24e82`.
Preserve both `distribution_version="0.4.0"` and its separate saved
`version="0.3.0"` field with their original meanings; do not rewrite one to match
the other. The outer receipt identifies its historical producer as
`35e43f23408c04a47b314fa66bc706b049b24db0`. Current reader and original producer
identities remain separate.

Saved source context declares 745 atoms, structure index 1 and reference chemical
state 0 (one saved state). The original request lists 744 source atoms in reverse
order, excluding source atom 248. The saved sorted source selection has 744
entries and the filtered provider input has 638 entries. These are different
axes, not interchangeable counts or an inferred correspondence. Preserve the
original selection, its supplied order, owner input maps and selected frame.
Source frame 1 may become problem-local frame 0 only through an explicit map.

The raw original output reports pocket IDs 0 and 1. Propose the first reported
pocket as a deterministic software control, subject to TopoMT confirming the
exact recovered public `record.source_id`, role and mapping. This choice does
not rank biological relevance. Obtain membership and geometry from public
recovery; the manifest counts and raw sphere indices cannot supply a BindingMT
atom map.

During the initial read-only inspection on 2026-11-14, all outer manifest member
lengths/hashes and inner artifact hashes were checked
by ZIP/JSON inspection. No TopoMT import, molecular conversion, source recovery
or BindingMT adapter was executed during that initial inspection. Hashes bind bytes;
they do not establish the meaning or applicability of a finding. The tracked
generated-version edit in the TopoMT development checkout is preserved. Any
future runtime receiving records that exact worktree and actual import origins.

## Executed public-reader readiness — 2026-11-15

The small client check now uses the actual public MolSysMT loader and
`topomt.load_provider_output`, followed by public `output.save` and recovery.
It does not introduce a BindingMT import adapter or claim owner profile admission.
The unchanged full source has 745 atoms and two structures; recovery restores
source frame 1/reference state 0 and the original run/selection context.

| Recovered public record | Role | Source members | Member points | Lining points | Alpha spheres |
| --- | --- | --- | --- | --- | --- |
| `pocketeer:0` | `pocket_residue_atoms` | 84 | 84 | 32 | 44 |
| `pocketeer:1` | `pocket_residue_atoms` | 109 | 109 | 32 | 46 |

The public member population includes atoms of associated residues. The lining
view separately contains the union of sphere-defining atoms. Neither is the
alpha-sphere support or an inferred closed volume. The proposed first member
enclosure must name its 84-member view and role; it cannot silently use the
32 lining points or promise to enclose every sphere. The public measurements
are `volume` and `score`, with their original definitions/units retained.

Both nm/ps and pm/fs policies pass positive load/save/load, rich record fields,
geometry, quantities, original identities and input-coordinate comparisons.
There are four successful public recoveries and two expected refusals, not a
new pytest suite. Run/artifact/metadata identity is exact. Physical floating
replay uses 8 binary64 epsilons relative, zero absolute tolerance, accounting
only for unit-conversion roundoff; it is not a method-accuracy tolerance.

Moving original selected detector-input atom 45 by 0.1 nm in an independent
source copy, outside the first record's exported membership, causes the public
reader to refuse with ValueError under both policies. Original source positions
and source/bundle file hashes remain unchanged. This exercises the owner's
full-input context check, not a claim that unchanged lining alone proves a
current finding. Source compatibility and permanent current applicability
remain separate obligations for any future BindingMT application.

Historical `output.attribution` is None in this archive and remains None after
resave/recovery. No bibliography is invented, and positive preserved-citation
coverage is not established by this case. TopoMT #199 asks whether this explicit
absence is admissible or an attributed saved case should be supplied.

### Exact source, guards and limits

Final executed TopoMT reader source is
`c83ab219c205a89a74ef8a8b6778fb10a89a371f`, publicly published. Another TopoMT
work session advanced its documentation between attempts; its runtime package
bytes and the frozen archive are unchanged from c972b13d. Every completed cell
has identical before/after source/head/worktree/native/input snapshots. The
tracked generated-version edit is preserved and hashed, with actual import
origins bound to the development checkout. BindingMT performs no provider edit,
merge, build, installation or native replacement.

MolSysMT remains `6deb0200cd7e384128eded46caaf751a57e65bc6`, PyUnitWizard
`cc51c8e600f78f108bc563edd3221bd8788ead86`, ArgDigest
`7be790b7334e16715fceebbfa8a2b400edf768bb`, SMonitor
`4a4cd1e415ea57e50c22925c36c2e727772fb3ab`, DepDigest
`5e93eaff314fe7cddf4eee89842977023b3a8b00` and the inspected Ackredit checkout
`0bb26ecb42edcfc0984469ff20c9ff93f8fd7f56`. Actual import paths and the mapped
MolSysMT native SHA-256
`2dbcf13d148cb3ae13a5908707fcdd61c75a2de5027e50c566ce07de81bdb038`
are recorded. This receives only the measured reader/source operations, not
native build authenticity, complete ABI or every current upstream revision.

The final check permits no originating-engine import, child process or network
call during recovery. Resolve Python's harmless platform probe before enabling
those guards. CPU6/nice10/one worker, CPU45/wall180, RSS512MiB/AS2GiB and monitored
reserve13GiB apply to the process and heavy descendants. Final sampled peak is
below 310 MiB; minimum available RAM stays above 45 GiB. No new resource
interruption occurs; earlier native498/evaluator holds remain unretried.
Slurm availability is not established by the sandbox socket probe, and no active
allocation/reservation is claimed. Temporary source/bundle/resave files retire.

Initial attempts preserve a premature process guard blocking Python platform
detection, a private comparison helper misreading a plain string as a quantity,
and overstrict equality for cross-policy binary64 geometry. They are harness
failures, not TopoMT defects or successful full receiving. A narrower earlier
load/resave check passes separately; only the complete final cell supplies the
rich replay, refusal and actual native-origin account above.

Six applicable workspace and fresh committed-export checks pass, and the export
is retired. Raw attempts/programs/source hashes/resource/origin/replay and
publication/readback receipts retain #50 custody in
`.molsyssuite/topomt-case-111-20261115/`, review 2026-11-22. No new scientific
execution, tag/package/release, hosted/full-suite/installed/public/platform or
shared/consumer adoption is inferred. Prior 1118-case receiving remains at its
original source and is not repeated or relabelled as this new reader check.

## Agreement requested from TopoMT

1. Confirm this exact source/bundle pair as an appropriate small source-bound
   case, or supply a replacement with exact published paths/hashes. Confirm its
   supported current reader source and the known saved structure/state context.
2. Name the recovered public record and its molecular membership role/index
   space; explain any unavailable or partial map/support. Confirm which owner
   artifacts establish source association and which remain caller preconditions.
3. Confirm the first application as **an enclosure of the explicitly reported
   fixed molecular members**, using `assurance_scope="placed_geometry"`. This
   does not claim to enclose every reported sphere or the entire pocket support.
   If that operation is inappropriate for the chosen record, agree another
   explicit view before implementing T0/T1.
4. Preserve original measurements/definitions, public geometries, fields and
   bibliography. Keep sphere geometry, membership and unavailable closed-volume
   semantics distinct. Import/replay adds no original-engine execution credit.

The existing public recovery boundary is
`topomt.load_provider_output(original_bundle, molecular_system=full_original_source)`.
It restores saved selection/frame/chemistry; the bridge must not override them.
Recovery needs the full molecular source, independently of the evidence ZIP.
Use public MolSysMT loading and mapping. Do not copy the owner test fixture,
parser, source matcher or numerical kernel into BindingMT.

## Bounded T0/T1 work after agreement

T0 implements an optional public adapter for the agreed ProviderRecord/run
profile, composing the existing immutable BindingSite data core. Keep TopoMT
loading behind the selected adapter route through DepDigest. Data-only saved
BindingSite reading needs neither TopoMT nor the originating engine. Retain
original identity/context, role/index spaces, measurements/units, geometry
availability, attribution and recovery refs; preserve the owner output unchanged.
The exact saved representation is an implementation admission decision, not a
new universal owner class or a frozen schema established by this document.

T1 exercises the admitted source-bound public route:

`saved source/bundle → TopoMT recovery → BindingMT import/inspection → save/read → explicit fixed-membership box → ordinary problem/result history`.

Use explicit length padding and the existing [native C2 application](binding_site_application.md).
Retain its source/normalized structure and atom correspondence, actual prepared
payload checks and caller consistency preconditions. No native docking is
required to qualify this software handoff. The original provider support and
the new member-enclosing box retain separate meanings and provenance.

Required positive/refusal controls:

- One complete public saved route under normal and pm/fs policies, followed by
  fresh-process replay without the original engine. Preserve owner/source bytes,
  record/run IDs, scientific quantities and bibliography.
- Incompatible source/structure/chemistry refuses owner recovery; check explicit
  source-to-normalized frame transport in the accepted problem. Never substitute
  a target-local frame for the original index.
- Unavailable mapping and evaluated-empty membership remain different states;
  center-only input never gains invented membership, radius, extent or volume.
- Changed detector-input atoms outside exported membership do not justify a
  current-finding claim from unchanged member coordinates. Historical data stay
  readable; the first profile explicitly refuses stronger assurance without
  rerunning the originating detector.
- Missing support, incompatible actual prepared input and inadequate padding
  refuse the affected box application. Public import or historical reading does
  not grant permanent applicability.

Run future receiving serially with the existing shared-workstation limits and
record actual sources/native artifacts before execution. It supplies its own
evidence and never resumes the earlier native498/evaluator holds. Shared
agreement, actual consumer adoption and installed/public/platform/native-engine
qualification retain their independent owners.

Next: obtain the four owner dispositions above in [TopoMT #199](https://github.com/uibcdf/topomt/issues/199),
linked to BindingMT #111, then admit and implement
the selected T0 slice. No provider defect or new missing kernel is alleged by
this input request. The original read-only inspection retains its prior #50 custody under
`.molsyssuite/selection-resource-756-20261114/`, review 2026-11-21. The subsequent
public-reader attempts have the separate custody/review above.
