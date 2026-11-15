# BindingSite declared-geometry placement: first native profile

Owner: [BindingMT #111](https://github.com/uibcdf/bindingmt/issues/111).
Admission: 2026-11-15, before implementation. **Implemented/source-received for
the bounded original-native sphere/axis-aligned-box profile; 1191 selected cases
pass per unit policy, plus independent actual saved-history receiving.**

## Chosen workflow and public request

Apply an existing sphere or axis-aligned box declaration to its original resident
MolSysMT MolSys, one source structure and its known associated chemical state.
Extend `to_search_domain` with optional `frame_id`: geometry placement requires
the exact nonempty ID of its declared Cartesian/PBC-false frame. Supplying that
ID explicitly associates the saved geometry with this source's stored frame;
the caller owns this initial association. Matching strings, counts or coordinates
do not authenticate an asset or prove an independent alignment. No transform,
independent target, oriented-box repair or stronger finding assurance is admitted.

Use existing `sphere_enclosure`/`axis_aligned_box` derivation and BoxRegion.
Preserve optional padding, source structure/state scope, units, exact original
declaration and the complete unchanged offline recipe. The new application is a
separate experimental `bindingmt.site_geometry_application` 0.1 record, with its
own UUID, content digest, producer, explicit request and assurance/preconditions.
Earlier unplaced recipes never become placed by being read.

The view is a spatial declaration, not molecular membership. Do not select atoms
inside it, require receptor atoms inside the box, or infer pocket volume/support.
An empty receptor-space region can be a valid declared search region. Molecular
selection, coordinates, states and maps remain public MolSysMT operations.

## Checked context and lifetime

Require bounded nonempty source atom inventory, structures and known associated
state. Select the sole declared structure, explicit requested structure, or the
sole source structure; no first-frame default for multiframe independent geometry.
An explicit state must agree with the structure and any declared state scope.
Record checked counts and the exact ordered binary64 nm position digest of the
whole selected source frame, without copying molecular coordinates into Site or
recipe. Signed-zero normalization and exact numerical meaning reuse the existing
native application witness; a digest is not source authentication.

Capture and actual use recheck the original source. DockingProblem must choose
that same source structure and native sorted-unique subset; source structure 2
can become normalized structure 0. Verify its actual native coordinates against
the corresponding original positions and its known associated state. Record a
detached historical construction receipt with the original/local axes and checked
position digest. Existing problem/result persistence preserves and rechecks it.
Backend export needs the actual problem and repeats the checks. A replacement
application requires a newly constructed problem.

This first profile conservatively requires the whole captured source frame to
remain unchanged, even for edits outside the declared shape. This does not prove
that every such edit invalidates a fixed Cartesian shape: it refuses reuse of this
particular captured placement. The saved declaration and historical application
remain readable; explicit recapture can reaffirm placement after an edit. Edits
to other unselected frames do not change this selected-frame witness. Caller
association and consistent inputs throughout check/use remain preconditions; no
mutation observer, atomic transaction or current-finding guarantee is promised.

## Actual prepared payload and effective engine box

Reuse the existing owner preparation/writer-map composition and staged-byte
admission before Vina's constructor. For this first profile every normalized
original receptor atom must remain directly mapped and keep its full-precision
position. Generated atoms can be explicitly present in the owner map; omitted
original atoms refuse this conservative frame-check profile. This restriction
is placement evidence, not mandatory contacts or chemical equivalence.

Read the actual written mapped positions through public MolSysMT PDBQTFileHandler.
They must agree with the checked prepared positions under the owner's declared
three-decimal angstrom serialization, allowing only binary64 representation
roundoff. No receptor containment condition applies. Check the actual six-decimal
angstrom engine box against the original unpadded declared sphere/box enclosure;
insufficient explicit padding or rounded-zero extents refuse before engine use.
Do not silently enlarge a box or claim that requested padding is preserved exactly
by backend quantization. Retain the effective box, writer payload hash/map, checked
atom axes and numerical convention in a separate historical payload receipt.

## Receiving controls and remaining scope

Exercise public sphere/box placement, independent/scoped declarations, explicit
frame association, source 2/local 0, nonidentity receptor subsets, unit forms and
policies, actual problem/result save/read, strict historical recipe integrity,
recapture after edits and changed source/normalized/state/scope refusals. Compare
actual prepared/exported bytes and maps, generated/omitted atoms, shifted payload
coordinates, inadequate padding and a receptor completely outside the box.
Use software engine sentinels; no docking or preparation-method accuracy follows.

Fresh readers need no molecular source or originating engine; historical replay
preserves producers and grants no numerical export without actual input checks.
Keep original molecular/quantity recipes and frozen fixtures unchanged. Run
affected checks serially with one CPU/worker, low priority, task-local monitored
reserve and bounded CPU/wall/RSS/address space. Earlier heavy interruptions are
not resumed. Workspace and fresh committed-export resource checks precede delivery.

TopoMT #199 and PeptideMT #60 keep their independent owner-input scopes. This
local placement profile adds no shared schema freeze, universal frame receipt,
source matcher, molecular geometry kernel or automatic detector rerun.

## Public example and historical record

```python
import bindingmt as bmt

frame = {"id": "original-cartesian", "coordinate_policy": "cartesian", "pbc": False}
site = bmt.BindingSite.from_sphere(
    source={"owner": "application", "id": "original", "revision": "r1"},
    center="[30, 40, 50] nm", radius="0.5 nm", frame=frame,
    structure_indices=[2], chemical_state_indices=[1],
)
box = site.to_search_domain(
    original_molsys, method="sphere_enclosure",
    frame_id="original-cartesian", padding="0.2 nm",
)
problem = bmt.DockingProblem(
    original_molsys, prepared_partner, box, receptor_structure_index=2,
)
numbers = box.to_backend_box(problem=problem)
```

These are authored software-example coordinates/state/padding, not scientific
recommendations. `axis_aligned_box` uses the existing box declaration instead.
An independent geometry needs an explicit structure when the source has several;
an unknown declared state remains unknown in the original Site while the new
application records its actually observed known association. `inspect` reports
conditional placement capability for spheres/identity-axis boxes; declaration
inspection itself still reports placement not assessed. Oriented boxes retain
readable geometry and an unavailable placement capability.

The experimental geometry application retains `geometry` (the whole old offline
recipe), `request`, checked `observation`, `positions_sha256`, `assurance`,
`consumer_binding`, producer, UUID and digest. Its typed historical reader checks
exact relationships, bounds, counts/state/frame and known experimental version.
Consistently reauthored history is not authenticated; actual use checks its real
source again. Whole-frame SHA-256 is ordered finite binary64 nm JSON with signed
zero normalized, not a tolerant matcher or molecular identity. No molecular
coordinates are duplicated into the Site/application.

Construction uses separate `bindingmt.site_geometry_use` 0.1, recording source/
normalized frame/state, checked original/local atom axes, full-source and subset
position witnesses and the exact application ref. Final staging uses separate
`bindingmt.site_geometry_payload_use` 0.1 in the existing receptor artifact:
actual payload SHA/map, checked axes, effective backend box, coordinate precision
and remaining caller assumptions. No support-membership field is invented.
Readers preserve historical producers/checks and require actual problem context
for numerical export. Replacing an application requires a new problem.

## Executed source qualification — 2026-11-15

BindingMT base is `fc4e8d9d6af531e959f5aa541211dbf14276da06`, plus this source
increment. `_private/site_geometry_application.py` owns the new consumer record,
original-frame checks and effective-box admission; existing `site_application.py`
shares the owner preparation/writer-map/byte checks. `core/site.py`, existing
BoxRegion, DockingProblem and Vina route both explicit application discriminators.
No molecular parser, selection/inclusion/alignment kernel, replacement writer,
new molecular snapshot class or optional scientific dependency is introduced.
A native-only public MolSys class gate refuses unsupported inputs before generic
form-adapter discovery. All actual molecular queries/maps remain MolSysMT-owned.

**1191 selected cases pass under each normal/pm-fs policy**:

| Population | Cases per policy |
| --- | --- |
| Existing/public-backend application controls | 63 |
| Declared-geometry application controls | 70 |
| Selection snapshots | 129 |
| Public region projection | 70 |
| Earlier declaration/membership/geometry | 195 |
| Resource/helper/serialization/internal geometry | 463 |
| Arguments/bootstrap | 37 |
| Legacy search-domain/box/protocol | 164 |

The 133-case application module contains 60 preceding cases and 73 added cases,
including three expanded public-backend cases. Two disjoint complete selectors
qualify that module with lower process footprints; their node union is verified.
Deselected nodes belong to the other declared selector, not skipped evidence.
All final cells have no failures/errors/skips. Before/after runtime/test/resource/
provider/native hashes match, including original frozen fixtures. Two existing
NumPy shape warnings in the 63/463 populations and three legacy-H5MSM warnings in
the 164 population remain reported. This is a selected consumer qualification,
not two independent whole suites or a new public/platform admission.

Controls include actual source 2/local 0, receptor subsets/reordering, frame/state/
budget contradictions, independent and unknown-state declarations, immutable
recapture, changed source/normalized/other-frame/chemical data, strict resealed
history, actual prepared/written bytes/maps, generated and missing atoms,
wrong written coordinates despite a consistent payload hash, quantized-box
refusal and a receptor wholly outside the declared shape. Existing ordinary
box/native-membership behavior is received alongside the new route. Actual public
Vina calls use software engine sentinels after successful staged admission;
stale actual inputs refuse earlier. Native docking/preparation-method accuracy
is not executed.

An additional eight authored same-problem controls use both shape kinds,
None/positive padding and both unit policies. They retain the actual staged
receptor artifact and application receipt in an ordinary empty DockingResult
built from the same problem. Exact unpadded geometries pass when representable.
These eight are separate from the 1191 pytest count. Actual JSON files read in
a fresh Pint/unyt process with all molecular source/query/build routes blocked,
zero molecular queries and unchanged saved bytes/producers. Numerical export
still refuses without a problem. The actual archived fc4e8d9 decoder rejects
the new recipe and retains plain-box export; its temporary runtime export retires.

The provider closure is MolSysMT `6deb0200cd7e384128eded46caaf751a57e65bc6`,
PyUnitWizard `cc51c8e600f78f108bc563edd3221bd8788ead86`, ArgDigest
`7be790b7334e16715fceebbfa8a2b400edf768bb`, SMonitor
`4a4cd1e415ea57e50c22925c36c2e727772fb3ab` and DepDigest
`5e93eaff314fe7cddf4eee89842977023b3a8b00`. Actual development origins and mapped
MolSysMT native SHA-256 remain
`2dbcf13d148cb3ae13a5908707fcdd61c75a2de5027e50c566ce07de81bdb038`.
No provider edit/merge/install/build/native replacement occurs. This extends only
measured operations, not native builder/source authenticity or the whole ABI.
Existing generated editable-version metadata is hashed separately from committed
source; exact Git commits identify this increment. The 164-case population also
binds the unchanged bundled alanine H5MSM asset; no new shipped helper/data input
or resource declaration is introduced.

Initial attempts preserve a wrong test containment method and a test factory's
duplicate structure keyword. The 133-case combined attempt then reaches complete
XML with one unsupported-object expectation failure but is **RSS-interrupted**
before completion: 562823168 sampled bytes against the 512 MiB group threshold.
Progress/XML do not qualify it. Available RAM remains above 43 GiB, but the local
limit is respected and only owned work is stopped. Generic form discovery on
that unsupported object loaded unnecessary adapters. The native-only consumer
now uses the owner's public class gate; it does not change MolSysMT detection or
unwrap provider errors. A tiny guarded selector receives that correction first.
After a fresh resource review, use disjoint 63/70-case application selectors,
CPU6/nice10/one worker and unchanged CPU45/wall180/RSS512MiB/AS2GiB/reserve13GiB.
The interrupted combined invocation is never restarted or counted. All final
smaller populations terminate without interruption. Earlier native498/evaluator
holds stay unretried. An initial supplemental launcher missing the explicit
unshipped helper path is retained separately; the final launcher records that
existing helper's exact hash and succeeds without installation.

Six applicable workspace and fresh committed-export checks cover Ruff/format,
report indexes, proposal syntax/references, shared exact resources and Conda
contract. The source export retires. Raw attempts/XML/programs/source/data hashes,
resource/origin/saved-record/old-reader/publication receipts retain #50 custody
under `.molsyssuite/geometry-placement-111-20261115/`, review 2026-11-22.
Hosted CI's separately pinned route, full-suite/installed/public/platform/native
engine/final-pair and shared/consumer gates remain unexecuted. No tag/package/
release or manual remote CI occurs; skipped CI is not a pass.

Next receive TopoMT #199's finite input/view/context/attribution decisions before
T0; PeptideMT #60 retains its independent P handoff. Independent-target transforms,
oriented geometry, partial prepared-frame retention, ensembles and stronger
finding assurance remain explicit future profile work, not implicit extensions
of this source/control receipt.
