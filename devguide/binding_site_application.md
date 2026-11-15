# BindingSite native enclosing-box application: first C2 profile

Owner: [BindingMT #111](https://github.com/uibcdf/bindingmt/issues/111).
Admission and source delivery: 2026-11-14.
Status: **bounded public native C2 application and construction/final-payload
admission source-received; 60 new controls and 1070 selected cases per policy**.

The implemented operation is `site.to_search_domain(original_molsys, *, method, view=None,
structure_index=None, chemical_state=None, expansion=None, padding=None,
assurance_scope="placed_geometry", context_policy="require_match",
max_source_atoms=100000, max_structure_occurrences=4096)` for native atom/group
membership only. Admit `method="enclosing_box"`, explicit length-valued padding
and complete C1 coverage. Use public MolSysMT region/coordinate queries and the
existing BoxRegion points/bounds construction; introduce no selector, molecular
geometry/contact algorithm, source matcher or mutation observer. Empty support,
degenerate unpadded extents, geometric-only views, independent target/transform,
trajectory aggregation and stronger `current_finding` assurance refuse.

The guarantee is an axis-aligned Cartesian enclosure of the captured fixed
member positions in the original source's stored coordinates, with padding.
It is a search box, not an inferred binding volume or a new binding finding.
PBC is false; no images, alignment or unwrapping are inferred. Caller source/
revision association and consistent inputs during checking/use remain required.

Persist a separate experimental `bindingmt.site_domain_application` 0.1 recipe
in existing BoxRegion derivation: full C1 projection/original Site, request,
sealed support bounds/output quantities, selected-coordinate SHA-256, assurance,
normalization/lifetime conditions and losses. No molecular coordinates are copied
into Site or recipe. The position digest uses finite canonical binary64 nm JSON
in projected atom order, with signed zeros normalized. It is an exact numerical
witness, not tolerant source identity, chemical equivalence or authentication.
Changing representation to different binary64 values requires recapture even
when a geometric tolerance would consider them close; unit policy changes alone
must not change the stored native data or their normalized extraction.

Data-only reading validates owned recipe relationships and reconstructs the box
from retained bounds/padding. It needs no molecular source or original engine;
it cannot prove that historical bounds/hash came from an authentic source.
It must preserve original producers and never upgrade a read into a new check.
Earlier plain/selection boxes and unplaced offline geometry retain their behavior.

Admit the actual use at DockingProblem construction and again at numerical
backend export. Original source structure 2 is normalized structure 0; preserve
that relation and transport ordered support atoms through the problem's explicit
sorted original-source atom axis. Every support atom must survive. Recheck C1
membership/state/count witnesses, the original source positions and actual
normalized receptor positions, and recompute the box from the latter before
allowing export. Require the declared native normalization route, one normalized
structure and matching known state position. Prepared/file-only inputs cannot
supply this proof in the first profile. Edits to either original support or the
normalized support refuse, including edits that leave the enclosing bounds equal.

Record a detached historical `search_domain_application` construction receipt
on the problem, with exact application ref, source/local structure/state and
ordered support correspondence. Existing problem/result records preserve it;
reconstruction with explicit inputs must recheck and match it. Historical
serialization remains readable after live mutation, while export reassesses.
No receipt grants permanent validity. New boxes require
`to_backend_box(problem=...)`; missing actual problem context refuses. Vina's
normal adapter supplies that context before preparation; plain boxes retain
their old call path. Changing atoms outside the fixed support does not establish
or invalidate the unassessed stronger finding claim.

## Final prepared-payload admission

Before constructing Vina's native engine, retain one actual receptor owner export
and its map while staging its exact bytes. Compose the atomic native preparation
projection with the writer map using public MolSysMT. Count-check the bounded
reports before normalization; the first profile refuses nested map stages.
Every fixed support atom needs a direct retained identity: a merged parent is
not a substitute for an omitted support atom. Actual prepared full-precision
positions must match the captured numerical witness in projected order.

Use MolSysMT's public PDBQTFileHandler to read the mapped support from those
actual bytes. Require writer payload SHA-256, one model, matching atom counts,
and inclusive containment of the written positions in the **actual six-decimal
angstrom box** that the existing adapter supplies. No hidden enlargement or
coordinate tolerance is invented; insufficient explicit padding refuses when
serialization/grid rounding puts a written member outside. The receipt records
preparation/writer correspondence, payload hash and effective box in the
receptor backend artifact. It asserts no chemical equivalence or current finding.

Construction receipts refer to an exact application ID/content; replacing the
application requires reconstructing the problem. Export repeats the actual
source/normalization checks, including after preparation. Existing scoring and
other engine routes without this actual application context remain unsupported
by the new recipe. Ordinary old box paths keep their established behavior.
Tests use authored prepared representations/software engine sentinels; their
success is not native docking or scientific preparation-method qualification.

## Public use and historical records

```python
import bindingmt as bmt

site = bmt.BindingSite.from_atoms(
    source={"owner": "application", "id": "original", "revision": "r1"},
    atom_indices=[4, 1, 3], structure_indices=[2], chemical_state_indices=[1],
)
box = site.to_search_domain(original_molsys, method="enclosing_box", padding="0.2 nm")
problem = bmt.DockingProblem(
    original_molsys, prepared_partner, box, receptor_structure_index=2,
)
history = problem.search_domain_application
numbers = box.to_backend_box(problem=problem)
```

These are authored software-example indices/state/padding, not scientific
recommendations. Group membership requires an explicit C1 expansion policy.
`inspect().capabilities.domain_placement` reports support for the molecular view;
source access and actual placement remain unchecked until the explicit operation.
The first profile accepts complete nonempty original-source atom/group views
from Site 0.2 and captured selection 0.4. There is no independent mapped target,
placed saved sphere/box or current-finding adapter in this increment.

`BoxRegion.to_dict`/`from_dict` preserve the full application recipe. Its typed
reader checks version/digest/refs, C1 region and sealed bounds/padding/output;
current geometry must match to 8 binary64 epsilons, zero absolute tolerance and
without shape broadcasting. Support minimum/maximum reconstruct output exactly.
Integrity is not authentication: consistently reauthored history can be read,
but actual use recomputes its positions and enclosure against real inputs.

The optional problem field/property `search_domain_application` is present only
for this recipe and carries experimental `bindingmt.site_domain_use` 0.1:
application ref, source/local structure/state, normalized-source atom axis,
ordered support source/local correspondence and position digest. It is labelled
as a historical construction check, never permanent validity. Normal problem
reconstruction with explicit molecular inputs must reproduce it exactly;
missing/changed or boolean-for-integer receipts refuse. Old plain problem
records keep their existing fields. Ordinary results retain the whole recipe
and receipt in their existing problem-info plane. Reading history performs no
current application and grants no numerical export without the actual problem.

Vina's final receptor artifact gains a `site_domain_application` payload receipt
with actual payload SHA-256, composed owner map, ordered support transport and
the actual rounded angstrom box. Existing `capture_backend_inputs` controls
whether complete backend input bytes are retained; it remains false by default.
Qualification retains hash/map receipts, and owned temporary payload files are
retired. No source coordinates, graph parser or new universal snapshot class is
introduced. Old offline geometry is still unplaced and refuses backend use.

## Executed source qualification — 2026-11-14

Base BindingMT: `468ed09859a78b30388941fc0e01d3b87231aa44`. The implementation
composes `core/site.py`, private `site_application.py`, existing BoxRegion and
its factored point-bounds code, normal DockingProblem and the Vina export/staging
adapter. The existing shared sealed-quantity helpers accept an owned field
namespace; old defaults/recipes remain unchanged. Product source introduces no
provider-private molecular operation, selector, preparation or contact kernel.

**1070 complete selected cases pass per normal/pm-fs policy**:

| Population | Cases per policy |
| --- | --- |
| New native application/actual staged-payload controls | 60 |
| Selection snapshots and public regions | 81 + 70 |
| Earlier Site declaration/membership/geometry | 195 |
| Resource/helper/serialization/internal geometry | 463 |
| Arguments/bootstrap | 37 |
| Existing search-domain/contracts/box-admission/protocol modules | 164 |

These are one selected software qualification under two policies, not two
independent full suites. Actual native source 2/local 0, a nonidentity support
map, all supported Site versions, original/normalized/interior/outside edits,
missing support, exact historical receipts, reauthored records, quantity forms,
equivalent coordinate units and unsupported assurance are exercised. Actual
owner PDBQT writing and reading validates final staged bytes/maps, removed
support and zero-padding/rounding refusal. An actual public backend call uses
authored preparation/software engine sentinels: valid admission reaches only the
sentinel constructor; stale source/prepared input refuses first. No native Vina
calculation or preparation-method accuracy is executed or inferred.

Actual saved domain, construction receipt, empty-result and staged-map history
are read in a fresh Pint/unyt process with all MolSys queries blocked. Original
recipe identities/producers and saved file bytes remain intact; outer `nm` versus
`nanometer` labels may follow the existing quantity form without changing values.
Numerical export still refuses without a real problem. The actual archived
468ed09 decoder refuses the new application recipe, retains plain-box export,
and its temporary runtime export is retired. No prefix install/build occurs.

The provider closure stays MolSysMT `27ba3d4c8d63a6a40a2d77f0190309ff15f4deef`,
PyUnitWizard `cc51c8e600f78f108bc563edd3221bd8788ead86`, ArgDigest
`7be790b7334e16715fceebbfa8a2b400edf768bb`, SMonitor
`4a4cd1e415ea57e50c22925c36c2e727772fb3ab`, and DepDigest
`5e93eaff314fe7cddf4eee89842977023b3a8b00`. Existing native SHA-256 remains
`2dbcf13d148cb3ae13a5908707fcdd61c75a2de5027e50c566ce07de81bdb038`; this extends
only the measured operations, not builder/source authenticity or the whole ABI.
The 164-case compatibility population additionally pins the provider's bundled
22-atom alanine H5MSM example before/after; its three legacy-H5MSM warnings per
policy are preserved. Two NumPy shape warnings per application/internal-geometry
population come from deliberate shape controls; they are not hidden failures.

Initial attempts retain a wrong test QuantityRecord field (`value` versus
`values`), an overstrict cross-form unit-label assertion, and a wrong exception
namespace in the first cold-reader launcher. A 53-case attempt reached progress
53/53 but exceeded the 120-second wall limit: **resource interruption, not pass**.
After a fresh resource check, keep CPU45/RSS512MiB/AS2GiB/reserve13GiB unchanged
and allow bounded wall180. CPU0 was measured at 100% while other CPUs were idle;
final cells use only CPU6/nice10/one worker, including heavy descendants. A free
CPU observation is not a reservation. All final source selections are rerun in
full; interrupted progress is not included. Existing native498/evaluator holds
remain unretried. Recorded affinity/limits/worker counts and source/provider/
native/frozen-input bytes accompany every final cell.

Six applicable workspace and fresh committed-export checks cover Ruff, format,
report indexes, proposal syntax/references, exact shared resources and Conda
contract; committed export is retired. Raw attempts, XML, source/data hashes,
resource and current/old-reader/publication receipts retain #50 custody in
`.molsyssuite/site-application-111-20261114/`, review 2026-11-21. CI's separately
pinned hosted dependency route is not qualified by this local source pair.
Remote/full-suite/installed/public/platform/native-engine/final-pair and shared
consumer acceptance remain unexecuted; no new tag/package/release is made.

This receives a first bounded native C2/M2 **source/control** profile. Full plan
completion, richer addressed/ensemble and placed declared geometry, stronger
finding assurance and T/H/P/shared/consumer gates remain separate. Next obtain
the offered exact saved TopoMT profile/input before implementing T0/T1, and
receive MolSysMT #756's resource-diagnostic correction when available. No owner
finding or originating scientific method is silently rerun.
