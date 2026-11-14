# BindingSite native region projection: first C1 profile

Owner: [BindingMT #111](https://github.com/uibcdf/bindingmt/issues/111).
Admission and source delivery: 2026-11-14.
Status: **bounded public native atom/group projection and historical reader implemented;
70 new controls and 765 complete selected cases per unit policy received**.
Prerequisites: [the actual #746/#750 client receiving](binding_site_molecular_contract_receiving.md).

## Public operations and caller association

The public method is `site.resolve_region(molecular_system, *, target=None, atom_map=None,
group_correspondence=None, structure_correspondence=None, view=None,
structure_index=None, target_structure_index=None, chemical_state=None,
expansion=None, coverage_policy="require_complete", context_policy="require_match",
max_source_atoms=100000, max_structure_occurrences=4096)`.

The positional molecular system is the resident **original native MolSys**
described by the site's source declaration. That association, unchanged relevant
source/revision/chemical contents and consistent inputs throughout the operation
are caller preconditions. Counts, labels or internally valid maps never prove
common origin, revision identity or a current binding finding.

`target=None` queries that source directly. A separately supplied native target
requires an explicit owner-supported atom-map report and target-to-source
structure-occurrence list. Missing maps are not identity maps. There is no
conversion, molecular reconstruction, coordinate copy, parser, matcher,
preparation or numerical geometry/contact operation. Use public MolSysMT
get/select/compose routes and preserve the existing immutable declaration.

`BindingSite.region_from_dict(record, *, experimental=False)` is the static
typed admission route returning a detached Mapping. It serves actual saved-result
replay without a new public molecular class or redundant file-IO aliases. Ordinary
JSON save/load plus this explicit reader suffices. Historical replay checks owned
relationships and public owner map data, without molecular source queries,
reconstruction or current-use reassessment.

## Admitted views, structure and chemical scope

Only the source-index atom/group membership view is eligible. It defaults to the
sole eligible `membership` view. Residue-address annotations, criteria and
geometric inclusion remain unsupported rather than being guessed.

`structure_index` names the original source structure. Omission can use the
single distinct declared source structure, or structure 0 on a one-structure
source; otherwise it refuses. Preserve the entire original declaration scope,
including repeated references. An explicit target's structure list has one source
index for every target-local occurrence. Build its inverse relation, preserving
repetitions. If the chosen source has several target occurrences, require explicit
`target_structure_index`; never silently choose the first. Record which declaration
occurrences share the chosen source and which target occurrence is actually used.

The first profile requires known selected source/target state associations.
`chemical_state=None` or `"structure"` uses those associations, never reference
state substitution. An explicit integer must agree. Selected known declaration
states must match. The admitted derived operation preserves state collection
positions/cardinality; check positions while retaining that caller precondition,
without calling it chemical equivalence. Mixed or missing required context refuses.

## Membership, maps and coverage

For atoms, omitted expansion means `listed_atoms`; no group enlargement occurs.
For groups, require explicit `source_members` or `target_group_members`:

- `source_members` transports the original group's atom members in entry-aligned
  source order through the checked atom map.
- `target_group_members` queries the corresponding current target group, retaining
  generated target atoms with absent original identity and no inherited evidence
  or confidence. It never hides removed original members.

MolSysMT owns group expansion. Query one atom list per requested group, preserving
group entry order; expressions/masks are not ordered expansion receipts. Unknown
parentage cannot certify complete/empty expansion and refuses this first profile.
Known empty groups remain distinguishable from omission and a negative finding.

Use `msm.compose_atom_maps` for typed directions, counts and inverse admission.
This profile admits retained/omitted/generated atoms, without merged atom/group
identity. Derive group correspondence only from aligned surviving atom parents.
An optional source-to-target `group_correspondence` integer vector uses -1 for an
explicit omission; it must agree with every surviving witness and be injective.
Requested groups lacking any surviving witness require that explicit vector,
including empty groups; do not infer their absence or identity from labels.

Default complete coverage refuses any known omitted requested source atom/group
before indexed target lookup. `allow_partial` retains every omission, entry and
weaker coverage state; it does not excuse contradictory/missing correspondence
or unknown parentage. -1 is never passed to a molecular lookup. Stable union of
entry outputs gives the exported target atom order. Exported here means returned
as a projection, not observed engine use or enforcement.

Read only relevant counts, source selection, atom parentage, selected state associations and group
membership, and recheck those facts before publishing. This detects observed
changes in those witnesses, not all coordinate/chemistry changes or concurrent
mutations. Counts/parents and map/occurrence metadata work is bounded by the two
explicit budgets, checked before map normalization or full-axis queries. The atom
budget applies to both atom/group axes; the occurrence budget bounds structure
and state-collection axes, declared references and correspondence vectors. Map
nesting is limited to 32; total stored map data is bounded before owner axis
normalization. Increase limits explicitly
for a separately resource-admitted larger workflow, never parallelize implicitly.

## Experimental result and replay

The method emits a detached finite-JSON `format="bindingmt.site_region_projection"`,
`format_version="0.1"`, `experimental=True` record with UUID/content digest,
complete unchanged site and exact definition/view reference, request/effective
choices, source/target observations, original owner map data/canonical map,
group/structure correspondence, entry-aligned requested/resolved/exported indices,
known omissions/generated identities, coverage and checks/caller preconditions.
Retain original claim/evidence refs without transferring their meaning to generated
members or reporting new binding support. Frame/coordinates, current finding,
chemical equivalence, binding support and consumer execution remain unassessed.

The typed reader requires explicit experimental opt-in and exact version/fields.
It verifies content integrity, original site refs and the owned projection
relationships. Public MolSysMT map validation is data-only and permitted during
replay; get/select/extract/convert/build/contacts and optional tools are not.
Retained observations are historical checked inputs, not authenticated molecules.
JSON object order cannot decide atom/group/occurrence order. Earlier Site records
and producers keep their versions/bytes; no shared `@1` freeze occurs.

## Ordinary use and saved results

The source descriptor is a caller association, not a loaded molecular object:

```python
import json
import bindingmt as bmt

site = bmt.BindingSite.from_atoms(
    source={"owner": "application", "id": "my-original-system", "revision": "r1"},
    atom_indices=[5, 1, 2],
)
region = site.resolve_region(original_molsys, structure_index=2)
indices_on_original = region["projection"]["exported_target_atom_indices"]
saved_json = json.dumps(region, allow_nan=False)
restored = bmt.BindingSite.region_from_dict(json.loads(saved_json), experimental=True)
```

For groups use `BindingSite.from_groups(...)` and explicitly choose
`expansion="source_members"` or `"target_group_members"`. A separately supplied
target uses `target=...`, an owner atom-map `atom_map=...`, and
`structure_correspondence=[...]` naming the original structure for every target
occurrence. `target_structure_index` disambiguates repeated occurrences.
`coverage_policy="allow_partial"` retains omissions; inspect `projection.coverage`
and every entry rather than treating an empty exported axis as no binding.
`observations.atom_group_indices` stores unknown parentage as -1; that sentinel
never becomes an index lookup. Generated output source identities are JSON null.

The result's exact root fields are `format`, `format_version`, `experimental`,
`record_id`, `site`, `definition_ref`, `view_ref`, `request`, `observations`,
`correspondence`, `projection`, `producer`, and `content_sha256`. Original owner
report and normalized map remain separate. The projection retains effective
axes, entry-aligned source/output identities, missing/generated members, exported
order, coverage and the explicit assurance limits. Only this discriminator/version
is accepted by the typed reader; ordinary Site IO keeps its existing profile.

## Executed qualification — 2026-11-14

The source increment starts at BindingMT `7b11bab5989a60f9bab530530826fa59f7acdfc8`.
Runtime is `core/site.py` plus the owned private `site_region.py` arrangement/
record checks and the existing ArgDigest registry. All selection, group expansion,
state/count inspection and typed atom maps use public MolSysMT operations.
No sibling source, native file, environment or scientific input is changed.

Actual provider source closure remains MolSysMT
`27ba3d4c8d63a6a40a2d77f0190309ff15f4deef`, PyUnitWizard
`cc51c8e600f78f108bc563edd3221bd8788ead86`, ArgDigest
`7be790b7334e16715fceebbfa8a2b400edf768bb`, SMonitor
`4a4cd1e415ea57e50c22925c36c2e727772fb3ab`, and DepDigest
`5e93eaff314fe7cddf4eee89842977023b3a8b00`. Loaded development native bytes remain
SHA-256 `2dbcf13d148cb3ae13a5908707fcdd61c75a2de5027e50c566ce07de81bdb038`,
whose actual-process origin was received in the preceding owner-fact population.
This increment uses that existing native context without claiming builder/source
authenticity, the whole ABI, or a different public artifact. Official suite-status
fetch observes newer provider commits but does not merge them into this tested
closure. Their receiving is separate work.

Complete selected modules pass **765 cases per policy**, normal and pm/fs:

| Population | Cases per policy | Scope |
| --- | --- | --- |
| `tests/test_site_region.py` | 70 | Actual public native projection, generation/omission, record admission and refusals |
| Existing site declaration/membership/geometry modules | 195 | Earlier factories, records, identities and quantities |
| Resource/helper/serialization/internal C2 modules | 463 | The complete previously admitted lower-resource selection |
| Argument/bootstrap modules | 37 | Public call contracts and fresh-process application units |

This is one selected source qualification repeated under two policies, not two
independent full suites. The 70-case module includes Pint/unyt application-policy
preservation and an unchanged original 0.2 specimen/producer. A separate fresh
process admits the **actual JSON saved by the public API test** under Pint and
unyt with molecular get/select/extract/convert/build/structure routes blocked:
zero molecular queries, four public data-only map admissions, unchanged input
bytes and no unrelated scientific tools imported. It requires no molecular input.
The first 40-case attempt had 39 passes and one incorrect expected exception
type; MolSysMT correctly refused an unknown chemical association with its own
StructuralInconsistencyError. The test now preserves that provider refusal.
Later complete 40- and 70-case populations pass. Initial partial populations
are not added to the final counts. Two pre-existing NumPy shape deprecation
warnings occur per 463-case selection; no new test failure or resource
interruption is suppressed.

All cells run serially on CPU0/nice10 with one native/build worker, monitored
13 GiB available reserve, 512 MiB group RSS threshold, 2 GiB address space and
45 CPU/120 wall seconds. Observed group RSS stays below 318 MiB and available
memory above 49 GiB. All owned groups are terminal. Earlier native498 and evaluator
interruptions remain held, unretried and unqualified. Source/provider/native and
frozen-input before/after bytes, XML, actual saved JSON, replay, resource outcomes
and applicable workspace/fresh committed-export checks retain #50 custody in
`.molsyssuite/site-region-111-20261114/`, review 2026-11-21. Checks cover Ruff,
format, report indexes, proposal syntax/references, shared resource validation
and Conda contract. Remote CI, installed/public delivery and cross-platform
qualification are not executed. Existing source tag 0.1.0 is unchanged.

## Qualification scope and next gates

This admits one selected known-state native membership projection at a time.
It rechecks observed witnesses; there is no universal source snapshot, mutation
observer, current-finding reassessment or inferred chemical equivalence. The
small controls exercise original coordinate/declaration nonmutation and owned
map/axis/parent/state observations, not every possible molecular attribute.
Historical valid data can be consistently reauthored and resealed; integrity
checks are not authentication or a fresh molecular assessment.

Next implement `from_selection` through MolSysMT and qualify inventory versus
spatial/state-dependent selection and retained declarations. C1's richer
addressed resolution, public C2/M2, placement/current-finding, T/H/P and shared/
consumer/full-suite/installed/public/platform/final-pair acceptance keep separate
gates. No heavy held population is automatically restarted. This bounded C1
source delivery does not close #111 or advance the source tag/package version.
