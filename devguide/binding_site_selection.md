# BindingSite native selection snapshot: first profile

Owner: [BindingMT #111](https://github.com/uibcdf/bindingmt/issues/111).
Admission and source delivery: 2026-11-14.
Status: **native selection snapshot/typed historical admission implemented;
81 new selection cases; provider resource-exception identity pending MolSysMT #756**.

The public operation is `BindingSite.from_selection(molecular_system, selection="all", *,
source, element="atom", syntax="MolSysMT", structure_indices=None,
chemical_state=None, name=None, role="selected_members", extensions=None,
max_source_atoms=100000, max_structure_occurrences=4096, max_dense_pairs=1000000)`.
The supplied source descriptor associates the original resident native MolSys;
it is not authenticated by counts, labels or the selection result.

Keep a **fixed evaluated membership**, never a live rule or binding finding.
The result uses a separate experimental Site 0.4 profile retaining original
selection/syntax/element, effective state and structural policy, checked axis
counts and actual ordered members. Earlier 0.1/0.2/0.3 bytes/producers remain
unchanged. The existing immutable lifecycle and historical reader admit this
profile without executing selection or accessing a molecular source.

With `structure_indices=None`, pass the source's public native Topology and
ChemicalStates domains to public `msm.select`, withholding Structures entirely.
That lets MolSysMT evaluate its coordinate-independent inventory/chemistry syntax
and reject coordinate-requiring queries itself. BindingMT implements no parser,
dependency classifier or geometry kernel. Structural scope remains independent;
the captured chemical-state dependency is explicit. Known explicit/reference
state positions are required. `chemical_state=None` means reference only in
this inventory mode, recorded explicitly rather than inherited from a provider.

With one explicit source structure, query the complete source using that single
occurrence. `chemical_state=None` means its actual association (`structure`);
an explicit state/reference must agree. Refuse unknown associations and multiple
occurrences rather than aggregating or choosing the first. Preserve expression
units/PBC syntax exactly; distances/contacts remain MolSysMT-owned. An opaque
structure-bound expression is conservatively budgeted for the full source
atom-pair square before execution. Numeric/all inputs need no pair allocation.

Admit the MolSysMT syntax, an ordinary nonempty string without external `@`
bindings, or unique nonnegative scalar/flat integer indices. Nested selections,
boolean inputs, implicit external variable capture, duplicate membership and
unsupported syntax/elements refuse; there is no hidden union or deduplication.
Only atom/group axes are admitted. Bound atoms/groups, structure/state counts,
request and result axes before substantial work. Recheck counts, effective state
and the selected result before publication; this is not a universal molecular
snapshot or a guarantee against concurrent coordinate/chemical changes.

The payload's typed selection receipt is historical checked input. It stores no
coordinates, native handle or optional engine output. `inspect()` distinguishes
captured selection from fresh molecular projection/current applicability.
`resolve_region` checks captured atom/group/state collection counts against its
original source and requires the captured chemical-state position to match its
selected source/target state, even on a structure-independent selection; it
does not rerun the predicate on a new structure or silently transfer its meaning.
Full source/revision and consistent check/use remain caller preconditions.

Qualify real inventory and structure/state/spatial owner operations, known empty
groups/empty outcomes, original request/order/identity, wrong scopes/states/budgets,
actual saved-Site and region replay, immutable derivation and existing records.
Use normal and pm/fs plus Pint/unyt contexts on small controls, one CPU/worker,
low priority and the monitored 13 GiB reserve/512 MiB owned RSS limit. Retain
source/native/XML/resource and workspace/fresh committed-export/publication
receipts under `.molsyssuite/site-selection-111-20261114/`, #50 custody, review
2026-11-21. Existing native498/evaluator holds are not restarted. Public C2/M2,
richer addressed/ensemble resolution, shared/consumer and full-suite/installed/
public/platform/final-pair gates remain separate; no new tag/package or manual
remote CI is authorized by this bounded increment.

## Public use and fixed meaning

```python
import bindingmt as bmt

source = {"owner": "application", "id": "original-system", "revision": "r1"}
inventory_site = bmt.BindingSite.from_selection(
    original_molsys, "group_name=='LIG'", source=source, element="group",
)
scoped_site = bmt.BindingSite.from_selection(
    original_molsys, "formal_charge>0", source=source, structure_indices=[2],
)
region = scoped_site.resolve_region(original_molsys)
selected_atoms = region["projection"]["exported_target_atom_indices"]
```

The strings are software examples, not recommended binding criteria.
The required source descriptor is the caller's association to its original
resident MolSys, not a computed/authenticated molecular identity. In inventory
mode the selected state is explicit/reference and structural scope is independent.
With a structure, the selected known association is mandatory; an explicitly
incompatible reference/integer refuses. A structure-dependent expression is
evaluated at most on that one supplied occurrence. No ensemble union, first-frame
default or implicit physical threshold is introduced.

The single membership view has `role="selected_members"` and a payload containing
the existing ordered `element`, `index_space`, `source_indices`, plus
`resolution={"state":"resolved","cause":"captured_selection"}` and `selection`.
The latter has exactly `request`, `effective`, `observations`, `provider`,
`membership_policy="fixed_snapshot"`, and `current_applicability="not_assessed"`.
Request retains the canonical query/syntax/element/context and limits; effective
records inventory versus single-structure mode, state policy/position, selected
occurrence and the conservative pair budget. Observations retain checked source
axis counts/reference/selected association; provider records MolSysMT/version.
The root producer is `from_selection`; immutable revision keeps the capture and
original evidence while updating interpretation through the existing lifecycle.

`inspect().selection_capture` exposes that history; current molecular projection
remains `not_performed` until requested. The result does not promise that a
predicate remains true after coordinate/chemical edits. Reading and projecting
fixed membership never reevaluate the original expression; current-finding,
chemical-equivalence and source-authenticity assertions remain unavailable.
The reader checks typed owned relationships and explicit numeric/all ordering;
it does not reparse or execute an opaque historical predicate. A consistently
reauthored and resealed record is not authenticated merely by its digest.

Use existing Site `save`/`read`/`from_dict(..., experimental=True)`; no new IO
aliases or source snapshot class are needed. Actual region output can nest this
0.4 Site through the existing experimental region result. Earlier profile bytes,
producer identities, source tag 0.1.0 and public-delivery status are unchanged.

## Resource-diagnostic provider gap: MolSysMT #756

[MolSysMT #756](https://github.com/uibcdf/molsysmt/issues/756) owns a reproduced
defect in public expression selection: `MemoryError` and `OSError(errno.ENOMEM)`
from the expression backend become `ArgumentError` with a parsing message; their
original object survives only as `__cause__`. `KeyboardInterrupt`/`SystemExit`
already retain identity. The consumer does not unwrap causes, replace the
selector or suppress the provider diagnosis. Its direct boundary preserves
whatever reaches it, but **end-to-end original resource-error identity through
expressions is unqualified pending the owner fix and client receiving**.

The finite one-atom probe injects exception objects; it does not allocate until
failure or stress the machine. Its completed reproduction is defect evidence,
not a passing resource-identity gate and not part of the software pass count.
Executed provider source is `27ba3d4c8d63a6a40a2d77f0190309ff15f4deef`; the same
public selector file remains at fetched origin
`baedf4c702e6b9438a63b097e1b048bf1ab3df6d`. No source merge, provider edit, native
replacement, environment installation or build was performed.

## Source qualification and limits — 2026-11-14

The implementation starts at BindingMT `d8143e205245c82e9fe5eff49d4cf8b338b57f28`.
`core/site.py` supplies the factory and existing immutable lifecycle;
`_private/site_selection.py` arranges public resident owner domains/queries and
data-only capture checks; `site_records.py` admits explicit 0.4 separately;
`site_region.py` enforces captured source counts/state without predicate replay.
No local molecular parser, manipulation, contact or geometry implementation is
added. Limits apply before owner-dependent work: atoms/groups default 100000,
structure/state axes 4096 and opaque structure-bound atom-pair square 1000000.
Capture repeats the query and observed context serially; these are consistency
witnesses, not an atomic snapshot or a bound on every possible parser operation.

Actual source closure remains MolSysMT `27ba3d4c8d63a6a40a2d77f0190309ff15f4deef`,
PyUnitWizard `cc51c8e600f78f108bc563edd3221bd8788ead86`, ArgDigest
`7be790b7334e16715fceebbfa8a2b400edf768bb`, SMonitor
`4a4cd1e415ea57e50c22925c36c2e727772fb3ab` and DepDigest
`5e93eaff314fe7cddf4eee89842977023b3a8b00`. Existing native SHA-256 remains
`2dbcf13d148cb3ae13a5908707fcdd61c75a2de5027e50c566ce07de81bdb038`, using the
preceding actual-process origin receiving. This does not authenticate the native
builder/source or whole ABI. Official suite status fetches newer commits without
merging them into the tested closure.

Final complete selected modules total **846 cases per normal/pm-fs policy**:
81 selection (including four direct owner facts), 70 public region, 195 earlier
site declaration/membership/geometry, 463 resource/helper/serialization/internal
C2, and 37 argument/bootstrap. Selection assertions cover real inventory chemistry,
known empty groups, no-structure input, all/numeric ordering, one-frame spatial
queries with equivalent nm/angstrom expressions, state/axis contradictions,
coverage, budgets, fixed meaning after coordinates change, immutable revision,
typed resealed data refusals and Pint/unyt application-policy preservation.
The two policy repetitions are not independent full suites. The normal and pm/fs
463-case groups retain their two existing NumPy shape deprecation warnings.

Actual public Site 0.4 and nested region JSON saved by the tests are admitted
in a separate fresh process under Pint/unyt with molecular source routes blocked:
zero molecular queries, no molecular input, four public data-only map admissions
and unchanged saved bytes. No optional scientific engine is imported. Source,
test, provider, native and original frozen input bytes are checked before/after
each final cell and against the committed payload. The early owner probe had a
fixture using an unsupported setter (corrected); a second five-assertion probe
is rejected by its source-drift check because runtime files changed during the
cell. Neither preliminary attempt qualifies a fixed source. Later complete
77- and final 81-case populations are separate, not added to the final counts.

All cells run serially, CPU0/nice10/one native/build worker, monitored 13 GiB
available reserve, 512 MiB group RSS threshold, 2 GiB address space and
45 CPU/120 wall seconds. No new resource interruption is admitted as a pass.
Earlier native498/evaluator holds remain unretried and unqualified. Retain raw
source/XML/resource/probe/saved-JSON/replay and workspace/fresh committed-export/
publication evidence under `.molsyssuite/site-selection-111-20261114/`, #50 custody,
review 2026-11-21. The six applicable checks cover Ruff, format, report indexes,
proposal syntax/references, shared resource declaration and Conda contract;
the committed export is retired. Hosted CI/full-suite/installed/public/platform/
final-pair and cross-component consumer acceptance remain unexecuted.

Next receive #756's owner correction when available and admit the public C2
placement/assurance profile before implementation, keeping the existing internal
geometry recipe distinct. Richer addressed membership, mixed/ensemble operations
and T/H/P/shared/consumer gates retain their owners. #111 stays active; this is
a useful bounded C1 selection increment, not full plan completion or a release.
