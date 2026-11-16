---
summary: Receive the latest-policy rule, update the remaining policy caller and synchronize the common guide through its owner tool.
issue: uibcdf/bindingmt#119
status: resolved
opened: 2026-11-15
closed: 2026-11-15
verification: inspected
area: [governance, ci]
guard:
normative: https://github.com/uibcdf/molsyssuite/blob/6e728809064e74860455963279dcc02327a63f5a/devguide/adoption_lifecycle.md
blocked_by: []
supersedes: []
---

# Latest published MolSysSuite policy receiving

## What

[MolSysSuite #201](https://github.com/uibcdf/molsyssuite/issues/201) requests
current-policy receiving through BindingMT #119. Technical compatibility with an
older immutable policy is distinct from adopting the latest published release.
The target is policy-v1.5.12, tag object
`7294f684dbaf3f0267c6e978937475d3deb8f633`, peeled source
`bb641488bc22b262549f795857d6ede529b14ca1`.

## How

Review both actual policy job calls, registered guides, local instructions and
independently qualified workflow/SDK pins. The ordinary policy caller already
uses policy-v1.5.12. The current-source qualification workflow's policy job still
uses policy-v1.5.4 and needs the same current target. Its historical executions,
wheel producers and scientific/installed/operator pins retain their identities.
No extra admission input is needed for BindingMT, present in the target registry.

Receive the published guide through the official registered synchronizer, using
an isolated central-source checkout and the real consumer checkout. Do not edit
the guide copy manually or update the human's central/sibling worktrees.

## Why and evidence boundary

The official before inventory freshly resolves the published tags and reports
`compatible`, observed policy-v1.5.12 plus policy-v1.5.4, with exit 1 under the
new rule. This is expected current-adoption evidence, not a retroactive failure
of the old policy's historical execution contract. Runtime/scientific and
installed/public/platform/hosted evidence remain separate.

## Acceptance criteria

Both actual callers use the verified latest immutable target, required inputs
are reviewed, official guide sync and fresh local/published-source inventories
agree, and relevant central/local administrative and exact-resource checks
pass. Return the exact owner disposition and still-pending evidence to #119/#201.
Preserve qualified pins, instructions, frozen resources, native bytes and the
0.1.0 tag. Use the authorized [skip ci] route without a manual CI dispatch.

## Review completed: delivered source and semantic scope

The remaining actual policy job in
`.github/workflows/qualify-current-0.1.0.yml` now calls policy-v1.5.12. The
ordinary `.github/workflows/molsyssuite-policy.yml` caller was already current.
No extra admission SHA is needed: BindingMT is registered in the immutable
target, and the target reusable workflow's only admission input is optional.
The change is exactly one reusable-job policy pin, not replacement of the
qualification workflow's scientific/build/installed recipes. Its source equality
request, version/wheel requirements and independently pinned source consumers
retain their existing contracts and histories.

The official `sync_vendored_guides.py` receives only
`MOLSYSSUITE_GUIDE.md` from an isolated, clean central checkout at published
`4ba79c9562716fc05e8040e1be8c387799a96e19`, after checking it against remote
main. The notice's originating source is
`6e728809064e74860455963279dcc02327a63f5a`; the guide bytes agree with both.
Old guide SHA-256:
`0188fee9330856bb1d40bdd5739c275a73471b602df20ba3a3f0ac8476677678`;
received SHA-256:
`2e13d5659ea1ce0f499c0935f5d6df25b1915781ef69ac52787120c351ac6546`.
No consumer guide was edited manually. The semantic diff adds the latest-policy
obligation, distinguishes historical compatibility from current adoption and
routes allostery-specific work to AlloMT, whose initial source foundation
keeps scientific methods pending qualification. This
requires no new BindingMT runtime or scientific-ownership change.

Root `AGENTS.md` and `devguide/AGENTS.md` are reviewed and unchanged: they already
require common-guide precedence and owner reporting. The maintainer's direct
commit/[skip ci]/no automatic manual-CI route remains explicit human authority;
missing remote evidence remains missing. No new exception or instruction edit
is needed for that authorized route.

Other local root-guide hashes, workflows, the current qualification operator,
runtime code and exact-resource declaration are inventoried before/after. The
SDK's official sync changes only the common guide; a separate inspection proves
that the workflow changes only the policy ref. In particular, SDK pins
`738fe8dd731abb99dac421b0bb2acc8564110180`,
`95f693be10c9205ab9b912cb53453710a207bf45` and the scientific/provider checkout
pins remain unchanged. Local review is not a new freshness certification for
all other canonical owner guides or completed support-library adoption (#19).
Frozen fixture guides/producers, scientific/native bytes, installed receipts,
human sibling checkouts and the immutable 0.1.0 tag retain their exact identities.

## Executed checks and receiving limits

The current official central inventory at source `4ba79c9` freshly resolves all
published official policy tags. Before correction, its selected BindingMT row
is `compatible`, observed policy-v1.5.12 plus policy-v1.5.4, exit 1. After
correction it is `current`, observed policy-v1.5.12 only, exit 0. The remote
immutable target is exactly the tag object/peeled source named above; no tag is
created or moved and no offline target override is used. The post-publication
receipt repeats this check on the committed consumer. Registered guide checking
also agrees with the owner bytes. This is one member's receiving, not closure
of the central all-member rollout.

The current source `check_repository.py` passes with no findings, using the
same `--skip-guide-content` division as the declared workflow; official guide
sync/check separately validates the changed common guide. This local live-source
checker run is **not execution of the immutable hosted policy job**. Its exact
published consumer commit is linked in
[#119](https://github.com/uibcdf/bindingmt/issues/119).

The existing reporting/CI-source/history/backlog administrative selection has
**20 passes**. Six normal workspace and fresh committed-export checks cover
Ruff/format, generated report indexes, document syntax/targets, exact shared
resources and the Conda contract. The fresh export is retired. The initial
administrative attempt preserves 19 passes and one stale-index failure before
the new report's final regeneration; it is not a scientific/resource failure.

No scientific calculation/test or held work is restarted, no native rebuild,
package install/tag/release/publication occurs and no hosted CI is dispatched.
The existing 268-case normal/pm-fs #120 and earlier scientific receipts retain
their original source scope and are not relabelled as executions on this guide
commit. #30/#47/#19/#100 retain unexecuted full-suite/installed/public/platform/
exact-source-hosted and broader adoption evidence. Skipped CI is not a pass.

## Custody and next action

Administrative cells use one logical CPU, nice 10, one native/build worker,
monitored 13 GiB available reserve, RSS 512 MiB/AS 2 GiB/CPU 45 s/wall 180 s.
Their source preflight/registration/resource receipts show no new interruption.
The controller remains unknown after the earlier sandbox socket denial; no
Slurm allocation or whole-workstation reservation is inferred.

Private source/tag/guide/workflow/admin/publication receipts remain under
`.molsyssuite/latest-policy-119-20261115/`, BindingMT #50 custody, review
**2026-11-22**. The task-owned central checkout and temporary workspace symlink
are retired after final published-source/guide verification; they are not a new
scientific/installed SDK pin. Qualified existing SDKs and shared environments
are retained.

Return **Review completed** with the exact public consumer/source/caller/guide
identities and executed/pending evidence to #119 and MolSysSuite #201. Central
owns its all-member observation/acceptance; #30 owns future exact-candidate
hosted receiving at an explicitly authorized checkpoint. Resume the independent
BindingSite/BindingInterface plan through the first ready finite owner case;
this administrative closeout does not complete #111 or #112.
