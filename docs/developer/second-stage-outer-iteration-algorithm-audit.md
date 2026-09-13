# Second-stage decisions and evidence

## Status and authority

Current decisions below are checked against production source **C**:
`c174458294c1a328058f96bf37a2a21ae16c06a4` (`develop`). This documentation-only
cleanup changes no numerical policy, public interface, build option, test,
quality baseline, or tolerance and performs no new numerical experiment.

- [Second-stage local fitting](second-stage-local-fitting.md) specifies current
  execution, ownership, acceptance, recovery, and persistence contracts.
- [Second-stage phase audit](second-stage-phase-audit.md) explains optional
  observations, schemas, and diagnostic interpretation.
- [fold-168 iteration regression](fold-168-iteration-regression.md) remains an
  unresolved issue with its existing investigation and future repair proposal.
- This page owns current decisions, evidence limits, and historical references.
  Retired reports are available at fixed Git revisions, outside the current
  reading path; they do not override the current specification.

## Decision table

Every current decision applies to **C**. Evidence IDs link to the inventory below,
which records each investigation's separate baseline. Reopening conditions are
requirements for future work, not experiments or tests added by this cleanup.

| Issue | Current decision | Applicable version | Reason | Evidence limits | Reopen when | Evidence |
| --- | --- | --- | --- | --- | --- | --- |
| Convergence certificate | Retain qualification, accepted active p99, complete nominal operator p99, and blockers. | C | Small accepted steps do not certify a fixed point; unavailable endpoints fail closed. | No retained predicate is proved implied by the others. Maximum is not a certificate predicate, but has other algorithm uses. | A proposed simplification addresses each protected failure mode and establishes redundancy. | [E0](#evidence-index) |
| P0 and component ownership | Retain builder/consuming transaction and shared component helpers with separate outer/final removal policies. | C | Rejection preserves required lifecycle updates; shared assembly does not unify acceptance policy. | Focused comparisons do not cover every fallback, tie, or salvage branch or establish general numerical equivalence. | Ownership or publication changes require revisiting these contracts. | [E1, E4](#evidence-index) |
| Cluster-history observer | Retain observation-only history; keep `best_audit_state` in production. | C | Diagnostic availability must not affect acceptance, radius, quarantine, stop, or persistence. | Extraction neutrality does not establish neutrality of earlier member-best gate removal. | A decision starts depending on observer data, or observer ownership changes. | [E5, E13](#evidence-index) |
| Objective/recovery revisions | Keep independent ownership and existing partition/background triggers. | C | Objective reevaluation is not a retry event. | Captures lack applied partition changes and successful full-loop recovery; current trigger counts being equal does not make ownership interchangeable. | Revision ownership or trigger semantics change. | [E6](#evidence-index) |
| P1 finalization and radius | Keep converged-only final polish at independent radius 1.0; production Shrink with implicit Keep, no Grow. | C | Independent Grow removal passed existing tests with rescue retained. | Combined rescue/Grow removal failed intensity scaling and was reverted; it is not independent rescue evidence. | A separate policy proposal supplies independent evidence for the changed mechanism. | [E2a, E2b](#evidence-index) |
| Objective exhaustion | Exclude it from quarantine failure evidence; retain search rejection and other hard/invalid/guard evidence. | C | Search non-improvement alone is not freeze or recovery-failure evidence. | This was a numerical policy change: activity, patience timing, and later trajectories can change. Mixed-failure coverage is bounded. | New evidence warrants changing failure classification or recovery requirements. | [E7](#evidence-index) |
| P2 acceptance changes | Fixed-order atomic component acceptance remains withdrawn; greedy salvage stays. Current member-best removal remains under investigation. | C | Unavailable global baselines caused healthy remote updates to be rejected in two existing regressions. | P2 was a numerical ablation. Its objective-revision retry wording and fixed-false rejection fields are superseded; fold-168 remains unresolved. | A replacement defines acceptance with unrelated missing objective evidence, or a separate member-best repair is validated. | [E3, E13](#evidence-index) |
| Global-best gate | Production ON, including cooperative protection. | C | Historical-best and previous references need not imply the same decision. | Zero best-only rejection in 312 complete-state and 12 cooperative comparisons; OFF did not disable cooperative best. No released best-only candidate was observed. | Actual best-only release and retention consequences are evaluated; cooperative removal needs separate evidence. | [E10](#evidence-index) |
| Final dependency polish | Retain enabled polish and strict operator recertification. | C | Applied-path value and removal safety remain unmeasured. | First round: five converged executions, zero applied. Follow-up: 128 inputs, 15 entries, zero applied; the specified three-dataset matrix stopped at fold-168's failing baseline. | The baseline is repaired and the specified dataset evidence/removal conditions are completed, with actual application distinguished from provenance. | [E8, E9](#evidence-index) |
| Cooperative rescue | Retain enabled capability and shared component evaluation. | C | ON improves response MSE; OFF improves audit objective and uses one fewer attempt in the exposed fixture. | Accepted and retained rescue updates are observed; neither policy uniformly dominates. Passing tests do not prove redundancy. | A separate policy decision resolves the quality/objective/cost tradeoff with relevant evidence. | [E11](#evidence-index) |
| Background-trigger policy | Retain any-change Frozen recovery. | C | Stronger thresholds changed termination and persisted peeling despite equal Gaussian parameters. | No successful production release or applied-partition coverage; reducing failed retries does not establish safe removal. | Evidence covers those missing paths and persisted-output quality. | [E12](#evidence-index) |
| fold-168 | Unresolved; preserve the original <=25 iteration gate and all quality gates. | C; recorded run at E13 baseline | 100 accepted iterations, `maximum-iterations`; quality/atom-cutoff gates pass, iteration gate fails. | Historical member-best intervention is not a repair implemented in C. | A separately implemented repair passes both original iteration and quality gates with the validation described in the open report. | [E13](#evidence-index) |

## Unresolved fold-168 regression

The recorded boundary is adjacent commits `f50a742c` (11 accepted iterations)
and `49d3516a` (100, best iteration 27). Restoring only the historical-best
acceptance gates for local candidates and ordinary boundary members on the
latter version, retaining its new quarantine policy and original stopping
conditions, returns the run to 11 iterations with every original gate passing.
This intervention supports a trajectory effect from member-best removal;
it does not justify treating later observer extraction as the cause.

The recorded current-source run at `6587542638d570685a11ef281aebd015d5eb8c06`
also accepts 100 iterations. The CLI succeeds; the existing runner fails only
`accepted_iterations <= 25`. No run at a newer source is claimed here. Patience
can reset on improvement over the previous state without improving historical
best; selecting iteration 27 earlier or relaxing the threshold is not a
validated repair. Final-polish workload comparisons remain incomplete at this
baseline failure. See the [open investigation](fold-168-iteration-regression.md)
for the proposed independent repair, exact inputs, and subsequent validation.

## Evidence index

Inventory date: **2026-09-13**. All historical report links below are pinned to
**C**, which preserves their text before retirement. The baseline column refers
to the investigation, not the revision used to retrieve the report. Numerical
results remain those reported at their original baselines.

Local paths are relative to the repository root and intentionally use code
format. **Local present / untracked** means the directory and listed key records
were found in this workspace; it does not certify completeness, hashes, or
reproducibility. **Local not found** means absent at the reported path in this
workspace, not proved lost everywhere. No hashes or experiments were regenerated.
Existing manifests are listed only as pointers.

Git preserves report text, not ignored raw evidence. This index is not a backup.
Retain the local directories and external inputs for reproduction or transfer;
this cleanup neither copies/uploads evidence nor deletes local artifacts.

| ID / fixed report | Investigation baseline | Local artifact location | Key records / existing hash manifests | Availability on inventory date |
| --- | --- | --- | --- | --- |
| [E0 / Consolidated audit](https://github.com/yslianTim/rhbm_gem/blob/c174458294c1a328058f96bf37a2a21ae16c06a4/docs/developer/second-stage-outer-iteration-algorithm-audit.md) | a4354e698e77398154009d231907ebf3ed4b1d52 (earlier per-atom-offset review; later consolidation at C) | Not specified | Failure-mode/safeguard table and historical provenance in the fixed report. | No separate local artifact location specified. |
| [E1 / P0](https://github.com/yslianTim/rhbm_gem/blob/c174458294c1a328058f96bf37a2a21ae16c06a4/docs/developer/second-stage-p0-structure.md) | c5417717154896e53d31a1d43910cde15e3f80f6; interface follow-up 1041f13b6eefcc81ed8f02f41f9e918a5acb9f85 | `build/p0-structure/` | Build/test logs referenced by the report; no individual filename specified. | Local not found; Git report remains available. |
| [E2a / P1 combined removal](https://github.com/yslianTim/rhbm_gem/blob/c174458294c1a328058f96bf37a2a21ae16c06a4/docs/developer/second-stage-p1-ablation.md) | 6e7ac3c5f452812da8aeb1b1ae2497b5cd0143ec (report baseline) | `build/p1-ablation/` | Reverted combined-ablation evidence referenced by the report. | Local not found; Git report remains available. |
| [E2b / P1 independent Grow](https://github.com/yslianTim/rhbm_gem/blob/c174458294c1a328058f96bf37a2a21ae16c06a4/docs/developer/second-stage-p1-ablation.md) | 6e7ac3c5f452812da8aeb1b1ae2497b5cd0143ec (report baseline) | `build/p1-grow-only/` | Independent Grow validation referenced by the report. | Local not found; Git report remains available. |
| [E3 / P2](https://github.com/yslianTim/rhbm_gem/blob/c174458294c1a328058f96bf37a2a21ae16c06a4/docs/developer/second-stage-p2-structure.md) | f50a742cd5d8f3c64fa70cbc11436c7677ff0370 | `build/p2-structure/` | `report.md`, independent patches and failed-step evidence (reported, unavailable here). | Local not found; Git report remains available. |
| [E4 / Component assembly](https://github.com/yslianTim/rhbm_gem/blob/c174458294c1a328058f96bf37a2a21ae16c06a4/docs/developer/second-stage-component-assembly.md) | ac9ca72ce0e5093abb8e2c97966c41644b494053 | `build/component-assembly/` | `comparison.json`, `source-verification.json` | Local present / untracked. |
| [E5 / Cluster history](https://github.com/yslianTim/rhbm_gem/blob/c174458294c1a328058f96bf37a2a21ae16c06a4/docs/developer/second-stage-cluster-history-observer.md) | a58713d1d37fe326e55c1ccc5bc29928a45aae96 | `build/cluster-history/` | `report.md`, `trace-comparison.json` | Local present / untracked. |
| [E6 / Revision separation](https://github.com/yslianTim/rhbm_gem/blob/c174458294c1a328058f96bf37a2a21ae16c06a4/docs/developer/second-stage-recovery-revision-decoupling.md) | e6296fdf604779f7259c99326de1258158df492a | `build/recovery-revision/` | `matrix.json`, `comparison.json`, `trace-comparison.json`, `source-verification.json` | Local present / untracked. |
| [E7 / Objective exhaustion](https://github.com/yslianTim/rhbm_gem/blob/c174458294c1a328058f96bf37a2a21ae16c06a4/docs/developer/second-stage-quarantine-objective-exhausted.md) | d511c88fe80e0665a92f84a6ac13371a0b6635c5 | `build/quarantine-objective-exhausted/` | `report.md`, `behavior-comparison.json`, `trace-comparison.json` | Local present / untracked. |
| [E8 / Final polish, first round](https://github.com/yslianTim/rhbm_gem/blob/c174458294c1a328058f96bf37a2a21ae16c06a4/docs/developer/second-stage-final-polish-only-ablation.md) | bf8840b0ea9cbc6f199518783896fb09c8462104 | `build/final-polish-only/` | `matrix.json`, `final-comparison.json`, `source-isolation.json` | Local present / untracked. |
| [E9 / Final polish, follow-up](https://github.com/yslianTim/rhbm_gem/blob/c174458294c1a328058f96bf37a2a21ae16c06a4/docs/developer/second-stage-final-polish-production-evidence.md) | fb926a5be58f0ade7ffd855bba44cc84789e79f0 | `build/final-polish-production/` | `summary.json`, `baseline.json`, `case1-on-j4-trace-off/fold-report.json`, `artifact-manifest.json`, `object-hashes.json` | Local present / untracked. |
| [E10 / Global-best](https://github.com/yslianTim/rhbm_gem/blob/c174458294c1a328058f96bf37a2a21ae16c06a4/docs/developer/second-stage-global-best-only-ablation.md) | 5174961e2373f97738ec8201666b643531c1bbf1 | `build/global-best-only/` | `best-analysis.json`, `exposure.json`, `comparison.json`, `source-isolation.json`, `capture-source-hashes.json` | Local present / untracked. |
| [E11 / Rescue](https://github.com/yslianTim/rhbm_gem/blob/c174458294c1a328058f96bf37a2a21ae16c06a4/docs/developer/second-stage-rescue-only-ablation.md) | ddc165052a2bce5828829cfd04d43385c3671a50 | `build/rescue-only/` | `report.md`, `event-comparison.json`, `quality.py`, `applied-source-verification.json` | Local present / untracked. |
| [E12 / Background trigger](https://github.com/yslianTim/rhbm_gem/blob/c174458294c1a328058f96bf37a2a21ae16c06a4/docs/developer/second-stage-background-trigger-ablation.md) | 277ab397be3b788da04f9d808464424a5fff57e2 | `build/background-trigger/` | `analysis.json`, `peeling-comparison.json`, `artifact-manifest.json`, `source-verification.json` | Local present / untracked. |
| [E13 / fold-168 investigation](https://github.com/yslianTim/rhbm_gem/blob/c174458294c1a328058f96bf37a2a21ae16c06a4/docs/developer/fold-168-iteration-regression.md) | 6587542638d570685a11ef281aebd015d5eb8c06 | `build/fold-168-regression-investigation/` | `run-index.json`, `65875426-official/report.json`, `65875426-official/actual.json`, `instrumentation-manifest.json`, `binary-hashes.json` | Local present / untracked. |

E13 retains the original input SHA-256 identities and reproduction commands in
its report; `run-index.json` describes 21 full dataset runs, with failed builds
kept separately. Neither that count nor directory existence implies a new
validation of all raw data. E8 and E9 are separate rounds, not a pooled sample:
five converged executions in E8 are not five unique cases, and E9's interrupted
or unrun workloads have no application conclusion. E10's 12 cooperative
comparisons cover endpoints/corrections, not 12 independent attempts.

The two dated verification sections removed from the current specification
remain available as [Workspace verification (2026-09-04)](https://github.com/yslianTim/rhbm_gem/blob/c174458294c1a328058f96bf37a2a21ae16c06a4/docs/developer/second-stage-local-fitting.md#workspace-verification-2026-09-04) (baseline
`b2c0bc5c9a64afea94d131939ee2b0eb437a9613` plus the then-uncommitted cleanup)
and [Readability refactor verification (2026-09-06)](https://github.com/yslianTim/rhbm_gem/blob/c174458294c1a328058f96bf37a2a21ae16c06a4/docs/developer/second-stage-local-fitting.md#readability-refactor-verification-2026-09-06) (paired baseline
`a285a63a29ff134f346b3697743d5cd2ce0032a5`). Neither section specifies a local
artifact directory, so raw-evidence availability is not established by this
inventory. Their old test results are historical records, not current validation.

### Earlier convergence provenance

These three existing historical records remain in place and are not production
specifications:

1. [Convergence safeguard audit](audit-history/second-stage-convergence/convergence-safeguard-audit.md)
2. [Stationarity and active-coordinate population audit](audit-history/second-stage-convergence/stationarity-active-coordinate-audit.md)
3. [Counterfactual convergence continuation audit](audit-history/second-stage-convergence/counterfactual-convergence-continuation-audit.md)
