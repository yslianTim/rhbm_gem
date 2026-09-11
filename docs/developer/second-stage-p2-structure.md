# P2 structural simplification: retained steps and acceptance limit

Baseline: `f50a742cd5d8f3c64fa70cbc11436c7677ff0370`, with a clean worktree.
This is a numerical ablation, not a numerical-equivalence refactor.

## Retained

1. Active/Frozen quarantine, with one retry per changed objective-domain revision,
   including background response changes. Retry uses the existing minimum radius
   and ridge. Complete recovery evidence is required. Next-iteration lifecycle
   publication no longer changes the already audited model; fallback re-audit is
   removed. Legacy probation counter positions now report domain retries.
2. Cluster best is history/tie-break only. Previous gates, strict local polish,
   rescue tolerance/global strict improvement, global-best gate, patience and
   final-state selection remain. Only passing candidates evaluate historical
   references; unavailable history cannot reject a candidate. Local and rescue
   share history updates. Removed best-rejection fields remain false.

The [execution map](second-stage-p0-structure.md#ownership-and-execution) reflects
these retained changes. Production Keep/Shrink, certificates and converged-only
final polish at radius 1.0 remain unchanged.

## Withdrawn: fixed-order atomic dependency-component acceptance

The attempted implementation formed complete components from all partition
keys/shared samples, included singletons, and checked prefixes in key order
without greedy subset removal. It failed the existing
`SystemBuildFailureDoesNotBlockRemoteCluster` and
`PersistentEmptySystemDoesNotBlockRemoteCluster` regressions.

Both fixtures have two dependency components and locally accepted keys, but an
unavailable complete global baseline (NaN response or empty samples in the
failing cluster). Mandatory complete-global acceptance rejects the healthy
remote update too. Its MSE remains at 0.024571741450482135 and
0.029314173978686155 respectively, whereas both tests require improvement.
Serial/parallel and intensity-scale assertions passed during the experiment.

By the requested acceptance rule, step 3 was withdrawn. The existing global
greedy salvage remains. No missing-domain exemption, tolerance change or weakened
remote-improvement assertion was introduced. A future component transaction
policy must first specify acceptance when unrelated objective evidence is absent.

## Verification and coverage

No test cases/framework were added, and numerical tolerances were not loosened.
The existing lifecycle case was adapted to domain retry, and assertions directly
requiring cluster-best rejection were adapted to the new policy. Initial step-2
failures concerned obsolete rejection logs, not model equality or scale safety.

Baseline and retained-step build/test results, four diagnostic-flag configurations,
independent patches and failed-step evidence are recorded locally under
`build/p2-structure/`. See its `report.md` for final validation status.

Existing cases do not exhaustively cover overlapping Frozen/retry targets,
background revision sequences or every incomplete-history branch. These paths
also received source review; passing tests is not complete coverage. No corpus,
6Z6U or fold-168 run, real-data quality claim or performance claim is made.
