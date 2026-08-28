---
Status: Historical audit record.
Current policy: second-stage-outer-iteration-algorithm-audit.md
---

# Second-stage convergence safeguard audit

> Current baseline (2026-08-28): production requires solver qualification,
> accepted active-DOF p99 below `1e-4`, and a complete nominal-DOF fixed-point
> residual p99 below `1e-4`, with invariants and orthogonal blockers clear.
> These terms are evaluated by one internal `ConvergenceCertificate`. Maximum
> remains a tail diagnostic. The first-round measurements below remain
> historical evidence.

## Current resolution

- Retain accepted-state p99 because post-raw correction can move a committed
  state after a small raw endpoint.
- Use strict operator residual p99 as the production residual predicate. It is
  computed before guard feasibility, terminal freeze, trust limiting, objective
  gates, and polish.
  Unavailable nominal DOFs produce a restricted classification, not a zero.
- Evaluate accepted movement over active optimization DOFs. Fixed and
  quarantined coordinates do not dilute that population, and each active shared
  offset contributes one group-level sample. Evaluate the strict operator
  residual over the complete nominal-DOF population, including fixed and
  quarantined shapes and shared offsets; unavailable evidence makes the
  operator incomplete rather than substituting a zero residual.
- Remove accepted/operator maximum from the production stopping conjunction.
  The historical `1e-3` gate is retained only by the
  `historical-cluster-active-proposal-maximum` and `production-maximum`
  diagnostic comparators.
- Require solver qualification in production.
- Treat objective-accepted final dependency polish on a `converged` path as
  provisional. Persist it only after the polished state independently passes
  the same solver-qualified, complete nominal-DOF strict-operator p99
  certificate; otherwise retain the already converged base state.

## Historical purpose and scope

This is the first-round Failure Mode × Safeguard audit for
`RunSecondStageLocalFitting`. It separates the existing convergence gate into
five predicates without changing the stopping policy:

- active-block stationarity;
- accepted-state p99 change;
- accepted-state maximum change;
- guarded-proposal p99 change;
- guarded-proposal maximum change.

Historically, accepted change was `|z(S[k+1]) - z(S[k])|` and the guarded
proposal was described as `|z(F(S[k])) - z(S[k])|`. That second description
was imprecise because guard damping and fixed/quarantined rollback had already
been applied. Historical schema 6 called it `guarded-proposal` and reserved
`fixed-point-residual` for the undamped offset-to-shape operator; legacy
schema 7 reports the strict operator residual directly and is normalized by a
read-only adapter. Current schema 8 serializes the canonical certificate. The
transformed coordinates `z` are log peak
height, log width, and offset/peak. Every coordinate must have p99 below
`1e-4` and maximum below `1e-3`. Objective-domain changes, quarantine
transitions, suspicious updates, and rejected clusters remain orthogonal
blockers and are not removal candidates in this audit.

At the time, the audit added debug-only observation without changing the
production convergence expression.

## Static implication results

| Relationship | Result | Reason |
| --- | --- | --- |
| accepted change ⇒ guarded-proposal change | False | Trust clipping, objective backtracking, or rejection can make the committed step small while the guarded proposal remains large. |
| guarded-proposal change ⇒ accepted change | False | Joint polish, boundary reconciliation, and rescue can move the accepted state away from an otherwise small guarded endpoint. |
| stationarity ⇒ small change | False | Stationarity describes inner-solver/refit qualification, not the magnitude of the outer fixed-point step. |
| small change ⇒ stationarity | False | Fallback, damping, fixed blocks, or an unfinished inner solve can return a numerically small step without qualification. |
| p99 ⇒ maximum | False in general | p99 intentionally ignores a sufficiently sparse extreme tail. |
| maximum ⇒ p99 | False at the current thresholds | All values can be below `1e-3` while at least 1% remain above `1e-4`. |

For the repository's linearly interpolated p99 and nonnegative change values,
there is one conditional implication. For a common population of at most 91
members, `p99 < 1e-4` implies `maximum < 1e-3`. At 92 members a counterexample
already exists: 91 zeros plus one `1.05e-3` value yield p99 `9.45e-5` while
the maximum fails. This makes the maximum predicate conditionally redundant
for small common populations, not globally redundant.

## Failure Mode × Safeguard matrix

`U` means the predicate has unique coverage in a targeted case, `C` means it
is complementary, and `—` means it is not the responsible signal.

| Failure mode | Stationarity | Accepted p99/max | Raw p99/max | Audit result |
| --- | --- | --- | --- | --- |
| Accepted no-op or small guarded step with a large fixed-point residual | C | — | U | Retain raw change. |
| Small raw endpoint followed by material polish/reconciliation | — | U | — | Retain accepted change. |
| Inner solve/refit unfinished but numerical change is small | U | — | — | Retain stationarity, audit its definition. |
| Guard-damped material proposal | U | C | C | Stationarity is the direct qualification signal; changes remain complementary. |
| Coherent mid-tail drift between `1e-4` and `1e-3` | — | U (p99) | U (p99) | Retain p99. |
| Isolated extreme-tail drift above `1e-3` | — | U (maximum) | U (maximum) | Retain maximum for general population sizes. |
| Fixed/quarantined zeros dilute the all-selected population | — | p99 semantic risk | p99 semantic risk | Active-block population needs redesign review. |
| Invalid/non-finite transformed coordinates | — | C | C | Both summaries fail safely; no removal conclusion. |

## Reproduced synthetic evidence

The targeted tests exercise the predicate truth table directly:

| Case | Stationarity | Accepted p99/max | Raw p99/max | Finding |
| --- | ---: | ---: | ---: | --- |
| Accepted small, raw large | 1 | 1/1 | 0/0 | Accepted does not imply raw. |
| Accepted large, raw small | 1 | 0/0 | 1/1 | Raw does not imply accepted. |
| Both changes small, stationarity false | 0 | 1/1 | 1/1 | Small change does not imply stationarity. |
| 999 zeros plus one `2e-3` tail | — | 1/0 | 1/0 | Maximum has unique sparse-tail coverage. |
| 1000 values at `5e-4` | — | 0/1 | 0/1 | p99 has unique coherent-drift coverage. |
| 990 fixed zeros plus ten active values at `5e-4` | — | all-selected p99=1; shadow p99=0 | same | The historical all-selected population can hide active-block drift. |
| Active cluster with IRLS maximum-iteration status | historical=1, full=0 | — | — | The historical cluster rollup can hide soft joint non-convergence. |

The debug trace is also checked for deterministic equality between serial and
parallel selection and for behavioral neutrality between Info and Debug runs.

## Current Debug trace contract

At Debug verbosity (`-v 4`) every accepted attempt emits schema `8` on a line
beginning with `Convergence safeguard audit:`. The record contains:

- `certificate-definition=1` and `comparator-set=1`;
- attempt, accepted iteration, selected/quarantined population;
- accepted active-DOF, historical all-selected, and operator nominal-DOF
  population sizes for all three coordinates;
- the seven-term serialized certificate and the five ordered decisions:
  `production`, `historical-all-selected`,
  `historical-cluster-active-proposal-maximum`,
  `historical-active-proposal`, and `production-maximum`;
- accepted active and operator nominal residual p99/maximum values;
- strict-operator unavailable and sparse-tail counts plus residual-state
  classification;
- isolated stop-candidate flags for the five decisions and exposure flags for
  the four non-production comparators;
- whether accepted and operator states are equal;
- unified-search invalid, trust-skipped, guard-rejected, objective-rejected,
  accepted-factor, and terminal counts, plus polish, boundary, and rescue paths;
- production and solver qualification, exact local-refit status counts, and
  joint-solver status counts;
- fixed-block, suspicious, rejection, quarantine-transition, and
  objective-domain-change counts.

These stable key/value fields support pairwise 2×2 tables, implication
counterexample searches, blocker/unique-blocker counts, and stratification by
population size, active ratio, quarantine ratio, and proposal path. A lack of
observed counterexamples is empirical dominance only; it is not treated as a
mathematical implication.

The schema-8 consolidation was replayed against the same frozen 600-case
reference on 2026-08-28. All 600 normalized production semantic digests match,
with zero safety regression and zero comparator exposure; the stop distribution
remains 42/372/163/23 for converged/audit-patience/all-rejected/maximum-
iterations. The tracked compact result is
[`convergence_certificate_baseline.json`](../../../../tests/benchmarks/convergence_certificate_baseline.json).

## Historical first-round decisions

| Safeguard | Classification | Next action |
| --- | --- | --- |
| Accepted-state change | Retain | It has a constructed unique case when post-raw correction moves the accepted state. |
| Raw fixed-point change | Retain | It uniquely detects suppressed or rejected raw residuals. |
| Stationarity | Semantic redesign needed | Decide whether soft joint non-convergence must clear active-block eligibility while preserving quarantine semantics. |
| p99 | Retain | It uniquely rejects coherent drift below the maximum threshold. Recompute on truly active coordinates in a follow-up. |
| Maximum | Retain globally; conditional second-round ablation candidate | It is redundant under the p99 predicate only for a common population `N <= 91`, but uniquely protects larger populations. |

No safeguard is deleted in this round. The next review should first resolve
stationarity semantics, then decide whether convergence populations should be
coordinate-specific active sets. Only after those definitions are stable
should a trajectory-changing ablation compare removal of the maximum gate for
small populations or any empirically dominated predicate on fold-168 data.

That follow-up is implemented as the shadow-only
[stationarity semantics and active-coordinate population audit](stationarity-active-coordinate-audit.md).
It preserves the production stopping expression while comparing strict
block-level stationarity, active-member changes, and one-sample-per-shared-DOF
offset changes.

The subsequent audits completed this follow-up. Production now uses the
preferred active-DOF population for accepted movement and the complete
nominal-DOF population for strict operator residuals, without a maximum gate.
The accepted/guarded and qualification independence results remain valid. The
active-member comparison was retired after the shared-DOF population became
authoritative; its results above remain historical evidence only.

The external fold-168 regression remains optional for this round because its
model and map inputs are not stored in the repository. When hash-matching inputs
are available, its Debug trace can be aggregated with the same contract without
changing the fitting configuration.

## Historical refresh verification (2026-08-27)

- Audit-enabled CTest passes 21/21, including analyzer fixtures, runner
  smoke/determinism, and the external counterfactual fold-168 audit.
- Audit-disabled CTest passes 19/19, including the external fold-168 regression.
- Both fold-168 runs stop after seven accepted iterations with
  `audit-patience`; the audit report is `no_convergence_trigger` with zero
  legacy-population, maximum-gate, and solver-qualification exposures.
- Audit-enabled and audit-disabled fold-168 `actual.json` files are
  byte-identical, and repository lint passes.
- The 600-case exposure corpus was not rerun.

## Historical verification status (2026-08-26)

- All seven new convergence-audit tests pass, including serial/parallel trace
  equality and Info/Debug behavioral equivalence.
- The 110 second-stage defense tests outside the repository's existing
  log-format mismatch pass after rebuilding, and repository lint passes.
- Full `ctest` passes 14 of 15 grouped tests. The remaining core-estimator group
  has seven pre-existing assertions that expect the older compact log spelling
  (for example `accepted_iterations=`), while the tracked implementation emits
  the current spaced, multi-line spelling (`accepted_iterations =`). This audit
  does not change that unrelated output contract.
- The fold-168 runner self-tests pass. The external fold-168 simulation was not
  run because the hash-matching model and map inputs are unavailable locally.
