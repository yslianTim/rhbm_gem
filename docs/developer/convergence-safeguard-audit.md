# Second-stage convergence safeguard audit

> Current baseline (2026-08-27): production requires current active-block
> stationarity plus accepted/raw active-DOF p99 below `1e-4`. Maximum change
> remains diagnostic and drives topology-drift checks, but is not a convergence
> predicate. The first-round measurements below remain historical evidence.

## Current resolution

- Retain accepted-state p99 because post-raw correction can move a committed
  state after a small raw endpoint.
- Retain raw fixed-point p99 because clipping, backtracking, or rejection can
  hide a material raw residual from the accepted step.
- Evaluate both summaries over active optimization DOFs. Fixed and quarantined
  coordinates do not dilute the population, and each shared offset contributes
  one group-level sample.
- Remove accepted/raw maximum from the production stopping conjunction. The
  historical `1e-3` gate is now the isolated `legacy-maximum` comparator.
- Keep strict stationarity diagnostic-only until continuation demonstrates a
  material objective or truth-error benefit.

## Historical purpose and scope

This is the first-round Failure Mode × Safeguard audit for
`RunSecondStageLocalFitting`. It separates the existing convergence gate into
five predicates without changing the stopping policy:

- active-block stationarity;
- accepted-state p99 change;
- accepted-state maximum change;
- raw fixed-point p99 change;
- raw fixed-point maximum change.

Accepted change is `|z(S[k+1]) - z(S[k])|`; raw change is
`|z(F(S[k])) - z(S[k])|`. The transformed coordinates `z` are log peak
height, log width, and offset/peak. Every coordinate must have p99 below
`1e-4` and maximum below `1e-3`. Objective-domain changes, quarantine
transitions, suspicious updates, and rejected clusters remain orthogonal
blockers and are not removal candidates in this audit.

At the time, the audit added debug-only observation without changing the
production convergence expression.

## Static implication results

| Relationship | Result | Reason |
| --- | --- | --- |
| accepted change ⇒ raw change | False | Trust clipping, objective backtracking, or rejection can make the committed step small while the raw fixed-point residual remains large. |
| raw change ⇒ accepted change | False | Joint polish, boundary reconciliation, and rescue can move the accepted state away from an otherwise small raw endpoint. |
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
| 990 fixed zeros plus ten active values at `5e-4` | — | all-selected p99=1; shadow p99=0 | same | Current population can hide active-block drift. |
| Active cluster with IRLS maximum-iteration status | current=1, full=0 | — | — | Current active-block flag can hide soft joint non-convergence. |

The debug trace is also checked for deterministic equality between serial and
parallel selection and for behavioral neutrality between Info and Debug runs.

## Current Debug trace contract

At Debug verbosity (`-v 4`) every accepted attempt emits schema `4` on a line
beginning with `Convergence safeguard audit:`. The record contains:

- attempt, accepted iteration, selected/quarantined population;
- production active-DOF, legacy all-selected, and active-member population
  sizes for all three coordinates;
- production, `legacy-population`, `legacy-maximum`, and `strict-dof`
  predicate vectors;
- accepted/raw p99 and diagnostic maximum values for every population;
- isolated stop-candidate and exposure flags for the three comparators;
- whether accepted and raw states are equal;
- trust limiting, backtracking, polish, boundary, and rescue path counts;
- current/full stationarity, exact local-refit status counts, and soft/hard
  joint status counts;
- fixed-block, damping, suspicious, rejection, quarantine-transition, and
  objective-domain-change counts.

These stable key/value fields support pairwise 2×2 tables, implication
counterexample searches, blocker/unique-blocker counts, and stratification by
population size, active ratio, quarantine ratio, and proposal path. A lack of
observed counterexamples is empirical dominance only; it is not treated as a
mathematical implication.

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
preferred active-DOF population without a maximum gate; the accepted/raw and
stationarity independence results remain valid.

The external fold-168 regression remains optional for this round because its
model and map inputs are not stored in the repository. When hash-matching inputs
are available, its Debug trace can be aggregated with the same contract without
changing the fitting configuration.

## Current refresh verification (2026-08-27)

- Audit-enabled CTest passes 21/21, including analyzer fixtures, runner
  smoke/determinism, and the external counterfactual fold-168 audit.
- Audit-disabled CTest passes 19/19, including the external fold-168 regression.
- Both fold-168 runs stop after seven accepted iterations with
  `audit-patience`; the audit report is `no_convergence_trigger` with zero
  legacy-population, maximum-gate, and strict-stationarity exposures.
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
