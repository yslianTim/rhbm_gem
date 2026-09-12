# Frozen recovery background-trigger ablation

## Scope and baseline

Baseline: clean `277ab397be3b788da04f9d808464424a5fff57e2`.
This experiment compares recovery eligibility only. Production retains the
any-background-change policy. No test cases, fixtures, assertions, public options,
or production diagnostic fields are added or changed.

All experiment sources and instrumentation are isolated under
`build/background-trigger/`. The production shrink-level representation and
shared component assembly/salvage remain unchanged.

## Policies

| Label | Eligibility after a newer recovery environment revision |
| --- | --- |
| A | Any exact background response change or applied changed partition |
| B001 | Applied partition or target-local change greater than 0.001 |
| B01 (primary) | Applied partition or target-local change greater than 0.01 |
| B1 | Applied partition or target-local change greater than 0.1 |
| C | Applied changed partition only |

The B thresholds are fractions of the fixed objective residual scale, not
percentages of the background value. Each target retains its existing ShapeAtom,
OffsetAtom or HardFailureCluster identity.

For a target, the sample set contains its own samples plus samples on other
selected atoms whose physical neighbor contributors include a target atom.
Owner/sample traversal produces a sorted, unique SampleRef list. No pruned
topology edge, healthy whole-cluster expansion or objective-mask approximation
is used to determine this physical sample set.

The reference is the background at freezing or the last actual retry. The B
experiment stores each reference sample's background and owner's positive finite
`ObjectiveScale::fit`, floored at `1e-12`, and computes

```text
D(target) = max_sample abs(current_background - reference_background)
                       / reference_scale
retry when D(target) > threshold
```

Suppression does not replace the reference. Several subthreshold changes can
therefore accumulate across attempts. Empty samples, missing references or
unavailable scales fall back to A eligibility for that target and carry a reason
in the build-local capture. Existing invalid-background errors are retained.

## Wiring and ownership

All variants preserve objective and recovery-environment revision advancement.
A newer revision is necessary but B/C additionally filter eligibility before
`QuarantineState::BeginIteration` creates activity masks. The existing single-argument
entry remains intact for direct callers/tests; only the isolated iteration path
uses the added eligibility input. Last-recovery revision updates retain the
existing freeze/retry semantics.

The experiment reference map belongs to iteration state, not a diagnostic observer.
After commit it records new Frozen targets, refreshes targets that actually retried
and remain Frozen, and removes released targets. References for suppressed targets
are retained. A/C non-capture builds do not populate this B-only reference map;
capture builds also maintain references for A/C to observe local-change exposure.

Applied partition change authorizes all existing Frozen targets once, including
simultaneous background change. Merely queued topology does not authorize retry.
Overlapping masks, minimum radius, 10x ridge, failure streaks, release rules,
objective gates, convergence and stopping are not changed. In particular,
audit-patience and all-rejected are not extended to manufacture later retry
opportunities.

## Measurement and comparability

The adapted recovery-revision capture runs the complete existing C++ executable.
It records target eligibility and reference attempts, local-change values,
retry/failed-retry/release state, active masks, selected/rejected keys, Gaussian
trajectories, selection objectives, production convergence certificates,
finalization states, persisted models/peeling, response MSE and production work
counters. Scheduled retries still blocked by another Frozen target are counted
separately. Operator evidence is the production certificate where evaluated;
absence of a certificate on an all-rejected path is not a passing certificate.

Each terminal state is additionally evaluated using the same initial objective
domain, scales and weights, rebuilding its selected-median background by the same
rule. Initial reference-domain records are compared across variants. This common
objective is separate from the native final-context audit and stored-best value;
measurement does not feed decisions or production counters.

Non-finite MSE in deliberately extreme existing fixtures is retained in the raw
records and classified as unavailable for finite-quality comparisons, not silently
turned into zero or counted as a new policy regression.

## Reproduction and build isolation

`prepare.py` archives the baseline and creates A/B001/B01/B1/C source copies.
`build.py` uses the baseline Ninja compiler/linker commands for trace OFF/ON with
trust-model experiment OFF. It recompiles the changed Quarantine and IterationProcess
translation units against each isolated source and relinks its own library, CLI
and test executable, reusing the baseline's unchanged object files and dependencies.
The eligibility overload does not alter QuarantineState layout; all changed private
iteration state is confined to IterationProcess. This avoids rebuilding unrelated
translation units without mixing policy implementations.

Each variant receives the unchanged generated CTest definitions with executable
paths relocated to its isolated build. Variant rpaths select its own library;
`library-resolution.log` also verifies actual library loading. `build-commands.json`
records exact commands. Existing tests and their assertions are reused unchanged.

`capture.py` relinks separate instrumented libraries and executable copies;
`run-rest.py` gates the numerical policy runs on completed baseline/A equivalence.
`analyze.py` aligns runs within each existing case, checks reference domains and
reports first eligibility, model-state and selection-objective divergences.
`benchmark.py` waits for concurrent validation to finish, then runs the two existing
recovery-active command cases three consecutive times per trace-OFF variant,
without capture. Request job count is 1 and thread-limit environment variables
are fixed to 1. Wall times include each complete fixture's setup and command work;
they are not isolated solver timings or parallel CTest performance measurements.

## Results

All ten policy/configuration combinations passed the complete existing CTest
suite (16/16 each): A, B001, B01, B1 and C, each with audit-trace OFF and ON,
trust-model experiment OFF. Unmodified baseline OFF/ON also passed 16/16.
The full capture executable passed 861 tests in 76 suites for baseline and each
of the five policies. No assertions or tolerances were changed.

A matched the unmodified baseline exactly in 15,740 event records, 3,142 trajectory
records, 9,426 recovery records and 3,588 non-timing finalization records. This
checks selection, recovery, state, objectives, stopping and persistence; it does
not compare timing or unordered worker logging.

| Policy | Total attempts (75 runs) | Retry / failed retry | Directly suppressed | Frozen target-attempts |
| --- | ---: | ---: | ---: | ---: |
| A | 3,142 | 190 / 190 | 0 | 265 |
| B001 (0.1%) | 3,142 | 190 / 190 | 0 | 265 |
| B01 (1%) | 2,998 | 20 / 20 | 26 | 121 |
| B1 (10%) | 2,958 | 0 / 0 | 6 | 81 |
| C | 2,958 | 0 / 0 | 6 | 81 |

All policies had eight Frozen entries and zero successful releases. No scheduled
retry was blocked by an overlapping Frozen mask. The retry reductions are not
all direct suppression: A→B01 removes 26 retries within shared attempts and 144
retries after the new stop; A→B1/C removes six within shared attempts and 184
after the new stop.

Only two production runs changed: `PotentialAnalysisCliAcceptsOnlyBackbone`
(run 1) and `PotentialAnalysisOnlyBackboneFiltersOnlyBackboneAtoms` (run 2).
Each ran 100 accepted attempts and stopped at MaximumIterations under A/B001,
28 under B01 and eight under B1/C; the latter policies stopped at AuditPatience.
The separate all-atoms run in the second case remained eight attempts. Across
all 75 runs, Converged remained five and AllRejected remained 21; AuditPatience
increased from 25 to 27 and MaximumIterations fell from 24 to 22. The changed
runs therefore provide no newly successful convergence certificate.

First eligibility divergence from A is attempt 6 for B01/B1/C, at the shape
atom target. B01 versus B1/C first differs at attempt 7. Common-prefix Gaussian
states and selection objectives do not diverge; the trajectories instead end
earlier. The cumulative reference behaves as intended: B01 suppresses attempt 6
with reference attempt 5 and D=0.00722513281064138, then retries at attempt 7
against that same reference with D=0.01439638136525152. At attempt 28 its
reference is attempt 25 and D=0.0092963582629549263, still below the threshold.
The existing stopping rule ends the run before another material change occurs.

### Terminal quality and persistence

Initial objective domains/scales/weights matched for all 75 paired runs. All 75
persisted Gaussian parameter sets matched exactly. All 72 finite common-objective
records and all 68 finite response-MSE records matched exactly; seven existing
non-finite MSE records remained non-finite and are not quality passes. Missing
or non-comparable measurements are retained in the raw capture.

**Peeling did not match in the two affected runs.** Each has 200 samples, of which
33 changed relative to A. Results below apply separately to each affected run.

| Policy | Common objective | Selected-response MSE | Native final audit | Peeling max absolute difference | Peeling RMS difference |
| --- | ---: | ---: | ---: | ---: | ---: |
| A / B001 | 0.86227276075258086 | 3.0706893604672603 | 0.85505688724164552 | 0 | 0 |
| B01 | same | same | 0.85897859908238083 | 0.13610176090954518 | 0.05202132723683727 |
| B1 / C | same | same | 0.86141380940947399 | 0.2217934940625934 | 0.08477474396014613 |

Thus only 73/75 combined model-and-peeling records are identical for B01/B1/C.
B001 matches A throughout; B1 and C match the observed numerical results, which
does not make their policies equivalent under other trajectories.

The source explains the persistence difference: `ApplyFitState` creates its model
snapshot with `context.frozen_background`; `CalculateSecondStageAdjustedResponse`
subtracts that background when preparing persisted adjusted samples. Earlier
termination retains a different last-refreshed background even when the final
Gaussian state is identical. The existing selected-response MSE helper sums only
selected Gaussian responses, so it does not measure this background/peeling
change. The common objective deliberately rebuilds background from each terminal
Gaussian state and likewise cannot establish peeling equality. Native final audits
use differing terminal contexts and must not be ranked as a common objective.
This existing persistence behavior was not changed by the experiment.

The final operator certificate in each affected run remains solver-qualified and
operator-complete, but not production-converged. Its nominal p99 residual vectors
(in the existing three-channel order) are A/B001 `(0.11157413, 0.10640018,
0.03236284)`, B01 `(0.11083000, 0.10622034, 0.04972636)` and B1/C
`(0.11036727, 0.10610849, 0.06066288)`. The third residual is larger at the earlier
stops. Full per-attempt certificates are retained in the per-case comparison;
operator completeness alone is not convergence or equal residual quality.

### Work and timing

Actual linear-solve calls across the two existing command cases were 3,353 for
A/B001, 1,049 for B01 and 409 for B1/C. These include the unchanged all-atoms run.
Supplemental solve-call instrumentation reproduced the primary capture's terminal
model, common-objective and persistence records. Detailed materialization,
objective-sample, correction and operator records remain in the capture and
per-case comparison JSON.

Trace OFF, capture disabled, fixed thread count 1; three consecutive complete
executions of the same two command cases per variant:

| Policy | Seconds, repetitions 1 / 2 / 3 | Median seconds |
| --- | --- | ---: |
| A | 2.045572 / 1.978807 / 1.977274 | 1.978807 |
| B001 | 2.082842 / 1.994507 / 1.985381 | 1.994507 |
| B01 | 1.074065 / 1.022089 / 1.022000 | 1.022089 |
| B1 | 0.858717 / 0.811983 / 0.821046 | 0.821046 |
| C | 0.860857 / 0.808004 / 0.804677 | 0.808004 |

These measurements describe the exposed fixtures, including setup. Much of the
work reduction comes from earlier termination. Three repetitions do not establish
small timing differences or general production speedups. CTest wall times are
not used as performance evidence.

### Decision

Keep production A. B01 reduces failed work, but sensitivity results differ
substantially, termination changes and downstream peeling changes. The available
metrics do not establish that this output change preserves quality. With no
successful production release or applied-partition exposure, B01 does not meet
the agreed conditions for a validated follow-up candidate. C eliminating the 190
failed retries is not evidence that recovery is safe to delete.

### Reproduction artifacts

All raw files and tools are under `build/background-trigger/` (build-local,
not tracked production source). `source-verification.json` records source hashes
and the precise variant delta; `artifact-manifest.json` hashes scripts, captures,
commands and result records. `A-equivalence.json`, `analysis.json`,
`peeling-comparison.json`, `solver-work.json` and `benchmark.json` contain the
machine-readable comparisons. Each capture retains complete executable output,
event/trajectory/recovery/target/finalization TSV files and captured case stdout.

From the repository root, the retained scripts reproduce the isolated experiment:

```sh
python3 build/background-trigger/prepare.py
python3 build/background-trigger/build.py A off
python3 build/background-trigger/build.py A on
python3 build/background-trigger/build-rest.py
python3 build/background-trigger/capture.py baseline
python3 build/background-trigger/capture.py A
python3 build/background-trigger/run-rest.py
python3 build/background-trigger/analyze.py
python3 build/background-trigger/benchmark.py
python3 build/background-trigger/solver-work.py
python3 build/background-trigger/verify-artifacts.py
```

These commands depend on the baseline trace-OFF/ON Ninja builds identified in
`build.py`; exact compiler/linker commands are retained per isolated build.
The ignored artifact directory must be preserved with this report when moving
between machines; the tracked report alone is not the raw experimental archive.

## Coverage limits

The baseline exposes 8 entries into Frozen and 190 retries, all failed, across
75 full-loop runs. Its retry-active cases are the two existing backbone command
cases. The direct release tests are retained, but do not supply production-loop
successful release coverage for B/C. No new scenarios are introduced to fill gaps.

No 600-case corpus, fold-168 numerical dataset or independent real-data benchmark
is run. The existing fold-168 CTest exercises its runner, not that dataset. The
trust-model-ON matrix is outside this experiment. Reduced failed retries alone
cannot establish that future recoveries are safe to suppress.

No queued/applied partition transition, production-loop successful release,
overlap-blocked scheduled retry or eligibility fallback was exposed by these
runs. Exact-threshold behavior and wider target/physical-dependency combinations
are not dynamically certified by this dataset. Existing direct checks passing
cannot substitute for these missing production trajectories: evidence is
insufficient for changing the recovery policy.

Final documentation validation: the existing core-contract CTest passed (1/1),
`git diff --check` passed, and source verification confirmed that production
`src`, `include` and `tests` remain identical to the baseline.
