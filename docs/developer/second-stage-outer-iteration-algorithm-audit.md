# Second-stage outer-iteration algorithm audit

## Status and authority

This document is the current decision and evidence authority for the outer
iteration implemented by `detail::RunSecondStageIterations`. The normative execution
description remains in [Second-stage local fitting](second-stage-local-fitting.md).
The four earlier audits are immutable historical records and are linked under
[Historical provenance](#historical-provenance).

The reviewed production baseline is
`a4354e698e77398154009d231907ebf3ed4b1d52` (per-atom offsets). The audit consolidation and diagnostic cleanup do
not change `FitOptions`, command-line options, model persistence, convergence
thresholds, candidate selection, stop precedence, or the production
trajectory.

## Scope and canonical states

The review covers the second-stage outer loop from a validated accepted state
through proposal construction, candidate selection, post-processing,
convergence, final dependency polish, and persistence. The statistical
derivation of the local MDPDE estimator and the later group-fitting stage are
outside this audit.

Three states must remain distinct:

- `S(k)` is the previous validated accepted state.
- `F(S(k))` is the complete, undamped joint-offset-to-local-shape operator
  endpoint, including availability and solver evidence.
- `S(k+1)` is the candidate that is actually accepted after the geometric
  factor search, objective gates, joint polish, boundary reconciliation, and
  rescue.

Accepted movement is `T(S(k+1)) - T(S(k))`. The strict fixed-point residual is
`T(F(S(k))) - T(S(k))`. `T` uses three transformed coordinates: log peak,
log width, and per-atom physical offset normalized by peak.

Accepted movement samples active optimization DOFs. Every shape-active atom
contributes one log-peak and one log-width sample; every offset-active atom
contributes one absolute offset-to-peak-ratio change sample. Fixed and quarantined coordinates do not
dilute this population.

The operator residual instead samples the complete nominal-DOF population,
including fixed and quarantined shapes and per-atom offsets. Missing or
non-finite endpoint evidence makes the operator incomplete; it must never be
replaced by the previous state to manufacture a zero residual.

## Per-atom availability and conditioning

Each selected atom has independent shape and offset availability masks. Offset
availability is set only for a non-hard-failure joint-offset result whose model
passes validity checks. A soft failure may therefore supply an available endpoint
without being solver-qualified. Shape and offset availability are independent;
a missing shape sets both shape residual coordinates to infinity, and a missing
offset sets its own residual to infinity. `operator_complete` is the conjunction
of the masks for every nominal selected atom; transformed finiteness and solver
qualification are separate requirements. Inactive atoms are not removed from this
nominal population, and each coordinate has its own percentile population.

Joint-offset and joint-polish conditioning normalize each design-matrix column,
form the normalized Gram matrix, and use LDLT `min(D)/max(D)` as a conditioning
proxy before ridge. A ratio at or below `1e-8` triggers a ridge multiplier floor
of `10`. Empty/invalid columns, failed factorization, or nonpositive/nonfinite
pivots return the existing zero sentinel and require the guard. This ratio is
neither a singular-value condition number nor an effective-rank certificate.
Other existing ridge safeguards may also increase multipliers; diagnostics report
the actual minimum/maximum multipliers after all guards, not only the floor.

Schema-1 diagnostic records separately report conditioning, solve status, and
per-atom availability. Phases distinguish outer operator, candidate polish,
boundary reconciliation, final dependency polish, and final recertification.
Offset status codes follow `JointOffsetSolveStatus` (0 converged; 1 system build
failed; 2 empty; 3 initial solve failed; 4 IRLS solve failed; 5 objective
deteriorated; 6 iteration limit). Hard-failure classification is recorded
separately. Joint-polish solve records report `solved`/`failed`; this does not
imply objective acceptance or final persistence approval.

## End-to-end state machine

```text
validated S(k)
  -> complete undamped joint per-atom offset endpoint
  -> complete undamped local-shape endpoint
  -> strict operator evidence F(S(k))
  -> geometric candidate factors: validity -> trust -> guard -> objective
  -> active-column joint polish
  -> boundary reconciliation and cooperative rescue
  -> complete-state global previous/best audit
  -> trust-radius and quarantine/probation transition
  -> assembled validated S(k+1)
  -> production convergence certificate
  -> stop policy selects a base final state
  -> final uncut dependency polish candidate
  -> strict or residual-non-regression persistence safety check
  -> persist Gaussian and peeling state
```

Validity establishes that a candidate can be represented. Trust limits the
step tested in the current iteration and updates the next radius. Guard tests
domain feasibility. Objective gates accept or reject candidates. None of
these responsibilities substitutes for fixed-point evidence.

## Authoritative production certificate

`ConvergenceCertificate::ProductionConverged()` is the only production stop
decision and requires all of the following:

```text
solver qualified
&& accepted active-DOF p99 < 1e-4
&& complete nominal-DOF operator
&& nominal fixed-point residual p99 < 1e-4
&& orthogonal blockers clear
```

The percentile predicate is coordinate-wise: the p99 for each of log peak,
log width, and per-atom offset must pass independently. Solver qualification
requires full, undamped, non-fallback active endpoints. Operator completeness and non-finite residuals fail
closed. Orthogonal blockers cover
objective-domain changes, quarantine transitions, suspicious offset fallback,
and rejected clusters.

Maximum values remain diagnostic measurements and do not define a separate
production policy.

## Failure mode and safeguard coverage

| Failure mode | Accepted p99 | Strict operator p99 | Qualification | Invariants / blockers |
| --- | ---: | ---: | ---: | ---: |
| Trust clipping or objective backtracking makes the committed step small while the full endpoint remains material | Detects the small committed step | Blocks the false fixed point | Provides endpoint quality | Records the limiting state |
| Polish, reconciliation, or rescue moves the committed state after a small operator endpoint | Blocks convergence | Detects the small endpoint | Confirms the endpoint solve | Records post-processing blockers |
| Soft solver failure, damping, or fallback produces small numerical movement | Observes movement only | Observes residual only | Blocks convergence | Preserves failure classification |
| Fixed or quarantined coordinates hide an unavailable nominal endpoint | Excludes inactive DOFs by design | Fails closed on incomplete evidence | Reports restriction | Enforces population completeness |
| One atom has unavailable or non-finite offset evidence | Includes its offset only when active | Unavailable evidence makes the operator incomplete; non-finite residual fails the percentile test | Availability alone does not qualify a solver | Each nominal atom retains its own coordinate |
| Objective domain or quarantine changes during the iteration | May still be small | May still be small | May still pass | Orthogonal blocker prevents a premature stop |

No retained predicate is implied by the others. A zero-exposure corpus result
is empirical evidence, not a mathematical redundancy proof.

## Final dependency polish recertification

Final dependency polish is objective-accepted provisionally. On a `converged`
path, a changed polished state is persisted only when a new certificate built
at that state passes `StrictOperatorPassed()`: solver qualification, complete
nominal operator evidence, and residual p99 must all pass. Failure,
incomplete evidence, or evaluation error retains the already converged base
state.

Non-convergence stop reasons use a residual non-regression policy. A strict
candidate is always safe to apply. Otherwise, both the selected base state and
polished candidate must have solver-qualified, complete,
finite nominal operator evidence, and every candidate coordinate must satisfy
`candidate_p99 <= max(base_p99, 1e-4)`. This comparison is coordinate-wise for
log peak, log width, and per-atom offset. If the base is not comparable, only a
strict candidate can pass. Any evaluation error, unavailable evidence, or
residual regression retains the base state. Maximum residual stays diagnostic.
The fallback does not resume the outer loop, increment accepted iterations, or
change the original stop reason.

## Current diagnostic contract

The current Debug trajectory is schema 10 and serializes the production
certificate plus its active and nominal populations, p99 and maximum values,
operator completeness, and four orthogonal blockers. Earlier trajectory
schemas are not accepted by the current analyzer; frozen schema-9 baselines
remain historical data.

The P0 recertification below compares this exact HEAD with documentation and
diagnostic/tooling changes only. Historical schema-9 measurements do not certify
the current implementation.

Frozen-IRLS predicted-reduction and rho instrumentation is not part of a
normal or routine audit build. It is available only through the developer-only
`RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT` build option and never controls the
production trajectory. Its logs are consumed only by
`analyze_trust_model_experiment.py`; the production corpus analyzer does not
import or aggregate them.

## Historical evidence and retired experiments

The frozen baseline expands the checked manifest to 600 deterministic cases.
At `03cdb6ef`, all 600 cases completed with no safety regression, and the stop
distribution was:

| Stop reason | Cases |
| --- | ---: |
| `converged` | 42 |
| `audit-patience` | 372 |
| `all-rejected-backtracking-exhausted` | 163 |
| `maximum-iterations` | 23 |

The consolidation baseline records a manifest SHA-256 of
`2b0d74c249a4723575696d98f99df9d3c449a30837adf8117aaf744145d0a87a`,
case identity SHA-256 of
`63985f81c5b188cfc9992742cad1f4b9418cc36c88b12398e5775a324289bd98`,
and frozen truth SHA-256 of
`04e3cda3b49857b2d5e4f63e973b2392dfe3095360974db988770e7468edd628`.

Historical all-selected, active-proposal, cluster/maximum, and
production-maximum policies are retired. Accepted-only persistence is retained
only as a targeted negative unit scenario: small accepted movement cannot
declare convergence while the strict operator residual remains material.
Frozen-IRLS/rho remains a separate diagnostic experiment because the observed
coverage and action divergence did not justify a production controller.
Coarse-to-fine factor refinement remains rejected on objective and truth
outcomes. Fold-168 remains an optional quality regression, not convergence
policy evidence.

## P0 recertification acceptance

The paired run must complete 600/600 cases on each side with identical manifest,
case identities, seeds, and frozen truth. Safety, quality, and efficiency are
reported independently; missing evidence cannot pass. Comparison schema 5 replaces
the former combined blocking gate. Case summaries use schema 14, aggregates schema
9, compact baselines schema 4; trajectory 10 and terminal 2 are unchanged.

- **Safety:** zero failed cases, finite positive terminal shape parameters,
  finite offsets/objectives, production certificate evidence for convergence,
  and strict/non-regression final-polish persistence evidence.
- **Quality:** every production semantic and normalized terminal-state digest,
  stop reason, objective, transformed-truth RMSE, and accepted-iteration count
  matches. Median/p90/p99 deltas supplement, never replace, per-case equality.
- **Efficiency:** all elapsed times must be present and positive; median and p90
  must both decrease to pass strict-speedup. Per-case deltas, p99, improvement
  fraction, and family/topology strata are reported. Failure here does not revoke
  safety/quality certification and a single timing pair is not a reliable speedup
  study.

Timing and schema-1 diagnostic records are excluded from production digests.
No production convergence/solver/persistence policy is altered. The original
baseline has no schema-1 diagnostic events; its absence is reported, not filled
with zero conditioning failures. Existing baseline final-polish logs serialize
residuals at limited precision: status plus rounded residual comparisons are
checked, not represented as an independent full-precision re-solve.

### Frozen truth provenance and coverage

The checked manifest and case-identity hashes still match the historical compact
baseline. Exact HEAD `a4354e6` emits frozen truth hash
`4d67c57ed11ed3129b44ee275447b3c1b6693c03cfb7c5be6d2265608664933d`,
which differs from the historical `04e3cda3...` hash. Historical per-case truth
files are not present in this workspace; the old compact baseline is retained
unchanged. This is an explicit historical-truth continuity gap, not a newly
established equivalence with the old controller. Natural truth is itself produced
by the HEAD's reference-estimation path, which has changed since that evidence.
P0 freezes the exact HEAD outputs for both sides and verifies emitted candidate
truth against them before computing metrics. Legacy topology labels such as
`unbalanced-shared-groups` remain frozen case identities; they do not assert
shared-offset or mixed-group invariants in the current implementation.

Natural scenarios expose truth for the target atom only, while terminal output
may also contain neighboring atoms. Truth RMSE keeps that existing target-truth
population; terminal validity and terminal-state digests still cover all emitted
atoms. Missing/duplicate target truth evidence fails; neighbors without truth do
not silently acquire invented ground truth.

### Replaying the HEAD pair

Use a fresh `build/p0-replay` directory. Build both binaries before timing and run
no other builds or tests during either corpus. The baseline source is exported
from the exact commit; the candidate is the reviewed working tree. The default
checked manifest is used, without a case filter. Both builds keep the trust-model
experiment disabled.

```sh
mkdir -p build/p0-replay/baseline-source
git archive a4354e698e77398154009d231907ebf3ed4b1d52 | tar -x -C build/p0-replay/baseline-source
for p0_side in baseline candidate; do
    p0_source=.
    if [ "$p0_side" = baseline ]; then
        p0_source=build/p0-replay/baseline-source
    fi
    cmake -S "$p0_source" -B "build/p0-replay/$p0_side-build" -G Ninja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo -DBUILD_TESTING=ON \
        -DBUILD_PYTHON_BINDINGS=OFF -DRHBM_GEM_ENABLE_UMAP=OFF \
        -DRHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE=ON \
        -DRHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT=OFF \
        -DRHBM_GEM_DEP_PROVIDER=SYSTEM
    cmake --build "build/p0-replay/$p0_side-build" --target convergence_exposure_case_runner -j 4
done
python3 resources/tools/developer/run_convergence_exposure_corpus.py \
    --executable build/p0-replay/baseline-build/bin/RHBM-GEM-CONVERGENCE-EXPOSURE \
    --output-dir build/p0-replay/baseline --threads 1 --jobs 1
python3 resources/tools/developer/run_convergence_exposure_corpus.py \
    --executable build/p0-replay/candidate-build/bin/RHBM-GEM-CONVERGENCE-EXPOSURE \
    --output-dir build/p0-replay/candidate \
    --reference-truth-dir build/p0-replay/baseline --threads 1 --jobs 1
```

The comparison JSON contains `safety_gate`, `quality_gate`, and `efficiency_gate`.
The runner exit status reports case execution/parsing failure; gate decisions
must be read from that JSON. Per-case logs and schema-1 diagnostics, trajectory-10,
terminal-2, and summary-14 artifacts remain under each run's `cases` directory.

### P0 result (2026-09-07)

**HEAD safety and trajectory-neutral quality recertification passed. Efficiency
strict-speedup failed. Historical frozen-truth continuity is not certified.**
The candidate adds diagnostics and audit tooling only; no numerical policy was
changed. Both binaries used AppleClang 21, system dependencies, RelWithDebInfo,
OpenMP AUTO, ROOT AUTO, UMAP OFF, and the trust-model experiment OFF. Each formal
run used one estimator thread and one job, without overlapping builds/tests.

| Gate / measurement | Result |
| --- | --- |
| Completed baseline / candidate | 600/600 / 600/600 |
| Failed cases / safety regressions | 0 / 0 on both sides |
| Safety gate | passed |
| Quality gate | passed; 600/600 semantic and terminal digests match |
| Per-case stop reason, objective, truth RMSE, accepted iterations | identical for all 600 |
| Objective / truth RMSE / accepted-iteration delta median, p90, p99 | all zero |
| Stop distribution | converged 54; audit-patience 289; all-rejected 165; maximum-iterations 92 |
| Efficiency gate | failed |
| Elapsed median, baseline → candidate | 0.170482 s → 0.174479 s (+2.34%) |
| Elapsed p90, baseline → candidate | 0.733874 s → 0.756307 s (+3.06%) |
| Elapsed p99, baseline → candidate | 1.738902 s → 1.723046 s (-0.91%) |
| Cases with lower elapsed time | 259/600 (43.17%) |

A single timing pair measures this audit execution, including additional Debug
serialization; it does not establish a production speedup or isolate logging
cost causally. Efficiency failure is retained independently of the two passing
gates. Final-polish safety statuses were relative-passed 39, failed 297, and
not-evaluated 264. Exactly the 39 relative passes were applied; the other 561
cases retained their base state.

Conditioning diagnostics cover 600/600 candidate cases. Counts below are system
construction events, not distinct atoms or cases. Positive-ratio distributions
exclude the explicitly counted zero failure sentinel. Full family/topology
strata, p01/median/p90 distributions, actual ridge ranges, solve status counts,
and active/nominal populations are in the machine-readable report.

| Phase / solve | Events | Conditioning guard | Zero sentinel | Positive ratio minimum |
| --- | ---: | ---: | ---: | ---: |
| candidate-polish/joint-polish | 15523 | 4783 | 0 | 1.618397e-13 |
| final-dependency-polish/joint-polish | 3475 | 1695 | 1282 | 1.618397e-13 |
| final-recertification/joint-offset | 1733 | 2 | 0 | 3.697002e-09 |
| outer-operator/joint-offset | 59691 | 8 | 0 | 9.348240e-10 |

The outer operator reports 719,712 nominal atom observations across 19,685 events:
508,817 shape endpoints were unavailable and zero offset endpoints were
unavailable. These are repeated observations, not unique failing atoms. Offset
solves nevertheless include 1,121 IRLS iteration-limit statuses; available
finite offset evidence must not be confused with solver qualification. No
boundary-reconciliation conditioning events were observed in this corpus; its
absence is not proof that the path is unreachable.

Validation passed the five core CTest groups (including estimator defense),
corpus contract/smoke/determinism tests, repository lint, and whitespace checks.
The corpus contract now has 11 tests, covering independent gates, missing and
nonfinite evidence, duplicate/incomplete pairs, diagnostic neutrality, and
target-only natural truth. Existing per-atom activity and logging-neutrality
regressions were retained; conditioning zero-column/nonfinite/threshold-boundary
coverage was added.

Baseline logs were reanalyzed after correcting the zero-accepted-iteration
safety check, without changing timing or numeric output. An initial candidate
parser-debugging batch, including one interrupted case, is retained separately
as `build/p0-recertification/candidate-preflight` and excluded from every gate
and timing statistic. The formal candidate is a fresh uninterrupted 600-case run.

- [Machine-readable P0 report](../../tests/benchmarks/per_atom_offset_recertification.json)
  records implementation patch/source/binary/library hashes, build settings,
  environment, artifact hashes, three gates, and diagnostic distributions.
- [Per-atom certificate baseline](../../tests/benchmarks/per_atom_offset_certificate_baseline.json)
  retains all 600 case identities and semantic/terminal digests. It does not
  overwrite the historical certificate baseline.
- Full logs, frozen truth, case summaries, diagnostics and comparisons remain in
  `build/p0-recertification/baseline` and `build/p0-recertification/candidate`.
  `build/p0-recertification/implementation.patch` captures the reviewed code/tool
  changes against the exact source baseline; these build artifacts are local.

## Historical recertification results

The following measurements apply only to their explicitly named historical
baselines. Their stop distributions are not acceptance targets for P0.

### Cleanup recertification result (2026-08-28)

The paired run used AppleClang 21, `RelWithDebInfo`, the checked manifest and
truth, one estimator thread, sequential jobs, and an exported `03cdb6ef`
baseline executable. All blocking conditions passed:

| Gate | Result |
| --- | ---: |
| Baseline / candidate completed | `600/600` / `600/600` |
| Failed cases | `0` / `0` |
| Production semantic digest matches | `600/600` |
| Normalized terminal-state digest matches | `600/600` |
| Stop distribution match | exact |
| Objective delta median / p90 | `0 / 0` |
| Transformed-truth RMSE delta median / p90 | `0 / 0` |
| Accepted-iteration delta median / p90 | `0 / 0` |
| Safety regressions | `0` |
| Elapsed median, baseline → candidate | `0.101151 s → 0.088944 s` |
| Elapsed p90, baseline → candidate | `0.708533 s → 0.527565 s` |

Timing is excluded from both semantic digests. With the experiment flag off,
the trust-model data structures and calculations are not compiled. Schema 9
contains maximum and tail diagnostics but no rho, comparator, exposure, or
accepted-only persistence fields.

### Non-converged final-polish residual safety result (2026-08-29)

The paired comparison used `8e65fa41` as the baseline and the same checked
600-case manifest, frozen truth, AppleClang 21 `RelWithDebInfo` build, one
estimator thread, and four corpus jobs. Both sides completed 600/600 cases with
zero failed cases and zero safety regressions. Production semantic digests
matched 600/600, stop-reason distributions remained exactly `converged 42`,
`audit-patience 372`, `all-rejected-backtracking-exhausted 163`, and
`maximum-iterations 23`, and accepted-iteration delta median/p90 was `0/0`.

The new safety result distribution was `relative-passed 11`, `failed 278`, and
`not-evaluated 311`; no corpus candidate required the absolute-pass fallback.
Every applied non-converged polish was one of the 11 relative passes, while all
278 failed candidates retained their base state. Of those failures, 134 had an
incomplete candidate operator and 144 had comparable evidence but regressed at
least one residual coordinate; one incomplete case was also solver-unqualified.

| Stop reason | Objective-accepted | Applied | Safety-rejected |
| --- | ---: | ---: | ---: |
| `all-rejected-backtracking-exhausted` | 104 | 8 | 96 |
| `audit-patience` | 167 | 3 | 164 |
| `maximum-iterations` | 18 | 0 | 18 |
| `converged` | 0 | 0 | 0 |

Terminal-state digests matched 322/600; the 278 intentional differences are
exactly the objective-accepted polishes rejected by the new safety gate.
Candidate-minus-baseline objective delta median/p90 was `0/0.160499662`, while
truth-RMSE delta median/p90 remained `0/0`: truth RMSE improved in 241 changed
cases and worsened in 37. The generic cleanup blocking gate reports failure
because it requires trajectory-neutral terminal states and zero outcome deltas;
those conditions do not apply to this intentional persistence-policy change.

## Historical provenance

The following records preserve the original measurements and decisions in
chronological order:

1. [Convergence safeguard audit](audit-history/second-stage-convergence/convergence-safeguard-audit.md)
2. [Stationarity and active-coordinate population audit](audit-history/second-stage-convergence/stationarity-active-coordinate-audit.md)
3. [Counterfactual convergence continuation audit](audit-history/second-stage-convergence/counterfactual-convergence-continuation-audit.md)
4. [Convergence exposure and counterfactual outcome audit](audit-history/second-stage-convergence/convergence-exposure-counterfactual-outcome-audit.md)

The historical records are provenance, not current production specification.
