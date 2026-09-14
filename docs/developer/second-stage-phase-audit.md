# Second-stage phase audit

This developer diagnostic does not change candidate acceptance, fitting objectives,
trust radii, quarantine, patience, or convergence. It owns a copy of the iteration
context and model snapshots. Local workers write separate cluster buffers; the
selection synchronization point merges them in cluster-key order. Isolated solver
runs take place after production decisions, with fresh solver workspaces.

Production phase notifications now delegate capture and event preparation to
this observer. Trust-model trials and funnels have a separate `TrustModelAudit`
collector and translation unit. Phase replay uses owned numerical snapshots
and does not retain cluster-history references. Neither collector is stored in
candidate results. See the [current state ownership](second-stage-local-fitting.md#implementation-responsibilities-and-state-ownership)
for the transaction boundary and certificate/diagnostic ownership.

## Observer and diagnostic contract

`RunSecondStageIterations` owns a non-copyable `SecondStageObservationSession`
separately from its numerical context. The session retains best trace, cluster
history, phase/trust collectors, and iteration/final-polish sidecars. Each
iteration's output survives until logging completes, then resets before the
next attempt; cluster history retains its existing cross-attempt lifecycle.
Candidate selection receives only a borrowed session pointer. Phase replay
copies the numerical context and model snapshots without retaining the live
session.

`SecondStageDiagnostics` holds history payloads, formatted comparison lines,
source IDs, observation tokens and output-only trial/sample/scale details.
Decision evidence is copied into these records for output, never read back by
production. Local worker records are allocated by key before parallel search;
final accepted/rejected output follows the existing key/rejection-event order.
Boundary scopes retain each trial's record index and publish the observation
matching the selected source, including endpoint fallback after correction.
Final-polish diagnostics likewise remain separate from its accepted state and
objective. Production progress labels and numerical log schemas remain unchanged;
phase events and counters use schema 2.

`ClusterHistoryObserver` owns Debug history, tie-break records, provisional
publication/rollback, and source IDs. It is created only for non-quiet Debug
runs; Info/quiet runs allocate no history maps and do not score historical
patches. Local workers use preallocated per-key entries. Component observations
start from iteration-baseline history, and only selected observations publish.
History selects retained diagnostic records, never production candidates.

Production publishes state, provenance, quarantine, and radii before notifying
the observer. Observer entry points contain exceptions; a failure disables
further history work for that run and reports unavailable diagnostics without
rejecting a candidate, rolling back production, or stopping estimation.
History evaluations use observer-owned counters. A change in those counters is
not evidence of a numerical trajectory change or a performance improvement.
`best_reference_unavailable` is history-update evidence, not a rejection gate.

`best_audit_state` remains production state for global-best acceptance, audit
patience, and non-converged final-state selection. It is independent of optional
cluster history. The earlier removal of member-best acceptance gates is a
separate numerical policy change; observer isolation does not resolve the
[open fold-168 regression](fold-168-iteration-regression.md).

The current Debug trajectory uses schema 10. It serializes the production
certificate, active/nominal populations, p99/maxima, operator completeness, and
four orthogonal blockers. The current convergence analyzer rejects earlier
trajectory schemas. The schema-compatible `suspicious-offset` label represents
`suspicious_block_fallback`, which also includes shape and hard-failure evidence.

Trust-model shadow records use schema 3, without the former fixed-false
`rejected-by-best` field; `analyze_trust_model_experiment.py` also accepts
historical schema 2. Funnel and trust analyzer-summary schemas remain unchanged. Keep/Grow/Shrink are diagnostic actions; production observation
receives a boolean shrink request. Trust-model work requires its independent
build flag, and disabled phase/trust entry points neither capture snapshots nor
run diagnostic solvers.

Schema-1 conditioning/solve records distinguish conditioning, solve status, and
per-atom availability across outer operator, candidate polish, boundary
reconciliation, final dependency polish, and final recertification. Offset
status codes follow `JointOffsetSolveStatus`: 0 converged, 1 system build failed,
2 empty, 3 initial solve failed, 4 IRLS solve failed, 5 objective deteriorated,
6 iteration limit. Hard failure is recorded separately. Joint-polish
`solved`/`failed` does not imply objective acceptance or persistence approval.
Reported ridge multipliers include all guards, not just the conditioning floor.

`final_uses_polish` records the persisted state's polish provenance, including
outer accepted-state polish. It is not a final dependency-polish application
counter. A `maximum-iterations` run may report it as true while making zero
final dependency-polish calls; use the explicit final-polish application records.

## Enable and analyze

```sh
cmake -S . -B build/phase-audit -DRHBM_GEM_ENABLE_SECOND_STAGE_AUDIT_TRACE=ON
cmake --build build/phase-audit
```

The option defaults to OFF. Both Debug logging (`-v 4`) and non-quiet execution are
required. The CMake build type is independent of this logging requirement. No CLI,
FitOptions, or database fields were added. Use the project's existing dependency
configuration when configuring a new build directory.

```sh
python3 resources/tools/developer/second_stage_phase_audit.py build/phase-audit/run.log \
  --output-dir build/phase-audit/analysis
```

The analyzer writes `candidates.json`, `attempts.json`, `counters.json`,
and `report.md`, with attempt summaries ordered by attempt number. Keep external
inputs and generated artifacts outside tracked source files. The existing
`tests/integration/fold_168_regression.py` runner can exercise
the trace by wrapping `build_command` to set the returned command's verbosity to 4;
leave its fixed template, baseline validation, and quality gates unchanged. Compare the same source with the CMake
option OFF and ON.

## Schema 2

Each `Second-stage phase audit: schema=2, payload=...` line contains one JSON
object. Candidate IDs identify an attempt, stage, key, and per-key sequence;
`key` contains zero-based selected-atom indices. Records are ordered by key and ID,
not execution chronology. Follow `parent_id` to reconstruct the update chain.
A `parent` record materializes an evaluator environment that was not already
captured. No comparisons cross attempts, even if their domain IDs match.

The frozen domain ID advances when the partition/domain is rebuilt or the frozen
background responses change. The collector owns the domain, sample assignments,
scales, weights, and immutable background of its attempt. Every objective is a
full global evaluation, including contributions to other atoms' samples.

- `objective`: fit, weighted tail, offset, and total; null means unavailable.
- `delta_parent` and `delta_baseline`: total-objective differences in this domain.
- `disposition` and `reason`: results of existing production checks. An observed
  state has no extra acceptance decision. `unavailable` includes missing search
  or polish proposals.
- `final_retained`: the candidate's member parameters exactly match final selection
  after boundary/quarantine fallback. This is parameter retention, not a claim
  that the entire candidate was independently committed. Rejected candidates are
  never marked retained. Global best-audit fallback after the iteration loop is
  outside this final-selection stage and remains visible in the production log.
- `operator`: isolated `T(candidate) - candidate`, not candidate-minus-parent.
  Coordinates are absolute changes in log peak height, log width, and
  offset-to-peak ratio, using the production percentile calculation. Population
  includes every selected atom for all three coordinates. Top atoms are the five
  largest residuals per coordinate, with deterministic index tie-breaking.
- `operator.status`: available only for complete, solver-qualified evidence.
  Unqualified complete solves retain raw residual measurements, explicitly marked
  unavailable for inference. Missing unrestricted shape-solver qualification on a
  protected-offset path is also unavailable; it is not inferred from a protected
  solve. Missing operator coordinates never become zero residuals.
- `operator_reproduced`: baseline residuals match production evidence atom by atom
  within `1e-10 + 1e-8 * abs(reference)` and both solves are qualified. The baseline
  also includes `production_operator`. No consistency conclusion is made for a
  round lacking reproducible, qualified evidence.

The producer and analyzer support schema 2 only. Earlier phase logs require the
analyzer from the matching Git revision. The current format records actual
candidates and operator evidence; it does not generate additional sampled states
or attempt-specific solver, member-gate, direction, or gradient investigations.

`Second-stage phase audit counters: schema=2` records `attempt`,
`objective_evaluations`, `objective_sample_evaluations`, `operator_evaluations`,
`failures`, and `elapsed_ms`. The sample count is nominal full-domain work, an
upper bound when evaluation exits early. These diagnostic counts do not increment
production objective/solver counters. Production wall time
necessarily includes diagnostic overhead; capture/serialization and evaluation
costs also affect observed timings. Compare operation counts and decisions, not
identical timing values.

Baseline and all requested main candidates receive isolated operator evaluations;
actual backtracking candidates receive objectives only (`operator=null`).
Each diagnostic operator uses the fixed input activity and ridge multipliers,
a copied context, fresh workspaces, one solver thread, and disabled diagnostic
trace. Debug logging already serializes production cluster execution; the option
does not change that policy. Tests exercise worker collection independently in
serial and concurrent execution.

## Interpretation

Read attempts in numerical order. Compare proposal to baseline, each polish to its search
parent, assembled polish to assembled search, and correction to its endpoint.
Do not sum local improvements to predict assembly: samples and neighbors couple
clusters. A correction improving its parent may still be worse than the baseline;
check both deltas and the production rejection reason.

The analyzer reports mixed objective/residual changes, assembly deterioration,
and globally improving corrections rejected by gates. It cannot establish which objective should replace
an existing one, or whether untested intermediate candidates pass production
gates. Those require a separately reviewed optimization change.

## Objective and operator contracts

Let `r` denote raw map residual, `X,y` a linear system, `D` the diagonal ridge,
`c₀` the previous offset, and `H` peak height. Definitions follow the current
implementations, not an assumed common optimizer objective.

| Stage | Parameters and samples | Equation / scale / reference |
|---|---|---|
| Joint offset | Physical per-atom offsets; raw samples of active cluster targets; outside-cluster selected neighbors and frozen background enter RHS | Solve `(XᵀWX+D)c=XᵀWy+Dc₀`; each IRLS step uses Cauchy weights from the preceding residual and its MAD scale, floored at `1e-12`; ridge starts at `1e-3` times column squared norm with conditioning multipliers |
| Shape MDPDE | Per-atom log response after subtracting neighbor and own offset response; retain positive responses inside fit range; `X=[1,-r_distance²/2]` | Beta solves `XᵀW(Xβ-y)=0`; weights and variance update together as below; trained alpha is fixed per atom |
| Local polish | Log-height/log-width and physical-offset increments for a cluster, on the supplied affected samples; other parameters/background frozen | One linearized weighted-ridge direction `(AᵀWA+D)δ=AᵀWr`; anchor is zero increment; current residual MAD/Cauchy weights are frozen for this solve |
| Boundary correction | Same increment coordinates on active interface/halo parameters and affected samples | Same frozen linearized surrogate; trust, suspicious-update checks and objective gates follow the direction solve |
| Audit objective | Full candidate responses; frozen owner fit/tail masks, sample counts and scales | Sum of normalized Cauchy fit loss + `0.25` times tail loss + offset-excess penalty; no solver ridge term |
| Certificate | Nominal selected-atom population in log-height/log-width/offset-to-peak ratio | Qualified complete `T(x)-x` and accepted movement must both have all coordinate p99 below `1e-4`, with no production blockers |

For shape MDPDE, with residual `e=y-Xβ` and variance `v`, the implemented weights
are `w=max(weight_min, exp(-alpha*e²/(2v)))` (alpha zero gives unit weights).
Invalid variance uses the weight floor. Given these weights, beta is weighted
least squares and `v_new=sum(w*e_new²)/(sum(w)-n*alpha*(1+alpha)^(-3/2))`.
Nonpositive denominator and nonfinite variance follow the existing numerical
fallback rules. Stopping requires both **squared beta change** and
`abs(v_new-v_old)/max(abs(v_new),abs(v_old),data_weight_min)` below the configured
tolerance. Invalid variances and the maximum-double fallback sentinel cannot
satisfy this check. Each call starts from OLS and its sample residual variance.

Audit uses `rho(u)=0.5*k²*log(1+(u/k)²)` with `k=1.345`, applied to residual divided
by the frozen owner scale. Fit/tail weights additionally use owner atom fraction
and respective sample count. Its offset penalty is `0.01/N` times squared excess
of the peak offset/signal ratio beyond the existing bound. This is not the
moving-anchor ridge penalty in a solver surrogate.

The complete operator composes joint offsets and local shape refits. Polish,
backtracking, member gates and final selection are not part of that operator.
The code therefore does not establish that `T(x)=x` implies audit stationarity,
nor that lowering audit loss lowers the operator residual. Dynamic scales,
log-response weighting, distinct samples/normalization, ridge anchoring and finite
inner stopping must be separated before attributing a trajectory conflict to a
single definition difference. A solved frozen surrogate is not proof that a
fixed nonlinear objective was minimized.

## Frozen-IRLS trust-model observation

The independent `RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT` option enables this
observation. It remains diagnostic-only and does not alter acceptance, radius
actions, stopping, or persistence. Logs are read with
`resources/tools/developer/analyze_trust_model_experiment.py`.

The audit-enabled build records a developer-only frozen-IRLS directional
prediction for every material base or polished trial that reaches the objective
gate. For the complete outer-previous-to-trial step `p`, it computes:

```text
ared = J(previous) - J(candidate)
r_lin = r_previous + J_r p
pred = sum(0.5 * sample_coefficient * frozen_Cauchy_weight
           * ((r_previous / scale)^2 - (r_lin / scale)^2))
       + exact_offset_penalty_reduction
rho = ared / pred
```

Fit samples use weight `1.0`, tail samples use `0.25`, and the production owner
cluster normalization and fixed fit/tail scales are reused. `J_r p` includes
selected targets and neighbours plus unselected contributors derived from
selected-group medians. The ratio is unavailable unless `pred` is finite,
positive, and larger than `1e-8 + 1e-3 * abs(J(previous))`.

The record status is one of `available`, `nonmaterial-step`,
`objective-unavailable`, `model-unavailable`, `residual-unavailable`,
`nonfinite`, `nonpositive-prediction`, or `nonmaterial-prediction`. A reported
counterfactual action uses rho bands at `0.25` and `0.75`, with `0.8` boundary
utilization required for growth. Objective backtracking remains the first
shrink rule, and unusable prediction falls back to the current actual-only
action. Only the final locally accepted candidate may be action-ready.
Boundary-reconciled, rescued, globally rejected, and non-final local records
are suppressed from action comparison while remaining in coverage and
calibration statistics. A separate funnel records generated, invalid,
trust-skipped, guard-rejected, nonmaterial, base-objective, and polish-objective
counts without evaluating a model before the objective gate.
