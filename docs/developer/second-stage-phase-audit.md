# Second-stage phase audit

This developer diagnostic does not change candidate acceptance, fitting objectives,
trust radii, quarantine, patience, or convergence. It owns a copy of the iteration
context and model snapshots. Local workers write separate cluster buffers; the
selection synchronization point merges them in cluster-key order. Isolated solver
runs take place after production decisions, with fresh solver workspaces.

Production phase notifications now delegate capture and event preparation to
this observer. Trust-model trials and funnels have a separate `TrustModelAudit`
collector and translation unit. Historical member references are read-only
iteration-baseline snapshots from `ClusterHistoryObserver`; replay owns its copies
and does not hold a live history observer. Missing history remains unavailable
diagnostic evidence. Neither collector is stored in candidate
results. See [P0 structural refactoring](second-stage-p0-structure.md) for the
transaction boundary and certificate/diagnostic ownership.

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
python3 tools/second_stage_phase_audit.py build/phase-audit/run.log \
  --output-dir build/phase-audit/analysis
```

The analyzer writes `candidates.json`, `attempts.json`, `counters.json`,
`direction_samples.json`, and `report.md`. Keep external inputs and generated artifacts outside tracked source
files. The existing `tests/integration/fold_168_regression.py` runner can exercise
the trace by wrapping `build_command` to set the returned command's verbosity to 4;
leave its fixed template, baseline validation, and quality gates unchanged. Compare the same source with the CMake
option OFF and ON.

## Schema 1

Each `Second-stage phase audit: schema=1, payload=...` line contains one JSON
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
- `direction_samples`: production damping interpolation at 1, 0.5, 0.125, 0.03125,
  and 0.001. Boundary corrections on attempts 5–8 additionally include an isolated
  `operator`, `delta_baseline`, and `member_gates`; alpha 1 reuses the event's
  objective and operator. Other probes compute objectives only. These are
  observations, not a mathematical descent proof or a new acceptance search.
- `member_gates`: all component members are evaluated, including those after a
  failed member. Each record contains candidate, previous, stored-best, and
  candidate-neighbor reevaluated best objectives, plus previous/best gate status,
  tolerance, and margin. Previous objectives and best parameters are owned
  snapshots from the correction event. Positive margin is remaining slack;
  progress gates accept zero margin. Absent best history is `not-applicable`;
  missing required evidence is `unavailable`, never a pass. The global strict
  gate uses the correction's original improvement reference and requires positive
  margin. These objective gates do not certify full production acceptance.

The analyzer retains schema-1 compatibility: missing optional probe fields mean
not measured. Its direction-sample comparisons report objective and p99 deltas
against both parent and same-attempt baseline, using `R∞ = max(p99)` for the scalar
residual. Joint improvement flags require reproduced baseline evidence and
complete, solver-qualified operators at both comparison endpoints. The report
places attempts 7–8 at alpha 0.5 first and retains every measured point and blocker.

`Second-stage phase audit counters` records separate objective/sample evaluations,
operator evaluations, unavailable/error counts, and evaluator elapsed milliseconds.
The sample count is nominal full-domain work (an upper bound for member-only
evaluations or early exits). Member and best-reference evaluations contribute to
these diagnostic counts. These do not increment production objective/solver counters. Production wall time
necessarily includes diagnostic overhead; capture/serialization and evaluation
costs also affect observed timings. Compare operation counts and decisions, not
identical timing values.

Baseline and all requested main candidates receive isolated operator evaluations;
backtracking intermediate points and other direction samples receive objectives only.
Each diagnostic operator uses the fixed input activity and ridge multipliers,
a copied context, fresh workspaces, one solver thread, and disabled diagnostic
trace. Debug logging already serializes production cluster execution; the option
does not change that policy. Tests exercise worker collection independently in
serial and concurrent execution.

## Interpretation

Read attempts 5–8 first. Compare proposal to baseline, each polish to its search
parent, assembled polish to assembled search, and correction to its endpoint.
Do not sum local improvements to predict assembly: samples and neighbors couple
clusters. A correction improving its parent may still be worse than the baseline;
check both deltas and the production rejection reason.

The analyzer separates sampled oversized steps, no sampled descent, mixed
objective/residual changes, assembly deterioration, and globally improving
corrections rejected by gates. It cannot establish which objective should replace
an existing one, or whether untested intermediate candidates pass production
gates. Those require a separately reviewed optimization change.

## Objective / fixed-point compatibility observations

Attempts 5–8 also capture `post-joint-offset`; `unrestricted-proposal` remains the
complete post-shape operator endpoint. This is observation-only, including on
protected paths: an intermediate production state is not automatically a qualified
unrestricted endpoint.

`Second-stage solver audit: schema=1` records carry attempt, source, atom index
(for shape solves), and kind. `production` sources observe actual solver inputs and
outputs. Final-selection sources observe isolated operator evaluations at attempts
5 and 8. A private thread-local scope in `utils/hrl/EstimationAudit.hpp` is enabled
only by the trace; candidate-selection workers inherit owned scopes, and isolated
operator solves use one thread. No public execution options or result types change.

- Weighted-ridge records report the used MAD scale, weight summaries, ridge range,
  anchor norm, rows/columns, final update and normalized normal-equation residual.
  Their elapsed time measures the diagnostic calculation.
- MDPDE records distinguish the last beta/variance update, the normal residual
  with the actually used weights, and the residual after recomputing weights from
  returned beta/variance. Recomputed weights and variance are never committed.
  `elapsed_ms` includes the solve; `diagnostic_elapsed_ms` measures the additional
  equation evaluation. Shape-support records preserve retained sample indices.
- Normal residuals are divided by `max(1, ||right-hand side||₂)`. Variance equation
  deltas are divided by the maximum magnitude of current/refreshed variance and
  the configured weight floor. These are diagnostics, not new stopping criteria.

`Second-stage compatibility audit: schema=1` adds direction and gradient records.
Coordinates are log peak height, log width, and physical offset divided by the
**fixed parent** peak height. Linear paths in these coordinates agree with the
existing production damping path; they differ from the state-dependent offset
ratio used to report the fixed-point certificate. Directions are infinity-norm
normalized, and central differences use h=1e-3, 3e-4, 1e-4 without clipping invalid
perturbations. A derivative is stable only when all slopes have the same sign,
spread is at most 10% of maximum absolute slope, and each exceeds the roundoff
estimate `32*epsilon*max(1, |J+|, |J-|)/h`. Uncertain derivatives are not zeros.

Directions cover unrestricted proposals, local polish, assembled polish and
boundary corrections. Full coordinate gradients are limited to final-selection
on attempts 5 and 8. The analyzer writes `solver_audit.json`, `compatibility.json`
and `compatibility.md`, retaining legacy-log support. Gradients rank reliable
components alongside operator top atoms; they do not reuse the `1e-4` certificate
threshold. Compatibility objective counts/time, extra operator counts/time and
failures are separate counter fields; their work also contributes to total audit
counts. Solver equation failures remain explicit records when available, and must
not be interpreted as a successful nonlinear solve.

### Definitions being compared

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
satisfy this check. The optional `last_variance_relative_change` trace field
records this stopping metric (null when unavailable), distinct from the
refreshed-weight variance equation residual. A refreshed-weight normal residual
is not a stopping gate, so success alone does not certify that residual is small.
Each call starts from OLS and its sample residual variance.

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
