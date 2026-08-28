---
Status: Historical audit record.
Current policy: second-stage-outer-iteration-algorithm-audit.md
---

# Convergence exposure and counterfactual outcome audit

> Current result (2026-08-28): the actual-reduction-aware radius-growth update
> passed a paired 600-case comparison against the guard/radius-decoupled
> `2d9b878c` baseline. The sampled trajectories were neutral. Production now
> requires accepted active-DOF p99 plus complete nominal-DOF fixed-point
> residual p99, with solver qualification, clear invariants, and clear
> orthogonal blockers. One `ConvergenceCertificate` owns the production
> decision. Maximum remains an independent diagnostic comparator.

The unified production controller uses one geometric factor search in the
order validity, trust, guard, then objective. Guard is feasibility-only,
guard-only factor reduction does not request radius shrink, and radius growth
requires boundary utilization of at least `0.8` plus material actual objective
reduction. On a `converged` stop, objective-accepted final dependency polish is
persisted only if the polished state independently passes the same strict
operator certificate; otherwise the already converged base state is retained.

## Actual-reduction radius-growth refresh

Radius growth now requires a boundary-active accepted step and an actual
objective reduction larger than the existing progress tolerance
`1e-8 + 1e-3 * abs(previous)`. It does not construct a predicted reduction or
an actual/predicted reduction ratio. Candidate acceptance, objective-induced
shrink, terminal rejection, and radius bounds are unchanged.

The update was evaluated against `2d9b878c` with the same 600-case manifest,
seeds, frozen truth, and single-thread fitting configuration. Both paired runs
completed 600/600 cases, and the existing comparison-schema-2 blocking gate
passed:

| Measure | Baseline | Actual-reduction growth |
| --- | ---: | ---: |
| Production convergence | 42 | 42 |
| Accepted-only shadow checkpoint | 138 | 138 |
| Accepted-iteration median | 12 | 12 |
| `audit-patience` stops | 372 | 372 |
| `all-rejected-backtracking-exhausted` stops | 163 | 163 |
| `maximum-iterations` stops | 23 | 23 |
| Safety regressions | 0 | 0 |

Paired objective, transformed-truth RMSE, and accepted-iteration deltas all
had median and p90 `0`; there were no material objective-benefit,
objective-harm, truth-benefit, or truth-harm cases. The result verifies that
the stricter growth semantic is safe for the sampled trajectories, but it does
not demonstrate a quality or efficiency improvement. The existing trust
controller test directly covers weak-reduction keep, material-reduction grow,
objective-backtracking shrink precedence, and guard-only backtracking.

## Predicted-reduction and rho shadow audit

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

This instrumentation does not modify acceptance, radius action, trajectory,
stopping, public settings, or the production trajectory schema. The corpus
report summarizes model coverage, rho calibration, current/shadow action
confusion, elapsed audit cost, and family, topology, base/polish, trial
disposition, factor, prediction status, boundary, cluster-size, and unselected
dependency strata. A conservative paired-corpus evidence gate may set
`model_based_controller_experiment_recommended = true`; it always reports
`production_promotion_recommended = false`. Production adoption requires a
separate shadow-controller trajectory replay.

The experiment recommendation requires all 600 paired cases to preserve their
schema-7 trajectory and terminal artifacts without safety regressions; usable
rho coverage of at least 70% overall, 60% per family, and 50% per topology;
at least 100 low, mid, and high ratios with 20 of each band per family plus 50
high-boundary ratios; median/p90 absolute calibration error no greater than
0.5/2.0 and family medians no greater than 0.75; at least 100 action changes
covering 1% of action-ready records, 25 shrink opportunities, and 25
growth-related opportunities; and audit/candidate-phase cost ratios no greater
than 25% median and 40% p90. Every failed condition and up to 30 priority replay
cases are written to the aggregate report.

The latest paired shadow audit used the same production trajectory before and
after adding objective-trial instrumentation. The baseline and two shadow runs
completed 600/600 cases with frozen truth and one fitting thread per case.
There were no failed cases or safety regressions, the paired blocking gate
passed, and terminal-state and schema-7 trajectory artifacts were identical for
every before/after case. Objective, transformed-truth RMSE, and
accepted-iteration deltas had median and p90 `0`; both objective and truth had
zero material benefit and harm cases. The two current
`trust-model-shadow-schema-2.json` artifacts were byte-identical for all 600
cases.

The candidate funnel generated 369,204 trials: 118,787 were trust-skipped, one
was guard-rejected, 24,245 were nonmaterial, 226,171 base trials reached the
objective gate, and 11,824 polish trials reached it. The 237,995 trial records
therefore exactly match the objective-evaluated count. Of 237,060 material
objective trials, only 5,507 (`2.323%`) produced a usable finite rho. The main
statuses were 225,799 `nonpositive-prediction`, 5,754
`nonmaterial-prediction`, and 935 `nonmaterial-step`.

Rho diversity was also insufficient: low/mid/high-interior/high-boundary counts
were `8/5/5,427/67`. Only three of 6,787 action-ready records changed action
(`0.0442%`), all growth-related, with no shrink opportunity. Instrumentation
cost relative to candidate phase was `35.47%` median and `45.00%` p90, above
the `25%/40%` gate. Thus
`model_based_controller_experiment_recommended = false` and
`production_promotion_recommended = false`; the frozen-IRLS rho model remains
diagnostic-only. The earlier 39,634-record, `11.10%` final-local-patch result is
historical and is not a valid current coverage denominator.

## Guard/trust-radius decoupling refresh

The guard-to-radius-shrink decoupling was evaluated against the unchanged
`7fb994a3` baseline with the same 600-case manifest, seeds, frozen truth, and
single-thread fitting configuration. Both paired runs completed 600/600 cases.
The comparison-schema-2 blocking gate passed:

| Measure | Baseline | Decoupled |
| --- | ---: | ---: |
| Production convergence | 42 | 42 |
| Accepted-only shadow checkpoint | 138 | 138 |
| Accepted-iteration median | 12 | 12 |
| `audit-patience` stops | 372 | 372 |
| `all-rejected-backtracking-exhausted` stops | 163 | 163 |
| `maximum-iterations` stops | 23 | 23 |
| Safety regressions | 0 | 0 |

Paired terminal objective and transformed-truth RMSE deltas both had median
and p90 `0`; there were no material objective-benefit, objective-harm,
truth-benefit, or truth-harm cases. Accepted-iteration delta also had median
and p90 `0`, and the objective-harm ratio was `0`. The corpus therefore verifies
that the semantic separation is neutral for its sampled trajectories; the
existing trust-controller test supplies direct coverage that guard-only factor
reduction does not request shrink while objective backtracking still does.

## Historical purpose and scope

The audit originally searched for production convergence checkpoints at which
the legacy all-selected population, the removed maximum gate, full solver
qualification, or either strict-operator policy rejected the then-current stop.
It then used the isolated continuation from the
[third-round audit](counterfactual-convergence-continuation-audit.md) to compare
outcomes. A mismatch alone is an exposure, not evidence of quality loss.

Current production has since promoted full solver qualification and the
complete nominal-DOF strict-operator p99 predicate. Consequently the historical
`fixed_point_operator` comparator is now equivalent to production, while
`solver_qualification` below still names the older active
operator-proposal-residual policy rather than the current gate.

The implementation remains behind
`RHBM_GEM_ENABLE_COUNTERFACTUAL_CONVERGENCE_AUDIT=ON`. It adds a test-only case
runner and developer scripts, but no production CLI, `FitOptions`, public API,
threshold, stopping expression, or model field. The search is not registered as
a normal CTest because its fixed budget is 600 complete fits.

## Historical reproducible corpus

[`convergence_exposure_manifest.json`](../../../../tests/benchmarks/convergence_exposure_manifest.json)
expands to three families, 25 variants per family, and eight replicas per
variant. Seeds are fixed by

\[
410000 + 10000\,\text{family} + 100\,\text{variant} + \text{replica}.
\]

- `natural` reuses `PotentialModelScenario`, its five `rhbm_test -t 5` atomic
  topologies, the full potential-fitting workflow, and five noise levels. Its
  element reference truth now uses the same zero-noise C/N/O reference method
  as `RunSimulationTestOnAtomicModel`; the pre-fix truth is frozen and reused by
  the post-fix run.
- `stationarity` varies near-collinearity, boundary conflict, weak peaks,
  noise, separation, and transformed initial perturbation.
- `population` creates fixed/quarantined blocks through finite estimator and
  fallback inputs, and varies active counts and shared-group imbalance. It does
  not write convergence masks directly.

The runner accepts one manifest-expanded `case_id`, seed, and thread count, so
every discovery can be replayed independently. The orchestrator runs all cases,
isolates a failed case, resumes only schema-compatible complete artifacts, and
does not stop when an exposure quota is reached. `--jobs` controls independent
case concurrency; `--threads` controls the fitting thread count inside a case.

```bash
python3 resources/tools/developer/run_convergence_exposure_corpus.py \
  --executable build-counterfactual/bin/RHBM-GEM-CONVERGENCE-EXPOSURE \
  --output-dir build-counterfactual/benchmark-results/convergence-exposure-corpus \
  --threads 1 --jobs 4
```

Historical case directories used `trajectory-schema-2.json` and
`counterfactual-schema-1.json`. The current runner writes `run.log`,
`scenario-truth.json`, `trajectory-schema-8.json`,
`counterfactual-schema-4.json`, `shadow-continuation-schema-2.json`, the
legacy shadow/terminal audit artifact, `trust-model-shadow-schema-2.json`, and
schema-11 `case-summary.json`. The trust-model artifact records the candidate
funnel plus trial identity, status, source, disposition, rejection causes,
prediction components, reductions, rho, boundary utilization, current/shadow
actions, objective-backtracking, and unselected dependency count. Per-record measured
audit time remains in `run.log` and case/aggregate analysis but is
excluded from the canonical trust-model artifact, so identical calculations
can be verified byte-for-byte without wall-clock noise. The aggregate report is
schema 7. A paired run additionally writes schema-3
`comparison.json`, including terminal accepted-iteration deltas and the
guard/trust decoupling blocking gate. Old summaries
are intentionally not resumed across the baseline change. A
`compact-baseline.json` freezes manifest, case/seed, truth, production semantic
trajectory, terminal-state digests, stop reasons, evaluation aggregates, and
the versioned certificate/comparator definitions without elapsed time.

## Historical exposure and outcome definitions

At the production checkpoint `T0`:

| Exposure | Comparator at `T0` | Isolated question |
| --- | --- | --- |
| `historical-all-selected` | Historical cluster qualification + all-selected p99 | Would the former population delay the stop? |
| `historical-cluster-active-proposal-maximum` | Historical cluster qualification + active-DOF proposal p99/max | Would the former proposal/maximum conjunction delay the stop? |
| `historical-active-proposal` | Solver qualification + accepted active-DOF and active operator-proposal p99 | Would the historical active-proposal policy delay the stop? |
| `production-maximum` | Current production predicates + accepted/operator maximum | Does a sparse residual tail justify restoring the maximum gate? |
| policy agreement | All comparators accept | No continuation is required. |

No production convergence trigger is a separate negative control. Outcomes are
evaluated on the objective domain frozen at `T0` and report transformed truth
RMSE for log peak, log width, offset/peak, and their aggregate. Objective and
truth improvements use the existing absolute `1e-8` plus relative `1e-3`
materiality rule. The report separates material benefit, neutral numerical
churn, material harm, unresolved budget exhaustion, and termination by an
existing safeguard.

Safety regressions are counted independently: new hard failures, non-finite
states, mixed shared-group activity, or new quarantine after `T0`. Maximum
statistics remain diagnostic and support the isolated rollback comparison.

The accepted-only shadow checkpoint is separate from those four comparators.
It requires the historical cluster qualification, accepted active-DOF p99, a stable
objective domain, no quarantine transition, no genuine suspicious/hard failure,
and no rejection, but ignores fixed-point residual p99. Only the
first shadow checkpoint is retained and final polish runs on isolated state, so
it cannot alter the production trajectory or output.

## Historical replay selection and production decision

Replay cases are deterministic: take at most ten sorted case IDs from each
isolated exposure, then fill to 30 from the remaining sorted exposures. All
exposures, not only this replay subset, feed the decision statistics.

Each comparator requires at least 15 exposures and at least five from every
corpus family. At least 70% must show material objective or truth benefit,
material harm must be at most 10%, aggregate raw-residual median must not
worsen, and no safety regression may appear. Historical population/maximum
results are reported as `rollback_candidate`; `historical-active-proposal` is
reported as `redesign_candidate`; and `production-maximum` is reported as
`promotion_candidate`. The strict-operator and solver-qualification redesign
was subsequently promoted into production. An incomplete run or unmet quota
still forces a diagnostic-only conclusion for any remaining comparator.

## Verification and observed evidence

The 600-case corpus was run before and after the benign fixed-endpoint and
audit-patience/quarantine lifecycle fix, using the same manifest, seeds,
samples, and frozen reference truth. Both runs completed 600/600 cases without
runner failure.

| Measure | Before | After |
| --- | ---: | ---: |
| Production convergence | 0 | 42 |
| Accepted-only shadow checkpoint | 29 | 137 |
| `audit-patience` stops | 406 | 364 |
| `all-rejected-backtracking-exhausted` stops | 126 | 126 |
| `maximum-iterations` stops | 68 | 68 |
| Safety regressions | 0 | 0 |

The 42 post-fix convergence cases comprise 41 `natural` cases and one
`stationarity` case. All C/N/O zero-noise positive controls converge, and all
eight `natural-v00` single-carbon zero-noise replicas converge. By natural
topology the convergence counts are C-C 8, CA-C 8, N-N 8, O-O 9, and UNK-C 8.

The accepted-only shadow checkpoint appears in 137 post-fix cases: 125
`natural` and 12 `stationarity`. Relative to the terminal outcome it could have
saved a median of 3 attempts, p90 22.4, and at most 57 attempts. The isolated
continuation follow-up below evaluates whether persistence can turn that
opportunity into a safe stopping rule.

On the frozen objective domain, paired objective delta has median/p90 0/0 with
zero material benefit and zero material harm. Transformed aggregate truth RMSE
delta also has median/p90 0/0; one case (`natural-v07-r4`) is classified as
material truth harm (`+1.7153815119742974e-05`) while its objective delta is not
material. No non-finite, hard-failure, suspicious/quarantine, or other safety
regression is observed.

All four non-production comparators have no genuine production-checkpoint
exposure in this corpus, so their minimum exposure/family quotas remain unmet.
The 42 production convergence checkpoints agree with every comparator.
Consequently `rollback_candidate=false`, `redesign_candidate=false`,
and `promotion_candidate=false` mean only that no comparator has enough
exposure for a policy change; they do not prove the production policy is
optimal. The accepted-only shadow result is reported independently and is not
counted as comparator exposure.

## Accepted-only isolated continuation outcome

The 137 shadow cases were replayed twice with identical production terminal
states and frozen truth. Each candidate checkpoint ran final dependency polish
in isolated workspaces and the production trajectory continued unchanged. All
827 non-time JSON artifacts and both decision reports are deterministic across
the two runs.

| Audit-only policy | Effective exposures (natural/stationarity) | Attempts saved total | Objective harm | Truth-harm cases | Safety | Decision |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| accepted-only K=2 | 70 (61/9) | 954 | 0 | 34 | 0 | retain shadow-only |
| accepted-only K=3 | 63 (55/8) | 847 | 0 | 22 | 0 | retain shadow-only |
| accepted-only K=5 | 59 (51/8) | 713 | 0 | 12 | 0 | retain shadow-only |
| dynamic raw | 3 (2/1) | 24 | 0 | 2 | 0 | insufficient exposure and truth harm |

The persistence policies meet the exposure quota and save substantial work,
but none satisfies the zero-material-harm truth rule. K=5 reduces the number of
truth-harm cases relative to K=2/K=3 but does not eliminate them. Dynamic raw
uses `min(1e-3, 1e-4 * 2^(streak-1))`; it reaches only three effective
exposures and two of those have material truth harm. Every candidate has finite
complete comparisons, zero endpoint safety violations, zero continuation
safety events, and zero material objective harm.

The earlier accepted-only persistence and dynamic-raw experiments remain
historical controls. They are not the production policy after the unified
stabilisation change.

## Rejected coarse-to-fine candidate-factor refinement

A production candidate evaluated a main-cluster-only step-size refinement on
2026-08-28. The existing `1, 1/2, 1/4, ...` search first found a passing lower
bound, then attempted at most two arithmetic midpoints against the immutable
pre-search objective state. The largest passing factor was retained. Boundary
reconciliation and rescue backtracking were unchanged.

The same 600-case manifest, seeds, frozen truth, and single-thread estimator
configuration were run against `32e484a7` and the candidate. Both sides
completed 600/600 cases with zero runner failures and zero safety regressions.
Refinement was materially exercised: 8,853 accepted local records in 521 cases
used non-dyadic factors.

| Measure | Baseline | Candidate |
| --- | ---: | ---: |
| Production convergence | 42 | 42 |
| Accepted-iteration median | 12 | 11 |
| Accepted-iteration delta p90 | - | +3.1 |
| Objective delta median / p90 | - | 0 / +0.000370533 |
| Objective material benefit / harm | - | 55 / 62 |
| Truth RMSE delta median / p90 | - | +2.54642e-7 / +0.000630830 |
| Truth material benefit / harm | - | 138 / 277 |
| Unified search trials | 361,133 | 364,169 |
| Objective-rejected trials | 201,068 | 207,249 |
| Candidate-phase time sum (ms) | 83,330.522 | 92,886.559 |
| Candidate-phase time median / p90 (ms) | 41.1025 / 460.048 | 45.3795 / 469.713 |

The existing blocking gate failed because objective material harm reached
`62/600 = 10.33%`, above the permitted 10%. Accepted-iteration p90 and truth
RMSE p90 also worsened, while unified trials grew 0.84%, objective-rejected
trials grew 3.07%, and aggregate candidate-phase time grew 11.47%. The strategy
was therefore not promoted; production code and the main algorithm description
remain on the original geometric first-passing policy.

Current verification on 2026-08-28:

- the converged final-dependency-polish path now re-evaluates a strict operator
  certificate on the polished candidate; the existing core-estimator suite
  passes with the certificate and base-state fallback behavior;
- the certificate-consolidated audit build and complete C++/Python CTest suite
  pass 20/20;
- fold-168 still stops after seven accepted iterations with `audit-patience` in
  both builds and produces byte-identical `actual.json` files;
- C/N/O positive controls pass 3/3 and `natural-v00` passes 8/8;
- the schema-8 certificate corpus completes 600/600 with zero failed cases,
  zero safety regressions, 600/600 normalized production semantic matches, and
  zero exposure for every renamed comparator. Stop counts remain 42
  `converged`, 372 `audit-patience`, 163
  `all-rejected-backtracking-exhausted`, and 23 `maximum-iterations`;
- the paired objective, transformed-truth RMSE, and accepted-iteration deltas
  all have median/p90 `0`; the 1.3 GB case artifacts remain under the audit
  build, while the tracked compact baseline is
  [`convergence_certificate_baseline.json`](../../../../tests/benchmarks/convergence_certificate_baseline.json);
- the frozen manifest, case identity, and truth hashes are respectively
  `2b0d74c249a4723575696d98f99df9d3c449a30837adf8117aaf744145d0a87a`,
  `63985f81c5b188cfc9992742cad1f4b9418cc36c88b12398e5775a324289bd98`,
  and `04e3cda3b49857b2d5e4f63e973b2392dfe3095360974db988770e7468edd628`;
- in the earlier unified-stabilisation comparison against its frozen schema-6
  baseline, transformed truth RMSE median improved
  from `0.0067853` to `0.0054550` and accepted-iteration p90 improved from
  `100` to `37`, but objective median worsens from `0.16498` to `0.23613`,
  accepted-iteration median worsens from `11` to `12`, and only `68.7%` of the
  150 overlap cases reduce objective evaluations. The planned quality and 70%
  case-level efficiency acceptance gates therefore remain unmet;
- the objective-trial rho audit records 237,995 objective-evaluated trials,
  only 5,507 usable ratios (`2.323%`), just 8 low and 5 mid ratios, three
  action divergences, and `35.47%/45.00%` median/p90 instrumentation cost;
  its evidence gate is no-go and diagnostic-only;
- the 137-case isolated continuation replay completes twice with no failures,
  no non-time artifact mismatch, and no production-policy recommendation;
- repository lint passes.
