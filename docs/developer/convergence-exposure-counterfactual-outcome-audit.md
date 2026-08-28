# Convergence exposure and counterfactual outcome audit

> Current result (2026-08-28): the actual-reduction-aware radius-growth update
> passed a paired 600-case comparison against the guard/radius-decoupled
> `2d9b878c` baseline. The sampled trajectories were neutral. Production now
> requires accepted active-DOF p99 plus complete nominal-DOF fixed-point
> residual p99, with solver qualification. Maximum remains an independent
> diagnostic comparator.

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

## Purpose and scope

The refreshed audit searches for production convergence checkpoints at which
the legacy all-selected population, the removed maximum gate, solver
qualification, or either strict-operator policy rejects the stop. It then uses the isolated continuation from the
[third-round audit](counterfactual-convergence-continuation-audit.md) to compare
outcomes. A mismatch alone is an exposure, not evidence of quality loss.

The implementation remains behind
`RHBM_GEM_ENABLE_COUNTERFACTUAL_CONVERGENCE_AUDIT=ON`. It adds a test-only case
runner and developer scripts, but no production CLI, `FitOptions`, public API,
threshold, stopping expression, or model field. The search is not registered as
a normal CTest because its fixed budget is 600 complete fits.

## Historical reproducible corpus

[`convergence_exposure_manifest.json`](../../tests/benchmarks/convergence_exposure_manifest.json)
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
`scenario-truth.json`, `trajectory-schema-7.json`,
`counterfactual-schema-3.json`, `shadow-continuation-schema-2.json`, the
legacy shadow/terminal audit artifact, and schema-8 `case-summary.json`. The
aggregate report is schema 5. A paired run additionally writes schema-2
`comparison.json`, including terminal accepted-iteration deltas and the
guard/trust decoupling blocking gate. Old summaries
are intentionally not resumed across the baseline change.

## Exposure and outcome definitions

At the production checkpoint `T0`:

| Exposure | Comparator at `T0` | Isolated question |
| --- | --- | --- |
| `legacy_population` | Production qualification + all-selected p99 | Would the former population delay the current stop? |
| `maximum_gate` | Production qualification + active-DOF p99/max | Would restoring maximum delay the current stop? |
| `solver_qualification` | Solver qualification + active-DOF p99 | Would solver qualification delay the current stop? |
| `fixed_point_operator` | Production qualification + accepted active-DOF p99 + complete nominal-DOF operator p99 | Was the guarded proposal hiding a material fixed-point residual? |
| `fixed_point_operator_maximum` | Fixed-point operator policy + accepted/operator maximum | Does a sparse residual tail justify the maximum gate? |
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

The accepted-only shadow checkpoint is separate from those five comparators.
It requires production qualification, accepted active-DOF p99, a stable
objective domain, no quarantine transition, no genuine suspicious/hard failure,
and no rejection, but ignores fixed-point residual p99. Only the
first shadow checkpoint is retained and final polish runs on isolated state, so
it cannot alter the production trajectory or output.

## Replay selection and production decision

Replay cases are deterministic: take at most ten sorted case IDs from each
isolated exposure, then fill to 30 from the remaining sorted exposures. All
exposures, not only this replay subset, feed the decision statistics.

Each comparator requires at least 15 exposures and at least five from every
corpus family. At least 70% must show material objective or truth benefit,
material harm must be at most 10%, aggregate raw-residual median must not
worsen, and no safety regression may appear. Legacy population/maximum results
are reported as `rollback_candidate`; `solver-qualified` is reported as
`redesign_candidate`; the strict-operator policies are reported as
`promotion_candidate`. An incomplete run or unmet quota forces a
diagnostic-only conclusion.

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

All five counterfactual policies have no genuine production-checkpoint
exposure in this corpus, so their minimum exposure/family quotas remain unmet.
The 42 production convergence checkpoints agree with both strict-operator
comparators. Consequently `rollback_candidate=false`, `redesign_candidate=false`,
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

Current verification on 2026-08-28:

- audit-enabled CTest passes 21/21 and audit-disabled CTest passes 19/19;
- fold-168 still stops after seven accepted iterations with `audit-patience` in
  both builds and produces byte-identical `actual.json` files;
- C/N/O positive controls pass 3/3 and `natural-v00` passes 8/8;
- the schema-7 unified-stabilisation corpus completes 600/600 with 42 production
  convergence cases, zero failed cases, and zero safety regressions;
- against the frozen schema-6 baseline, transformed truth RMSE median improves
  from `0.0067853` to `0.0054550` and accepted-iteration p90 improves from
  `100` to `37`, but objective median worsens from `0.16498` to `0.23613`,
  accepted-iteration median worsens from `11` to `12`, and only `68.7%` of the
  150 overlap cases reduce objective evaluations. The planned quality and 70%
  case-level efficiency acceptance gates therefore remain unmet;
- the 137-case isolated continuation replay completes twice with no failures,
  no non-time artifact mismatch, and no production-policy recommendation;
- repository lint passes.
