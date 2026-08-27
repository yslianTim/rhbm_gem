# Convergence exposure and counterfactual outcome audit

> Current result (2026-08-27): the post-fix (`after`) 600-case corpus is the
> canonical result for the current implementation. The instrumentation-only
> pre-fix (`before`) corpus is retained as the paired historical baseline.
> Production still requires its existing convergence qualification and
> accepted/raw active-DOF p99 without a maximum gate; the accepted-only
> checkpoint remains audit-only evidence and does not change that predicate.

## Purpose and scope

The refreshed audit searches for production convergence checkpoints at which
the legacy all-selected population, the removed maximum gate, or solver
qualification rejects the stop. It then uses the isolated continuation from the
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
`scenario-truth.json`, `trajectory-schema-5.json`,
`counterfactual-schema-3.json`, the shadow/terminal audit artifacts, and
schema-7 `case-summary.json`. The aggregate report is schema 4. Old summaries
are intentionally not resumed across the baseline change.

## Exposure and outcome definitions

At the production checkpoint `T0`:

| Exposure | Comparator at `T0` | Isolated question |
| --- | --- | --- |
| `legacy_population` | Production qualification + all-selected p99 | Would the former population delay the current stop? |
| `maximum_gate` | Production qualification + active-DOF p99/max | Would restoring maximum delay the current stop? |
| `solver_qualification` | Solver qualification + active-DOF p99 | Would solver qualification delay the current stop? |
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

The accepted-only shadow checkpoint is separate from those three comparators.
It requires production qualification, accepted active-DOF p99, a stable
objective domain, no quarantine transition, no genuine suspicious/hard failure,
and no rejection, but ignores raw p99 while recording it. Only the first shadow
checkpoint is retained and final polish runs on isolated state, so it cannot
alter the production trajectory or output.

## Replay selection and production decision

Replay cases are deterministic: take at most ten sorted case IDs from each
isolated exposure, then fill to 30 from the remaining sorted exposures. All
exposures, not only this replay subset, feed the decision statistics.

Each comparator requires at least 15 exposures and at least five from every
corpus family. At least 70% must show material objective or truth benefit,
material harm must be at most 10%, aggregate raw-residual median must not
worsen, and no safety regression may appear. Legacy population/maximum results
are reported as `rollback_candidate`; `solver-qualified` is reported as
`redesign_candidate`. An incomplete run or unmet quota forces a diagnostic-only
conclusion.

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
saved a median of 3 attempts, p90 22.4, and at most 57 attempts. This is evidence
for a future raw-p99 study, not permission to remove or relax the production raw
p99 `< 1e-4` gate.

On the frozen objective domain, paired objective delta has median/p90 0/0 with
zero material benefit and zero material harm. Transformed aggregate truth RMSE
delta also has median/p90 0/0; one case (`natural-v07-r4`) is classified as
material truth harm (`+1.7153815119742974e-05`) while its objective delta is not
material. No non-finite, hard-failure, suspicious/quarantine, or other safety
regression is observed.

The three counterfactual policies still have no genuine production-checkpoint
exposure in this corpus, so their minimum exposure/family quotas remain unmet.
Consequently `rollback_candidate=false` and `redesign_candidate=false` still
mean only that no evidence currently requires a rollback or redesign; they do
not prove the production policy is optimal. The accepted-only shadow result is
reported independently and is not counted as comparator exposure.

Current verification on 2026-08-27:

- audit-enabled CTest passes 21/21 and audit-disabled CTest passes 19/19;
- fold-168 still stops after seven accepted iterations with `audit-patience` in
  both builds and produces byte-identical `actual.json` files;
- C/N/O positive controls pass 3/3 and `natural-v00` passes 8/8;
- both 600-case runs are complete and deterministic, with zero safety
  regression;
- repository lint passes.
