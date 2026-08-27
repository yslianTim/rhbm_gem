# Convergence exposure and counterfactual outcome audit

> Historical audit status: the production baseline is now active-DOF p99
> without a maximum-change gate. The outcome criteria below remain a record of
> the experiment that informed that decision.

## Purpose and scope

This fourth-round audit searches deliberately for production convergence
checkpoints at which strict stationarity, the active shared-DOF population, or
both reject the stop. It then uses the isolated continuation from the
[third-round audit](counterfactual-convergence-continuation-audit.md) to compare
outcomes. A mismatch alone is an exposure, not evidence of quality loss.

The implementation remains behind
`RHBM_GEM_ENABLE_COUNTERFACTUAL_CONVERGENCE_AUDIT=ON`. It adds a test-only case
runner and developer scripts, but no production CLI, `FitOptions`, public API,
threshold, stopping expression, or model field. The search is not registered as
a normal CTest because its fixed budget is 600 complete fits.

## Reproducible corpus

[`convergence_exposure_manifest.json`](../../tests/benchmarks/convergence_exposure_manifest.json)
expands to three families, 25 variants per family, and eight replicas per
variant. Seeds are fixed by

\[
410000 + 10000\,\text{family} + 100\,\text{variant} + \text{replica}.
\]

- `natural` reuses `PotentialModelScenario`, its five `rhbm_test -t 5` atomic
  topologies, the full potential-fitting workflow, and five noise levels.
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

Each case directory contains `run.log`, `scenario-truth.json`,
`trajectory-schema-2.json`, `counterfactual-schema-1.json`, and
`case-summary.json`. Full logs and generated state stay in build artifacts.
Only the declarative manifest and compact aggregate evidence are repository
artifacts.

## Exposure and outcome definitions

At the production checkpoint `T0`:

| Class | Strict/current (`P1`) | Current/active-DOF (`P2`) |
| --- | --- | --- |
| stationarity-only | rejects | accepts |
| active-DOF-only | accepts | rejects |
| combined | rejects | rejects |
| policy agreement | accepts | accepts |

No production convergence trigger is a separate negative control. Outcomes are
evaluated on the objective domain frozen at `T0` and report transformed truth
RMSE for log peak, log width, offset/peak, and their aggregate. Objective and
truth improvements use the existing absolute `1e-8` plus relative `1e-3`
materiality rule. The report separates material benefit, neutral numerical
churn, material harm, unresolved budget exhaustion, and termination by an
existing safeguard.

Safety regressions are counted independently: new hard failures, non-finite
states, mixed shared-group activity, or new quarantine after `T0`. The report
also aggregates active-population `N<=91` evidence and maximum-only catches;
the maximum gate is not changed in this round.

## Replay selection and production decision

Replay cases are deterministic: take at most ten sorted case IDs from each
exposure class, then fill to 30 from the remaining sorted exposures. All
exposures, not only this replay subset, feed the decision statistics.

`P1`, `P2`, or `P3` is a fifth-round production candidate only when its
applicable set contains at least 15 exposures and at least five from every
required class, at least 70% show material objective or truth benefit, material
harm is at most 10%, the aggregate raw-residual median does not worsen, and no
safety regression appears. An incomplete run or unmet exposure quota forces a
shadow-only conclusion.

## Verification and observed evidence

The checked-in compact outcome records the completed fixed-budget result and
any corpus shortfall. Fold-168 remains the external `audit-patience / no
convergence trigger` negative control; continuation must never suppress that
stop. Audit-OFF builds retain the production trajectory and output behavior.
