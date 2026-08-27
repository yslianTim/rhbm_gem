# Convergence exposure and counterfactual outcome audit

> Current baseline (2026-08-27): production uses current stationarity and
> accepted/raw active-DOF p99 without a maximum gate. The original 600-case
> search remains historical evidence for that decision and was not rerun during
> the current schema refresh.

## Purpose and scope

The refreshed audit searches for production convergence checkpoints at which
the legacy all-selected population, the removed maximum gate, or strict
stationarity rejects the stop. It then uses the isolated continuation from the
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

Historical case directories used `trajectory-schema-2.json` and
`counterfactual-schema-1.json`. The refreshed runner writes `run.log`,
`scenario-truth.json`, `trajectory-schema-4.json`,
`counterfactual-schema-2.json`, and schema-5 `case-summary.json`. Old summaries
are intentionally not resumed across the baseline change.

## Exposure and outcome definitions

At the production checkpoint `T0`:

| Exposure | Comparator at `T0` | Isolated question |
| --- | --- | --- |
| `legacy_population` | Current stationarity + all-selected p99 | Would the former population delay the current stop? |
| `maximum_gate` | Current stationarity + active-DOF p99/max | Would restoring maximum delay the current stop? |
| `strict_stationarity` | Strict stationarity + active-DOF p99 | Would stricter qualification delay the current stop? |
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

## Replay selection and production decision

Replay cases are deterministic: take at most ten sorted case IDs from each
isolated exposure, then fill to 30 from the remaining sorted exposures. All
exposures, not only this replay subset, feed the decision statistics.

Each comparator requires at least 15 exposures and at least five from every
corpus family. At least 70% must show material objective or truth benefit,
material harm must be at most 10%, aggregate raw-residual median must not
worsen, and no safety regression may appear. Legacy population/maximum results
are reported as `rollback_candidate`; strict-DOF is reported as
`redesign_candidate`. An incomplete run or unmet quota forces a diagnostic-only
conclusion.

## Verification and observed evidence

The 600-case corpus was not rerun for the schema-4/schema-2 refresh, so this
document makes no new corpus outcome claim. Existing analyzer fixtures and
small runner tests validate the new classifications and artifact contract.
Fold-168 remains the external `audit-patience / no_convergence_trigger`
negative control; continuation must never suppress that stop. Audit-OFF builds
retain the production trajectory and output behavior.

Current refresh verification on 2026-08-27:

- audit-enabled CTest passes 21/21 and audit-disabled CTest passes 19/19;
- fold-168 stops after seven accepted iterations with `audit-patience` in both
  builds and produces byte-identical `actual.json` files;
- schema-4 aggregation reports seven records, no convergence trigger, and zero
  legacy-population, maximum-gate, or strict-stationarity exposures;
- repository lint passes;
- the 600-case corpus was intentionally not run, so no new rollback or redesign
  decision is claimed.
