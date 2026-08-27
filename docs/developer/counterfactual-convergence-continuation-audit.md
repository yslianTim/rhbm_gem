# Counterfactual convergence continuation audit

> Historical audit status: production now uses current stationarity plus
> accepted/raw active-DOF p99 and no maximum-change gate. The legacy policy
> comparisons below remain developer diagnostics.

## Purpose and isolation boundary

This third-round audit asks whether a production convergence stop rejected by
strict stationarity or the active shared-DOF population would yield material
progress if iteration continued. It follows the
[stationarity and active-coordinate audit](stationarity-active-coordinate-audit.md).

The experiment is compiled only when
`RHBM_GEM_ENABLE_COUNTERFACTUAL_CONVERGENCE_AUDIT=ON`. The option defaults to
`OFF`, requires `BUILD_TESTING=ON`, and adds no CLI, `FitOptions`, public API,
threshold, or model-format field. A normal build retains the original stopping
expression and trajectory.

The audit suppresses only a production `converged` return when at least one
candidate policy rejects that checkpoint. Audit patience, all-rejected
resolution, empty active sets, quarantine, topology changes, objective gates,
and the original 100-attempt limit keep their production behavior. A
continuation additionally stops after 10 accepted iterations or 25 attempts
beyond the production checkpoint.

## Policies and trajectory

One continuation trajectory evaluates the first stop checkpoint for:

| Policy | Stationarity | Change population |
| --- | --- | --- |
| `production` | Current cluster rollup | All selected atoms |
| `strict-current` | Strict active-block qualification | All selected atoms |
| `current-dof` | Current cluster rollup | Active shape members and shared-offset DOFs |
| `strict-dof` | Strict active-block qualification | Active shape members and shared-offset DOFs |
| `strict-member` | Strict active-block qualification | Active members; diagnostic only |

If all production candidates agree at the original stop, the run ends normally
as `policy-agreement`. Otherwise the state immediately before the production
stop remains the common trajectory origin, and no estimator or safeguard is
changed while continuing.

Every first-stop checkpoint runs final dependency polish in isolated solver
workspaces. The simulated finalized state is logged for comparison but never
feeds back into the continuation state. This preserves the pre-checkpoint
trajectory and makes each checkpoint equivalent to stopping and finalizing at
that iteration. Final polish uses the checkpoint's current domain, while every
reported objective is reevaluated on the objective domain frozen at `T0`, so
cross-policy `delta J` does not confound fit progress with a changing domain.

## Records and analysis

Debug output uses an independent schema so the second-round schema `2` remains
stable:

- `Counterfactual convergence checkpoint: schema=1` records policy, attempt,
  elapsed continuation time, activity/domain sizes, latest and best-audit
  objectives, final-polish result, population, and accepted/raw median, p99,
  and maximum transformed changes.
- `Counterfactual convergence atom: schema=1` records the finalized
  amplitude, width, and offset by atom serial ID together with group,
  shape/offset activity, and quarantine state.
- `Counterfactual convergence termination: schema=1` records policy
  completion, budget exhaustion, or the unchanged production safeguard that
  ended continuation.

The checkpoint also reports multi-member shared groups, the distribution of
`Hmin/Hmedian`, and how often the weakest-peak member supplies the largest raw
offset/peak change. The analyzer joins these records with the existing schema-2
trajectory aggregation, which retains solver/path strata, active population
sizes, `N<=91` p99-to-maximum evidence, and unique blockers.

```bash
python3 resources/tools/developer/analyze_counterfactual_convergence.py \
    run.log --format json --output counterfactual-audit.json
```

An optional truth JSON may be supplied with `--truth`. It contains either an
`atoms` list or a serial-ID mapping with `amplitude`, `width`, and `offset`.
Optimization outcomes and truth-based estimation errors are reported
separately. The aggregate report includes exposure counts and overlaps,
termination categories, unresolved policies, and material-improvement counts.
Material objective improvement uses the existing progress tolerance: absolute
`1e-8` plus relative `1e-3`.

## Fold-168 evidence

The external fold-168 CIF and map hashes match the fixed repository fixture.
The audit runner uses Debug verbosity, keeps all artifacts in the specialized
build directory, and never rewrites the checked-in quality baseline. The
regular runner now parses both the legacy one-line and current multiline
second-stage summaries without changing its command template.

The observed fold-168 trajectory has seven accepted iterations and stops with
`audit-patience`. All seven schema-2 records have converged joint solves, no
quarantine, and identical current, active-member, and active-DOF predicate
sequences. Therefore the build-gated run correctly reports
`no_convergence_trigger`; it does not suppress audit patience or manufacture a
counterfactual continuation. Its 168-atom quality gate remains within the
checked-in schema-5 baseline. The normal and audit builds produced byte-identical
`actual.json` state, summary, and quality records for this negative control.

Configure the external negative control without committing local input paths:

```bash
cmake -S . -B build/counterfactual-audit \
  -DBUILD_TESTING=ON \
  -DBUILD_PYTHON_BINDINGS=OFF \
  -DRHBM_GEM_ENABLE_COUNTERFACTUAL_CONVERGENCE_AUDIT=ON \
  -DRHBM_GEM_ENABLE_FOLD_168_REGRESSION=ON \
  -DRHBM_GEM_FOLD_168_MODEL=/path/to/fold_test_model_0.cif \
  -DRHBM_GEM_FOLD_168_MAP=/path/to/sim_map_gaus_grid0.10_charge1_bw0.50.map
cmake --build build/counterfactual-audit --target tests_all -j
ctest --test-dir build/counterfactual-audit \
  -R counterfactual_fold_168_audit --output-on-failure
```

## Decision rules

- A strict or active-DOF production redesign requires a real convergence
  exposure followed by material objective improvement, or improved truth-based
  parameter error without a new failure.
- A predicate mismatch without continuation benefit remains semantic evidence,
  not observed quality loss.
- Policies that only add attempts, numerical churn, or a later failure remain
  shadow-only.
- The maximum gate is retained. Its next ablation decision uses the active
  `N<=91` frequency and unique catches from the joined schema-2 report.
- An empty active set remains restricted/all-fixed termination, not full
  convergence.

The current fold-168 result is a valid negative control but supplies no evidence
for changing production convergence. A representative trajectory with an
actual production convergence exposure is still required before a fourth-round
redesign or ablation.

The targeted discovery and cross-case outcome rules are continued in the
[convergence exposure and counterfactual outcome audit](convergence-exposure-counterfactual-outcome-audit.md).

## Verification status

- Normal and audit-gated builds compile successfully.
- The four pure continuation-controller cases pass in both builds; the audit
  build also passes quiet execution and audit-patience precedence checks.
- The analyzer and fold-runner fixtures pass, including legacy and multiline
  summary parsing.
- Fold-168 passes in both builds with seven accepted iterations,
  `audit-patience`, seven schema-2 records, and no convergence trigger.
- Repository lint and diff hygiene checks pass.
- The full normal CTest run passes 16 of 17 registered test groups. The sole
  failing group contains the same seven pre-existing log-format assertions as
  the baseline; the defense suite now passes 121 of 128 cases with no new
  failure.
