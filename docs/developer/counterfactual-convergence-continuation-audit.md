# Counterfactual convergence continuation audit

> Current baseline (2026-08-28): production uses solver qualification plus
> accepted active-DOF and complete nominal-DOF fixed-point residual p99, with
> invariants and orthogonal blockers clear. The single
> `ConvergenceCertificate` is authoritative.

The production-equivalent `fixed-point-operator` duplicate has been removed.
The older active-operator-proposal residual comparator is now named
`historical-active-proposal`; it does not mean that solver qualification is
absent from current production.

## Purpose and isolation boundary

This historical audit asks whether continuing past a production convergence
stop under an isolated alternative population or residual policy yields
material progress. It follows the
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
| `production` | Full solver qualification | Accepted active-DOF p99 plus complete nominal-DOF strict operator p99 |
| `historical-all-selected` | Historical cluster qualification | All-selected accepted/active-proposal p99 |
| `historical-cluster-active-proposal-maximum` | Historical cluster qualification | Active accepted/active-proposal p99 plus maximum `< 1e-3` |
| `historical-active-proposal` | Full solver qualification | Accepted active-DOF plus active operator-proposal p99; historical comparator |
| `production-maximum` | Full solver qualification | Current production predicates plus accepted/operator maximum `< 1e-3` |

The historical-all-selected comparator includes fixed and quarantined zero
changes. The historical cluster/maximum comparator stays on active DOFs so it
isolates the older qualification/proposal/maximum conjunction. Active-member
summaries remain diagnostic but are not a continuation policy.

If all policy candidates agree at the original stop, the run ends normally
as `policy-agreement`. Otherwise the state immediately before the production
stop remains the common trajectory origin, and no estimator or safeguard is
changed while continuing.

Every first-stop checkpoint runs objective-only final dependency polish in
isolated solver workspaces. The simulated finalized state is logged for
comparison but never feeds back into the continuation state. This preserves the
pre-checkpoint trajectory and makes each checkpoint comparable to stopping and
finalizing at that iteration. Production `converged` finalization now adds a
separate persisted-state strict-operator certificate; that production-only
certificate is not retroactively part of these historical counterfactual
checkpoint outcomes. Final polish uses the checkpoint's current domain, while
every reported objective is reevaluated on the objective domain frozen at `T0`,
so cross-policy `delta J` does not confound fit progress with a changing domain.

## Records and analysis

Debug output uses an independent schema:

- `Counterfactual convergence checkpoint: schema=4, comparator-set=1` records policy, attempt,
  elapsed continuation time, activity/domain sizes, latest and best-audit
  objectives, final-polish result, population, and accepted/raw median, p99,
  and maximum transformed changes.
- `Counterfactual convergence atom: schema=4, comparator-set=1` records the finalized
  amplitude, width, and offset by atom serial ID together with group,
  shape/offset activity, and quarantine state.
- `Counterfactual convergence termination: schema=4, comparator-set=1` records policy
  completion, budget exhaustion, or the unchanged production safeguard that
  ended continuation.

The schema-3 reader is a normalization-only adapter for frozen baselines; a
record may not mix schema-3 names with schema-4 fields. The checkpoint also
reports multi-member shared groups, the distribution of
`Hmin/Hmedian`, and how often the weakest-peak member supplies the largest raw
offset/peak change. The analyzer joins these records with the schema-8
trajectory aggregation, which retains solver/path strata, active population
sizes, strict-operator availability and tail evidence, `N<=91`
p99-to-maximum evidence, and unique blockers.

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
`audit-patience`. The refreshed schema-7 records have converged joint solves,
no quarantine, and no production convergence checkpoint. Therefore the
build-gated run correctly reports
`no_convergence_trigger`; it does not suppress audit patience or manufacture a
counterfactual continuation. Its 168-atom quality gate remains within the
checked-in quality baseline. The normal and audit builds produced byte-identical
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

## Historical decision rules

- Restoring the legacy population or maximum gate requires a real isolated
  exposure followed by repeatable material benefit without a new failure.
- At the time of this audit, promoting solver qualification used the same
  requirement and was reported as a redesign candidate. Solver qualification
  has since become part of production together with the strict operator gate.
- A predicate mismatch without continuation benefit remains semantic evidence,
  not observed quality loss.
- Policies that only add attempts, numerical churn, or a later failure remain
  shadow-only.
- Maximum remains diagnostic; the historical cluster/maximum and
  `production-maximum` comparators may delay an audit build after a production
  convergence checkpoint.
- An empty active set remains labelled restricted/all-fixed; this diagnostic
  classification does not add a production stop condition.

The fold-168 result was a valid negative control but supplied no evidence for
the then-pending production redesign. It remains historical evidence rather
than a description of the current gate.

The targeted discovery and cross-case outcome rules are continued in the
[convergence exposure and counterfactual outcome audit](convergence-exposure-counterfactual-outcome-audit.md).

## Historical refresh verification (2026-08-27)

- Audit-enabled CTest passes 21/21; audit-disabled CTest passes 19/19.
- Fold-168 stops after seven accepted iterations with `audit-patience` in both
  builds. Schema-5 aggregation reports seven records, no convergence trigger,
  and zero exposures for all three comparators.
- Audit-enabled and audit-disabled fold-168 `actual.json` files are
  byte-identical; repository lint passes.
- The 600-case exposure corpus was not rerun.

## Historical verification status

- Normal and audit-gated builds compile successfully.
- The four pure continuation-controller cases pass in both builds; the audit
  build also passes quiet execution and audit-patience precedence checks.
- The analyzer and fold-runner fixtures pass, including legacy and multiline
  summary parsing.
- Fold-168 passed in both builds with seven accepted iterations,
  `audit-patience`, seven historical schema-2 records, and no convergence trigger.
- Repository lint and diff hygiene checks pass.
- The full normal CTest run passes 16 of 17 registered test groups. The sole
  failing group contains the same seven pre-existing log-format assertions as
  the baseline; the defense suite now passes 121 of 128 cases with no new
  failure.
