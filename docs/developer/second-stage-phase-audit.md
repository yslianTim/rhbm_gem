# Second-stage phase audit

This developer diagnostic does not change candidate acceptance, fitting objectives,
trust radii, quarantine, patience, or convergence. It owns a copy of the iteration
context and model snapshots. Local workers write separate cluster buffers; the
selection synchronization point merges them in cluster-key order. Isolated solver
runs take place after production decisions, with fresh solver workspaces.

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

The analyzer writes `candidates.json`, `attempts.json`, `counters.json`, and
`report.md`. Keep external inputs and generated artifacts outside tracked source
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
  and 0.001. Only objectives are computed at additional sample points. These are
  observations, not a mathematical descent proof or a new acceptance search.

`Second-stage phase audit counters` records separate objective/sample evaluations,
operator evaluations, unavailable/error counts, and evaluator elapsed milliseconds.
The sample count is nominal full-domain work (an upper bound if an objective
evaluation exits early). These do not increment production objective/solver counters. Production wall time
necessarily includes diagnostic overhead; capture/serialization and evaluation
costs also affect observed timings. Compare operation counts and decisions, not
identical timing values.

Baseline and all requested main candidates receive isolated operator evaluations;
backtracking intermediate points and direction samples receive objectives only.
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
