# Joint component active-set terminology audit

## Scope and baseline

The audit starts at `149481138e939c0f20ac906dc53f60e1870ff3ab` on
`develop`, with a clean working tree. Searches cover `src`, `include`,
`tests`, `docs`, resources, build registration and Python bindings. The
case-insensitive keyword scan of the Joint runtime and its developer guide
finds 638 lines (`active`, `inactive`, `free`, `working`, `fixed`, `blocked`,
`blocking`, `release`, `face`, `bound`, `rank`, `kkt`, `linear`, `qr`).
Definitions, callers, test references and serialization consumers were traced
before editing. The table groups occurrences of the same symbol and meaning.

Classes: **AC** actual active constraint; **WS** solver working set; **FV**
free variable; **ND** numerical frozen direction; **ALG** algorithm name;
**EN** ordinary English; **LOG** diagnostics/logger; **DATA** serialized field.
All changed C++ declarations are private under `src/core/detail`; none is an
installed public API. `Numerics.hpp`, `LinearSolve.hpp`,
`ActiveSetBoundedLM.hpp` and `GaussSeidelBcd.hpp` have cross-translation-unit
callers. Damping traces are compiled only with test instrumentation.

## P0 decisions and complete change list

| File | Symbol | Current meaning / class | Recommended name | Action / consumers | Contract risk |
| --- | --- | --- | --- | --- | --- |
| LinearSolve.hpp | `LinearSolveWorkspace` | A/C constrained QR scratch, ALG | `ConstrainedAcWorkspace` | P1; Numerics, Search, numerical tests | Private type only |
| Numerics.hpp, LinearSolve.cpp | `SolveLinear` | Exact fixed-width constrained A/C least squares, ALG | `SolveConstrainedAcLeastSquares` | P1; all callers and numerical tests | Private function only |
| LinearSolve.cpp | `UnweightedSolveImpl` | Implementation of that same A/C solve, ALG | `SolveConstrainedAcLeastSquaresImpl` | P1 | Translation-unit local |
| Numerics.hpp/.cpp | `CertifyLinear` | Feasibility/KKT/residual certificate for A/C, ALG | `CertifyConstrainedAcSolution` | P1; refinement, state evaluation, tests | Private function only |
| LinearSolve.hpp/.cpp | `free` | QR eligibility mask, FV; all C plus A outside WS | `is_free_variable` | P1; workspace and local alias | No extra mask |
| LinearSolve.hpp/.cpp | `columns`, `previous_columns` | Ordered QR columns and previous iteration's columns, FV | `free_column_indices`, `previous_free_column_indices` | P1; capacity assertions | Preserve column order |
| LinearSolve.hpp/.cpp | `selected` | Normalized QR design over free columns, FV | `free_design` | P1 | Preserve sparse assembly |
| LinearSolve.cpp | `has_previous_face` | Previous WS/free-column comparison available, WS | `has_previous_working_face` | P1 | Local state only |
| LinearSolve.cpp | `blocking` | First blocking amplitude column, WS | `blocking_amplitude_column` | P1 | Same tie/order policy |
| LinearSolve.cpp | `release` | Amplitude released from WS, WS | `released_amplitude_column` | P1 | Same gradient/tie policy |
| Numerics.hpp, LinearSolve.cpp | `free_face_changes`, `free_face_unchanged` | Free-variable face comparison counters, FV | `free_variable_face_changes`, `free_variable_face_unchanged` | P1; numerical tests | Not serialized |
| Numerics.hpp, LinearSolve.cpp | `LinearResult::releases` | A/C WS release count, WS | `working_set_releases` | P1; EvaluateProfile, tests | Not serialized |
| Numerics.hpp/.cpp | `Certificate::free_rank` | Rank of solved A/C free design, FV | `ac_free_design_rank` | P1; BCD diagnostics, test JSON adapter | Keep JSON `free_rank` |
| Numerics.hpp/.cpp | `Certificate::active_atoms`, certification local `active` | Atom indices with exactly zero A, AC | `active_amplitude_atom_indices` | P1; streaming certificate, tests, test JSON adapter | Keep JSON `active_atoms` |
| Search.cpp | `linear_workspace` | Reused A/C scratch, ALG | `ac_workspace` | P1 | Same lifetime/storage |
| ActiveSetBoundedLM.hpp/.cpp | Workspace, step, trace `active` | Width bounds held equal during subproblem, WS | `working_bound_mask` | P2; Search, InstrumentedLM, LM tests | No new mask; not persisted |
| ActiveSetBoundedLM.cpp | `free_indices` | QR step coordinates excluding WS and frozen directions, FV | `free_step_indices` | P2 | Preserve order |
| ActiveSetBoundedLM.cpp | `hit_indices` | Newly blocking upper bounds, WS | `newly_blocked_indices` | P2 | Same blocking order |
| ActiveSetBoundedLM.cpp | `stationary_null` | Unresolved sensitivity with roundoff gradient, ND | `stationary_null_direction_mask` | P2; comments distinguish from WS | Remains separate storage |
| ActiveSetBoundedLM.cpp | local `release` | Width coordinate released from WS, WS | `released_bound_index` | P2 | Same multiplier checks |
| ActiveSetBoundedLM.hpp/.cpp | `active_set_sequence` | Per-iteration WS masks, WS | `working_set_sequence` | P2; existing bit-pattern assertions | Test-only trace |
| ActiveSetBoundedLM.hpp/.cpp | `free_ranks` | Free augmented QR subproblem ranks, FV | `free_subproblem_ranks` | P2; LM tests | Test-only trace |
| ActiveSetBoundedLM.hpp/.cpp | `damping_releases` | WS releases across damping attempts, WS | `damping_working_set_releases` | P2; LM tests | Not serialized |
| ActiveSetBoundedLM, Numerics, Search | Step/trace/work `releases`, `LmTrial::active_set_releases` | WS release counts, WS | `working_set_releases` | P2; merge, Components evidence formatter | Preserve evidence strings |
| Search.cpp | `active_eta_max` | Configured upper bound, EN | `eta_upper_bound` | P2 | Same optional bound |
| Numerics.hpp/.cpp | `Certificate::linear_solves`, `linear_releases` | A/C QR solves and WS releases, LOG | `ac_linear_solve_count`, `ac_working_set_releases` | P3; BCD diagnostics, JSON adapter | Keep JSON `linear_solves` |
| GaussSeidelBcd.hpp/.cpp | `linear_rank` | Initial profile's A/C free-design rank, LOG | `ac_free_design_rank` | P3; FailureDiagnostics/formatter | Logger key changes; no JSON/SQLite path |
| GaussSeidelBcd.hpp/.cpp | `linear_solve_count` | Initial profile's A/C solves, LOG | `ac_linear_solve_count` | P3; same source certificate | Logger key only |
| GaussSeidelBcd.hpp/.cpp | `active_set_releases` | Initial profile's A/C WS releases, LOG | `ac_working_set_releases` | P3; same source certificate | Logger key only |
| GaussSeidelBcd.hpp/.cpp | `scaled_projected_kkt` | Initial profile's normalized A/C projected KKT, LOG | `ac_scaled_projected_kkt` | P3; same source certificate | Logger key only |
| JointCliProgressReporter.cpp | Progress `KKT=` | `JointBcdProgress::global_ac_kkt`, LOG | `ac-kkt=` | P3; progress assertions | Text label only |

## Names retained after semantic review

| File / consumers | Symbols | Meaning / class | Action and contract |
| --- | --- | --- | --- |
| Numerics.hpp; tests/support/JointTestContext.hpp | `LinearPolicy`, `LinearResult`, `EvaluationContext::linear` | Linear least-squares policy/result, ALG | Retain correct algorithm names; no need to rename every linear operation |
| Numerics, LinearSolve, ActiveSetBoundedLM, Components | `active_set_iterations`, `active_set_iteration_factor`, `maximum_active_iterations`, `damping_active_set_iterations`, `ActiveSetStep`, `ActiveSetBoundedLM` | Active-set algorithm, ALG | Retain; does not assert geometric constraint status |
| Numerics.hpp/.cpp | `WidthBoundKkt::active`, `active_width_bound_count`, `active_width_bounds` | Bound status determined from eta/slack with existing tolerance, AC | Retain exact checks and consumers |
| Search/Numerics | `active_bound_indices`, `active_bound_multipliers` | Returned working bounds are attained bounds of the valid step, AC | Retain; a subset, not an exhaustive geometric-bound mask (stationary-null directions can also lie on a bound) |
| GaussSeidelBcd | `width_bound_active` | Trial endpoint reaches width upper bound, AC/LOG | Retain endpoint comparison |
| LocalConditioning; public JointBlockConditioning; JointResultJson | `active_amplitude_count`, `active_width_bound_targets`, `free_ac_columns`, `free_ac_rank`, `free_ac_conditioning` | Endpoint active face and its free A/C design, AC/FV/DATA | Retain public types, JSON reader/writer and validation |
| LocalConditioning.cpp, TiledDerivative.cpp | `free_columns`, `free`, `free_design`, `free_design_rank` | Endpoint-face A/C columns (all C and nonzero/positive A), FV | Retain; add comment that this is not mutable solver WS; do not modify Profile Jacobian |
| LinearSolve.cpp | refinement `face`, `qr`, `rank`, `reduced`, `scales` | QR matrix/factorization/rank/normalization, ALG | Retain correct mathematical names |
| Numerics.hpp/.cpp | `linear_solver_valid`, `linear_solver_reason`, `projected_kkt` | A/C linear solver/certificate result, LOG | Retain; source is explicit in enclosing certificate and diagnostic adapter |
| Search, Components | `active_bound_occurrences` | Count of attained working bounds in successful returned steps, AC/LOG | Retain existing evidence key and meaning; not an exhaustive geometric-bound census |
| Components.cpp | `active-set-lm-subproblems`, `active-set-iterations`, `releases`, `active-bound-occurrences` | Algorithm/work evidence text, DATA | Retain exact evidence names and reason text |
| BcdPolicy, JointResultJson, public results | `active-set-lm`, solver provenance, convergence identifiers | Algorithm identity, ALG/DATA | Retain schema 23, SQLite 20, Stage JSON 4, provenance v6 |
| LinearSolve/LM/Numerics/LocalConditioning | `rank-deficient-free-subproblem`, `active-set-budget-exhausted`, `active-set-iteration-budget`, other status/reason strings | Failure/qualification contracts, ALG/DATA | Retain every literal |
| JointProgress.hpp | `global_ac_kkt`, `projected_width_kkt`, `width_gradient_inf_norm` | Already identifies distinct quantities, LOG | Retain observer shape, cadence and event order |
| ResourceWork, CompactSvd, StructuralPartition, TiledQR | free dimension estimates, rank boundaries, QR and topology names | FV/ALG/EN | Retain memory formulas, rank policies and traversal |
| PostFitPeeling, JointGroupFitting, WidthCalibration, JointComponentEstimator | fixed charge/neighbor/support, working initialization, bounded widths | Physical fixed contributions or ordinary English, EN | Retain unrelated domain meanings |
| tests/support/JointRuntimeJson.hpp and JSON-based fixture checks | `active_atoms`, `free_rank`, `linear_solves` | Test-facade serialized certificate, DATA | Retain literal keys; update only C++ member access |

No retired `FixedNeighborBlockCoordinate.cpp` link occurs in the current
developer documentation. No compatibility alias, persistent format change,
schema bump, extra dependency, new branch/worktree or Notion update is needed.
