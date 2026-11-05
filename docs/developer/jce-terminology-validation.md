# Joint component terminology validation

## Configuration and scope

- Baseline: `149481138e939c0f20ac906dc53f60e1870ff3ab`, clean `develop`.
- Existing `build/release`: Release, AppleClang, C++20, system dependencies,
  testing and Python bindings enabled, UMAP/ROOT/OpenMP enabled.
- Builds: `cmake --build build/release --parallel 4`.
- Tests: sequential CTest with `OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4
  OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`. Existing per-group
  environment overrides remain unchanged. No simultaneous build/test runs.
- No branch, worktree, dependency, public API, persistence format, numerical
  policy or Notion change; no push.

The [P0 audit](jce-terminology-audit.md) records the complete rename list,
retained names, mathematical classifications and serialization consumers.
The [runtime guide](joint-component-runtime.md#constraint-status-working-sets-and-free-variables)
defines the production terminology.

## Tests and contract checks

| Stage | Release build | CTest | Exact numerical capture |
| --- | --- | --- | --- |
| P0 baseline | Passed | 23/23 | Captured before source edits |
| P1 constrained A/C | Passed | 3/3 Joint groups | Identical to P0 |
| P2 bounded LM | Passed | 3/3 Joint groups | Identical to P0 |
| P3 diagnostics/logger | Passed | 5/5 Joint, persistence and CLI groups | Identical to P0 |
| P4 final | Passed | 23/23, zero failures, 39.74 seconds | Identical to P0 |

Final registration coverage: **1062 registered, 1062 enabled, 1062 uniquely
covered, zero missing**. The intentional overlap between the general
persistence group and the focused Joint persistence filter is retained.

The final run includes numerical, conditioning, Gauss–Seidel BCD, Progress,
Runtime, Snapshot/contract, JSON/SQLite round-trip and CLI integration tests,
plus command, quickstart, Python binding and Python example smoke tests.
Existing expected values and tolerances were not weakened or removed.

Targeted semantic coverage:

- The existing A/C interior/boundary/release test also checks an exactly zero
  amplitude that remains free outside the working set. It verifies both signed
  C columns stay free, a blocked A is not free, a released A is free, and the
  final QR column order is `{0,1,3}` for the release fixture.
- The additional LM test places a stationary-null direction exactly on its
  bound beside a blocked informative direction. Both attain equality, but
  only the informative direction enters `working_bound_mask`; the trace's
  free-subproblem rank excludes the frozen direction.
- Existing LM tests retain their exact step/gradient/multiplier bit assertions,
  zero-valued released-coordinate trace, rank deficiency, damping bracket and
  bisection, trust-radius, weak-width and failure-reason coverage.
- The existing failed-block recovery test checks all four new A/C diagnostic
  keys and rejects the old keys. Existing Progress tests check `ac-kkt=` and
  retain their exact output-count/throttling assertions.

The initial P0 attempt used `OMP_NUM_THREADS=1`; 22/23 groups passed and the
unrelated Eigen thread-restoration test failed when it requested two threads.
Before any source edit, the run was repeated with the configuration above:
23/23 passed in 37.03 seconds. No test assertion, production setting or
numerical tolerance was changed. All subsequent stage runs passed.

## Numerical equivalence

An ignored test probe was compiled with the existing Release test flags and
linked to the existing instrumented test library. It is not part of production
or the registered test suite. Floating-point values are emitted as `uint64_t`
bit patterns, with sizes, masks, indices, counters, booleans and strings.
Timing and measured RSS are excluded; deterministic memory estimates remain.

| Capture | Cases | Lines | Bytes | SHA-256, identical before and after |
| --- | --- | ---: | ---: | --- |
| A/C and LM | 10 A/C + 10 LM | 2648 | 76579 | `9ee93751d848d4adcdf91374461af58f8df8e125cd3ece1127eb4d85ab7c3f7e` |
| BCD and component assembly | 8 scenarios, each direct and orchestrated | 764 | 858826 | `aa45a5e4256345c1410b087fd44a9844140e0629497ac5c4a4e103e498c973dc` |

`cmp` succeeds for both captures after every implementation stage and after
the final test/documentation changes. There are no numerical exceptions.

A/C cases cover positive/zero/blocked amplitudes, signed C, release, rank
deficiency, nonuniform column scales and four near-cancellation refinement
fixtures (including the existing unsuccessful certification outcome). Captured
values include beta, normalized free-design entries, column indices/masks,
rank, solves, face/iteration/release/refinement counts, feasibility, residual
statistics, objective, projected KKT and active amplitude indices.

LM cases cover interior/trust-radius boundary, new blocking, initially bound
release, constrained damping, rank-deficiency recovery, invalid dimensions,
infeasible bounds, stationary-null at a bound, small nonstationary sensitivity
and tiny-diagonal radius control. Captures include every damping attempt's
step, gradient, multipliers, working-set sequence, free-subproblem ranks,
KKT tolerances, phase, acceptance, status and counts, plus the returned step,
damping, KKT, reason/status and work counters.

BCD capture includes eta/B, beta/A/C, prediction, residual, objective, sweep
A/C and Width KKT, block acceptance/recovery, qualification, Independent
Replay, component assembly, stop reasons, progress events and work counts.

| BCD scenario | Objective, identical before/after | A/C KKT | Width KKT | Sweeps | Profiles |
| --- | ---: | ---: | ---: | ---: | ---: |
| Single block | 2.0199078132865863e-10 | 3.885780586188048e-16 | 4.0945025148175773e-13 | 3 | 6 |
| Three forward blocks | 2.0199079006559097e-10 | 3.1542578549093037e-09 | 1.1782208542143735e-10 | 6 | 36 |
| Active Width upper bound | 2.025400264756001e-10 | 3.3306690738754696e-16 | 1.1102230246251565e-16 | 2 | 4 |
| Nominal sweep limit | 2.0958526701391608e-10 | 2.9403302566066714e-06 | 1.0939575745272379e-07 | 4 | 24 |
| Failed-block recovery | 1.7934999880322473e-10 | 1.7828503395822537e-08 | 8.188935085584603e-10 | 5 | 20 |
| Adaptive extension | 8.938632823304983e-11 | 7.259212841947971e-09 | 6.936572427562737e-10 | 34 | 136 |
| Memory rejection | Not evaluated | Unavailable | Unavailable | 0 | 0 |
| Hard sweep cap | 8.938746938111263e-11 | 3.0175387077946425e-08 | 2.883422611787978e-09 | 31 | 124 |

The active-bound case retains one active Width bound; the recovery case
retains one successful recovery out of one attempt. Memory rejection preserves
the original resource-limit status without treating its default objective
storage as an evaluated objective.

## Source, names and persistence

A lexical comparison against the initial HEAD checks every changed production
C++ file after only the declared per-file identifier replacements, the five
logger-label replacements and comment removal. **All 14 changed production
C++ files have identical executable tokens under that mapping.** Thus the
arithmetic, constants, conditions, iteration limits, traversal, assembly,
refinement, damping, acceptance and certification expressions are preserved.
The Local Conditioning and Profile Jacobian changes are comments only.

Source/reference and exported-symbol checks find the renamed A/C entrypoints
and no retired `SolveLinear`, `CertifyLinear` or `LinearSolveWorkspace`
references/aliases in executable code. Old diagnostic keys occur only in
explicit negative test assertions and the audit table. Legitimate active
constraint, free-design, fixed-domain and algorithm names remain.
All 942 original test macro declarations are preserved under four declared
test-name changes; exactly one additional test declaration was added.

Public headers, persistence implementations, schema constants, Python source,
dependency and CMake files are unchanged. Every string literal in the
test-certificate JSON adapter is identical to P0, including `active_atoms`,
`free_rank` and `linear_solves`. SQLite v20, Stage JSON v4, Joint JSON v23,
solver provenance v6, convergence identifiers, evidence/status/reason strings
and the current strict readers remain unchanged. `git diff --check` passes.

Local build evidence is under `build/release/jce-terminology/`: stage logs,
baseline/current TSV captures, probe source/build runner, BCD case summary and
source-equivalence output. Validation was performed on the existing macOS
Release configuration and deterministic fixtures; no new platform or external
real-map run is claimed.
