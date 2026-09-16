# Production failed-only refinement

Second-stage shape solves now enable failed-only refinement by default. The native
OLS / MDPDE fixed-point solve still runs first, with its existing iteration limit,
tolerance and weight floor. A successful native result is returned unchanged.
The switch is explicit in `FitOptions` and `PotentialAnalysisRequest`, and exposed
as `--second-stage-failed-only-refinement false` to recover native-only behavior.

Only the second-stage iteration engine passes this switch into prepared shape
fitting. First-stage fitting, alpha training, group solvers and direct generic
local/MDPDE calls do not apply refinement. Later group fitting may see different
second-stage models, as it did in the original closed-loop experiment.

## Numerical contract

The private core is shared by production and the existing endpoint experiment.
It refines the native endpoint with the same Powell hybrid method in
`(beta0, log(beta1), log(variance))`, central-difference Jacobian, `xtol=1e-12`
and initial trust-region factor `1`. Candidate equation work, including endpoint
verification, is bounded by 128 evaluations. Acceptance requires a valid Gaussian,
positive finite variance, full weighted rank, positive denominator and fresh
floor-inclusive scaled equations with maximum norm at most `1e-8`.

Every candidate is compared with a continuation of the original fixed-point
endpoint. The reference must reach residual `1e-10` within 10,000 total updates,
including native iterations. Each transformed-coordinate relative difference is
bounded by `1e-6` using `max(1,abs(candidate),abs(reference))`; the weight maximum
absolute difference is at most `1e-6`, and the floor masks must match exactly.
This is numerical branch evidence, not a mathematical proof of branch identity.

Accepted beta, variance, fresh weights and covariance are published together.
Covariance uses the existing formula and must remain finite and positive. Rejection
preserves the original endpoint, weights, covariance, native status and native
iteration diagnostics. Exact-fit/variance boundaries, invalid rank or denominator,
budget exhaustion, reference failure and branch mismatch remain unqualified.
The reference is never a fallback. Successful native solves are not newly
certified to `1e-8`, and a failed endpoint is not promoted merely for already
having a small fresh residual.

## Qualification and observation

`RHBMBetaEstimateResult.status` remains the native status. `Qualification()` derives
`NativeSuccess`, `RefinedSuccess` or `Unqualified` from it and the optional
refinement evidence. Proposal health, active-coordinate qualification, nominal
operators, recovery and final certification consume this effective qualification.
Missing attempted-solve evidence stays unqualified even if an offset-only fallback
retains a previous fitted model.

Audit shape records retain effective `status` for existing readers and add
`native_status`, `qualification` and optional `refinement`. Evidence includes
residuals, branch differences, floor agreement, candidate equation work and separate
reference updates/stop reason. Audit settings record the switch. These transient
fields require no database schema migration. Public C++ struct additions require
consumers of the shared library to rebuild.

Testing-only wrappers keep legacy/failed-only/fresh-residual comparisons, using the
same numerical core. An explicitly active experiment owns the policy to prevent
double refinement. Production has no environment policy or shared mutable switch.
The historical capture runner explicitly disables production refinement when the
binary supports the switch. A passive testing capture records each actual solve
without changing its parameters or running extra equations.

## Validation

Fresh Debug/SYSTEM/OpenMP builds and hash-verified fold-168 inputs were used.
Reproducible output is retained under `build/failed-only/`; historical
`build/endpoint-refinement/` evidence is not overwritten.

- Offline replay: all 337 captured shapes plus four finite-variance fixtures passed
  equation and branch checks. The frozen budget remains 128; equation work is
  exactly the historical 6,522 evaluations (minimum/median/p99/maximum 18/19/28/35),
  with 7,306 reference updates and maximum 287 per case.
- Native-only CLI: all 168 complete persisted records, quality metrics and final
  certificate match the historical baseline exactly. There are 32 attempts, 31
  accepted iterations, best iteration 28, and an unqualified final operator.
  All 5,544 shape calls perform zero refinement work.
- Default-enabled CLI: audit ON j4/v4, ON j4/v3 and OFF j1/v3 exactly match the
  historical failed-only saved records, quality and final certificate. Each run
  has 3,192 shape calls and 20 accepted refinements, with 519 candidate equations
  (maximum 33) and 4,439 reference updates. Successful native calls do not perform
  the 3,172 fresh diagnostic evaluations used by the earlier experiment.
- Recovery accepts four steps. Its residual mean square follows
  `211.3306154800788 → 172.42359421900008 → 164.4348236898434 →
  160.62693493957704 → 158.7712681370311`. The final eight trials are still rejected
  by the best-objective bound.

The default-enabled run stops after 14 attempts / 13 accepted iterations and saves
best iteration 6, with final complete/qualified both true and no final polish.
Its saved-state nominal p99 is
`[0.0062227567253929665, 0.008322470379508595, 0.0010021782426235788]`.
The stop reason remains `recovery-failed`; this does not satisfy outer convergence.
Reference trajectories and saved-state residuals are reported separately.

Regression tests cover native-success bypass, real failed endpoints, covariance
consistency, native/effective status separation, missing evidence, offset failures,
explicit prepared-fit opt-in, generic-fit isolation, insufficient data, exact-fit
and near-zero noise, rank/denominator rejection, the known other root, floor-mask
mismatch, candidate/reference budgets and reference stagnation. C++ command,
Python binding, audit/parser, replay and install-consumer checks accompany the
numerical runs. All 21 CTest groups passed with the captured-data replay enabled.
The audit-OFF build passed 30 focused solver/production tests and both numerical
probe tests; all 12 probe records exactly match the audit-ON build. Both enabled
and disabled fold-168 runs also match their historical references in the
`BUILD_TESTING=OFF` build. Its library contains production refinement and no test
instrumentation symbols. `lint_repo` and the installed C++ consumer smoke passed.
[Compact validation evidence](figures/failed-only-refinement/results.json)
records the final artifacts and their checks.

The production comparison runner is `tests/integration/failed_only_refinement.py`.
It accepts `--executable`, `--model`, `--map`, `--output`, `--reference`, optional
`--legacy`, `--jobs`, `--verbosity` and `--capture` (testing builds only). Enabled
runs omit the new CLI option so they exercise the default. It clears experimental
policy variables, verifies input/source/binary provenance, compares complete saved
records and final certificates, and validates every recorded refinement.

The [original experiment](endpoint-refinement-experiment.md) remains historical
evidence. Improved inner qualification does not change the outer `1e-4` threshold,
offset IRLS or best-objective bound, and does not establish outer convergence.
