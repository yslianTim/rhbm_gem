# Joint-ABC certification contract, version 2

Retained mathematical and historical scope contract. The multi-version/matrix
procedures below describe the frozen baseline, not the current daily workflow.
Use the [runtime guide](joint-component-runtime.md) for current commands and the
[evidence index](joint-component-evidence.md) for archived evidence.

This testing-only experiment extends the frozen `ecf55f45` coverage results. It does not replace production fitting, change the statistical objective, or reinterpret the earlier uniform-width fixture. The original coverage reports remain immutable.

## Inputs and objective

The nine datasets, 72 precision/start combinations, first-stage initialization, four main-data fixed-B controls and truth definitions are inherited from the [fixture catalog](../../tests/fixtures/joint_component/catalog.json) (the original v1 hash/width contract is preserved in `tests/fixtures/joint_component/simulation-contract.json`). Each dataset runs initialization once, using its float32 map and the original deterministic sampling, local alpha training and fixed-offset first-stage fitting. All versions use the resulting MDPDE B vector. A/C are profiled again. Truth is used only for simulation scoring and separately labelled generating-model diagnostics.

The objective is equal-weight least squares with alpha zero, nonnegative A, positive B and signed C. No regularization, additional noise, width bounds or step clipping is introduced. Search retains 200 profile trial evaluations and 100 accepted updates, factor 0.1, ftol 1e-14, xtol 1e-12 and gtol 1e-12. Reference solves and endpoint audits are charged separately.

## Simulation manifest v2

Production simulation emits `schema_version: 2`. Each atom includes `effective_gaussian_width` and `effective_charge_width`; inapplicable quantities are null. `ElectricPotential::GetEffectiveWidths` supplies both single-Gaussian evaluation and metadata. It applies O × 0.8, N × 0.9 and other elements × 1 to the requested blurring width. The five-Gaussian model has no single effective Gaussian width; its charge width follows the existing minimum-width rule when the charge kernel is enabled. The user model has neither of these widths.

`kernel.version` identifies the existing potential model formula. `kernel.width_policy` is `element-scaled-v1` with the three factors for `single_gaus`, or `model-specific-v1` otherwise. The kernel retains near-zero distance, charge cutoff and minimum charge width fields. There is no global `effective_charge_width` in v2. Model/input hashes, source/configuration/build hashes, preparation order and accumulation order remain recorded.

`support.version` is `sphere-fma-v1`:

- coordinates: `fma(index, spacing, origin)` on each axis;
- squared distance: `fma(dz, dz, fma(dy, dy, dx*dx))`;
- membership: the strict comparison `squared_distance <= cutoff*cutoff`;
- charge-kernel cutoff and outer support cutoff remain separate policies.

No radius tolerance is permitted. Full-map float32 equivalence is checked before fitting the main dataset. Main-data support must remain 139,551 rows and 407,237 memberships; near-0.10 must retain its frozen support. The generator's scalar kernel still receives the square root of the contracted squared distance.

The new Python reader accepts v1 only for the exact model/map/manifest hash triples in the frozen coverage or uniform-width baseline. It does not infer a width policy for unknown v1 files. The coverage C++ adapter likewise requires the frozen coverage v1 manifest hash, or a validated v2 contract. Historical experiment entry points and their existing Python reader reject v2. A historical uniform-width map must be regenerated with its historical generator, never by rewriting its truth to fit the current generator.

## Immutable observation snapshot

`voxels.csv` preserves voxel IDs, explicit coordinates, double/float32 observations, and zero/negative values. `contributors.csv` stores `(row, atom, square)` in row/atom order. `snapshot.json` contains CSR row offsets, dimensions, policy version and hashes of both tables. Every fit records the snapshot hash. Altered tables or snapshot metadata are rejected. The experiment checks equality of the legacy coordinates with explicit FMA coordinates before creating the snapshot.

Search, numerical replay, derivative/high-precision audit and boundary scans consume the saved support and squared distances. Python separately verifies the snapshot geometry and replays predictions, residuals, KKT and width gradients directly from the saved memberships. No solver iteration updates the observation domain.

## Search variants

`legacy` uses the instrumented adapter with Eigen's original arithmetic and acceptance behavior. `guarded` retains the original diagonal metric. `guarded-log` sets D = I in log-B space. Both guarded variants validate the initial state and each candidate before acceptance using independent inner reference solves and scalar prediction/gradient replay. The copied testing adapter is derived from Eigen 5.0.1 `NonLinearOptimization/LevenbergMarquardt.h`, with the original MPL-2.0 notice retained. System Eigen is unchanged.

Validity gates use the original tolerances: primary/reference KKT 1e-10, scaled coefficient difference 1e-10, prediction absolute/relative tolerances 2e-12/2e-13, KKT replay difference 1e-13, and gradient absolute/relative tolerances 1e-13/2e-9. Rank, conditioning and cancellation are also retained as evidence. Identifiability is not a hidden constraint on B.

An invalid trial leaves accepted parameters, residual and Jacobian unchanged. The LM radius is multiplied by 1/4 and damping by 4 before recomputing a candidate. Exhaustion, an unrepresentable step or stagnation after invalid candidates returns the last trusted state with an explicit failure, not a convergence claim. Every trial retains accepted/candidate parameters, step, diagonal, radius, damping, predicted/actual reduction and ratio. Audit values do not alter Legacy calculations.

## Certificates and representative solutions

The raw `joint_qualified`, `qualification_failure` and `qualification_checks` are preserved. The new certificate is a separate artifact, not a rewrite of the legacy report. It records all failures, unavailable checks and evidence paths, numerical values and thresholds. Unavailable never means passed.

Numerical optimality covers feasibility, both inner KKT checks, coefficient agreement, both B gradients and local correction. Identification includes design, projected-width and profile-Jacobian rank evidence. Derivative evidence is labelled verified, transition-unverified, resolution-unverified or unavailable. An exactly active A with dual gradient no larger than the existing 1e-10 KKT resolution is conservatively labelled weakly active/nonregular. A small positive A is never rounded to zero.

The regular certificate requires every applicable check, independent endpoint replay, trusted endpoint evaluation and completed search. When an independently solved 50/100-digit fixed-face reference is feasible and agrees, its local correction is assessed against the unchanged 1e-10 gate, with the source explicitly recorded alongside the legacy double value. Its KKT and gradient gates also remain applicable. Derivative success alone cannot pass local correction.

Legacy representatives remain the minimum-RSS legacy-qualified branches. New representatives are the minimum-RSS regular-qualified branches, or null if none qualify. Truth error is never used in either selection or in the certificate. Generating-model rank/sensitivity evidence is stored separately, so a full numerical endpoint rank cannot erase the zero-signal or duplicate model's structural nonidentifiability.

## Derivative resolution audit

The two legacy step sizes remain in each fit. All registered precision targets (weak-1e-4, near-0.02, active-a and baseline first-stage controls), including targets already verified by the legacy check, use the full step ladder. Other unverified small endpoints also retain a ladder as diagnostic evidence. The three frozen directions are all-positive, alternating and the weakest projected-width direction. Samples use h = 0.01 × 2^-k, k = 0…16, retaining both primary/reference solves, active sets, KKT, coefficient disagreement, sensitivity, finite-difference errors and cancellation.

For adjacent central differences F_k, form R_k = (4 F_(k+1) - F_k)/3. A candidate interval uses R_k and R_(k+1). Its uncertainty is their norm difference plus the propagated primary/reference residual disagreement, normalized by max(1e-12, the two extrapolation norms). Only unchanged-face intervals are eligible. Select the least estimated uncertainty; exact ties retain the larger step. Neither truth nor the analytic Jacobian discrepancy enters selection.

The selected interval must have estimated relative error at most 1e-7, and both extrapolations must agree with the analytic derivative within 1e-6. The independent high-precision derivative must also pass and the 50/100-digit solutions must agree to scaled 1e-20. If no trustworthy interval exists, the derivative remains unverified.

High precision uses Boost decimal floats, independent Gaussian/charge kernel evaluation, dense Householder QR, differentiated fixed-face KKT and local correction. It does not promote a rounded double coefficient solution into a reference answer. Exact repeated squared distances may reuse kernel values within one column. An identical endpoint may reuse an already completed audit within the same run; the original artifact is explicitly named. A fresh run recomputes the evidence.

## Boundary and negative controls

Active-A scans retain exact active sets, scaled slack, dual and complementarity. At generating truth and actual endpoints they evaluate the positive path `delta eta > 0`, `delta A = 4*pi*B^2*C*delta eta`, `delta C = 0`, for atoms 2, 6 and 10 and steps 2^-k, k = 4…24. Both fixed-face and re-constrained profiles are evaluated at 50/100 digits. Prediction changes, objectives, slopes and the actual near-zero branch's first-order contribution are recorded.

To avoid repeatedly factoring the unchanged tall matrix, an independent dense QR of `[X0, y]` retains its constant subspace. At each scan point the two changed columns are orthogonally reduced, followed by small dense QR solves for each active face. This is equivalent least squares; it does not form normal equations. A dedicated test checks high-precision agreement, feasibility, compensation and the constrained objective.

Low residual or a one-sided match does not establish regular parameter recovery. Zero-signal and duplicate remain negative controls, with independent generating-model rank and width sensitivity diagnostics.

## Reproduction and accounting

Build the testing executable, then run:

```sh
python3 tests/integration/joint_abc_certification.py run \
  --executable BUILD/bin/mdpde_experiment \
  --model MODEL --map MAP --manifest MANIFEST --output FRESH_DIRECTORY
python3 tests/integration/joint_abc_certification.py summarize FRESH_DIRECTORY
python3 tests/integration/joint_abc_certification.py audit FRESH_DIRECTORY \
  --executable BUILD/bin/mdpde_experiment
python3 tests/integration/joint_abc_certification.py compare \
  --left FIRST_DIRECTORY --right SECOND_DIRECTORY --output COMPARISON.json
python3 tests/integration/joint_abc_certification.py plots FIRST_DIRECTORY --output FIGURES
```

The C++ entry points are `joint-abc-certification MODEL MAP MANIFEST OUTPUT` and `joint-abc-certification-audit OUTPUT`. The runner verifies source/build/input stability and requires fresh run directories. All 216 searches, including failures, are included in exact scientific JSON/CSV comparison after excluding time, process peak RSS and execution paths. Identical-endpoint audit reuse follows a sorted deterministic order.

Initialization, search, search reference calls, legacy endpoint audit, new derivative/precision audit and boundary scans record elapsed time. Memory fields are process high-water RSS, not isolated allocations attributable solely to a phase; reports must preserve that distinction. Production builds contain manifest v2 but no experiment solver or multiprecision audit.

References: [Boost decimal floating-point types](https://www.boost.org/latest/libs/multiprecision/doc/html/boost_multiprecision/tut/floats/cpp_dec_float.html), [Ceres derivative and extrapolation discussion](https://ceres-solver.googlesource.com/ceres-solver/+/987d3b6b370ab65205a97cf8377d8848483458a4/docs/source/derivatives.rst).
