# Target estimability and 6Z6U experiment results

This implements the target contract on top of direct A/C and the existing
Guarded LM / LegacyCompact search. It does not promote OperatorPcg or SpqrBounds,
add iterative A/C, alter other experimental tools, or establish large-case performance.
The complete 6Z6U solve has **not run**; its budget remains unset.

## Case preparation and truth results

The experiment checked input identities, selection, grid/support arithmetic,
normalization, component layout and actual First B0 values. Its manifests, truth
files, receipts and reproduction tools have been removed; only the results are
retained here. Truth was used for validation only, never by the estimator.

Production preparation on 2026-09-24 independently reproduced 1,559 targets,
2,192 contributors, 262,801 rows and 866,216 memberships. Structural component
sizes are 2,167 and 25 atoms. There are 27 singleton halos on 26 profiled rows,
2,165 FullABC atoms, 262,775 informative rows and 17 two-row FullABC halos.
No active-face classification is claimed for this preparation-only run.

Truth uses complete sidecar-v2 atom identities, A=atomic number,
B=effective_gaussian_width and C=charge_used. B contains 0.40, 0.45 and 0.50
Angstrom. A and C follow the supplied-map kernel units (U*Angstrom^3 and
U*Angstrom). Independent scalar forward replay over the frozen memberships
passed the float32-ULP plus existing double-forward allowance: maximum absolute
error 2.382330182015835e-7, RMSE 1.8633138857787364e-8, normalized objective
3.185058819046246e-16. The quantization residual is retained, not subtracted.

## Target evidence at the returned state

`TargetRuntimeConvergence()` is separate from the unchanged all-parameter
`RuntimeConvergence()`. Inputs without selection metadata get `NotRun`, not an
inferred target interpretation. Component-local and assembled-global evidence
use the saved coefficients and log widths. They do not reprofile to replace the
state under examination. Existing feasibility, inner/KKT, full width stationarity,
scalar replay and assembly checks remain required.

At an interior state, build the complete observation Jacobian in `(A,C,log B)`
order after existing singleton profiling. Normalize each column by its saved
2-norm (a zero column stays zero with scale 1), reduce with tiled QR, and SVD the
compact factor. Retain the original row count in the shared rank context.
For `Z=[Z_T Z_H]`, testing whether target selector rows annihilate `null(Z)` is
algebraically equivalent to testing target rank after `(I-Z_H Z_H^+) Z_T`.
This preserves halo coupling without constructing an observation-space projector
or a normal matrix. A selector leakage above 1e-10 fails identifiability. Singular
values in `[0.5*tau,2*tau]` are conservatively threshold-sensitive/Unavailable.
The saved evidence includes scales, spectra, thresholds, rank and per-atom null
leakage. This is a **local** interior certificate, not a global uniqueness proof.

The undamped full profile Jacobian includes the residual correction. Its
pseudoinverse local step is mapped into target A/C/log-B; A/C corrections are
scaled by `1+abs(coefficient)`. Target correction must be at most 1e-10, and mapped
profile-null directions must vanish to the same tolerance. Rank-sensitive profile
steps remain unavailable even when the observation Jacobian is well resolved.
All width gradients still satisfy the existing 1e-12 stationarity condition.
Any FullABC A=0 in the component makes its target geometry unavailable.

This implementation retains dense **parameter-space compact factors and SVDs**;
it does not claim scalable 6Z6U rank/endpoint memory. A separately budgeted run is
required before asserting large-case feasibility or success.

## Output, uncertainty and persistence

JSON v5 adds `target_evidence` and `target_runtime_convergence` to components and
the assembled result. Readers still accept v3/v4 without inventing target evidence.
SQLite uses its existing JSON payload storage. CSV appends
`TargetRuntimeConvergence,ParameterIdentifiability`. Unidentified halo A/B/C are
blank in CSV and unavailable as atom-wise Second point/covariance. Internal
representatives remain in JSON for full prediction and peeling reconstruction.
Targets retain returned point states with their separate convergence status;
a point alone is not a passed certificate.

For an identifiable interior target, the full coupled Jacobian's identifiable
subspace supplies its covariance block in `(A,C,log B)`. Variance uses informative
rows minus effective rank; singleton nuisance rows/rank cancel together.
Boundary states, insufficient degrees of freedom, zero residual variance and
incomplete evidence retain explicit unavailable reasons. Covariance availability
is reported separately from point/truth success. Peeling still requires all
interpolation-stencil voxels to lie inside the saved observation domain.

The full-case success criterion remains all 1,559 target A/B/C scaled truth errors
at most 1e-3, target convergence, objective quality, complete prediction and
persistence/export validation. Covariance availability is reported separately.
A partial result or resource stop does not establish full-case success.

## Spectrum regression and validation

The reproduced `active-a / first-stage-double` endpoint was checked against the
independent reference for projected, normalized projected and full profile
spectra, column norms, active face, gradient and local correction. The experiment
endpoint file and its replay entry point have been removed. Near cancellation,
separately rounded triplet products were magnified
by column normalization. Only the cancellation fallback now uses reference-order
multiply-add accumulation, out-of-place QR/RHS and column-major projection. Normal
fast paths and acceptance tolerances are unchanged. Both EIGEN and SPQR yield
2.6179363200699017e-15 spectrum error at this endpoint (previously
1.2990450737934428e-9 on SPQR), below the original 1e-10 limit.

Small tests cover genuine halo-only nonuniqueness, the same geometry selected as a
target, nonzero residuals, singleton/two-row mixtures, independent augmented-SVD
covariance, equivalent-representative peeling, multiple components, unobserved
targets, A=0, rank-sensitive evidence, saved-state mismatch, JSON v3/v4/v5,
SQLite/CSV, CLI and public APIs. The historical dense covariance test explicitly
exercises its legacy contract; new target tests use valid noisy endpoints.

The following results were recorded before removal of experiment artifacts and
tools; they are historical validation results, not a current test inventory:

| Verification | Result |
|---|---|
| SPQR selected CTest groups | 15/15 passed; joint core 127 passed |
| EIGEN selected CTest groups | 13/13 passed; joint core 123 passed, 4 SPQR-only skipped |
| Acceptance tool unit tests | 6 passed |
| Existing validation runner | 18 passed after allowing process-memory inspection |
| Installed C++ consumer | Built, executed, JSON/SQLite round-trip passed |
| Python binding and examples | Both passed (included in SPQR groups) |

The broader documentation-link group still fails on three archives already absent
at the baseline: `joint-sparse-acceptance/measurements.json.gz`,
`joint-sparse-acceptance/command-exports.tar.gz` and
`joint-operator-search/receipts.tar.gz`. No archive was fabricated and no check was
relaxed. This is not an all-repository-tests-passed claim.

One diagnostic multi-component EIGEN start with a different halo representative
placed the local profile step at its numerical rank threshold; it correctly
remained unavailable. The passing multi-component control uses a well-resolved
start, and separate tests explicitly require unavailable threshold-sensitive
results. Functional small-case completion must not be presented as 6Z6U
full-case success.
