# Joint component estimator runtime

This guide describes the production target-only Pure Halo-Support Exclusion
estimator. The persisted parameterization is
`target-only-pure-halo-support-exclusion-v1`.

## Production flow

`RunJointComponentWorkflow` builds the Pure Exclusion problem and obtains
qualified FirstStage target widths. A width calibration is available only with
at least 20 qualified target donors and 95% target coverage. The sole global
width bound is `tau_max = 4 × Q0.95(qualified FirstStage target widths)`.
Nonfinite or nonpositive donor widths, median-fallback donors, initialization
failures, and overflowed bounds do not produce a usable calibration. Initial
widths are projected into this finite bound; unavailable calibration prevents
Joint search from starting.

Targets are the only free Joint components. Non-target Halo support is used
only for the Pure Halo-Support Exclusion row selection, with the production
support radius of 2.5 Å. The model is
`f_i(r) = A_i G(r; tau_i) + C_i K(r; tau_i)`, with `A_i ≥ 0`,
`0 < tau_i ≤ tau_max`, and signed `C_i`. C multiplies the charge kernel; it is not
a constant background offset.

The atomic parameter vector is `theta_i = (A_i, tau_i, C_i)`. Tau is the
physical width in angstrom; the optimization coordinate is
`eta_i = log(tau_i)`, with `eta_max = log(tau_max)`. The kernels are

```text
G(r; tau) = (2*pi*tau^2)^(-3/2) * exp(-r^2/(2*tau^2))
K(r; tau) = erf(r/(sqrt(2)*tau))/r
K(0; tau) = sqrt(2/pi)/tau
dG/deta = G(r; tau) * (r^2/tau^2 - 3)
dK/deta = -sqrt(2/pi)/tau * exp(-r^2/(2*tau^2))
```

The runtime retains its production support cutoff and near-origin handling.
`BasisValues::d_gaussian_d_log_tau` and `d_charge_d_log_tau` store these
log-coordinate derivatives. `JointState::tau` stores physical values,
`log_tau` stores eta, and `eta_gradient` stores `dF/deta`, where
`dF/deta_i = tau_i * dF/dtau_i`. The LM step is an eta step. Calibration uses
Hyndman–Fan Type 7 quantiles of qualified FirstStage Target tau values.

The Joint Component Estimator uses Gauss–Seidel Block Coordinate Descent
(BCD) to visit component blocks in forward order. Each visit forms the
Effective Response using the latest accepted state of the other Targets and
holds those Targets fixed throughout the local optimization. Only an accepted
candidate updates the shared state before the next block visit.

```text
Joint Component Estimator
  Gauss–Seidel BCD
    Variable Projection
      Constrained A/C Active-set QR
      Bounded Active-set LM
```

Variable Projection re-solves the constrained A/C problem at every Width
trial. Bounded Active-set LM solves the resulting profile log-width problem;
the Profile Jacobian retains its nonzero-residual derivative correction.

Production uses Gauss–Seidel BCD blocks with core size 12, forward serial
Gauss-Seidel order, with a nominal budget of 200 sweeps and a hard cap of
300. A block visit may accept at most
one trusted update. At fixed widths, the local A/C problem uses column
normalization, Eigen SparseQR, a nonnegative-amplitude active set, strict
`1e-10` local KKT certification, and the existing local QR iterative
refinement when needed. C remains signed. Width updates use bounded log-width
Active-set LM with a bounded trust-region subproblem, damping adjustment, and
the weak-width safeguard. The full profile derivative includes its
nonzero-residual correction and cancellation protection. Compact SVD supports
rank and weak-direction checks; it is not an alternate A/C solver.

The adaptive budget keeps the last 20 normalized projected-KKT values:
`q = max(A/C KKT / 1e-7, projected log-width KKT / 1e-9)`. At the nominal
boundary, all values must be finite and nonnegative, with no pending block
recovery. The best q in the recent ten sweeps must be at most 10 and improve
at least 5% over the previous ten-sweep best, which must exceed 1. An extension
continues the same accepted state up to 300 sweeps. A qualified endpoint stops
earlier. Denied extensions report `block-sweep-budget`; an exhausted extension
reports `block-hard-sweep-budget`. Pending recovery retains its distinct
`block-recovery-sweep-budget` reason. There is no public budget selector.

## Constraint status, working sets and free variables

For fixed widths `tau_i = exp(eta_i)`, the A/C model is exactly linear:
`min_beta 1/2 ||X(eta) beta - y_eff||^2`, subject to `A_i >= 0` and signed,
unconstrained `C_i`. This A/C problem is not a Taylor linearization.
Its constrained minimizer is `beta_star(eta)`, and the normalized profile
objective is `phi(eta) = ||X(eta) beta_star(eta) - y_eff||^2 / (2*s^2)`,
where `s` is the parent observation scale.

| Term | Meaning in the runtime |
| --- | --- |
| Active constraint | An inequality attains equality: `A_i = 0` or `eta_i = eta_max`. Bound diagnostics retain their existing floating-point tolerance. |
| Inactive constraint | An inequality has strict slack: `A_i > 0` or `eta_i < eta_max`. |
| Working set | Constraints temporarily held as equalities in an active-set subproblem; solver iteration state, not a geometric constraint census. |
| Free variables | Variables allowed in the current QR solve. A/C includes every C and every A outside the amplitude working set. |
| Frozen direction | A stationary-null Width direction excluded as a numerical safeguard, independently of the inequality working set. |

`SolveConstrainedAcLeastSquares` uses `ConstrainedAcWorkspace::is_free_variable`,
ordered `free_column_indices` and the normalized `free_design`. Blocking fixes
only an amplitude at zero. Releasing it immediately restores free-variable
eligibility even while its value is still zero. A zero amplitude can therefore
be an active constraint outside the working set. The solution certificate's
`active_amplitude_atom_indices` records geometric zero amplitudes, while
`ac_free_design_rank` records the QR solve's free-design rank.

Bounded Active-set LM uses `working_bound_mask`, `free_step_indices` and
`newly_blocked_indices`. `stationary_null_direction_mask` remains separate:
neither working-bound coordinates nor frozen directions enter the free step
QR. A frozen direction can itself lie at the upper bound without entering the
working set. Returned `active_bound_indices` and `active_bound_multipliers`
describe the attained working bounds of a valid step; they need not enumerate
every geometric active constraint. `working_set_releases` counts releases,
and `active_set_iterations` retains its algorithmic meaning.

Local Conditioning and the Profile Jacobian derive their free A/C design from
the endpoint's constraint face. Those columns are not a saved copy of the
iterative solver's working set. Their `free_ac_columns`, `free_ac_conditioning`
and `free_design` names retain that meaning.

The progress label `ac-kkt=` reports the accepted sweep's global A/C projected
KKT. Block failure diagnostics report the initial local A/C profile through
`ac_free_design_rank`, `ac_linear_solve_count`, `ac_working_set_releases` and
`ac_scaled_projected_kkt`. Width LM status and Width KKT remain separate fields.
Logger levels, throttling and event frequency are unchanged.

Progress reports `max_tau` for the maximum Target tau, `eta-grad` for the
raw eta gradient norm, and `eta-projected-kkt` for the bound-projected eta
KKT. Block failure diagnostics use `tau_upper_bound` and `eta_step_inf`.

## Qualification

The global search preserves objective monotonicity and practical-v2
convergence: projected A/C KKT ≤ `1e-7` and bound-aware projected Width KKT ≤
`1e-9`. Width stationarity uses the projected gradient on the feasible
upper-bound face, with feasibility and complementarity checks; an ordinary
gradient cannot veto a valid active-bound solution. Confirmation requires at least two completed sweeps and no block
failure on the confirming sweep. Global feasibility, the full sweep replay
audit, recovery, and candidate rollback remain part of search. Endpoint
Certification and Independent Replay verify the accepted state. Component
assembly consistency, Local Conditioning, target point
quality, and resource preflight remain explicit evidence; local LM termination
alone does not qualify a component. Failed initialization, unavailable width
calibration, incomplete search, unqualified endpoint, missing component state,
and nullable diagnostics remain representable saved outcomes.

Local Conditioning assesses each block's rank, spectra and active-face
behavior without changing its numerical classification or the accepted state.
Target point validation uses `target-point-validation-v2`: failed/unavailable
component convergence or parameters produce Unavailable; invalid conditioning
mapping, IllConditioned, Singular, ActiveFace or Unavailable conditioning
produce Warning; NotAssessed conditioning produces NotAssessed. Passed
numerical checks, valid parameters and WellConditioned evidence produce
Available. ActiveFace remains conditional evidence, not IllConditioned.

Point quality does not establish global identifiability or precision. Removing
a diagnostic can change quality without changing an estimate. Target runtime
qualification still requires numerical convergence, valid selected parameters,
complete row availability and assembly mapping/consistency. See the
[current point contract](joint-target-estimability.md).

## Halo reconstruction and downstream work

Qualified Joint Target points are the donors for geometry-local Halo Median
Reconstruction. Each Halo uses up to its four nearest qualified Target points
by Euclidean distance in Å; equal-distance ordering uses Atom ID. A result
records the donor IDs and the geometry-local kNN-4 contract. No retained Halo
or unavailable donors is a valid result state.

The reconstructed Halo points support Original-domain PostFitPeeling. The
workflow then runs Post-Joint Local MDPDE for qualified selected targets and
Group Fitting from those local results. Peeling collects
Original support rows at r² <= 4 Å², uses raw map voxel values and direct 2.5 Å
neighbor kernels, then directly subtracts the Target's own C/tau charge.
Post-Joint Local fitting uses the positive corrected grid observations within
1 Å, with one shared dataset for alpha, estimation and information. FirstStage
Fibonacci sampling is unchanged. These downstream stages retain their own
availability and qualification evidence; see the
[stage integration guide](joint-stage-integration.md).

Halo medians retain the `(A, tau, C)` order:
`theta_h = (median A_i, median tau_i, median C_i)` over each Halo's qualified
geometry-local donors. Charge correction remains
`y_i_Gaussian(r) = y_i_peel(r) - C_hat_i * K(r; tau_hat_i)`.
The generic `GaussianModel3D::GetWidth()` API and the Halo reconstruction
`width` field describe physical Width and retain their names.

Memory preflight and memory-bounded evaluation remain part of problem
construction, Gauss–Seidel BCD search, endpoint verification, local conditioning,
Halo reconstruction, and peeling.

## Persistence contract

| Payload | Accepted format |
| --- | --- |
| SQLite database | Schema v21 |
| Stage Result JSON | Version 5 |
| Joint Result JSON | Schema v24 |

SQLite opens a new empty database by initializing schema v21 and opens an
existing database only when its exact schema is v21. Stage JSON accepts only
version 5. Joint JSON accepts only schema 24 and the `strict-v2` or
`practical-v2` convergence identifiers. It stores target-only Pure Exclusion
results. No format performs migrations or compatibility decoding. Unsupported
inputs are rejected without repair or overwrite.

Joint state objects contain exactly `ac`, `tau`, `log_tau`, `eta_gradient`
and `objective`. Initialization uses `tau`, with `original_tau` and
`used_tau` in each captured atom record. Candidate diagnostics use
`candidate_eta_gradient_inf_norm` and `global_eta_gradient_threshold`;
the latter is the bound-projected eta stationarity threshold.
Units metadata uses `joint-kernel-map-units-v2` and `tau: angstrom`.
Joint CSV parameter columns are `A,tau,C`, in that order. Schema 24 accepts
only these current keys, including when the Joint payload is stored in a
SQLite v21 database. Stage JSON v5 and native Gaussian Width fields keep
their independent contracts.

Current documents can store qualified and unqualified results, nullable
diagnostics, missing component states, unavailable width calibration, and
unavailable Halo reconstruction. The writer does not select a schema from
result contents.

Every saved result, including an initialization or calibration failure before
search begins, requires the complete solver provenance object:

```json
{
  "search_method": "gauss-seidel-bcd",
  "contract_version": "joint-solver-provenance-v6",
  "sparse_backend": "EIGEN",
  "bcd_core_atoms": 12,
  "bcd_policy_version": "gauss-seidel-bcd-production-v1",
  "bcd_maximum_sweeps": 300,
  "bcd_order": "forward",
  "bcd_local_search": "active-set-lm"
}
```

Missing/null solver provenance, unknown fields, prior method/policy versions,
and mixed contracts are rejected. The provenance revision changes the naming
and persisted structure; numerical policy remains practical-v2 with the same
adaptive 200/300 sweep budget. Pre-search failures carry the configured solver
identity without claiming that search ran or that a width bound was available.

## Source index

- [ConstrainedAcLeastSquares.cpp](https://github.com/yslianTim/rhbm_gem/blob/develop/src/core/detail/joint_component/ConstrainedAcLeastSquares.cpp)
  implements Constrained A/C Active-set QR and its working-set transitions.
- [ActiveSetBoundedLM.cpp](https://github.com/yslianTim/rhbm_gem/blob/develop/src/core/detail/joint_component/ActiveSetBoundedLM.cpp)
  implements the bounded LM trust-region subproblem and damping search.
- [GaussSeidelBcd.cpp](https://github.com/yslianTim/rhbm_gem/blob/develop/src/core/detail/joint_component/GaussSeidelBcd.cpp)
  implements block search, acceptance and endpoint certification.
- [BcdPolicy.hpp](https://github.com/yslianTim/rhbm_gem/blob/develop/src/core/detail/joint_component/BcdPolicy.hpp)
  defines `BcdExecutionPolicy` for component orchestration; `BcdSearchPolicy`
  in `GaussSeidelBcd.hpp` owns the internal search settings.
- [ComponentSweepBudget.hpp](https://github.com/yslianTim/rhbm_gem/blob/develop/src/core/detail/joint_component/ComponentSweepBudget.hpp)
  owns the unchanged adaptive sweep decisions.
- Refactor validation (`bcd-refactor-validation.md`, absent from this checkout) records the baseline and
  numerical/contract checks.
- Tau naming validation (`tau-naming-validation.md`, absent from this checkout) records the symbol inventory,
  numerical equivalence and logger/persistence checks.

The [data I/O architecture](architecture/dataobject-io-architecture.md)
documents the database and JSON boundaries. The
[stage integration guide](joint-stage-integration.md) describes downstream
ownership and invalidation.
