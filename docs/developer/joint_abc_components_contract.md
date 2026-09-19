# Exact-component equivalence contract, version 1

Retained mathematical and historical scope contract. The multi-version/matrix
procedures below describe the frozen baseline, not the current daily workflow.
Use the [runtime guide](joint-component-runtime.md) for current commands and the
[evidence index](joint-component-evidence.md) for archived evidence.

This testing-only experiment extends the frozen `d9f22026` certification results.
It has three sequential gates: unchanged monolithic evaluation/certification and
structural partitioning; same-state numerical equivalence; independent Guarded
search followed by global endpoint certification. A failed gate is retained as
evidence and does not authorize weakening a tolerance or relabelling a failure.

## Observation and component identity

The original nine immutable voxel/contributor snapshots are the inputs. Each
source snapshot, voxel table, contributor table, starting-width vector and source
fit is hashed. The observation domain, memberships, squared distances, double and
float32 observations are never reconstructed during optimization. Four starts
(`first-stage`, `narrower`, `wider`, `mixed`) come directly from the frozen source
fits; remapping an atom ID cannot change the mixed-start parity rule.

Connectivity is structural. A union-find joins the first contributor of each
voxel with every other contributor on that row. No edge can be removed because
of a zero amplitude/charge, underflow, weak sensitivity or a coupling threshold.
Components are ordered by their minimum stable atom identity. Local-to-parent
and parent-to-local atom/row mappings preserve the original accumulation order.
Rows with no contributor have a constant objective `sum(y*y)/2`; atoms with no
support remain explicit unobserved components. A constant row occurs exactly
once in the assembled loss and never in a component's search residual.

The source census has one connected component per dataset. The main dataset has
168 atoms, 139,551 rows and 407,237 memberships. These are regression controls,
not evidence that the original problem can be split into smaller solves.

Seven additional datasets are disjoint unions of frozen snapshots:

| Dataset | Source blocks |
|---|---|
| regular-two | baseline, weak-1e-2 |
| regular-three | baseline, baseline, near-0.10 |
| regular-near | baseline, near-0.02 |
| regular-weak | baseline, weak-1e-4 |
| regular-active | baseline, active-a |
| regular-zero | baseline, zero-signal |
| regular-duplicate | baseline, duplicate |

Each source occurrence has its own identity namespace. A disjoint union is an
algebraic observation fixture, not a newly generated physical map. Original
coordinates are retained as provenance, not used to infer cross-block edges.
Small geometry tests separately verify sphere support and shared-voxel edges.

## Evaluation, rank and audit context

The objective is equal-weight least squares, alpha zero, A nonnegative, B
positive through log-B coordinates, C signed. Support, kernels and the full
variable-projection Jacobian are unchanged. No variance fit, regularization,
new B bounds, step clipping or robust loss is introduced.

The parent observation normalization is `max(1, norm(y))`. Double and 50/100-digit
calculations compute this norm independently in their own arithmetic, always
from the full parent observations. Component evaluation, KKT, replay, width
gradients and audits must not replace it with a norm of the selected rows.

The rank policy identifies each matrix family separately: column-normalized
design, normalized free design, projected-width matrix, column-normalized
projected-width matrix, and the full profile Jacobian. At a shared global state,
the SVD absolute threshold is

```
epsilon * max(global_row_count, global_family_column_count) * global_sigma_max
```

The maximum singular value is calculated at the current state. A union of sorted
component spectra provides the block matrix's global spectrum; the same absolute
threshold applies to every block. A sum of ranks obtained with different local
thresholds is not global rank evidence. SparseQR retains its separate absolute
pivot threshold based on the largest normalized column norm. Active-set release
and iteration policies retain factors 128 and 20, respectively, and same-state
solves use the shared global driver and its complete coefficient vector.

An explicit audit registry replaces numerical kernels' case-name/atom-count
decisions. Historical registrations preserve trial diagnostics, expansion after
an unverified legacy derivative, registered precision targets and boundary atoms.
Composite registrations inherit the source requirements even when the combined
atom count exceeds twelve. Shared-state directions are generated once, saved,
and restricted to components without renormalization. Final certification uses
all-ones, alternating and the weakest projected-width direction of the global
endpoint being certified. Search-time diagnostics may use the saved initial
reference directions, but cannot substitute that evidence for the endpoint audit.
Each monolithic and assembled endpoint has its own saved global direction set;
supplemental component audits restrict the assembled endpoint's set. Their
legacy two-step checks are refreshed before Richardson/high-precision checks.
If no global width spectrum exists, the saved unit-axis fallback is labelled
as missing weak-direction evidence and cannot establish a global certificate. Repeated singular values
require spectrum/subspace evidence rather than equality of arbitrary SVD bases.

Precision targets retain the complete `0.01 * 2^-k`, k=0..16 ladder. Richardson
interval selection uses successive finite differences and inner-solve noise,
never truth or analytic agreement. The uncertainty gate is 1e-7, both derivative
agreement gates are 1e-6, and 50/100-digit agreement is scaled 1e-20. KKT and
local correction remain 1e-10; width stationarity remains 1e-12. Boundary evidence
does not turn a weakly active or cross-face state into a regular certificate.

## Same-state acceptance

First compare raw predictions at supplied A/C/log-B, independently of whether
profiling succeeds. Then compare constrained profiles at the same B and valid
fixed-face derivatives. Both linear backends share blocking/release decisions;
the partitioned backend factors the disconnected free-face blocks independently.
No normal equations or numerical edge pruning are permitted.

| Quantity | Required agreement |
|---|---|
| IDs, memberships, stored squared distances | Exact |
| Raw prediction, each row | 2e-12 + 2e-13 times the larger absolute prediction |
| Regular profiled A/C | Maximum scaled difference at most 1e-10 |
| Objective | Absolute difference divided by parent scale squared at most 1e-12 |
| Width gradient | Absolute/relative tolerances 1e-13 / 2e-9 |
| Projected-width and full Jacobian | Relative Frobenius difference at most 1e-8, denominator floor 1e-12 |
| Singular values | Difference divided by global maximum at most 1e-10; equal rank at the shared threshold |

The scaled parameter difference is `max(abs(a-b)/(1+max(abs(a),abs(b))))`.
Original RSS and all numerical discrepancies are retained even when normalized
gates pass. Rank-deficient, infeasible, untrusted or changing-face derivatives
are explicitly unavailable or failed, never represented by a fabricated zero
Jacobian. Threshold-sensitive rank disagreements block the corresponding claim.

The derivative layer uses one canonical profiled A/C/B state after comparing
the independently solved coefficients. Its declared face probes use the two
legacy step sizes, 1e-4 and 5e-5, in shared normalized all-ones and alternating
directions (or explicitly supplied global directions). Every primary/reference
probe must preserve the canonical face. A face disagreement or transition is
recorded as unavailable for fixed-face derivatives, while raw predictions,
objective, profile values and design spectrum remain independently testable.
The failed exploratory run is retained. No spectrum tolerance is changed.
Single-component snapshots dispatch directly to the existing monolithic backend;
their identity comparisons are labelled and are not evidence of a genuine split.

## Independent search and assembled certificates

Guarded retains its original LM diagonal metric: factor 0.1, ftol 1e-14, xtol
1e-12, gtol 1e-12. Each component has 200 profile evaluations and 100 accepted
updates; reference and audit work is counted separately. These are per-component
budgets, not a claim of equal aggregate cost to a 200-evaluation monolithic fit.

Search uses parent observation normalization and explicitly recorded local
linear rank/active-set rules. It reads no other component's evolving parameters,
results or cache. This search policy is distinct from same-state/global
certification, whose rank dimensions and spectral thresholds remain global.
Search completion, component qualification and assembled qualification are
separate fields. After reassembly, evaluate and audit the actual global state.
The active-set response norm is also inherited from the parent. Search rank
dimensions remain local. Global certification uses the actual assembled A/C,
including its KKT, residual, gradient, Jacobian and rank evidence. A separately
reprofiled control must agree within 1e-10; it never replaces that state.

Compare same-named starts in the two precisions. Never assemble a solution by
choosing a different best start for each component. Regular paired endpoints
must both qualify, with scaled A/B/C and maximum log-B differences at most 1e-8,
prediction infinity difference divided by parent scale at most 1e-8, and scaled
objective difference at most 1e-12. Iteration paths need not agree. A discrepancy
triggers same-state evaluation at both endpoints to distinguish decomposition
errors from search differences.

Failed components retain their last trusted state and reason without cancelling
other components. If a component has no usable state, the assembled prediction
and objective are unavailable and a partial-row mask identifies available data.
Missing blocks must never be replaced by zeros. Zero-signal and duplicate
generating-model evidence remains separate from numerical endpoint rank; truth
does not enter search, interval selection or qualification.

An isolated rerun consumes the frozen component input, parent context, matching
start and audit registration. It must reproduce that component's scientific
record without requiring sibling executions.

## Reproduction and accounting

The Python runner provides regression, run, audit, summarize, compare and
rerun-component operations. All outputs require fresh run directories. Each
complete experiment includes nine originals and seven composites, both
precisions and all four starts, and is repeated independently in another fresh
directory. No audit cache crosses run boundaries. Numeric kernels are single
threaded; independent process orchestration is recorded separately.

Retain census, state parity, endpoint parity, certificates, every failure,
phase costs, input/source/build hashes and a scientific artifact inventory.
Repeat comparison excludes elapsed time, process high-water RSS and execution
paths. Search reference time and precision time are nested costs and must not
be added twice. Process peak RSS is not an isolated component allocation peak.

Production integration, a memory-layout rewrite and robustness experiments are
outside this contract. Testing-only linkage is verified in a production build.

High-precision global references may use the proven block structure, with a
single full-parent normalization in each precision. Global boundary objectives
include the unchanged blocks' path/fixed-face/profile losses and constant rows;
global KKT is the maximum across blocks. Small controls compare this calculation
against dense global QR, including a constrained release and nonzero constants.
Exact-input high-precision results may be memoized within one audit case. Keys
include frozen support, observations, eta, active face (or full beta for boundary
scans), directions and the immutable parent observations. Caches are cleared
before every case and isolated rerun. Reference computations and reuses are
recorded separately; no cache survives a fresh run or replaces a global test.

## Runner operations

The testing executable is `build/observation-matching/on/bin/mdpde_experiment`.
The Python entry point is `tests/integration/joint_abc_components.py`.

```text
regression --baseline FROZEN_CERTIFICATION --executable EXECUTABLE --output GATE
same-state --baseline FROZEN_CERTIFICATION --executable EXECUTABLE --output STATE_GATE
run --baseline FROZEN_CERTIFICATION --executable EXECUTABLE --regression-gate GATE/regression.json --output RUN_A
run --baseline FROZEN_CERTIFICATION --executable EXECUTABLE --regression-gate GATE/regression.json --output RUN_B
audit --run RUN_A --executable EXECUTABLE --output FRESH_AUDITS
summarize --run RUN_A [--audits FRESH_AUDITS]
compare --left RUN_A --right RUN_B [--left-audits AUDIT_A --right-audits AUDIT_B] --output COMPARISON.json
rerun-component --run RUN_A --executable EXECUTABLE --dataset DATASET --case CASE --component STABLE_ID [--audits FRESH_AUDITS] --output RERUN
```

`run` validates the prior regression's snapshot hashes, repeats all same-state
controls, performs paired searches, and runs fresh audits. `audit` independently
recomputes audits without repeating searches. An audit wrapper update is allowed
only when every numerical search kernel matches the run's saved source hashes;
the fresh audit executable, sources and provenance are recorded separately.
No previous audit or boundary cache is reused. Fresh audits use three independent
case processes, each with one numerical thread and its own empty cache. `rerun-component` reads only frozen
input, its frozen search context, final global audit context and paired source start;
the original component result is read afterwards solely for comparison.
Inputs, implementation sources and hashes are copied into each run. Reports are
derived into `summary/` (under the selected audit root for a fresh audit); complete fits, trials, spectra, probes, audits, boundary
scans, failures and contexts remain in `datasets/` and `audits/`.
