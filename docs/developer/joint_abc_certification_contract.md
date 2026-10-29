# Joint solver numerical contract

The production Joint estimator solves a bounded nonlinear least-squares problem
for Target amplitudes, signed charge coefficients, and widths. Amplitudes are
nonnegative; widths remain positive and bounded by the calibrated global upper
limit. The objective, coefficient units, and public workflow are defined in the
[Joint runtime guide](joint-component-runtime.md) and
[Joint command contract](commands/potential_analysis.md).

## Local solve and width search

At fixed widths, the local A/C problem uses Eigen SparseQR and an active set for
`A >= 0`; `C` remains signed. Width search uses bounded Active-set LM. Each
trial remains feasible, and an accepted local update must pass its trust checks
and the roundoff-safe global objective acceptance. A block visit accepts at
most one trusted update. Failed blocks preserve the last accepted component
state and are retried on the next sweep.

The full derivative includes the coefficient correction for nonzero residuals.
An active upper-width face uses the one-sided KKT sign, and the active set may
release a boundary coordinate when the objective points inward. Local KKT and
trust checks do not substitute for global component qualification.

## Component qualification

The default global criteria are practical-v2: projected A/C KKT at most
`1e-7` and bound-aware projected log-width KKT at most `1e-9`. Stationarity
confirmation requires at least two completed sweeps and a current sweep with no
block failure. Global feasibility, the full sweep replay audit, recovery, and
candidate rollback remain part of search. Endpoint Certification and
Independent Replay check the accepted state. Sweep-budget exhaustion is not
convergence; an unqualified candidate is not published as a Target estimate.

Local conditioning and parameter stability are read-only point-quality
assessments. They do not alter accepted parameters or revoke a qualified global
endpoint. Perturbation-based Parameter Stability applies deliberate local
parameter perturbations and records the resulting parameter, prediction, and
objective sensitivity. Current Joint Result JSON uses schema v21 and accepts
only `strict-v2` and `practical-v2` convergence identifiers.

## Focused regression coverage

The maintained tests use small deterministic inputs. `ActiveSetBoundedLM_test`
covers interior and upper-bound solutions, multiplier signs, release,
trust-region damping recovery from rank deficiency, and invalid problems.
`Numerics_test` checks derivatives,
active-face rebuilding, feasible width replay, and projected KKT behavior.
`FixedNeighborBlockCoordinate_test` covers retries, accepted-state safety,
complete-sweep confirmation, endpoint qualification, and budget exhaustion.
