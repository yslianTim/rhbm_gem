# Joint target point validation

Current target evidence uses `target-point-validation-v1` and describes selected
targets only. It reports whether each target has an available point estimate
from a converged Joint component, whether its parameters are valid, and what
target-only local conditioning and parameter stability diagnostics observed.
Supported halos use validated FirstStage `(A, B, 0)` snapshots and do not have
Joint target point evidence.

## Current point contract

The required statuses are component convergence, parameter validity, and target
point validation. A target can receive a converged Joint point estimate when
its component has a valid state under the fixed FirstStage halo contribution.
Poor local conditioning, an unstable bounded probe, or an unassessed probe is
carried as point-quality evidence; it does not silently exclude the target from
the component. This qualification does not certify FirstStage halo accuracy.

Per-atom quality records contain the component status, parameter validity,
local-conditioning status and block IDs, parameter-stability status and block
IDs, point-quality status, and a reason. Conditioning and stability summarize
FixedNeighbor blocks at the accepted endpoint. They do not establish global
identifiability, estimate marginal parameter uncertainty, or replace a global
multi-start analysis.

`JointStageAdapter` reads these records for targets and separately publishes
fixed halo values with FirstStage source provenance. `PostJointLocalValidation`
still applies the existing local qualification before publishing post-Joint
local fits. Local MDPDE covariance and group information remain based on those
target-local fits; Joint does not produce a global covariance.

## Historical target-estimability records

`target-estimability-v1` through `target-estimability-v4` retain the exact
meaning saved by their original schema. In particular, historical rank pass,
rank deficiency, threshold sensitivity, and unavailable rank remain historical
rank evidence. Loading does not reinterpret them as local conditioning or
parameter stability, and it does not recompute point evidence or relabel old
convergence states.

Historical schemas with global rank records remain readable through the result
decoder. New point-validation evidence has no global rank records, null-space
leakage, target Jacobian, or rank-derived qualification. The exact JSON schema
and compatibility range are documented in the
[Joint command contract](commands/potential_analysis.md).
