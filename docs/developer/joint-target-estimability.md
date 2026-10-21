# Joint target point validation

Current target evidence uses `target-point-validation-v1` and describes selected
targets only. It reports whether each target has an available point estimate
from a qualified Pure Exclusion Joint component, whether its parameters are
valid, and what target-only local conditioning and parameter stability
diagnostics observed. Excluded non-target atoms do not have Joint target point
evidence or fitted FirstStage halo parameters in the production path.

## Current point contract

The required statuses are component convergence, parameter validity, and target
point validation. A target can receive a formal point estimate only when its
Pure Exclusion component is qualified. Poor local conditioning, an unstable
bounded probe, or an unassessed probe is carried as point-quality evidence; it
does not silently exclude the target from the component. Halo reconstruction
availability is a separate downstream status.

Per-atom quality records contain the component status, parameter validity,
local-conditioning status and block IDs, parameter-stability status and block
IDs, point-quality status, and a reason. Conditioning and stability summarize
FixedNeighbor blocks at the accepted endpoint. They do not establish global
identifiability, estimate marginal parameter uncertainty, or replace a global
multi-start analysis.

`JointStageAdapter` reads these records for targets. In the production Pure
Exclusion path, `PostJointLocalValidation` still applies the existing local
qualification before publishing post-Joint local fits; it does not infer halo
estimates from the excluded contributors. Local MDPDE covariance and group
information remain based on those target-local fits; Joint does not produce a
global covariance. Historical FixedHalo results retain their separate
FirstStage-source records in persistence; they are not produced by current
estimation.

## Historical compatibility

The result decoder preserves the exact historical meaning of saved rank,
convergence, FixedHalo, and FullABC-pruning evidence. Loading does not
reinterpret retired rank records as local conditioning or stability, recompute
point evidence, or enable retired estimators for new fits. Current schemas and
payload compatibility are documented in the
[Joint command contract](commands/potential_analysis.md); the current fitting
and reconstruction behavior is described in the
[Joint runtime guide](joint-component-runtime.md).

Joint result JSON schemas 3 through 18 remain readable; earlier Joint JSON is
rejected without migration or database writes. New Pure Exclusion results use
schema 17, or schema 18 when they include geometry-local Halo reconstruction
provenance. Both preserve original/fitting row mappings and target-only point
evidence.
Schemas 3 through 13 preserve their saved evidence and convergence statuses.
Schema 10's rank records, schema 11's FullABC-pruning records, and schemas 13
and 14's FixedHalo provenance/snapshots remain decode-only. Schema 17's
GroupKey median snapshots remain readable with their original contract. Schema 16 records
global convergence criteria; schema 15 remains readable with its historical
strict-criteria defaults. Schemas 8 through 11 keep correction-free numerical
qualification separate from endpoint verification. Legacy reference agreement
is not reinterpreted as conditional KKT or Independent Replay.

SQLite remains schema v20. Estimator-neutral stage format 4 stores
`PostJointLocalResult`, paired peeling, compact local information, posterior
evidence, and sample geometry. Stage format 3 is read and converted using each
available local result's fixed Joint charge profile. Legacy sample BLOBs with
unavailable geometry and historical Group evidence remain unavailable rather
than being assigned new meanings.
