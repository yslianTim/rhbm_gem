# Joint target point validation

The current contract is `target-point-validation-v2`. Each selected Target
point records component convergence, parameter validity, Local Conditioning,
its block identities, point quality and a reason.

Rules are evaluated in this order:

| Condition | Point quality |
|---|---|
| Component convergence failed, unavailable or not run | Unavailable |
| Parameters invalid or unavailable | Unavailable |
| Invalid target/conditioning block mapping | Warning |
| Conditioning Singular, IllConditioned, ActiveFace or Unavailable | Warning |
| Conditioning NotAssessed, including no blocks | NotAssessed |
| Passed numerical checks, valid parameters, WellConditioned | Available |

ActiveFace retains the existing conservative warning for possible active-set
switching. It is not reclassified as IllConditioned. Local Conditioning's
numerical calculations are unchanged. A missing membership in a nonempty block
collection, duplicate identity/membership or out-of-range atom mapping remains
explicitly invalid. Empty diagnostics represent NotAssessed.

This contract deliberately removes parameter perturbation sensitivity from
point quality. Warning or NotAssessed may become Available when the remaining
checks pass. This change does not imply improved estimation accuracy or global
identifiability. IllConditioned and Singular still warn. Numerical convergence
and target parameter validity remain mandatory, independently of quality.

TargetRuntimeConvergence still merges component convergence, parameter validity,
selected target checks, assembly mapping/consistency and available rows.
IsQualifiedJointAnalysisResult still requires valid initialization, completed
search and Passed runtime/target convergence. Conditioning remains diagnostic;
it does not revoke or replace endpoint qualification. Halo reconstruction,
Original-domain peeling and Group Fitting keep their separate contracts.

Joint JSON schema 22 writes this contract and omits the retired diagnostic
members completely. The reader accepts only schema 22 and rejects unknown
component/point fields and noncurrent target contracts. Schema 21 documents
must be read with their original Git version. No legacy reader is retained.
SQLite schema 20 and Stage Result JSON version 4 are unchanged: SQLite stores
the versioned Joint JSON payload and performs no migration. The C++ public
struct layout also changes; downstream C++ consumers must rebuild.

See the [runtime guide](joint-component-runtime.md),
[data I/O architecture](architecture/dataobject-io-architecture.md) and
[cleanup acceptance record](joint-diagnostics-cleanup.md).
