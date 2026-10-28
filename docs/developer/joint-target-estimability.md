# Joint target point validation

The current target evidence contract is `target-point-validation-v1`. It
describes selected Target points from the production target-only Pure
Halo-Support Exclusion fit.

Each point records component convergence, parameter validity, Local
Conditioning, Parameter Stability, point quality, and a reason. Conditioning
and stability also retain the associated FixedNeighbor block identities.
Unavailable checks and `NotAssessed` statuses are valid diagnostic results.

Target point validation does not assert global identifiability, estimate
marginal parameter uncertainty, or change the global solver's convergence
criteria. Global convergence, endpoint verification, and Independent Replay
remain separate evidence. A Halo has no free Joint point; geometry-local
kNN-4 Halo Median Reconstruction is a separate downstream result.

The Joint Result writer emits schema v20 and the reader accepts only schema
v20. The payload retains target point evidence, convergence, component and
row mappings, candidate diagnostics, Local Conditioning, Parameter Stability,
and optional reconstruction evidence. It permits unqualified results,
missing states, nullable diagnostics, and unavailable reconstruction.

SQLite schema v20 stores Joint schema v20 and Stage Result JSON version 4.
Persistence performs no migrations or historical decoding. See the
[Joint runtime guide](joint-component-runtime.md) and
[data I/O architecture](architecture/dataobject-io-architecture.md).
