# Joint stage integration

This guide describes data ownership from the production Joint fit through
Original-domain peeling, Post-Joint Local MDPDE, and Group Fitting. Solver
details are in the [Joint runtime guide](joint-component-runtime.md).

Joint Components is the default potential estimator; `two-stage` remains an
explicit choice. Potential fitting excludes hydrogen from targets and eligible
contributors, including direct C++ workflow entrypoints. Sampling uses the
Sampler's deterministic Fibonacci default. These input policies do not change
Joint problem construction or its Pure Halo-Support Exclusion model.

## Result ownership

The Joint result stores the target-only Pure Exclusion state, component
identities, target point validation, convergence evidence, candidate
diagnostics, local conditioning, parameter stability, and optional
Halo-reconstruction evidence.

The estimator-neutral Second-stage records expose available Joint Target
A/B/C values, target/Halo roles, component and run identity, and convergence
status to existing model consumers. A missing or unqualified target remains
unavailable with its reason.

Post-Joint Local results are separate target-keyed records. They contain local
alpha, OLS and MDPDE results, solver diagnostics, sample counts, failure
reasons, and compact information used by Group Fitting. They do not modify the
Joint state or its objective. Halo atoms contribute to reconstruction and
peeling; they do not receive target-local or group results.

## Data flow

1. The workflow builds a Pure Halo-Support Exclusion problem with support
   radius 2.5 Å. Non-target Halo support selects rows for exclusion; Halo atoms
   are not free Joint parameters.
2. Qualified FirstStage Target widths initialize the Joint fit and provide the
   global width calibration. A result with unavailable calibration keeps its
   diagnostics and does not proceed as a successful Joint endpoint.
3. FixedNeighbor and its bounded Active-set LM / constrained A/C SparseQR
   solve produce the target-only state. Global convergence, endpoint
   verification, Independent Replay, target qualification, Local Conditioning,
   and Parameter Stability are recorded.
4. Qualified Target A/B/C points provide geometry-local kNN-4 Halo medians.
   Each reconstructed Halo records its qualified donor Target IDs. An
   unavailable reconstruction remains a separate status from Joint numerical
   convergence.
5. Original-domain PostFitPeeling uses available Target points and reconstructed
   Halo points to subtract other contributors from each target's raw samples.
   It keeps one result per raw sample, including unavailable reasons.
6. Post-Joint Local MDPDE fits each qualified target using the peeled samples
   and the target's fitted Coulomb charge profile. Group alpha training and
   Group Fitting consume qualified local estimates and their compact
   MDPDE-derived information.

Only selected targets enter local fitting, group alpha training, or group
statistics. Local results never rewrite Joint parameters. Group posterior
Gaussian values remain separate from each member's Joint charge profile.

## Invalidation and availability

Derived records retain their Joint source and run identity. A Joint endpoint
change invalidates the run's peeling, Post-Joint Local results, group results,
and Joint snapshot. A raw sample change invalidates peeling and its dependent
local/group results while retaining a still-valid Joint snapshot. A peeling
change invalidates local and group results. A local-result change invalidates
the affected group result.

Initialization failure, unavailable width calibration, failed search,
unqualified endpoint, missing component state, unavailable reconstruction, or
incomplete downstream results are legal current Joint outcomes. Downstream
consumers preserve unavailable states and reasons; they do not treat absent
data as zero.

## Persistence

SQLite schema v20 stores the current model structure, Joint Result JSON schema
21, and Stage Result JSON version 4. Joint JSON is target-only Pure Exclusion.
Stage JSON retains current FirstStage and SecondStage estimates, raw and
peeling samples, Post-Joint Local results, group posterior, fixed-charge
profiles, summaries, and source provenance.

Only an empty database is initialized. An existing database must have the
exact v20 schema. Stage and Joint readers accept only their current versions.
There are no migrations or historical decoders; unsupported inputs are
rejected without changing the database.

See the [data I/O architecture](architecture/dataobject-io-architecture.md)
for schema checks and the [Joint runtime guide](joint-component-runtime.md)
for solver qualification.
