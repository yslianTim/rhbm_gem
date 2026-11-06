# Joint stage integration

This guide describes data ownership from the production Joint fit through
Original-domain peeling, Post-Joint Local MDPDE, and Group Fitting. Solver
details are in the [Joint runtime guide](joint-component-runtime.md).

Joint Components is the default potential estimator; `two-stage` remains an
explicit choice. Potential fitting excludes hydrogen from targets and eligible
contributors, including direct C++ workflow entrypoints. FirstStage sampling
uses the Sampler's deterministic Fibonacci default. Post-Joint observations
come directly from original map voxels. These input policies do not change
Joint problem construction or its Pure Halo-Support Exclusion model.

## Result ownership

The Joint result stores the target-only Pure Exclusion state, component
identities, target point validation, convergence evidence, candidate
diagnostics, local conditioning, and optional
Halo-reconstruction evidence.

The estimator-neutral Second-stage records expose available Joint Target
A/tau/C values, target/Halo roles, component and run identity, and convergence
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
3. Gauss–Seidel BCD updates component blocks sequentially. Each block uses
   Variable Projection with Constrained A/C Active-set QR and Bounded Active-set
   LM to produce the target-only state. Global convergence, Endpoint
   Certification, Independent Replay, target qualification and Local
   Conditioning are recorded.
4. Qualified Target A/tau/C points provide geometry-local kNN-4 Halo medians.
   Each reconstructed Halo records its qualified donor Target IDs. An
   unavailable reconstruction remains a separate status from Joint numerical
   convergence.
5. Original-domain PostFitPeeling uses available Target points and reconstructed
   Halo points to subtract other contributors from original voxel responses.
   Candidate observations come from OriginalInput support with r² <= 4 Å²,
   including the 2 Å boundary, in increasing Original-row order. The Fitting
   mask never removes reconstruction observations. Contributions use their
   actual 2.5 Å kernel support, exclude hydrogen and exclude the Target itself.
   Each sparse support contribution is evaluated once and shared across Targets;
   direct non-self sums avoid cancellation against the own Gaussian signal.
   Missing qualified contributor state makes the affected response unavailable.
6. The Target's own Joint C/tau charge is evaluated directly at each candidate
   voxel and subtracted from its peeled response. No Post-Joint step interpolates
   map responses or predictions, or depends on Fibonacci sample count/selection.
   Qualified Local fitting uses the positive finite corrected observations in
   0 <= r <= 1 Å. Alpha training, OLS/MDPDE, QR qualification, working covariance
   and information all consume the same selected design. The candidate domain
   remains 2 Å; the production fitting range is a separate Post-Joint contract.
   Group alpha training and
   Group Fitting consume qualified local estimates and their compact
   MDPDE-derived information.

Only selected targets enter local fitting, group alpha training, or group
statistics. Local results never rewrite Joint parameters. Group posterior
Gaussian values remain separate from each member's Joint charge profile.

## Invalidation and availability

Derived records retain their Joint source and run identity. Grid observations
store Original row index, global grid index, squared distance, physical voxel
position and raw map response, alongside optional peeled response and its
unavailable reason. Local failures use the same Original row/global grid pair.
Finite nonpositive corrected responses are valid corrections, then excluded by
the Gaussian estimator's positive-response policy. Local counts distinguish all
2 Å candidates, corrected responses, positive responses and radius-retained
solver rows.

A Joint endpoint change invalidates the run's peeling, Local and Group results
and snapshot. Replacing or clearing the snapshot invalidates downstream grid
evidence. Peeling changes invalidate Local and Group results; Local changes
invalidate the affected Group. Unrelated FirstStage raw sample edits retain
independent grid results. Atom position edits clear analysis and spatial caches.

A live Joint analysis binds to its source map's observation revision. Replacing,
clearing or normalizing that map invalidates the bound model analysis before
its next access, including model copies. Map copies have independent revisions.
Map geometry is immutable. Destroying a map does not invalidate stored evidence;
loading a persisted model restores historical evidence without a live map binding.
A new workflow builds fresh Original observations and a new binding. Contexts
borrow immutable inputs and must not outlive or mutate them.

Initialization failure, unavailable width calibration, failed search,
unqualified endpoint, missing component state, unavailable reconstruction, or
incomplete downstream results are legal current Joint outcomes. Downstream
consumers preserve unavailable states and reasons; they do not treat absent
data as zero.

## Persistence

SQLite schema v21 stores the current model structure, Joint Result JSON schema
24, and Stage Result JSON version 5. Joint JSON is target-only Pure Exclusion.
Joint results require the complete BCD solver provenance v6 described in the
[Joint runtime guide](joint-component-runtime.md), including pre-search failures.
Stage JSON retains current FirstStage and SecondStage estimates, FirstStage raw
samples, original-grid peeling observations, Local results, Group posterior,
fixed-charge profiles, summaries and provenance. Readers validate ordered unique
Original/global grid identities, response/reason pairing, snapshot row mapping,
voxel positions against recorded map geometry, source/run identity, correction
counts and fixed charge. Grid peeling requires its Joint snapshot; available
Halo-subtracted responses require reconstruction provenance. Old raw-index
Post-Joint records are rejected, with no compatibility reader.

The numerical and statistical range assessment is recorded in the
[production decision record](post-joint-grid-validation.md). Information and
Group covariance remain conditional working quantities; voxel count alone does not
establish independent observations or calibrated posterior coverage.

Only an empty database is initialized. An existing database must have the
exact v21 schema. Stage and Joint readers accept only their current versions.
There are no migrations or historical decoders; unsupported inputs are
rejected without changing the database.

See the [data I/O architecture](architecture/dataobject-io-architecture.md)
for schema checks and the [Joint runtime guide](joint-component-runtime.md)
for solver qualification.
