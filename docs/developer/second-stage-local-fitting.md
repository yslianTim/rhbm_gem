# Second-stage local fitting

## Scope

`RunSecondStageLocalFitting(ModelObject&, const FitOptions&)` refines the
Gaussian estimates of the selected atoms after the first-stage local and group
fits. It accounts for overlapping responses from neighboring selected atoms
and unselected model atoms while updating each selected atom's amplitude,
width, and offset. Unselected atoms contribute background responses but are
never added to the optimizer state.

The stage keeps candidate states in memory and writes one validated terminal
state to `ModelObject`. Individual outer iterations do not partially update the
stored atom estimates.

## Model context and initialization

The fitting context contains, for each selected atom:

- its raw local-potential sampling entries and trained `alpha_r`;
- every model atom that contributes to each sample.

Neighbor candidates are searched within `kNeighborAtomSearchRange`. A neighbor
contributes to a sample only when its distance from that sample does not exceed
`kNeighborContributionDistanceMax`. Unselected candidates are deduplicated
across targets and retained only when they affect at least one sample. When
`FitOptions::exclude_hydrogen` is true, hydrogen atoms are removed from this
contributor set; other selection exclusions do not remove background
contributors.

Before the second stage starts, the first group fit consumes the selected raw
sampling entries and supplies the per-atom posterior and group prior used by
seed selection. The first-stage local MDPDE and OLS models are not passed to
the second-stage seed selector.

The initial Gaussian seed is rebuilt for every selected atom. The current local
MDPDE and local OLS estimates are not seed sources. The first valid source is
selected in this order:

1. group posterior;
2. group prior;
3. same-group parameter median;
4. global parameter median.

The median sources are bootstrapped in two passes. Each atom first contributes
at most one direct model, preferring its valid group posterior over its valid
group prior. The component-wise group and global parameter medians are computed
only from those direct models; atoms filled from a median are not fed back into
either median. The selected seed replaces the complete local MDPDE model,
including its offset and uncertainty. Direct posterior and prior seeds retain
their uncertainty, while median seeds use zero uncertainty.

If a valid seed cannot be obtained for every selected atom, the stage exits
without changing the stored estimates, peeling sampling entries, or group
results.

Each effective unselected neighbor then receives a transient seed from the
same-group median or, when that group has no selected direct seed, the global
median. Unselected seeds never feed back into either median pool. Selected and
unselected seed summaries are logged separately. If an unselected seed cannot
be built, the stage exits with `no-valid-unselected-neighbor-seed` before
changing stored results.

The valid initial state is used to build a weighted coupling topology. The
topology records atoms that jointly affect objective samples. Each iteration
partitions the active portion of this topology into deterministic cluster keys,
so independent clusters can be evaluated and accepted separately.

The initial weighted topology retains the fixed minimum edge weight `0.05`.
After accepted iterations, the stage adaptively rebuilds the topology from the
latest validated atom models when either the maximum transformed-coordinate
drift from the last topology reference state reaches `0.10`, or three accepted
iterations have elapsed since the last rebuild. Rejected attempts and
trust-radius retries do not advance this interval. An iteration that already
changes the objective domain through terminal isolation skips the adaptive
rebuild without clearing the accumulated drift or interval.

Adaptive rebuilds use edge hysteresis. A previously absent edge must have
weight at least `0.06` to enter the graph, while an edge present in the previous
post-residue-cutoff adjacency remains until its weight falls below `0.04`.
The existing ten-residue cluster limit is then applied in descending current
edge-weight order. Binary fallback remains the conservative response to an
invalid Jacobian.

## Iteration flow

Each outer attempt performs the following sequence:

1. Build the active atom list by excluding terminal-fallback atoms.
2. Partition the active coupling topology and reconcile the per-cluster
   objective and trust-region states.
3. Jointly estimate one shared offset per represented group within each cluster
   using robust IRLS and the fixed `kJointOffsetRidgeRatio`.
4. Build component-wise group-median models from the post-solve snapshot. For
   each active atom, subtract its selected neighbors and all effective
   unselected contributors from the observed sample responses.
5. Refit the atom's local Gaussian with its trained `alpha_r`, using its
   group-median model as the fixed offset model. These refits form the raw
   fixed-point state.
6. Limit each cluster's raw proposal to its trust region and score the resulting
   endpoint candidate. If the endpoint fails the objective guard, backtrack in
   the three transformed coordinates within the same outer attempt.
7. For a stationarity-eligible cluster without suspicious atoms, attempt one
   joint amplitude/width/offset polish. Keep the polish only when it strictly
   improves the base candidate on the same objective scale.
8. Build the accepted-induced interaction graph from shared boundary samples.
   Revalidate each connected component. Every eligible component attempts one
   joint correction over the direct boundary contributors: only their
   transformed amplitude and width may change, while one component-local offset
   column is shared by every closure atom with the same group key. A valid
   endpoint remains the fallback and is replaced only when the correction fits
   every member trust radius and strictly improves that endpoint. For an invalid
   endpoint, the correction must strictly improve the previous component
   objective; an unavailable or rejected correction falls through to the
   unchanged common-factor component backtracking. A failed component rolls back
   only its members, so unrelated components and remote singleton clusters remain
   eligible for commit. Then build maximal eligible boundary components from the
   safe accepted state and the best finite objective-rejected proposal retained
   for each rejected cluster. Components containing at least one such proposal
   receive the same endpoint, routine joint-correction polish, and common-factor
   checks. A failed correction retains a valid rescue endpoint, and a failed
   rescue leaves the safe accepted state unchanged.
9. Run the unchanged global previous/best audit on the complete assembled state.
   If tolerance-level deterioration from independent reconciliation units causes
   aggregate rejection, remove non-improving units from worst to best until the
   first passing subset is found. If only strictly improving units remain but the
   historical best gate still fails, roll back the attempt without weakening the
   gate.
10. Update trust radii and persistent-failure state. After an accepted state,
   conditionally rebuild the adaptive topology before updating the global audit
   state and applying the stopping conditions.

The neighbor-adjusted response for atom `i` and one of its samples is:

```text
adjusted response = observed response
                  - sum(selected-neighbor fitted responses)
                  - sum(unselected-contributor responses)
```

For every evaluated selected snapshot, an unselected contributor uses the
latest median of the selected models with its `GroupKey`. If no selected atom
has that key, it keeps its initialization-time global-median seed. The same
resolver is used by the coupling topology, joint-offset system, local refit,
objective guards, audit, and final peeling calculation.

This offset solve and neighbor-adjusted local refit constitute one raw
fixed-point update. These group-median-adjusted entries are temporary inputs to
the current raw proposal; they are not persisted as peeling sampling entries.

The joint-offset parameterization is cluster-local. Atoms with the same group
key share one offset column when they are in the same coupling cluster. The same
group key represented in separate clusters has a separate offset variable in
each solve, preserving independent cluster failure and acceptance behavior.
Joint polish uses the same cluster-local group-sharing rule, while also updating
each atom's amplitude and width.

## Parameter coordinates

Trust-region steps and convergence checks use three transformed coordinates:

```text
log peak height
log width
offset-to-peak ratio
```

The logarithmic coordinates keep amplitude and width positive when a candidate
is decoded. A candidate is invalid when its amplitude or width is not finite
and positive, its offset is not finite, or its transformed coordinates cannot
be decoded to a valid Gaussian model.

## Cluster objective

Every raw sample belongs to exactly one owner cluster: the cluster containing
the sample's target atom. Samples whose distances are inside the inclusive
`distance_min <= distance <= distance_max` interval form the fit-range domain;
all other raw samples form the tail-validation domain. Selection flags and the
sign of the response do not remove samples from either domain.

The initial validated state supplies two independent, fixed robust scales for
each cluster:

```text
fit scale  = max(MAD(fit residual),
                 1e-6 * MAD(fit adjusted response), 1e-12)
tail scale = max(MAD(tail residual),
                 1e-6 * MAD(tail adjusted response), 1e-12)
```

A cluster must have valid fit-range samples. Its tail may be empty, in which
case both tail loss fields are zero and no tail scale is needed. Scales are not
warmed up or updated after candidates are accepted.

For cluster `c`, the fixed-scale objective is:

```text
fit-range loss       = mean Cauchy(fit residual / fit scale, cutoff=1.345)
tail validation loss = mean Cauchy(tail residual / tail scale, cutoff=1.345)
tail penalty         = 0.25 * tail validation loss
offset penalty       = 0.01 * mean offset-plausibility penalty
cluster total        = fit-range loss + tail penalty + offset penalty
```

There is no width-prior term. Group posterior and prior models participate in
seed selection only. The offset-plausibility residual floor uses the owning
cluster's fixed fit scale.

The global objective weights clusters by active atom count:

```text
global objective = sum((cluster atom count / active atom count) * cluster total)
```

Owner assignment makes every sample appear once in this global sum, including
boundary samples. Local scoring still includes every sample affected by the
candidate cluster. It applies that sample's owner scale and exact global
normalization coefficient, so the local candidate-minus-previous difference
matches the corresponding full-global difference when only that cluster
changes.

Candidate scoring uses a provisional copy of the cluster objective state. A
rejected base, polish, backtracking trial, or boundary-component candidate does not
advance the previous or best references. The best objective and maximum
transformed change are retained to break objective ties.

All second-stage audit tolerances use:

```text
tolerance(reference) = absolute tolerance + relative tolerance * abs(reference)
```

Progress and deterioration guards use `1e-8 + 1e-3 * abs(reference)`. Strict
best, tie, polish-improvement, and trust-growth comparisons use
`1e-10 + 1e-8 * abs(reference)`. Candidate comparisons against previous and
best compute separate tolerances from their respective references. The
joint-offset IRLS objective retains its independent tolerance.

## Trust region

For a non-suspicious cluster, previous and raw offsets are reduced to one
physical offset per `GroupKey` by deterministic component medians. Base
proposal trials use factors `1, 1/2, 1/4, ...`; each trial interpolates the
atom-level log-peak and log-width coordinates and the physical group offset
with the same factor:

```text
C_g(t) = C_previous,g + t * (C_raw,g - C_previous,g)
```

The realized atom models are then re-encoded and measured against the trust
radius. At `t = 0`, the shape is the previous shape and the offset is the
previous group median; at `t = 1`, the complete raw shared-offset state is
recovered. If projecting an inconsistent previous group onto its median already
exceeds the radius, the proposal fails before objective evaluation with an
explicit diagnostic reason.

If the endpoint is valid but fails its local objective guard, the same cluster
attempt evaluates factors `1/2, 1/4, 1/8, ...` between the previous and endpoint
states. Local and boundary-component backtracking preserve one physical offset per
`GroupKey`; shapes remain interpolated in transformed coordinates. Search stops
when the largest transformed change is below
`kTransformedChangeTolerance`. The first passing trial is committed
with endpoint uncertainty, records its factor, and does not grow the radius.
Rejected trials do not mutate objective state or polish provenance. A cluster
that exhausts objective backtracking is immediately excluded from trust-radius
shrink. When another, radius-retryable cluster causes an unchanged-state retry,
the exhausted cluster is skipped and its diagnostic retains its actual radius.
It becomes eligible again after another cluster commits a state change.

The polish step is limited by the radius remaining after the accepted base
movement. A rejected polish keeps the base candidate and is not backtracked.
A radius-retryable rejected cluster shrinks its own radius once. A
non-backtracked accepted cluster grows its radius only when the objective
strictly improves and its step is close to the current boundary. Trust-region
updates are isolated by cluster.

Accepted clusters are first connected only when they both affect the same boundary
sample. Each multi-cluster component
first revalidates the factor-`1.0` assembled endpoint against every member's local
criteria and a unique-owner component audit relative to the committed state. An
eligible component then attempts one joint correction. A valid endpoint is kept
unless the correction strictly improves it. When the endpoint fails, an
unavailable or rejected correction evaluates common factors `1/2, 1/4, 1/8, ...`
for that component alone. Exhaustion rolls back only its member clusters;
independent components and remote singleton clusters retain their accepted
endpoints. A component-backtracked state does not grow its members' trust radii,
and polish provenance is retained only for atoms with a material polished endpoint
change.

After this accepted-only pass, every rejected cluster that produced finite
objective values retains its lowest-objective trust-region proposal in memory.
Clusters with suspicious atoms, hard joint-offset failures, offset-only refit
fallbacks, invalid proposals, unavailable objectives, or unchanged-state
exhaustion are excluded. The remaining rejected proposals and safe accepted
clusters form maximal eligible boundary components. A component containing at
least one rejected proposal is re-evaluated as a combined endpoint, then receives
one routine joint-correction attempt. A valid endpoint is replaced only by a
strictly better correction. Common-factor backtracking is used only when the
endpoint is invalid and the correction is unavailable or rejected. Successful
rescue promotes the rejected members without trust-radius growth or shrink. Failed
rescue is transactional: previously accepted members retain their safe state.

Boundary correction eligibility is intentionally broader than local polish
eligibility. A valid guarded refit and a non-hard joint-offset result are enough;
IRLS objective deterioration, IRLS iteration exhaustion, and valid non-success
local estimation statuses may participate. Local joint polish still requires
full stationarity. Hard solver failures, suspicious atoms, and offset-only
fallbacks remain fixed contributors and never enter a correction patch.
Every eligible accepted-only or rescue boundary component receives at most one
correction attempt per outer iteration.

After component reconciliation, the complete assembled state must still pass the
unchanged global previous/best audit. On aggregate failure, independent components
and singleton clusters are scored by their exact global objective delta. Only
non-improving units are removed, from worst delta to best, with the full audit
recomputed after each removal. The first passing subset is committed atomically.
If all remaining units strictly improve the previous objective but cannot satisfy
the historical best gate, the complete remaining attempt is marked exhausted;
objective tolerances are never relaxed.

## Numerical defenses and terminal isolation

- Joint-offset conditioning and column-collinearity guards can increase the
  ridge multiplier for affected atoms.
- Suspicious evaluation has two paths. An offset-only update checks finite
  zero-offset responses, offset magnitude, center sign flip, and radial
  rebound. A post-refit update first requires a valid second-stage model, then
  applies the same guards plus width growth and amplitude-offset compensation.
  A local-refit candidate uses the post-refit path. If that candidate fails,
  the fallback preserves the previous amplitude and width, applies only the
  jointly estimated offset, and uses the offset-only path. Width and
  compensation are intentionally not reevaluated for this fallback.
- The previous suspicious baseline is built in one fit-range scan. It records
  the innermost response, per-radius response medians, distance range, maximum
  absolute response, and residual scale `1.4826 * MAD`. Candidate profiles do
  not calculate a residual MAD. The post-refit candidate and its offset-only
  fallback reuse the same previous baseline.
- The center sign-flip guard treats the smallest sampled radius as the
  innermost response and only rejects a statistically significant
  positive-to-negative change. It estimates the previous atom's fit-range
  residual scale as `1.4826 * MAD`, requires the previous innermost response to
  exceed `3 * scale`, and requires the candidate to be below the negative of
  both that noise threshold and `0.25 * previous innermost response`.
  Near-zero noise crossings and negative-to-positive changes do not trigger
  this guard.
- The radial rebound guard uses the same previous residual noise estimate. A
  radial magnitude must exceed both `1.5 * abs(candidate innermost response)`
  and `max(0.25 * abs(previous innermost response), 3 * scale, 1e-12)`.
  Upward excursions use
  `max(0.20 * abs(previous innermost response), 3 * scale, 1e-12)` and become
  suspicious only after more than one excursion.
- Sign flip and rebound require a trustworthy previous radial shape. Width
  growth remains active for valid post-refit models and uses the available
  fit-range distance span without depending on radial monotonicity.
  Amplitude-offset compensation likewise does not require a trustworthy
  radial shape; it uses the previous innermost response when available and the
  previous center signal as its reference. Its offset response delta is
  evaluated exactly as
  `candidate.offset * candidate.OffsetBasisAtDistance(0) -
  previous.offset * previous.OffsetBasisAtDistance(0)`.
- An offset-only suspicious seed rolls back only atoms in the same coupling
  cluster that share its `GroupKey`, matching the cluster-local shared offset
  parameterization. It no longer propagates through a depth-limited atom
  overlap graph. If both a post-refit candidate and its fallback fail, the
  complete acceptance cluster is still rolled back atomically. On the next
  attempt, affected active atoms receive
  `kSuspiciousJointOffsetRidgeMultiplier`.
- A failed or suspicious local refit falls back to the previous Gaussian
  parameters with the newly estimated offset when that offset-only update
  passes finite-response, offset-magnitude, sign-flip, and radial-rebound
  guards. The owning cluster is then stationarity-ineligible for that attempt.
  If the fallback fails one of those guards, the affected post-refit cluster is
  rolled back.
- Repeated stable hard solver failures or suspicious rollbacks become terminal
  cluster fallbacks. Terminal atoms retain their previous validated states and
  no longer participate in later solves.
- Terminal isolation removes only the affected cluster, allowing independent
  active clusters to continue fitting.
- Rejected-cluster debug output distinguishes failures before objective
  evaluation from objective rejection. Pre-objective failures report their
  proposal reason, radius, available step norm, and `objective =
  not-evaluated`; `objective-unavailable` is reserved for an objective that was
  actually attempted but could not be calculated.

When terminal isolation changes the active partition, the implementation uses
the current validated state to rebuild the new clusters' fit and tail scales
and all atom-count normalizations. Per-cluster previous/best references and the
global best-audit baseline are reset to that state; objective values from the
old and new domains are never compared.

An adaptive topology rebuild always becomes the next hysteresis reference,
even when its active partition is unchanged. In that case the objective scales,
best audit, trust radii, and solver workspaces remain intact. When the complete
cluster-key, sample, or boundary-dependency mapping changes, the current accepted
state initializes new fit/tail scales, cluster objectives, and the best-audit
baseline. Exact
cluster keys retain their trust radii; merged or split keys start at the initial
radius. Solver workspaces and exhausted-backtracking keys are reset, audit
patience is cleared, and convergence is disabled for that attempt so the new
partition executes at least one complete iteration. Objectives from before the
partition change are never compared with the new domain.

## Global audit and stopping

The global audit uses the fixed per-cluster fit/tail scales and retains the
earliest state that improves the best objective beyond the strict tolerance.

The stage stops on the first applicable condition:

- no valid initial seed is available for every selected atom;
- every selected atom has become terminal;
- the accepted and raw fixed-point transformed changes both converge, all
  active clusters are accepted, and no active cluster is suspicious or
  unhealthy;
- `kLocalFittingAuditPatience` accepted iterations produce no strict global
  audit improvement;
- an all-rejected attempt reaches one of the terminal resolutions below;
- `kLocalFittingMaximumIterations` outer attempts are reached.

On an all-rejected attempt, the iteration limit has the highest resolution
priority. Otherwise, rejected clusters are partitioned into objective-
backtracking-exhausted and radius-retryable sets. If any retryable radius
shrinks, the stage continues. When none can shrink, the terminal reason is:

- `all-rejected-backtracking-exhausted` when every rejected cluster is
  exhausted;
- `all-rejected-minimum-radius` when there are no exhausted clusters and every
  rejected cluster is radius-retryable and saturated at the minimum radius;
- `all-rejected-no-retry-progress` when exhausted and saturated retryable
  clusters are both present.

Numerical and invalid-model rejections that cannot perform objective
backtracking retain the minimum-radius retry behavior. Trust-radius shrink
retries do not consume audit patience while a rejected cluster's radius is
still changing.

Convergence writes the current accepted state. Audit-patience, minimum-radius
all-reject, backtracking-exhaustion, no-retry-progress, and iteration-limit
stops always write the best validated audit state when one is available;
otherwise they write the latest validated state. Best tracking, best-relative
guards, audit patience, iteration history, and stop reasons are independent of
which validated state is ultimately written. Convergence and terminal
isolation always write the latest validated state. Terminal reconciliation
preserves validated progress from non-terminal clusters.

## Final state application and group fitting

After the stopping policy selects the final validated state, the stage builds
one atom-level snapshot from the MDPDE model stored in each selected result.
This is the actual best-audit or latest-validated state chosen for application
according to the stopping condition, not the last raw proposal and not a
group-median snapshot.

For every selected atom, the stage then rebuilds its persistent peeling
sampling entries from the raw entries:

```text
peeling response = raw response
                 - sum(final selected-neighbor MDPDE responses)
                 - sum(final unselected-contributor responses)
```

Selected neighbors use their final atom-level MDPDE models. Unselected
contributors use group medians derived from the final selected state, or their
initial global seed when no matching selected group exists. The calculation
preserves the original sampling-point order and metadata. Only selected local
Gaussian results are written; unselected seeds remain transient. The selected
results and rebuilt entries are persisted with `SetGaussianResult` and
`SetPeelingSamplingEntries` before group fitting begins.

`RunSecondStageLocalFitting` then calls
`RunGroupPotentialFitting(model_object, options, true)`. The group fit reads the
newly persisted complete entries through `GetPeelingSamplingEntries(false)` and
reuses the `alpha_g` trained before the first group fit.

After the second stage returns, the workflow passes the same persisted peeling
entries to `RunLocalAlphaTraining(..., true)` and
`RunFixedOffsetLocalFitting(..., true)`, retrains `alpha_g`, and calls
`RunGroupPotentialFitting(..., true)` for the final group fit. These later local
fits may update local Gaussian results, but they do not rebuild or overwrite
the peeling entries; the final group fit therefore consumes the same
atom-level peeling snapshot written during second-stage finalization.

## Performance architecture

The fitting context stores numeric group IDs, selected group membership,
unselected-contributor group IDs, flattened sample-neighbor edges, prepared
refit design templates, and profile-radius ordering. These structures are
rebuilt only when terminal isolation changes the partition.

Each outer iteration builds one selected/unselected Gaussian snapshot, one
adjusted-response cache for refits, and one residual/objective baseline.
Cluster candidates are represented by atom-local state patches. Candidate
evaluation overlays a patch on the previous state, recomputes medians only for
groups touched by the patch, and evaluates the objective as baseline plus the
changed sample and offset delta. For a boundary-component guard, affected sample
IDs are sorted and deduplicated so a boundary sample is recomputed once. The
   boundary dependency and deterministic accepted/rescue-induced components are derived
from the current partition and rebuilt with adaptive topology or terminal
   isolation. Without a multi-cluster accepted or rescue component, candidate
   selection stays on the existing fast path.

Joint-offset, joint-polish, and boundary-correction solvers retain their sparse
pattern analysis for the lifetime of a partition and refresh only numeric
values, weights, responses, and ridge terms. Boundary-correction workspaces use
the deterministic member/interface/closure/sample signature as their key and
are cleared with the other solver workspaces when the partition changes.
Independent boundary components are corrected serially in cluster-key order.
Cluster candidate workers own independent solver state and patches. With
OpenMP, more than one cluster, and `thread_size > 1`, base proposal, local
backtracking, polish, and local objective evaluation run per cluster in
parallel while Eigen uses one thread. Results are committed in the partition's
fixed cluster order; the serial path calls the same worker and merge code.

Non-quiet runs also emit non-blocking performance counters for complete-state
materializations, Gaussian cache hits/misses, recomputed/reused objective
samples, symbolic solver analyses, adaptive topology rebuilds and partition
changes, boundary reconciliation attempts/backtracks/rejections, and
boundary joint-correction attempts/acceptances/fallbacks, symbolic analyses,
boundary rescue attempts/acceptances/fallbacks/rejections, rescue exclusion
reasons, and elapsed time. The existing
iteration/candidate/topology/total elapsed-time field remains unchanged.
These counters are diagnostic evidence rather than acceptance thresholds.

## Logging

After valid seeds are available, non-quiet runs print a compact header and
update one progress row per outer attempt with `Logger::ProgressLine`. An outer
attempt may be an accepted iteration or an all-rejected trust-region retry, so
`Try` can advance without `Acc`. The header and progress rows use the same
fixed column widths, including enough space for both scientific-notation
values in `dMax A/R`.

| Column | Meaning after the current outer attempt |
|---|---|
| `Try/Acc` | One-based outer attempt / cumulative accepted iterations |
| `Atom A/T` | Remaining active / cumulative terminal atoms |
| `Cluster A/R` | Accepted / rejected candidate clusters |
| `Polish E/A/R/S` | Eligible / accepted / rejected / skipped polish clusters after component reconciliation and the final global guard; `E = A + R + S` |
| `Suspicious` | Atoms rolled back by the suspicious-offset checks in this attempt |
| `dMax A/R` | Maximum transformed change in the accepted/raw state; accepted is `-` on an all-rejected attempt |

Objective-domain startup diagnostics report the weights, cluster and unique
fit/tail sample counts, and fixed-scale median/p99/maximum. Debug rejection
diagnostics use `fit/tail-weighted/offset/total` order and also report raw tail
loss, weights, sample counts, fixed scales, and backtracking
trials/factor/exhaustion. Accepted local factors are logged at debug level.
Each multi-cluster unit emits a distinct `Boundary-component reconciliation`
record with its cluster, atom, and boundary-sample counts plus trials, factor,
accepted/rejected, exhausted status, accepted source, and previous/endpoint/final
component objectives. An attempted correction also emits
`Boundary-interface joint correction` with its interface/closure/parameter
counts, solver status,
damping, maximum normalized trust step, strict-improvement reference/candidate
objectives, acceptance result, and endpoint-fallback outcome. An all-rejected
debug record reports
`exhausted/retryable/radius-changed/radius-saturated` counts so its retry or
terminal classification can be audited directly. Terminal, convergence, and
summary messages finish the active progress line before normal line output.

Adaptive rebuild diagnostics use a distinct
`Adaptive local-fitting topology rebuild` record so the one-time initial
coupling and residue-cutoff summaries remain stable. Each record reports the
accepted iteration, drift or interval trigger, maximum drift, old/new cluster
and boundary-sample counts, added/removed adjacency edges, and whether the
partition and objective domain changed. Non-quiet runs also show a
`Rebuild local-fitting coupling topology` percentage progress bar with the
same sample-based work accounting as the initial topology build.

Non-quiet runs end with this summary format:

```text
Second-stage local fitting summary: accepted_iterations=<N>, best_iteration=<initial|N|unavailable>, stop_reason=<reason>, best_audit_objective=<value|unavailable>, final_uses_polish=<yes|no|unavailable>, final_state_source=<best-audit|latest-validated|unavailable>.
```

`final_uses_polish` describes the state actually written to `ModelObject`. It is
`yes` when at least one atom's most recent transformed-parameter update in that
state came from an accepted polish, `no` when none did, and `unavailable` when
the stage exits before a valid state can be formed. `best_iteration` describes
the audit result and `final_state_source` identifies the state actually
written.

`quiet_mode` suppresses the second-stage informational logging.
