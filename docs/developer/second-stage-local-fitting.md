# Second-stage local fitting

## Scope

`RunSecondStageLocalFitting(ModelObject&, const FitOptions&)` refines the
Gaussian estimates of the selected atoms after the first-stage local and group
fits. It accounts for overlapping responses from neighboring selected atoms
and unselected model atoms while updating each selected atom's amplitude,
width, and offset. Unselected atoms contribute background responses but are
never added to the optimizer state.

The implication and redundancy review of the stopping predicates is recorded
in [Second-stage convergence safeguard audit](convergence-safeguard-audit.md).
Its shadow-only follow-up is the
[solver-qualification and active-coordinate population audit](stationarity-active-coordinate-audit.md).
The third-round causal experiment is the
[counterfactual convergence continuation audit](counterfactual-convergence-continuation-audit.md).

The stage keeps candidate states in memory and writes one validated final
state to `ModelObject`. Individual outer iterations do not partially update the
stored atom estimates.

`FitOptions::second_stage_boundary_halo_depth` controls boundary-correction
shape expansion and defaults to one physical-dependency hop. A value of zero
keeps the direct-interface behavior. Final uncut-component polish is enabled by
`FitOptions::enable_second_stage_dependency_polish`; its nonlinear round limit
is `FitOptions::second_stage_dependency_polish_max_iterations`, which defaults
to ten. An enabled polish with a zero round limit is rejected before any model
write. These settings intentionally have no command-line flags.

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
iterations have elapsed since the last rebuild. Rejected attempts do not
advance this interval. Quarantined parameter blocks
remain in the graph and objective domain, so quarantine does not itself rebuild
or renormalize either one.

Adaptive rebuilds use edge hysteresis. A previously absent edge must have
weight at least `0.06` to enter the graph, while an edge present in the previous
post-residue-cutoff adjacency remains until its weight falls below `0.04`.
The existing ten-residue cluster limit is then applied in descending current
edge-weight order. Binary fallback remains the conservative response to an
invalid Jacobian.

## Iteration flow

Each outer attempt performs the following sequence:

1. Build the stage-local activity masks for shape blocks, shared-offset groups,
   and hard-failure clusters. Quarantined blocks remain in the atom list as
   fixed background.
2. Partition the complete selected-atom coupling topology and reconcile the per-cluster
   objective and trust-region states.
3. Jointly estimate one shared offset per represented group within each cluster
   using robust IRLS and the fixed `kJointOffsetRidgeRatio`. Preserve the finite,
   undamped endpoint as the offset part of the fixed-point operator. A hard
   joint solve marks the cluster endpoint unavailable.
4. Build component-wise group-median models from the post-solve snapshot. For
   each selected atom, subtract its selected neighbors and all effective
   unselected contributors from the observed sample responses.
5. Refit the atom's local Gaussian with its trained `alpha_r`, using its
   group-median model as the fixed offset model. The undamped offset-to-shape
   result is the complete fixed-point operator endpoint.
6. Search one geometric factor sequence `1, 1/2, 1/4, ...` for each cluster.
   Each factor constructs log-shape and shared-physical-offset coordinates,
   skips candidates outside the trust radius, then applies offset-only and
   post-refit feasibility guards, and finally the previous/best objective gate.
   Guard never damps, trust never accepts, and the objective gate never chooses
   a second independent factor. If every material factor is guard-infeasible,
   deactivate only the terminal shape or offset-group block and repeat the same
   search for the remaining blocks. Objective exhaustion rejects the cluster.
7. For a solver-qualified cluster, attempt one joint
   amplitude/width/offset polish over its active columns. Inactive shapes and
   offset groups decode to the endpoint values and remain fixed background.
   Keep the polish only when it strictly improves the base candidate on the
   same objective scale.
8. Build the accepted-induced interaction graph from shared boundary samples.
   Revalidate each connected component. Every eligible component attempts one
   joint correction over the boundary shape-active set. This starts with the
   direct boundary contributors and expands by the configured number of hops.
   A hop follows only the raw sample target and directly selected neighbors and
   remains inside the boundary component; it does not follow virtual edges
   induced by an unselected contributor's group median. Shape-active atoms may
   change transformed amplitude and width. The offset closure is then rebuilt
   from their group keys, with one component-local offset column shared by all
   same-group closure atoms; closure-only shapes remain fixed. A valid
   endpoint remains the fallback and is replaced only when the correction fits
   every member trust radius and strictly improves that endpoint. For an invalid
   endpoint, the correction must strictly improve the previous component
   objective; an unavailable or rejected correction falls through to the
   unchanged common-factor component backtracking. A failed component rolls back
   only its members, so unrelated components and remote singleton clusters remain
   eligible for commit. Then build maximal eligible boundary components from the
   safe accepted state and the best finite objective-rejected proposal retained
   for each rejected cluster. Components containing at least one such proposal
   receive the same endpoint, active-column joint correction, and common-factor
   checks. Cooperative rescue permits a member to deteriorate only within
   `1e-8 + 1e-3 * abs(previous member objective)`; the component and tentative
   assembled global objectives must both improve strictly, and the global
   historical-best tolerance still applies. A failed correction retains a valid
   rescue endpoint, and any failed rescue atomically leaves the safe accepted
   state unchanged.
9. Run the unchanged global previous/best audit on the complete assembled state.
   If tolerance-level deterioration from independent reconciliation units causes
   aggregate rejection, remove non-improving units from worst to best until the
   first passing subset is found. If only strictly improving units remain but the
   historical best gate still fails, roll back the attempt without weakening the
   gate.
10. Update trust radii and the reversible quarantine/probation state. After an accepted state,
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

This offset solve and neighbor-adjusted local refit constitute one undamped
fixed-point operator update. These group-median-adjusted entries are temporary
operator inputs; they are not persisted as peeling sampling entries.

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

Convergence statistics are evaluated over active optimization degrees of
freedom, not every selected atom. Shape-active atoms contribute one log-peak
and one log-width sample each. Every active shared offset within one
`(cluster, group_id)` contributes one offset-to-peak sample equal to the
maximum absolute transformed change among its members. Fixed and quarantined
coordinates are excluded; mixed-activity groups and non-finite member changes
fail convergence.

Production convergence requires accepted active-DOF p99 and complete nominal-
DOF fixed-point residual p99 below `1e-4`, with solver qualification and all
orthogonal blockers clear. Maximum transformed change remains a tail diagnostic
and topology-drift metric, but is not a convergence predicate.

The strict fixed-point operator
`F(S[k])`: one undamped joint-offset solve followed by undamped local shape
refits using those offsets. If production fixed/quarantine handling changes an
offset, all selected shapes are refit once in an isolated workspace. If an
unrestricted offset is unavailable, shape residuals are marked unavailable
rather than replaced with zero. This strict operator uses the nominal selected
shape and shared-offset DOFs, including fixed and quarantined blocks, and does
is the production residual source.

Production convergence qualification is the existing active-block cluster
rollup. Solver qualification is a separate developer comparator: active local
shape refits must report `SUCCESS`, active shared offsets must report
`Converged`, and both require a full undamped, non-fallback endpoint. A usable
soft endpoint may continue through candidate selection without being solver
qualified. This distinction does not change the production stopping policy.

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

The global objective weights clusters by selected atom count. Quarantined
blocks remain included in this normalization:

```text
global objective = sum((cluster atom count / selected atom count) * cluster total)
```

Owner assignment makes every sample appear once in this global sum, including
boundary samples. Local scoring still includes every sample affected by the
candidate cluster. It applies that sample's owner scale and exact global
normalization coefficient, so the local candidate-minus-previous difference
matches the corresponding full-global difference when only that cluster
changes.

Candidate scoring uses a provisional copy of the cluster objective state. A
rejected base, polish, candidate-search trial, or boundary-component candidate does not
advance the previous or best references. The best objective and maximum
transformed change are retained to break objective ties.

All second-stage audit tolerances use:

```text
tolerance(reference) = absolute tolerance + relative tolerance * abs(reference)
```

Progress, deterioration, and actual-reduction trust-growth comparisons use
`1e-8 + 1e-3 * abs(reference)`. Strict best, tie, and polish-improvement
comparisons use `1e-10 + 1e-8 * abs(reference)`. Candidate comparisons against
previous and best compute separate tolerances from their respective references.
The joint-offset IRLS objective retains its independent tolerance.

## Trust region

For a non-suspicious cluster, previous and operator offsets are reduced to one
physical offset per `GroupKey` by deterministic component medians. Base
proposal trials use factors `1, 1/2, 1/4, ...`; each trial interpolates the
atom-level log-peak and log-width coordinates and the physical group offset
with the same factor:

```text
C_g(t) = C_previous,g + t * (C_raw,g - C_previous,g)
```

The realized atom models are then re-encoded and measured against the trust
radius. At `t = 0`, the shape is the previous shape and the offset is the
previous group median; at `t = 1`, the complete operator shared-offset state is
recovered. If projecting an inconsistent previous group onto its median already
exceeds the radius, the proposal fails before objective evaluation with an
explicit diagnostic reason.

Each cluster owns one factor sequence. Every trial first constructs the
log-shape/shared-physical-offset candidate, then checks validity, trust
admissibility, guard feasibility, and the previous/best objective gates in that
order. Trust-inadmissible trials do not run guard or objective evaluation.
Search stops when the largest transformed change is below
`kTransformedChangeTolerance`; the first passing material trial is committed
with endpoint uncertainty and its factor is recorded. Rejected trials do not
mutate objective state or polish provenance.

When every material factor is guard-infeasible for one shape or shared-offset
group, that block is made locally inactive and the same function restarts the
factor sequence for the remaining blocks. This is an iterative block-isolation
loop, not a recursive candidate selection or a new outer attempt.

The polish step is limited by the radius remaining after the accepted base
movement. A rejected polish keeps the base candidate and is not backtracked.
Radius updates use one controller entry point while preserving the validated
baseline order: accepted objective-backtracking shrink, accepted growth, then
retryable rejection shrink. The first guard-feasible factor that reaches the
objective gate is the accepted-shrink reference. Guard-only factor reduction
therefore does not shrink the radius; a later objective rejection followed by
acceptance at a smaller factor does. Local terminal rejection remains
retryable; an exhausted boundary or final-audit search keeps its radius because
another shrink cannot produce a material candidate. An accepted cluster grows
its radius only when the actual objective reduction exceeds the existing
progress materiality tolerance (`1e-8 + 1e-3 * abs(previous)`), its step is
within the outer 20% of the current radius, and no objective-backtracking
shrink is already required. Guard-only backtracking does not suppress this
independent growth rule. No radius update reruns the same validated state;
updates remain isolated by cluster and clamp to `0.0625...4.0`.

This remains an actual-reduction-aware transformed-step cap, not a model-based
trust region: the production controller does not consume a predicted reduction
or actual/predicted reduction ratio. Candidate acceptance and historical-best
gates are unchanged; actual reduction only controls whether the next radius may
grow.

An audit build computes a developer-only frozen-IRLS directional model for
every material base or joint-polish trial that reaches the objective gate,
including locally accepted and objective-rejected candidates. It freezes Cauchy
weights and objective scales at the outer previous state, applies the existing
transformed response Jacobian to the complete previous-to-candidate step, and
includes selected targets, selected neighbours, and group-median-derived
unselected contributors. The residual surrogate uses the production fit/tail
sample coefficients and the exact previous-to-candidate offset-plausibility
penalty change. A ratio is emitted only when predicted reduction exceeds the
same progress materiality tolerance and both reductions are finite. Polish rho
uses the complete outer-previous-to-polished step; improvement relative to the
accepted base candidate is recorded separately.

The shadow status distinguishes `available`, `nonmaterial-step`,
`objective-unavailable`, `model-unavailable`, `residual-unavailable`,
`nonfinite`, `nonpositive-prediction`, and `nonmaterial-prediction`. Its
counterfactual action preserves objective-backtracking shrink precedence,
falls back to the current actual-only action when prediction is unusable, uses
rho bands `0.25` and `0.75`, and requires at least `0.8` boundary utilization
for counterfactual growth. Only the final locally accepted candidate can become
action-ready. Boundary-reconciled, rescued, globally rejected, and non-final
trial records remain available for coverage and calibration but are suppressed
from action comparison. Pre-objective validity, trust, guard, and nonmaterial
outcomes are counted in a separate candidate funnel without running the model.
None of these calculations run unless
`RHBM_GEM_ENABLE_COUNTERFACTUAL_CONVERGENCE_AUDIT=ON`, and no shadow result is
applied to candidate acceptance, radius updates, convergence, or output state.

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

After this accepted-only pass, every rejected cluster that produced a finite
objective proposal retains its lowest-objective trust-region patch in memory.
The rejected proposals and safe accepted clusters form maximal eligible
boundary components. Suspicious and hard-failure atoms are not removed from a
component; their inactive parameter blocks remain fixed while active neighbors
can supply correction columns. A cooperative endpoint may contain a member
whose objective is slightly worse than its previous value, but only within the
normal progress tolerance. The endpoint or its joint/backtracked replacement
must strictly improve the component audit, and the tentatively assembled state
must then strictly improve the previous global audit without violating the
historical-best tolerance. Member historical bests are updated only for actual
member improvements and do not independently veto rescue. Successful rescue
promotes the rejected members without trust-radius growth or shrink. Failed
rescue is transactional: every component member retains its safe state.

Boundary correction eligibility is intentionally broader than local polish
eligibility. A component is eligible whenever at least one shape or shared-offset
column remains active. Inactive shape columns use endpoint amplitude and width;
an inactive offset group has no offset column and all group members use the
endpoint offset. Candidate guards inspect only materially changed active blocks,
while fixed-block invariants require every inactive value to remain bitwise at
its endpoint. IRLS objective deterioration, IRLS iteration exhaustion, and valid
non-success local estimation statuses may participate. Local joint polish still
requires full solver qualification. Every eligible accepted-only or rescue boundary
component receives at most one correction attempt per outer iteration.

Before a joint-correction candidate is accepted, neighbor-adjusted profiles are
rebuilt from its formally resolved selected and unselected snapshots. Every
materially changed atom must pass the post-refit suspicious-profile guards.
Failure keeps the exact endpoint. Accepted-only reconciliation and rejected
cluster rescue use the same halo construction and correction path.

After component reconciliation, the complete assembled state must still pass the
unchanged global previous/best audit. On aggregate failure, independent components
and singleton clusters are scored by their exact global objective delta. Only
non-improving units are removed, from worst delta to best, with the full audit
recomputed after each removal. The first passing subset is committed atomically.
If all remaining units strictly improve the previous objective but cannot satisfy
the historical best gate, the complete remaining attempt is marked exhausted;
objective tolerances are never relaxed.

## Final uncut dependency polish

After the existing stop policy selects the best-audit or latest-validated base
state, but before peeling entries or atom models are written, the stage attempts
one final dependency polish. Current clusters are first treated as indivisible
DSU units. They are then merged using the complete
`GraphTopology::sample_dependency_list`, without the weighted-edge threshold or
ten-residue cutoff. This deliberately retains virtual dependencies introduced
when an unselected contributor resolves through a selected-group median.
Quarantined shape blocks and offset groups are not variables, but their fixed
models remain in every sample response and in the objective domain. A component
is skipped only when it has no active shape or offset column.

Each component is solved serially in fixed key order. Shape-active member atoms
have `log(amplitude)` and `log(width)` variables, and each offset-active group
represented in that component has one shared physical offset. Inactive blocks
decode from the endpoint. The same group in two independent components is not
tied across components. Sparse weighted-ridge directions,
robust weights, conditioning guards, and each original cluster's trust radius
are reused from boundary correction. Up to the configured number of nonlinear
rounds is attempted. A round linearizes at its latest endpoint, while every
cumulative trust step is measured from the pre-polish base state.

Direction construction fixes the current unselected-contributor snapshot.
Candidate validation uses the formal resolver again, so selected-group medians,
unselected responses, residuals, and the nonlinear objective reflect the
candidate. A component patch requires valid Gaussian parameters, post-refit
suspicious guards, every member trust radius, all member-cluster objective
guards, and strict component-objective improvement. Solver or validation failure
falls back only that component.

Accepted component patches are assembled and subjected to a complete global
audit. If it fails, exact full-audit deltas are used to remove the worst
non-improving component patches. Unless the surviving state strictly improves
the pre-polish global objective, the complete polish is discarded and the base
state is written unchanged. When the selected base state stopped by production
convergence, an objective-accepted polish remains provisional until the strict
fixed-point operator is evaluated again on the polished state. The polish is
applied only when every solver is qualified, the nominal-DOF operator is
complete, and every operator-residual p99 is below `1e-4`. An unavailable,
failed, or above-threshold certificate discards the polish and writes the
already converged base state unchanged. Other stop reasons retain the
objective-only final-polish policy. An applied final polish updates the final
audit and provenance but does not increment the outer accepted-iteration count
or change the stop reason. The final-polish diagnostic distinguishes objective
acceptance, certificate status, and actual application.

## Numerical defenses, partial active set, and quarantine

- Suspicious evaluation has offset-only and post-refit modes. Offset-only checks
  finite zero-offset responses, offset magnitude, center sign flip, and radial
  rebound. Post-refit first requires a valid model and additionally checks width
  growth and amplitude-offset compensation. Guard precedence is unchanged.
- Every assessment records the selected reason, guard mode, and signed
  normalized margin. A one-threshold margin is
  `observed / limit - 1`; AND predicates use the minimum constituent margin and
  OR predicates use the maximum. Positive means violated, zero is the boundary,
  negative is safe, and invalid or non-finite inputs use positive infinity.
- The previous suspicious baseline is built in one fit-range scan. It records
  the innermost response, per-radius response medians, distance range, maximum
  absolute response, and residual scale `1.4826 * MAD`. Candidate profiles do
  not calculate a residual MAD. The post-refit candidate and its offset-only
  fallback reuse the same previous baseline.
- The center sign-flip guard treats the smallest sampled radius as the
  innermost response and only rejects a statistically significant
  positive-to-negative change. It requires the previous response to exceed
  `3 * scale`, and the candidate to be below the negative of both that noise
  threshold and `0.25 * previous innermost response`. The radial rebound guard
  uses the same noise estimate, its existing magnitude thresholds, and more
  than one upward excursion. Sign flip and rebound require a trustworthy
  previous radial shape; width growth and amplitude-offset compensation do not.
- Guard is feasibility-only. The cluster controller owns the sole factor list
  and applies one common factor to log-shape and shared physical offsets. Trust-
  inadmissible candidates skip guard and objective evaluation. A full,
  solver-qualified endpoint below `kTransformedChangeTolerance` is stationary;
  reaching that tolerance only after factor reduction is step-limited.
- Failure is block-local. An unsafe offset update fixes the complete shared
  offset group; an unsafe shape update fixes only that atom's amplitude/width;
  and a hard joint-offset solve fixes the cluster's offset blocks while allowing
  safe shape refits. No post-refit guard causes a complete-cluster rollback.
- The same shape/offset/hard masks parameterize local polish, boundary
  correction, rescue, and final dependency polish. Suspicious or quarantined
  atoms and their samples remain in coupling components, residual evaluation,
  and the objective domain. Only inactive columns are omitted, and their decoded
  values must equal the endpoint.
- A fixed shape may retain a guard-safe jointly estimated offset. Conversely, a
  fixed offset group does not prevent safe shapes in that group or cluster from
  changing. The next attempt applies the `10x` suspicious ridge multiplier to
  affected atoms.
- Stable near-convergence failures enter stage-local quarantine only after the
  same target and reason occurs in five accepted iterations. Targets are an atom
  shape block, a shared-offset group, or a hard-failure cluster. A changed reason
  or missing observation resets the pre-quarantine count.
- Quarantine never removes atoms or rebuilds the objective domain. After two
  accepted iterations a target receives probation; a topology partition change
  may trigger it early. Each target gets at most three probes. Probes use the
  minimum trust radius `0.0625` and `10x` ridge. Overlapping probes are selected
  in cluster, group, then atom priority.
- A material probation proposal must pass guards, trust/fixed-block invariants,
  member/component/global objective gates, and the historical-best gate. A
  guard-safe non-material stationary result at the minimum radius may also
  release the target. Failure restores the previous fixed block and increments
  the probe count; after three failures the target remains fixed only until the
  current second-stage call ends.
- Rejected-cluster debug output distinguishes failures before objective
  evaluation from objective rejection. Pre-objective failures report their
  proposal reason, radius, available step norm, and `objective =
  not-evaluated`; `objective-unavailable` is reserved for an objective that was
  actually attempted but could not be calculated.

An adaptive topology rebuild always becomes the next hysteresis reference,
even when its active partition is unchanged. In that case the objective scales,
best audit, trust radii, and solver workspaces remain intact. When the complete
cluster-key, sample, or boundary-dependency mapping changes, the current accepted
state initializes new fit/tail scales, cluster objectives, and the best-audit
baseline. Exact
cluster keys retain their trust radii; merged or split keys start at the initial
radius. Solver workspaces are reset, audit
patience is cleared, and convergence is disabled for that attempt so the new
partition executes at least one complete iteration. Objectives from before the
partition change are never compared with the new domain.

## Global audit and stopping

The global audit uses the fixed per-cluster fit/tail scales and retains the
earliest state that improves the best objective beyond the strict tolerance.

The stage stops on the first applicable condition:

- no valid initial seed is available for every selected atom;
- accepted active-DOF p99 and complete nominal-DOF fixed-point residual p99 are
  both below `1e-4`, every solver is qualified, all clusters are accepted, and
  no orthogonal blocker is present;
- `kLocalFittingAuditPatience` accepted iterations produce no strict global
  audit improvement;
- an all-rejected attempt terminates after applying its per-cluster radius
  actions once;
- `kLocalFittingMaximumIterations` outer attempts are reached.

An all-rejected attempt does not rerun the unchanged `S[k]`. Retryable
rejections shrink their stored radius once, while exhausted terminal searches
keep it, then the attempt stops with `all-rejected-backtracking-exhausted`
unless the outer iteration limit has priority.

Convergence writes the current accepted state. Audit-patience, all-rejected,
and iteration-limit stops always write the best validated audit state when one
is available;
otherwise they write the latest validated state. Best tracking, best-relative
guards, audit patience, iteration history, and stop reasons are independent of
which validated state is ultimately written. Unresolved quarantine targets are
kept at their latest validated fixed values and do not prevent unrelated active
blocks from converging or being written.

## Final state application and group fitting

After the stopping policy selects the final validated state, the stage builds
one atom-level snapshot from the MDPDE model stored in each selected result.
This is the actual best-audit or latest-validated state chosen for application
according to the stopping condition, not the operator endpoint and not a
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
refit design templates, and profile-radius ordering. Quarantine changes only
activity masks; these structures are rebuilt only when adaptive topology changes
the partition.

Each outer iteration builds one selected/unselected Gaussian snapshot, one
adjusted-response cache for refits, and one residual/objective baseline.
Cluster candidates are represented by atom-local state patches. Candidate
evaluation overlays a patch on the previous state, recomputes medians only for
groups touched by the patch, and evaluates the objective as baseline plus the
changed sample and offset delta. For a boundary-component guard, affected sample
IDs are sorted and deduplicated so a boundary sample is recomputed once. The
The boundary dependency and deterministic accepted/rescue-induced components
are derived from the current partition and rebuilt with adaptive topology.
Without a multi-cluster accepted or rescue component, candidate selection stays
on the existing fast path.

Joint-offset, joint-polish, and boundary-correction solvers retain their sparse
pattern analysis for the lifetime of a partition and refresh only numeric
values, weights, responses, and ridge terms. Boundary-correction workspaces use
the deterministic member/interface/closure/sample signature as their key and
are cleared with the other solver workspaces when the partition changes.
Independent boundary components are corrected serially in cluster-key order.
Cluster candidate workers own independent solver state and patches. With
OpenMP, more than one cluster, and `thread_size > 1`, base proposal, local
candidate search, polish, and local objective evaluation run per cluster in
parallel while Eigen uses one thread. Results are committed in the partition's
fixed cluster order; the serial path calls the same worker and merge code.

Non-quiet runs also emit non-blocking performance counters for complete-state
materializations, Gaussian cache hits/misses, recomputed/reused objective
samples, symbolic solver analyses, adaptive topology rebuilds and partition
changes, boundary reconciliation attempts/backtracks/rejections, and
boundary joint-correction attempts/acceptances/fallbacks, symbolic analyses,
boundary rescue attempts/acceptances/fallbacks/rejections, rescue exclusion
reasons, final dependency-polish component/atom/parameter/round counts,
acceptances/fallbacks, and elapsed time. The existing
iteration/candidate/topology/total elapsed-time field remains unchanged.
These counters are diagnostic evidence rather than acceptance thresholds.

## Logging

After valid seeds are available, non-quiet runs print a compact header and
update one progress row per outer attempt with `Logger::ProgressLine`. An outer
attempt may be accepted or terminally rejected, so `Try` can advance without
`Acc`. The header and progress rows use the same
fixed column widths, including enough space for both scientific-notation
values in `dMax A/F`.

| Column | Meaning after the current outer attempt |
|---|---|
| `Try/Acc` | One-based outer attempt / cumulative accepted iterations |
| `Atom A/Q` | Atoms without quarantine targets / atoms covered by current quarantine targets |
| `Cluster A/R` | Accepted / rejected candidate clusters |
| `Polish E/A/R/S` | Eligible / accepted / rejected / skipped polish clusters after component reconciliation and the final global guard; `E = A + R + S` |
| `Suspicious` | Atoms with at least one shape, offset, or hard-failure block fixed in this attempt |
| `dMax A/F` | Maximum transformed change in the accepted/fixed-point operator state; accepted is `-` on an all-rejected attempt |

Objective-domain startup diagnostics report the weights, cluster and unique
fit/tail sample counts, and fixed-scale median/p99/maximum. Debug rejection
diagnostics use `fit/tail-weighted/offset/total` order and also report raw tail
loss, weights, sample counts, fixed scales, unified trial dispositions, and the
accepted factor. Accepted local factors are logged at debug level.
Each multi-cluster unit emits a distinct `Boundary-component reconciliation`
record with its cluster, atom, and boundary-sample counts plus trials, factor,
accepted/rejected, exhausted status, accepted source, previous/endpoint/final
component objectives, locally deteriorated member count, and maximum local
deterioration. Rescue records also include component and final assembled-global
improvements. An attempted correction also emits
`Boundary-interface joint correction` with its
direct-interface/shape-active/offset-active/offset-closure/parameter counts, suspicious count,
solver status,
damping, maximum normalized trust step, strict-improvement reference/candidate
objectives, acceptance result, and endpoint-fallback outcome. An all-rejected
debug record reports unified trial dispositions, terminal category, radius
action, and stop classification. Operator-assessment debug records report each
affected atom's reason, margin, and fixed shape/offset/hard block. Completion
warnings report cumulative quarantine
entries, releases, failed probation probes, and unresolved targets. Convergence
and summary messages finish the active progress line before normal line output.

Adaptive rebuild diagnostics use a distinct
`Adaptive local-fitting topology rebuild` record so the one-time initial
coupling and residue-cutoff summaries remain stable. Each record reports the
accepted iteration, drift or interval trigger, maximum drift, old/new cluster
and boundary-sample counts, added/removed adjacency edges, and whether the
partition and objective domain changed. Non-quiet runs also show a
`Rebuild local-fitting coupling topology` percentage progress bar with the
same sample-based work accounting as the initial topology build.

Finalization emits a separate `Final dependency polish` record with aggregate
component, atom, parameter, round, acceptance/fallback, objective-before/after,
and elapsed-time values. Debug logging adds one record per component, including
symbolic-analysis and suspicious-candidate counts. These records are emitted
before the existing stable final summary.

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
