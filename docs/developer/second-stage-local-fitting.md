# Second-stage local fitting

## Scope

`detail::RunSecondStageIterations(ModelObject&, const FitOptions&)` refines the
Gaussian estimates of the selected atoms after the first-stage local fit. It
accounts for overlapping responses from neighboring selected atoms
and unselected model atoms while updating each selected atom's amplitude,
width, and offset independently. Chemical `GroupKey`/`group_id` is neither
queried nor stored anywhere in this stage; there are no synthetic groups,
shared-offset columns, or selected group-median refits. Residue identity is also
neither queried nor stored: chain IDs and sequence IDs do not constrain
components. Optimization clusters still organize selected sample coupling,
numerical solves, and acceptance gates, with at most 100 selected atoms each.
Chemical group estimation runs once after the second local fitting stage.

Effective unselected contributors have no optimizer state or
estimation degrees of freedom. All three parameters of their fixed background
model are the component-wise MDPDE medians of all selected atoms, independent
of the sample target's optimization cluster. This background is refreshed from
the latest accepted selected state at each outer-iteration boundary and frozen throughout that
attempt, including candidate selection, polish, audit, and operator evaluation.
Unselected responses still contribute to fitting, residuals, and final peeling,
but no unselected model is persisted.

This page specifies the current production algorithm. Earlier decisions and
evidence are recorded in the
[Second-stage outer-iteration algorithm audit](second-stage-outer-iteration-algorithm-audit.md).
Earlier safeguard, population, and continuation reviews are
historical records linked only from that audit; they do not override this
page's current frozen-background contract.

Implementation ownership and the unchanged per-path validation sequence are
specified in [P0 structural refactoring](second-stage-p0-structure.md).
Candidate selection now uses a private builder and one transaction publication;
convergence decisions receive p99 evidence separately from diagnostic maxima.

The stage keeps candidate states in memory and writes one validated final
state to `ModelObject`. Individual outer iterations do not partially update the
stored atom estimates.

`FitOptions::second_stage_boundary_halo_depth` controls boundary-correction
parameter halo expansion and defaults to one physical-dependency hop. A value of zero
keeps the direct-interface behavior. Final uncut-component polish is enabled by
`FitOptions::enable_second_stage_dependency_polish`; its nonlinear round limit
is `FitOptions::second_stage_dependency_polish_max_iterations`, which defaults
to ten. An enabled polish with a zero round limit is rejected before any model
write. These settings intentionally have no command-line flags.

## Implementation responsibilities and state ownership

The internal implementation follows the iteration sequence rather than exposing
all second-stage services through candidate selection:

| Module in `src/core/detail` | Responsibility |
| --- | --- |
| `IterationProcess` | Initialization, frozen-background and pending-partition boundaries, convergence and stop decisions, final certification, and persistence |
| `IterationProposal` | Joint offsets, local shape refits, fallback, and unrestricted fixed-point operator evidence |
| `CandidateSelection` | Builder-owned per-cluster candidate search, local joint polish, and trust-radius control |
| `CandidateEvaluation` | Typed references with scopes only where policy differs; separate local/boundary results, original gate ordering, and proposed history updates |
| `CandidateTransaction` | Private selection builder, staged quarantine, and consuming publication of validated results |
| `BoundaryReconciliation` | Boundary correction, backtracking, rescue, complete-selection audit/salvage |
| `DependencyPolish` | Final uncut-component candidate generation, assembly and salvage; validation delegates to `CandidateEvaluation` |
| `ObjectiveEvaluation` | Objective domains, full and incremental evaluation, tolerances, and previous/best objective history |
| `SuspiciousUpdate` | Profile baselines, suspicious assessments, coordinate activity, and candidate/polish guards |
| `Quarantine` | Active/Frozen failure tracking, domain retry, and next-iteration activity |
| `Diagnosis` | Progress and certificate output, graph/objective diagnostics, and performance counters |
| `PhaseAudit` / `TrustModelAudit` | Observation-only snapshots, isolated probes, frozen-IRLS/rho trials and serialization |

`GaussianModelOperations`, `PreparedLocalGaussianFit`, `SecondStageFitting`,
`CouplingGraph`, and `JointFitting` retain the underlying model operations,
prepared design, state/residual representation, graph construction, and solvers.
The public Gaussian estimator workflow uses internal fitting-range constants;
`FitOptions` does not expose radial bounds.

`CandidateSelectionInputs` contains read-only algorithm inputs. The private
`CandidateTransactionBuilder` owns working activity, state, provenance and
objective history; boundary operations are its private methods. Evaluators
return values and proposed history updates without changing caller history.
Solver workspaces, counters and observers remain mutable working resources.
Quarantine is staged on a copy and changes only next-iteration activity,
without modifying the audited state or requiring fallback re-audit. The builder then freezes a read-only `CandidateTransaction`.
Its consuming commit publishes validated state, provenance, history, quarantine
and radius actions, including the required rejection updates on all-rejected
attempts. Convergence uses the committed candidate and the retained previous
state. The unrestricted operator evidence remains separate from these
production restrictions.

`SecondStageContext::atom_list` is accessed directly. Model snapshots capture
both the selected Gaussian models and the immutable background in effect when
they are built. `BuildSecondStageModelSnapshot` accepts a full fit state, a
patch-backed state view, or existing model snapshots; callers do not need a
separate model-projection step. The topology drift reference retains only model
snapshots, while accepted and best-audit states retain complete local results.
Patch/view and residual-overlay representations continue to avoid full-state
materialization during candidate evaluation.

The `IterationDiagnostics::proposal_maximum_transformed_change` field measures the constrained
proposal's maximum movement. It is distinct from the nominal operator residual
in `ConvergenceCertificate`. Maximum/population measurements reside in
`ConvergenceDiagnostics`; the existing progress label and audit schemas are
unchanged.

## Model context and initialization

The fitting context and optimizer state contain selected atoms only. Each
selected target stores its raw local-potential sampling entries, trained
`alpha_r`, prepared local design, and selected sample-neighbor edges.
Unselected background geometry stores only per-sample distance lists. Atom
identity is used only during initialization for neighbor search, deduplication,
and neighbor counting; it is not retained in the context. There is no
unselected support-row design, local result, `alpha_r`, or independent model.
Neither selected nor unselected chemical `GroupKey`/`group_id` is queried or stored.

Neighbor candidates are searched within `kNeighborAtomSearchRange`. A neighbor
contributes to a sample only when its distance from that sample does not exceed
`kNeighborContributionDistanceMax`. Each physical atom contributes at most
once to a sample. When `FitOptions::exclude_hydrogen` is true, hydrogen atoms
are removed from this contributor set; other selection exclusions do not remove
background contributors. One unselected atom can affect targets in different
clusters without connecting those clusters, and uses each target's own frozen
cluster model.

Before the second stage starts, the workflow independently trains one
first-stage local `alpha_r` for each selected atom from that atom's raw sampling
entries. Atoms below the ten-entry training threshold use the training minimum
of zero instead of inheriting another atom's value. The workflow then performs
fixed-offset local fitting and copies the complete local result from
`FittingStage::First` to `FittingStage::Second`. The first stage does not train
`alpha_g` or run group fitting.

The initial Gaussian seed is rebuilt for every selected atom. The first valid
source is selected in this order:

1. the atom's copied first-stage local MDPDE;
2. the global local-MDPDE parameter median.

The component-wise global amplitude, width, and offset medians are computed
once from all valid selected local MDPDE models. Invalid local models are
excluded, and atoms filled from the median are not fed back into its source
pool. A valid direct local seed retains its complete model and uncertainty. A
global-median fallback replaces the complete local MDPDE and uses zero
uncertainty. Local OLS, group posterior, group prior, and same-group medians are
not initialization sources.

If a valid seed cannot be obtained for every selected atom, the stage exits
without changing the stored estimates, peeling sampling entries, or group
results.

The selected initial state supplies the weighted coupling topology and
deterministic cluster keys. Unselected contributors are not topology nodes and
create no hard dependency edges or cluster ownership constraints. After the
selected-only partition is built, the stage builds a complete, finite frozen
background before constructing the initial objective baseline. The first attempt
reuses this cache and baseline; boundary refresh starts with the second attempt.

The global median includes every selected atom in the current state with equal
weight, including fixed and quarantined atoms. Amplitude, width, and offset
are reduced independently by the existing Gaussian parameter median helper:
odd populations use the central value; even populations average the two central
values. Repeated selected groups have no special weighting. Unselected atoms
never enter the pool, and background models have no uncertainty estimate.
Responses from all effective unselected contributors are summed into one
immutable cache per selected sample. A missing, incomplete, or non-finite cache
aborts the stage without writing Gaussian or peeling results; no partial cache
is committed. There is no separate unselected seed or refit failure path.

The initial weighted topology retains the fixed minimum edge weight `0.05`.
After accepted iterations, the stage adaptively rebuilds the topology from the
latest validated atom models only when the maximum transformed-coordinate
drift from the last topology reference state reaches `0.10`. Accepted iteration
count does not trigger a rebuild. The reference snapshot is updated only at
initialization and after an actual rebuild, so drift accumulates relative to
the last rebuild rather than the previous iteration. Quarantined parameter blocks
remain in the graph and objective domain, so quarantine does not itself rebuild
or renormalize either one.

Adaptive rebuilds use edge hysteresis. A previously absent edge must have
weight at least `0.06` to enter the graph, while an edge present in the previous
post-atom-cutoff adjacency remains until its weight falls below `0.04`.
Initial, adaptive, and binary-fallback topologies all apply the same internal
100-selected-atom component limit. Threshold-retained atom edges are processed
in descending current weight order, with equal weights ordered by canonical
atom-index pair. Union-find merges two components only when their combined
atom count does not exceed the limit. All retained edges internal to each final
component remain in its adjacency, not only the edges used for merging.
Disconnected atoms remain singletons, even when their residue labels match.
Fixed and quarantined selected atoms count toward the limit; unselected atoms
do not. Coupling Jacobians and dependency components include only selected
participants. Binary fallback remains the conservative response to an invalid
selected Jacobian.

The limit constrains topology partitions, not boundary or final-polish solve
sizes. The full physical sample dependencies and pre-cutoff retained edge list
are preserved. Boundary reconciliation and final dependency polish may connect
multiple optimization clusters and exceed 100 atoms, without changing the
global selected-atom background median pool. A changed
partition still takes effect together with its background at the next iteration
boundary, using the existing domain reset, previous/best rescoring, and
domain-change convergence blocker.

## Iteration flow

Each outer attempt performs the following sequence:

1. Apply any pending selected-only partition at the iteration boundary. Build
   and freeze the global-median background from the latest accepted selected
   state, then re-evaluate previous and retained best objectives under this same
   background. Rebuild the sampling domain and solver workspaces only if the
   partition changes. The first attempt uses the initialization background and
   baseline without rebuilding or logging the same background twice.
2. Build the selected per-atom shape, offset, and hard-failure activity masks.
   Fixed and quarantined selected atoms remain in the graph, median pool,
   samples, and objective domain. Reconcile per-cluster objective and trust
   states; unselected geometry has no activity mask.
3. Jointly estimate selected offsets with one column per atom, ordered by atom
   index within each cluster. Each seed and ridge anchor is that atom's previous
   offset. Robust IRLS, conditioning/collinearity guards, and
   `kJointFittingRidgeRatio` are unchanged. Preserve the complete undamped
   nominal solve as fixed-point evidence; production masks, validity, and
   fallback apply per atom. Outside-cluster selected neighbors and the frozen
   unselected response enter the fixed RHS, not parameter columns.
4. Freeze the post-offset atom-level snapshot. Subtract each selected neighbor's
   own complete Gaussian response and the unchanged unselected background cache.
   The target's offset response uses its own width and offset; no selected
   median replaces target or neighbor models.
5. Refit selected atoms with their own trained `alpha_r` and existing
   fallback/quarantine handling. There are no unselected refits. The unrestricted
   operator uses the same frozen background and does not manufacture an endpoint
   from production fallback.
6. Search one geometric factor sequence `1, 1/2, 1/4, ...` for each cluster.
   Each factor constructs selected log-shape and per-atom physical-offset
   coordinates, skips candidates outside the trust radius, applies offset-only
   and post-refit feasibility guards, and then the previous objective gate. Cluster best is history/tie-break only.
   The frozen background is unchanged at every factor. Guard never damps, trust
   never accepts, and the objective gate never chooses a second independent
   factor. If every material factor is guard-infeasible, deactivate the terminal
   independent atom shape or offset block and repeat the search.
   Objective exhaustion rejects the cluster.
7. For a solver-qualified, fully active cluster, attempt one joint
   amplitude/width/offset polish. Local polish retains its existing eligibility
   restrictions. In the shared polish parameterization, each inactive
   coordinate decodes to that atom's own endpoint value. Keep the polish only
   when it strictly improves the base candidate on the same objective scale.
8. Build the accepted-induced interaction graph from shared boundary samples.
   Revalidate each connected component. Every eligible component attempts one
   joint correction over active coordinates in the physical halo, starting with
   direct boundary contributors and expanding by the configured number of hops.
   A hop follows the raw sample target and every direct selected neighbor
   inside the boundary component. Shape and offset active lists independently
   filter this physical halo using their atom-level masks; their union is the
   parameter atom set. There is no chemical-group closure. Each active shape has
   two columns and each active offset has its own column. Merging boundary
   components does not change the frozen global background. A valid
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
10. Stage Active/Frozen quarantine for the next iteration without changing
    the audited state. Publish quarantine and trust-radius updates together. An accepted state's adaptive rebuild can queue a
    new partition for the next attempt; pending changes block convergence.
    Update the global best and audit patience using same-background scores.

Every response path uses the same frozen cache:

```text
adjusted response = raw response
                  - frozen unselected background
                  - sum(selected-neighbor fitted responses)
residual          = adjusted response - selected-target fitted response
```

The model snapshot holds selected models and an immutable reference to the
background cache. Capturing a later iteration's cache does not mutate an
earlier snapshot. Candidate overlays change selected model responses only.

One undamped selected joint-offset solve followed by selected local refits
constitutes the fixed-point operator for that frozen background. Adjusted
entries are transient operator inputs; only final selected peeling entries are
persisted.

Joint offsets are independent per selected atom, even for identical chemical
keys. Joint, boundary, and final polish use `2 * active shape atoms + active
offset atoms` columns. An offset derivative enters only its atom's column;
the shape Jacobian retains the width derivative of that atom's offset response.
Inactive coordinates retain their own endpoint, without averaging. Unselected
responses are fixed residual/RHS terms with no Jacobian column or median chain
derivative. Candidate decoding, damping, and solver fallback do
not refresh any background parameter.

## Parameter coordinates

Trust-region steps and convergence checks use three transformed coordinates:

```text
log peak height
log width
offset-to-peak ratio
```

Accepted-movement statistics are evaluated over active optimization degrees of
freedom, not every stored atom. Every shape-active selected atom contributes
one log-peak and one log-width sample; every offset-active atom contributes
its own absolute offset-to-peak change. Fixed and quarantined coordinates are
excluded from accepted movement only. No group maximum or mixed-group mask is
computed. Unselected contributors have no active or nominal coordinates,
solver qualification, operator endpoint, or quarantine target. Background
refresh does not create a latent-movement blocker or dilute the population.

`ActiveCoordinatePopulation` stores a shape-atom index list shared by the two
shape summaries and an independent offset-atom index list. Empty populations
retain zero p99 and maximum; non-finite evidence fails closed.

Production convergence requires accepted active-DOF p99 and complete nominal-
DOF fixed-point residual p99 below `1e-4`, with solver qualification and all
orthogonal blockers clear. Maximum transformed change remains a tail diagnostic
and topology-drift metric, but is not a convergence predicate.

The strict fixed-point operator `F(S[k])` is one undamped per-atom joint-offset
solve followed by undamped own-model selected refits under this attempt's
frozen background. If production fixed/quarantine handling changes an offset,
the operator is evaluated separately without production fallback. An unavailable
unrestricted offset makes dependent shape evidence unavailable rather than
zero. The nominal population includes all three coordinates of every selected
atom, including fixed and quarantined atoms, and does not reuse the accepted
active population.

Production uses full solver qualification: active local shape refits must
report `SUCCESS`; each active atom offset requires its owning cluster's solve
to report `Converged`. Both require a full undamped, non-fallback endpoint.
A usable soft endpoint may continue through candidate selection without being
solver qualified. Historical cluster rollups and active proposal residuals are
not evaluated by the current runtime and do not define production.

The logarithmic coordinates keep amplitude and width positive when a candidate
is decoded. A candidate is invalid when its amplitude or width is not finite
and positive, its offset is not finite, or its transformed coordinates cannot
be decoded to a valid Gaussian model.

## Cluster objective

Every selected raw sample belongs to exactly one owner cluster: the cluster
containing the sample's selected target atom. Unselected contributors never own
objective rows. `src/core/detail/FittingRanges.hpp` defines three internal
constants: `kSignalDistanceMax = 1.0`, `kTailDistanceMin = 1.2`, and
`kTailDistanceMax = 2.0`, in angstroms. The signal/fit domain is the inclusive
`[0, kSignalDistanceMax]` interval; the tail objective domain independently uses
`[kTailDistanceMin, kTailDistanceMax]`. Selection flags and the sign of the
response do not remove samples from either objective domain. Direct signal
refits and alpha training use the signal interval, with the refit's existing
positive adjusted-response filter.

The ranges may overlap or leave a gap. An overlapping sample contributes to
both objective terms using their respective scales and sample counts; its
residual is evaluated once per contribution evaluation. Samples in neither
region are skipped before residual and scale validation. Full and incremental
objectives use the same membership rules, and unique-sample performance counts
use the union rather than the sum of region counts.

The default gap `(1.0, 1.2)` is excluded from both residual objective terms.
Offset estimation still uses all raw samples, and joint polish retains its
existing sample sources, so this gap is not excluded from every solver. Tail
is an objective constraint region, not a strictly held-out validation set.
`PotentialAnalysisRequest` and `RHBMTestRequest` no longer expose fitting
bounds, and neither command accepts `--fit-min` or `--fit-max`.

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

There is no width-prior term. Group posterior and prior models do not
participate in initialization or the objective. The offset-plausibility
residual floor uses the owning cluster's fixed fit scale.

The global objective weights clusters only by selected target count.
Unselected responses affect selected residuals, but add no atom weight or
independent offset-plausibility penalty. Quarantined selected blocks remain
included in this normalization:

```text
global objective = sum((cluster selected-target count / total selected-target count)
                       * cluster total)
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
transformed change are retained to break objective ties. Cluster best is not an
acceptance gate: after the previous gate (and strict local-polish improvement),
reevaluate stored parameters in the candidate neighbor environment for history.
Unavailable historical evidence preserves the old history and disables its tie
break, without rejecting the candidate. Local and rescue share this updater.

At each background refresh, re-evaluate both the previous selected state and
the retained global best under the new cache before comparing or choosing
them. Reset each affected cluster's local best history to this attempt's
previous baseline; do not compare historical numbers computed under another
background. An unavailable retained-best objective discards that best entry.
Background-only refresh does not rebuild the sampling domain or its fixed
robust scales. Refresh itself is not an improvement: audit patience uses the
candidate's strict improvement over the recomputed previous baseline, alongside
the existing domain/quarantine/radius reset conditions.

All second-stage audit tolerances use:

```text
tolerance(reference) = absolute tolerance + relative tolerance * abs(reference)
```

Progress and deterioration comparisons use
`1e-8 + 1e-3 * abs(reference)`. Strict best, tie, and polish-improvement
comparisons use `1e-10 + 1e-8 * abs(reference)`. Global gate comparisons against previous and best retain separate tolerances;
cluster best comparisons only update history.
The joint-offset IRLS objective retains its independent tolerance.

## Trust region

Base proposal trials use factors `1, 1/2, 1/4, ...`; each trial interpolates
each atom's log-peak and log-width coordinates and its own physical offset
with the same factor:

```text
C_i(t) = C_previous,i + t * (C_raw,i - C_previous,i)
```

The realized atom models are then re-encoded and measured against the trust
radius. At `t = 0`, every atom retains its own previous model; at `t = 1`,
the complete per-atom operator endpoint is recovered, subject to production
activity masks. The unselected amplitude, width, and offset remain the same
frozen background at every factor, even when selected median ordering changes.
There is no median projection before trust evaluation.

Each cluster owns one factor sequence. Every trial first constructs the
log-shape/individual-physical-offset candidate, then checks validity, trust
admissibility, guard feasibility, and the previous objective gate in that
order. Cluster-best history is evaluated only after acceptance gates pass. Trust-inadmissible trials do not run guard or objective evaluation.
Search stops when the largest transformed change is below
`kTransformedChangeTolerance`; the first passing material trial is committed
with endpoint uncertainty and its factor is recorded. Rejected trials do not
mutate objective state or polish provenance.

When every material factor is guard-infeasible for one atom's shape or offset,
that block is made locally inactive and the same function restarts the
factor sequence for the remaining blocks. This is an iterative block-isolation
loop, not a recursive candidate selection or a new outer attempt.

The polish step is limited by the radius remaining after the accepted base
movement. A rejected polish keeps the base candidate and is not backtracked.
Radius updates apply accepted objective-backtracking shrink, then retryable
rejection shrink. The first guard-feasible factor reaching the objective gate
is the accepted-shrink reference. Guard-only factor reduction does not shrink
the radius; later objective rejection followed by acceptance at a smaller factor
does. Exhausted boundary or final-audit searches keep their radius. Existing
keys never grow: accepted steps otherwise keep their radius. Shrink multiplies
by `0.5` down to `0.0625`; newly introduced topology keys start at `1.0`.
Cooperative rescue retains its existing acceptance and lifecycle policy.
The production controller uses neither actual-reduction growth nor rho.

A trust-model experiment build computes a developer-only frozen-IRLS
directional model for every material base or joint-polish trial that reaches
the objective gate,
including locally accepted and objective-rejected candidates. It freezes Cauchy
weights and objective scales at the outer previous state, applies the existing
transformed response Jacobian to the complete previous-to-candidate step, and
includes selected target and selected-neighbor deltas. The frozen unselected
background is already included in the baseline residual and has no derivative
or candidate delta. The residual surrogate uses the production fit/tail
sample coefficients and the exact previous-to-candidate offset-plausibility
penalty change. A ratio is emitted only when predicted reduction exceeds the
same progress materiality tolerance and both reductions are finite. Polish rho
uses the complete outer-previous-to-polished step; improvement relative to the
accepted base candidate is recorded separately.

The shadow status distinguishes `available`, `nonmaterial-step`,
`objective-unavailable`, `model-unavailable`, `residual-unavailable`,
`nonfinite`, `nonpositive-prediction`, and `nonmaterial-prediction`. Its
shadow action preserves objective-backtracking shrink precedence,
falls back to the current actual-only action when prediction is unusable, uses
rho bands `0.25` and `0.75`, and requires at least `0.8` boundary utilization
for shadow growth. Only the final locally accepted candidate can become
action-ready. Boundary-reconciled, rescued, globally rejected, and non-final
trial records remain available for coverage and calibration but are suppressed
from action comparison. Pre-objective validity, trust, guard, and nonmaterial
outcomes are counted in a separate candidate funnel without running the model.
None of these calculations run unless
`RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT=ON`, and no shadow result is
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
global historical-best tolerance. Member historical bests use the same
history/tie-break updater as local candidates and do not independently veto rescue. Successful rescue
promotes the rejected members without trust-radius growth or shrink. Failed
rescue is transactional: every component member retains its safe state.

Boundary correction eligibility is intentionally broader than local polish
eligibility. A component is eligible whenever at least one atom shape or offset
column remains active. Each inactive coordinate uses its atom's endpoint value.
The physical halo determines both active lists, and their union determines the
parameter atom set; no same-key atoms are added. The frozen unselected
background does not change. IRLS
objective deterioration, IRLS iteration exhaustion, and valid
non-success local estimation statuses may participate. Local joint polish still
requires full solver qualification and all cluster coordinates active. Every eligible accepted-only or rescue boundary
component receives at most one correction attempt per outer iteration.

Before a joint-correction candidate is accepted, neighbor-adjusted profiles are
rebuilt from its selected model snapshot and fixed background. Every materially
changed selected atom must pass the post-refit suspicious-profile guards.
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

The attempted fixed-order atomic component replacement was withdrawn after two
remote-cluster regressions with unavailable global baselines; see
[P2 status and evidence](second-stage-p2-structure.md). Greedy salvage remains.

## Final uncut dependency polish

After the existing stop policy selects the best-audit or latest-validated base
state, but before peeling entries or atom models are written, only a `converged`
stop with polish enabled attempts final dependency polish. Current clusters are first treated as indivisible
DSU units. They are then merged using the complete
`GraphTopology::sample_dependency_list`, without the weighted-edge threshold or
100-atom cutoff. This retains every direct selected-target/selected-neighbor
dependency; unselected contributors never connect these components.
Quarantined atom shape and offset blocks are not variables, but their fixed
models remain in every sample response and in the objective domain. A component
is skipped only when it has no active shape or offset column.

Each component is solved serially in fixed key order. Shape-active member atoms
have log-peak and log-width variables, and every offset-active atom has one
independent physical-offset variable. The last validated frozen background is
reused unchanged; merged components do not reassign samples to new background
medians. Inactive coordinates decode from their own atom endpoint. Sparse
weighted-ridge directions,
robust weights and conditioning guards are reused from boundary correction.
Each member uses a final-polish-only trust radius of `1.0` on every round,
independent of outer radius history. Up to the configured number of nonlinear
rounds is attempted. A round linearizes at its latest endpoint, and each round
measures its trust step from that endpoint. There is no additional cumulative
trust cap from the pre-polish base state.

Direction construction and candidate validation use the same selected snapshot
and frozen background, so selected deltas, residuals, and the nonlinear
objective reflect the same candidate. The final operator certificate also
uses this background; no final median refresh is performed. A component patch requires valid Gaussian parameters, post-refit
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
already converged base state unchanged.

For all non-convergence stops, the chosen base state is persisted directly,
without polish or recertification. Maximum residual remains diagnostic. Applied
polish updates audit and provenance without changing accepted iterations or the
stop reason. Diagnostics report strict-fixed-point policy and absolute-passed,
failed, error, or not-evaluated status, candidate evidence and actual application.

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
  and applies one common factor to log-shape and individual physical offsets. Trust-
  inadmissible candidates skip guard and objective evaluation. A full,
  solver-qualified endpoint below `kTransformedChangeTolerance` is stationary;
  reaching that tolerance only after factor reduction is step-limited.
- Failure is atom-block-local. An unsafe selected offset update fixes only that
  atom's offset; an unsafe shape update fixes only that atom's shape. A hard
  joint-offset failure retains that cluster's previous state while independent
  clusters can continue.
- The same selected shape/offset/hard masks parameterize local polish, boundary
  correction, rescue, and final dependency polish. Suspicious or quarantined
  selected atoms and their samples remain in coupling components, residual
  evaluation, and the objective domain. Only inactive columns are omitted,
  and their decoded values must equal the endpoint. Unselected geometry has
  no parameter blocks to freeze.
- A fixed shape may retain a guard-safe jointly estimated offset. Conversely, a
  fixed atom offset does not prevent that atom's safe shape or other atoms'
  offsets from changing. The next attempt applies the `10x`
  suspicious ridge multiplier to affected nodes.
- Active targets freeze after five consecutive observations of the same terminal
  reason. Targets remain shape blocks, offset blocks or hard-failure clusters.
  A changed reason or absent observation resets active tracking; a background
  refresh alone does not clear the consecutive-failure count.
- Frozen targets record the last attempted objective-domain revision. The
  revision advances on an applied partition/domain change or changed background
  response, not on an unchanged rebuild or merely queued topology. Each new
  revision permits one retry; retries are transient, not a third lifecycle.
  Background updates can consequently permit a retry every iteration.
- Retry keys use minimum radius `0.0625` and `10x` ridge. Non-retrying Frozen
  masks retain priority; overlapping retry targets do not unlock coordinates
  still frozen by another target. No cooldown, attempt limit or Exhausted state
  remains. Frozen targets awaiting a domain change do not hold audit patience
  open as a scheduled recovery; active failure tracking and transitions still do.
- Recovery needs no affecting failure and all required target coordinates active.
  It requires a finally accepted material change, or a complete, guard-safe,
  solver-qualified nonmaterial unrestricted endpoint. Missing evidence cannot
  release a target. A failed retry records the revision and remains Frozen.
- Newly Frozen targets affect only the next proposal. Lifecycle updates never
  rewrite the audited model/provenance, so quarantine fallback re-audit is removed.
  Immediate block isolation and solver fallback remain. Final activity still
  fixes unresolved targets; transition convergence blockers remain.
- Rejected-cluster debug output distinguishes failures before objective
  evaluation from objective rejection. Pre-objective failures report their
  proposal reason, radius, available step norm, and `objective =
  not-evaluated`; `objective-unavailable` is reserved for an objective that was
  actually attempted but could not be calculated.

An accepted state's adaptive topology rebuild becomes the next hysteresis
reference when applied. An unchanged partition retains objective scales,
trust radii, and solver workspaces. A changed cluster/sample/boundary mapping is
queued without switching the running attempt's partition or background. At the
next iteration boundary, build the new partition's complete finite cache and
apply them together. A failed cache cannot commit a partial background or write
unvalidated results.

A partition change initializes new fit/tail scales, cluster baselines, and
solver workspaces from the latest accepted selected state with the new cache.
Both that state and any retained global best are scored on the new domain.
Exact cluster keys retain trust radii; merged or split keys start at the initial
radius. Audit patience is cleared, and the domain-change blocker covers both
pending and newly applied partitions. A pending change therefore cannot be
skipped by convergence; if the outer limit has been reached, finalization uses
the last validated partition and frozen background instead. Background-only
refresh leaves sampling-domain/scales/workspaces intact and adds no separate
movement blocker. Scores from different backgrounds or domains are never
compared directly.

## Global audit and stopping

The global audit uses the fixed per-cluster fit/tail scales and retains the
earliest state that improves the best objective beyond the strict tolerance.

One internal `ConvergenceCertificate` is the sole source of convergence truth.
`ProductionConverged()` requires solver qualification, accepted active-DOF p99
below `1e-4`, a complete nominal-DOF operator, nominal fixed-point residual p99
below `1e-4`, and clear orthogonal blockers. Fixed and
quarantined coordinates are excluded only from the accepted population; they
remain in the nominal operator population. An empty accepted population passes
its percentile check vacuously, but an all-fixed state still needs qualified,
complete, sufficiently small nominal operator evidence. Offset qualification
is checked per active atom against its owning cluster's solve status. An
unavailable endpoint makes the operator incomplete instead of substituting
the previous state as a zero residual. `StrictOperatorPassed()` reuses the same certificate for converged
final-polish certification without the accepted-movement or orthogonal-blocker
terms.

The stage stops on the first applicable condition:

- no valid initial seed is available for every selected atom;
- accepted active-DOF p99 and complete nominal-DOF fixed-point residual p99 are
  both below `1e-4`, every active coordinate is solver-qualified, all clusters
  are accepted, and no orthogonal blocker is present;
- `kLocalFittingAuditPatience` accepted iterations produce no strict candidate
  audit improvement over their same-background previous baselines;
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

## Final state application and final group fitting

After the stopping policy selects the final validated state and any certified
final polish, the stage captures its selected MDPDE models together with the
last validated immutable background. This is the actual best-audit or latest
validated state chosen for application, not an operator endpoint. Each selected
atom keeps its own final offset. Final polish, certificate, audit, and peeling
share this exact background; persistence neither averages selected offsets nor
recalculates a background median. OLS/MDPDE uncertainty and polish provenance
follow the existing application policy.

For every selected atom, the stage then rebuilds its persistent peeling
sampling entries from the raw entries:

```text
peeling response = raw response
                 - sum(final selected-neighbor MDPDE responses)
                 - last validated frozen unselected background
```

Selected neighbors use their final atom-level MDPDE models. Every unselected
contribution is already included in the frozen per-sample cache, with its last
validated global-median background model. The calculation preserves original
sampling-point order and metadata. Only selected local Gaussian results are
written; no atom-local analysis entry is created for an unselected contributor. The selected
results and rebuilt entries are persisted together with
`ApplyAtomLocalSecondStageResult`. `detail::RunSecondStageIterations` then returns
without training `alpha_g` or running group fitting.

After the second stage returns, the workflow directly runs
`RunGroupAlphaTraining(model, options)` and
`RunGroupPotentialFitting(model, options)` without further local fitting.
The group fit consumes the `FittingStage::Second` local models, alpha-r values,
and atom-level peeling snapshot written during second-stage finalization.
It preserves both local fitting results and the peeling entries, and produces
one group result and one optional posterior/outlier result per atom, stored
independently of local fitting stages. Group fitting and group-result access
do not accept a fitting-stage parameter.

## Performance architecture

The fitting context stores selected atom indexes, flattened selected
sample-neighbor edges, unselected per-sample distance lists, prepared selected
local designs, and profile-radius ordering. Geometry and
designs are fixed for the stage. Quarantine changes only selected activity masks.

Each outer iteration builds one immutable unselected background response cache,
selected model snapshots, one selected adjusted-response cache, and one
selected-row residual/objective baseline. Model snapshots hold immutable shared
references to their background so subsequent refreshes cannot alter them.
Cluster candidates are represented by atom-local state patches. Candidate
evaluation overlays direct selected model deltas on the previous
state and evaluates the objective as baseline plus the changed sample and
offset delta. Selected local refits use one frozen post-offset own-model
snapshot. For a boundary-component guard, affected sample
IDs are sorted and deduplicated so a boundary sample is recomputed once. The
boundary dependency and deterministic accepted/rescue-induced components
are derived from the current partition and rebuilt with adaptive topology.
Without a multi-cluster accepted or rescue component, candidate selection stays
on the existing fast path.

Joint-offset, joint-polish, and boundary-correction solvers retain their sparse
pattern analysis while its sparsity pattern is unchanged. Unselected background
responses have no sparse design columns or median chain derivatives.
Boundary-correction workspaces use
the deterministic shape-active/offset-active/sample signature as their key and
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
| `Atom A/Q` | Selected atoms without quarantine targets / selected atoms covered by current quarantine targets |
| `Cluster A/R` | Accepted / rejected candidate clusters |
| `Polish E/A/R/S` | Eligible / accepted / rejected / skipped polish clusters after component reconciliation and the final global guard; `E = A + R + S` |
| `Suspicious` | Atoms with at least one shape, offset, or hard-failure block fixed in this attempt |
| `dMax A/F` | Maximum transformed change in the accepted/fixed-point operator state; accepted is `-` on an all-rejected attempt |

Objective-domain startup diagnostics report selected-target weights, cluster and unique
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
direct-interface/shape-active/offset-active/parameter counts, suspicious count,
solver status,
damping, maximum normalized trust step, strict-improvement reference/candidate
objectives, acceptance result, and endpoint-fallback outcome. An all-rejected
debug record reports unified trial dispositions, terminal category, radius
action, and stop classification. Operator-assessment debug records report each
affected atom's reason, margin, and fixed shape/offset/hard block. Completion
warnings report cumulative quarantine
entries, releases, failed domain retries, and unresolved targets. The existing
probation counter positions retain their schema and now count domain retries. Convergence
and summary messages finish the active progress line before normal line output.
Joint-candidate member guards emit `Joint candidate objective rejection: schema=1`
Debug records alongside their component summary. Each record identifies the
candidate source (endpoint, joint correction, backtracking, rescue variants, or
final polish), component keys, evaluation sequence, factor, and polish round.
Only the first failing member is recorded; evaluation retains its existing
short-circuit order. Full member keys use internal selected-atom indexes.
Previous/best/candidate objectives use fit/tail-weighted/offset/total order,
with `max_digits10` precision. Enabled member gates report the reference,
candidate-minus-reference delta, absolute/relative tolerance, calculated
`tolerance`, and inclusive upper bound `reference + tolerance`. The tolerance
is `absolute + relative * abs(reference)` using `kObjectiveProgressTolerance`.
Best gates are `not-checked` in cooperative rescue and final polish; absent
values are `unavailable`. Member failure outcomes distinguish previous, best,
both, missing evidence, and nonfinite evidence. Records with `member=none`
distinguish a downstream global failure after member checks from a final-polish
global failure before member checks. Local rejection summaries remain local:
`rejected-by=none` does not imply that a later joint candidate passed.
Member rejection records reuse the evaluated objectives. Debug mode additionally
tracks cluster-best provenance and performs diagnostic-only reference comparisons;
Info and quiet modes create no provenance snapshots or extra evaluations.

`Cluster best source: schema=1` records best initialization, partition/background
resets, candidate improvements, and step tie-break updates. IDs combine attempt,
full cluster key, and per-key update sequence, independent of worker completion
order. Each record retains the effective model snapshot (including candidate
overlays), shared immutable background/domain, and contribution sample refs.
`predecessor` links updates and same-key resets. `retained` distinguishes provisional
updates from the source finally published for the attempt; `Cluster best publication`
identifies that source even when it originated in an earlier attempt. Copying or
rolling back a cluster objective state also copies or restores its source.

When a joint member fails with a best reference, `Cluster best comparison: schema=1`
links that source and reports stored, historically reproduced, and re-evaluated best
objectives in fit/tail-weighted/offset/total order. Re-evaluation holds historical
member parameters fixed and substitutes either the current previous-state or actual
candidate's external models. A separate ordered decomposition substitutes current
domain/sample refs, then frozen background, then external models. These deltas are
order-dependent diagnostics, not independent causal contributions. Original and
re-evaluated gate results use the existing progress-tolerance helper; only the
original gate participates in acceptance. Each record counts its additional
objective and residual-sample evaluations separately from production counters.

`Cluster best environment: schema=1` compares actual sample refs, ownership,
fit/tail masks, scales, normalization counts, and background responses. It lists
historical/previous/candidate parameters for the member and changed contributors
of those samples. Historical model/domain/background snapshots are immutable;
re-evaluation uses an independent sparse residual baseline and does not change the
live context, solver caches, best values, or convergence decisions. All numeric
records use `max_digits10`. Missing provenance is `unavailable`; incompatible keys
or snapshots are `not-comparable`.

Frozen-background debug diagnostics report the selected target serial,
global-median amplitude/width/offset, and effective background sample count.
No second-stage diagnostic reports a chemical group. Contributor-refit, hard-edge,
and hard-closure-overflow diagnostics have been removed. The developer-only
trust shadow retains its legacy `unselected-dependencies=0` trace field.

The current convergence trace is schema 10. It serializes only try/accepted
iteration and atom/quarantine counts, active and nominal populations,
accepted/operator p99 and maximum, the six-bit production certificate, and
the four orthogonal blockers. The current analyzer accepts only schema 10;
frozen schema-9 baselines remain historical data. Atom audit records separately
use schema 2 with serial/amplitude/width/offset and no `group` field. The analyzer
also reads legacy atom schema 1, where `group` remains required; it never inserts
a synthetic group into schema-2 records.

Adaptive rebuild diagnostics use a distinct
`Adaptive local-fitting topology rebuild` record so the one-time initial
coupling and atom-cutoff summaries remain one-time initial records. The cutoff
record uses `Local-fitting atom cutoff: atoms=N, limit=100, clusters=C,
max-atoms=M, cutoff-edges=E.` and counts selected atoms only. Threshold-sensitivity
summaries describe pre-cutoff connectivity; the formal component summary
describes post-cutoff connectivity. Each adaptive rebuild record reports the
accepted iteration, drift or interval trigger, maximum drift, old/new cluster
and boundary-sample counts, added/removed adjacency edges, and whether the
partition changed and is pending. The rebuild record reports no immediate
objective-domain reset; a distinct objective-domain log records application
at the next iteration boundary. Non-quiet runs also show a
`Rebuild local-fitting coupling topology` percentage progress bar with the
same sample-based work accounting as the initial topology build.

Finalization emits a separate `Final dependency polish` record with aggregate
component, atom, parameter, round, acceptance/fallback, objective-before/after,
residual-safety policy/status, application, base/candidate solver and operator
evidence, residual p99/maximum, and elapsed-time values. Debug logging adds one
record per component, including symbolic-analysis and suspicious-candidate
counts. These records are emitted before the existing stable final summary.

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

## Workspace verification (2026-09-04)

The implementation baseline was the clean revision
`b2c0bc5c9a64afea94d131939ee2b0eb437a9613`. The verified workspace is that
revision plus the current uncommitted atom-cutoff statistics cleanup. The
cutoff's completed atom union-find now supplies the component summary directly,
without rebuilding connectivity from adjacency. The cutoff summary retains
only its limit and cut-edge count; atom count and component statistics are read
from the existing topology and formal summary. The 100-selected-atom limit,
edge ordering, threshold/hysteresis behavior, frozen backgrounds, acceptance,
convergence, persistence, and diagnostic formats are unchanged. Verification
used AppleClang 21, RelWithDebInfo, system dependencies, OpenMP 5.1 AUTO, and
disabled UMAP/ROOT:

- `tests_all`, `rhbm_tests_core_estimator`, and `rhbm_tests_data_runtime` passed.
- Full CTest passed all 18 entries, including convergence analyzer/runner,
  fold-168 runner contracts, smoke, and serial/parallel determinism.
- `lint_all`, including repository lint and install-consumer smoke, passed.
- The developer-only trust-model experiment ON build passed the focused tests,
  all 18 CTest entries, and lint/install-consumer smoke. Its defense suite
  contains 101 cases. The experiment was restored to OFF, followed by a normal
  rebuild, focused tests, full CTest, lint, and another numerical capture; the
  normal defense suite contains 100 cases.
- Source counts remain 729 `TEST`/`TEST_F` and 17 `TEST_P`. No repository test
  file or case was added; the fold-168 runner retains its eight Python cases.
- Existing graph tests cover 100/101-atom boundaries, isolated atoms, small
  custom limits, strong-edge priority, equal-weight canonical tie-breaking,
  edge/active-index permutations, complete internal edges, empty/invalid input,
  capped binary fallback, post-cutoff hysteresis, and pre-cutoff sensitivity
  versus post-cutoff component summaries. Existing cases now compare the single
  component summary against actual partitions, including component count,
  maximum size, ratio, and empty-graph zero values, and assert the exact stable
  atom-cutoff log record.
- A connected 101-atom chain plus two remote selected atoms verifies actual
  cutoff and 101-atom boundary reconciliation, serial/parallel agreement,
  unchanged intensity-scale tolerances, and remote-cluster improvement. The
  chain uses slightly nonuniform spacing to avoid repeated-weight rounding ties;
  the production weight comparisons and numerical tolerances are unchanged.
  Uncut dependency components can exceed the 100-atom topology limit.
- The shared-background fixture connects selected atoms through physical
  sample coupling rather than matching residue labels. Relabeling selected and
  unselected chain IDs, sequence IDs, and residue names while holding selection,
  elements, geometry, samples, initial local state, and `alpha_r` fixed leaves
  second-stage models, peeling, stop reasons, topology records, and convergence
  evidence unchanged. The chemical-key independence coverage is retained.
- Frozen backgrounds, shared contributors, fixed/quarantined median pools,
  hydrogen exclusion, deduplication, best-audit rescoring, adaptive partitions,
  final polish/peeling, unselected non-persistence, and full workflow regressions
  passed.
- Eighteen existing fixture configurations were captured before the change.
  All 18 matched byte-for-byte afterward in OLS/MDPDE models, uncertainty,
  peeling, completion/stop records, and convergence evidence, using hexadecimal
  floating-point serialization. Initial coupling/cutoff, threshold-sensitivity,
  adaptive-topology, and frozen-background records also matched exactly, as did
  captured warnings. All 18 matched the clean baseline again after restoring
  OFF. No fixture geometry was changed or configuration excluded, including the
  shared-background and boundary fixtures changed in the preceding topology
  revision. Timing and performance counters were not compared.
- Fold-168 baseline/report schema remains 6, with 168 selected atoms, a 100-atom
  limit, and at least two topology clusters. Input hashes, reference quality
  metrics, quality tolerances, and the 25-accepted-iteration gate are unchanged.
  Runner tests passed, including rejection of legacy residue logs and schema 5.
  The external fold-168 benchmark was not run: its model/map inputs are not
  configured in this workspace.
- Reverse searches found no residue lookup, residue pre-merge, residue-count
  cutoff, chemical keys, shared-offset merge, or group-median refit in the
  second-stage production path. Public headers, fitting options, CLI, database
  schema, selection flags, stage flow, third-stage estimators, and standalone
  empty-selected behavior are unchanged. `git diff --check` passed.

The ROOT-disabled build retains existing unrelated painter warnings; those
files were not changed. Only this
workspace-verification record was updated in the normative document; its
algorithm description, the existing Notion algorithm page, and historical audit
documents are unchanged.


## Readability refactor verification (2026-09-06)

The paired baseline is revision `a285a63a29ff134f346b3697743d5cd2ce0032a5`.
Its executable and shared library were preserved before implementation. Both
build variants used Debug, system dependencies, OpenMP AUTO (5.1), UMAP, and
ROOT; the trust-model experiment was tested both OFF and ON.

- The six responsibility modules and diagnostic ownership described above are
  implemented. No public estimator options, numerical thresholds, stop rules,
  logging schemas, or persistence interfaces changed.
- Both variants built `rhbm_tests` and passed the estimator, algorithm, math,
  and HRL CTest entries.
- Repository lint and `git diff --check` passed. No test file or case was added:
  the modified defense source retains 101 declared cases, including its
  conditional experiment case.
