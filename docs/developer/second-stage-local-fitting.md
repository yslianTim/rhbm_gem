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

Implementation ownership and candidate acceptance references are specified below.
Candidate selection uses a private builder and one transaction publication;
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

Shared fitting support lives in `src/core/detail/gaussian_fit/`; second-stage
implementation lives in `src/core/detail/second_stage/`, with observation and
logging under its `observation/` subdirectory. All retain the existing
`rhbm_gem::core::detail` namespace.

The internal implementation follows the iteration sequence rather than exposing
all second-stage services through candidate selection:

| Module in `src/core/detail/second_stage/` | Responsibility |
| --- | --- |
| `IterationProcess` | Initialization, frozen-background and pending-partition boundaries, convergence and stop decisions, final certification, and persistence |
| `IterationResult.hpp` | Outer attempt summary and stop reasons |
| `ConvergenceCertificate` | Active-coordinate summaries, shared fixed-point operator evidence summaries, solver qualification, and convergence certificate |
| `CandidateEvidence.hpp` | Candidate decision evidence, pre-objective failure reasons, and boundary accepted sources |
| `CandidateState.hpp` | Selection data, boundary decisions, polish progress, and trust-radius update records |
| `IterationProposal` | Joint offsets, local shape refits, fallback, and unrestricted fixed-point operator evidence |
| `CandidateTransactionLocal.cpp` | Builder-owned per-cluster candidate search, local joint polish, and trust-radius control |
| `CandidateEvaluation` | Typed references with scopes only where policy differs; separate local/boundary results and complete global candidate acceptance with original gate ordering; no history inputs or results |
| `CandidateTransaction` | Per-key provisional selection, final classification, staged quarantine, and consuming publication of validated results |
| `CandidateTransactionBoundary.cpp` | Shared normal/cooperative component evaluation and application, complete-selection audit/salvage |
| `DependencyPolish` | Final uncut-component candidate generation and salvage policy; validation delegates to `CandidateEvaluation` |
| `ComponentAssembly` | Ordered patch application and excluded-component trial assembly |
| `ObjectiveEvaluation` | Objective domains, full and incremental evaluation, tolerances, previous objectives and the production global best |
| `SuspiciousUpdate` | Profile baselines, suspicious assessments, coordinate activity, failure masks, and candidate/polish guards |
| `Quarantine` | Active/Frozen failure tracking, domain retry, and next-iteration activity |
| `observation/SecondStageObservation` / `observation/SecondStageDiagnostics` | Passive session, bounded worker buffers, actual gate and lifecycle evidence, scalar progress diagnostics |
| `observation/SecondStageLogging` | Progress and certificate output, graph/objective diagnostics, and read-only performance formatting |
| `observation/PerformanceCounters` | Atomic counts, phase timings, and current/retired solver workspace totals; publishes once at scope exit |

`GaussianModelOperations` and `PreparedLocalGaussianFit` provide shared model
operations and prepared designs in `gaussian_fit/`, alongside `FittingRanges.hpp`.
`SecondStageState`, `CouplingGraph`, and `JointFitting` retain the second-stage
state/residual representation and seed selection, graph construction and topology
drift, and solvers. `IterationProcess.hpp` declares only `RunSecondStageIterations`;
convergence types and the outer attempt result have their own headers. Operator
evidence summarization serves production convergence and final certification.
Observation copies the resulting certificate without replaying the operator. Seed diagnostic records belong to
`observation/SecondStageDiagnostics.hpp`.
`CandidateEvidence.hpp` supplies diagnostics and candidate evaluation without
including the complete selection. `CandidateState.hpp` owns selection data.
`CandidateTransaction.hpp`
owns trust-radius/backtracking control declarations and `CandidateSelectionInputs`;
the transaction implementation is distributed across
`CandidateTransaction.cpp`, `CandidateTransactionLocal.cpp`, and
`CandidateTransactionBoundary.cpp`. The public Gaussian estimator workflow uses
internal fitting-range constants; `FitOptions` does not expose radial bounds.

`CandidateSelectionInputs` contains read-only algorithm inputs. The private
`CandidateTransactionBuilder` owns working activity, state and provenance;
boundary operations are its private methods. Evaluators return numerical decisions
without accepting or returning per-cluster history. Observation receives actual
decisions without a decision return value.
Solver workspaces, counters and observers remain mutable working resources.
The observation pointer in `CandidateSelectionInputs` is non-owning.
Boundary candidate and correction references borrow only the per-key samples,
objective domain, previous member objectives, optional best objective, component,
reference objectives, and performance counters. Correction also borrows its
endpoint and strict-improvement reference. Context comes from the candidate
overlay; evaluators do not include the transaction header or receive solver,
proposal, provenance, or observation-session state. An unavailable precomputed
correction objective is reused as unavailable evidence, never recomputed.

`BoundaryObservationScope` and `JointCandidateObservation` record actual endpoint,
correction, backtracking and rescue outcomes. Each event records the gate at its
original decision point; a later rejected correction cannot change an earlier
endpoint record. Workers keep at most five abnormal details and merge in fixed
stage/key/trial order. Counters retain all outcomes. Quiet mode and missing
sessions do not affect numerical decisions.

`RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT` applies the same private definition to the
library and tests through one helper. ON does not require testing. Extra payloads
exist only in enabled, non-quiet Debug sessions. Allocation, collection and writer
failures disable extra recording for that session without entering numerical
fallback branches. Debug's existing solver scheduling condition is independent
of this option. See the [audit contract](second-stage-audit.md).

`Finish` stages next-iteration quarantine without publishing it. Consuming `Commit`
publishes quarantine, applies trust updates and copies keys, then publishes the
accepted model/provenance (or restores the previous state when all candidates are
rejected), and finally records the committed outcome. It returns `CandidateCommitResult`;
the runner moves the returned keys and trust updates into `IterationResult`.
Transaction code does not depend on the outer result or entry header.
`RunSecondStageIterations` owns a non-copyable `SecondStageObservationSession`;
`SecondStageContext` contains only atom sampling/design data and the immutable
frozen background. The session owns a bounded audit payload and scalar
iteration/final-polish progress summaries. It owns no models or replay snapshots.

`CandidateDecisionEvidence` contains objective references, factors, failure
classification, and the invalid/guard/objective-rejection counts used by search
control. `StabilizationTerminalEvidence` retains ordered guard/invalid/exhausted
reasons and the affected atom/mode. Quarantine reads accepted evidence first,
then rejected evidence in rejection-event order; objective exhaustion still
does not constitute quarantine failure evidence. Accepted shrink compares the
accepted factor with the first objective-evaluated factor directly.

Candidate references may point to production-owned member-best patches, while
results contain no observer history tokens. `CandidateSelection` retains keys, state/provenance, decision evidence,
and lightweight boundary decisions. `IterationResult` retains accepted/rejected
keys, radius updates, progress state and stop information; audit-patience reset
uses rejected keys directly. Bounded output-only event records live in observation sidecars.
Optional boundary observation scopes keep endpoint/correction/backtracking
record associations independently of numerical candidate results.

Quarantine is staged on a copy and changes only next-iteration activity,
without modifying the audited state or requiring fallback re-audit. The builder then freezes a read-only `CandidateTransaction`.
Its consuming commit publishes validated state, provenance, quarantine
and radius updates, including the required rejection updates on all-rejected
attempts. Only after production publication does it record the commit outcome. Convergence uses the committed candidate and the retained previous
state. The unrestricted operator evidence remains separate from these
production restrictions.

Rejection rolls back affected models and provenance without discarding
rejection-driven radius shrink, fixed activity, failure evidence, or quarantine
transitions. An all-rejected attempt restores the previous model and still
publishes those lifecycle updates. Suspicious-failure evidence is sampled before
quarantine publication. Observer publication follows production commit and does
not supply a decision or undo production publication.

`SecondStageContext::atom_list` is accessed directly. Model snapshots capture
both the selected Gaussian models and the immutable background in effect when
they are built. `BuildSecondStageModelSnapshot` accepts a full fit state, a
patch-backed state view, or existing model snapshots; callers do not need a
separate model-projection step. The topology drift reference retains only model
snapshots, while accepted and best-audit states retain complete local results.
Patch/view and residual-overlay representations continue to avoid full-state
materialization during candidate evaluation.

The session's `IterationDiagnostics::proposal_maximum_transformed_change` field measures the constrained
proposal's maximum movement. It is distinct from the nominal operator residual
in `ConvergenceCertificate`. Maximum/population measurements reside in
`ConvergenceDiagnostics`; the existing progress label and audit schemas are
unchanged.

### Candidate acceptance references

`EvaluateCandidate` uses `(candidate, reference)`, with a separate optional
observation scope for boundary and final-polish records. Typed local and boundary
policies share numerical primitives while retaining their own acceptance rules.
Proposal construction, factor retry, correction generation, and salvage remain
orchestration responsibilities; failures retain their short-circuit order.

| Phase | Evaluation order and reference |
| --- | --- |
| Local search | Construction/validity and the nonmaterial shortcut precede trust, guard, and the previous-objective gate. Passive observation follows the production decision. |
| Local polish | Existing solver feasibility/trust checks, then strict objective improvement against the accepted local endpoint. |
| Ordinary boundary | Member previous gates in existing key order, followed by the combined objective. |
| Boundary correction | Suspicious-polish guard, raw objective evidence, member/combined acceptance, then strict improvement against the original correction reference. The raw evidence is reused, including when unavailable. |
| Cooperative rescue | Tolerated member deterioration, combined previous/best acceptance, and strict global improvement. |
| Global selection audit | When triggered by the [conditional selection audit](#conditional-selection-audit), affected-sample union and complete-state previous/best gates, followed by the outer salvage policy. |
| Final polish | Validity, suspicious-polish guard, strict global improvement, and member non-regression against the base. Converged finalization separately requires strict operator recertification. |

### Conditional selection audit

These call conditions were checked against `develop` commit `7500a45eb2fc20ff9d04768789883ef4b1d2eb65`
on 2026-09-14: [boundary reconciliation and audit/salvage](https://github.com/yslianTim/rhbm_gem/blob/7500a45eb2fc20ff9d04768789883ef4b1d2eb65/src/core/detail/second_stage/CandidateTransactionBoundary.cpp#L664)
and [runner publication and post-commit scoring](https://github.com/yslianTim/rhbm_gem/blob/7500a45eb2fc20ff9d04768789883ef4b1d2eb65/src/core/detail/second_stage/IterationProcess.cpp#L599).
This is production objective acceptance, separate from optional decision-recording
observation. Making these calls unconditional would change the algorithm.

`ReconcileSelectedBoundaries()` builds the ordinary boundary component list
from the currently selected keys before reconciliation. Each component has at
least two connected keys; remote singleton patches do not make the list
nonempty. It then evaluates `previous_audit_objective` from the current domain
and residual baseline, regardless of whether this component list is empty.

| Ordinary component list | Previous audit objective | Before cooperative rescue | Cooperative rescue |
| --- | --- | --- | --- |
| Nonempty | Available | Run ordinary reconciliation, then call `AuditAndSalvageFinalSelection()`. | Attempt eligible components. |
| Empty | Available | Skip ordinary reconciliation and the first audit/salvage call. | Still attempt eligible components. |
| Nonempty | Unavailable | Run ordinary reconciliation, then reject all remaining selections as exhausted; do not call audit/salvage. | Skip. |
| Empty | Unavailable | Skip ordinary reconciliation and the first audit/salvage call; this branch adds no global rejection. | Skip. |

After that branch, the post-rescue audit/salvage call occurs only when previous is
available and `ReconcileCooperativeComponents()` returns an accepted result.
Rescue eligibility is independent of whether the ordinary component list was
nonempty. An accepted rescue still has to survive this subsequent audit/salvage
to be retained at commit.

Whenever called, the selection audit covers the complete assembled state,
including remote singleton patches, using the affected-sample union and the
existing previous bound and retained global-best bound when a best exists.
With no ordinary components and no accepted rescue, local selections can
proceed to commit without a complete-selection acceptance gate.

Calling audit/salvage does not guarantee an objective result: an empty selection
returns no objective, including when salvage removes every selected unit. An
unsalvageable remaining selection is rejected as exhausted and clears the
result. Thus an absent `final_audit_objective` alone does not establish whether
audit/salvage ran. After these conditional paths, `MaterializeSelection()`
finalizes accepted/rejected lists once. The runner's subsequent objective
reuse or recomputation is described in [Global audit and stopping](#global-audit-and-stopping).

### Shared component infrastructure

Outer boundary and final uncut component builders share DSU participant merging
and root-to-key collection in `CouplingGraph`. Their participant sources,
minimum component sizes, halo expansion, selected-owner filtering and output
ordering remain distinct. Outer boundary components require two accepted keys;
final components may contain one multi-atom key.

`ComponentAssembly` applies borrowed patches in input order, skipping null entries
and an optional excluded position. Outer selection applies patches to its
builder-owned state; final polish constructs complete states from its base and
retained component patches. The module owns no selection, provenance, radius,
quarantine or diagnostic state.

Each stage owns its audit and removal loop directly, with outer calls governed by
the [conditional selection audit](#conditional-selection-audit). Only a failed
initial audit with a nonempty selection builds the outer rejection ranking, once.
Outer considers units that do not independently strictly improve previous, with
exact-delta scoring, unavailable evidence scored as infinity, worst-first ordering,
lexical key tie-breaking, previous/best gates, and exhausted fallback. It audits
after each removal and stops when accepted or empty. Final polish tries each accepted component for removal
using full-state audits and selects the best available single removal each round,
retaining first-position tie-breaking and strict improvement over the base. A selected
final removal reuses its already computed objective, without an extra audit.
With no accepted components, final polish reuses the base objective.
Ordinary/cooperative sweep ordering and polished-state recertification are
unchanged.

`FinalDependencyPolishResult` contains only state, objective and acceptance.
Component/trial details and timing live in the session's final-polish sidecar.
Acceptance compares the last successfully improved numerical objective with
the base objective; no successful update remains unavailable evidence. Neither
component acceptance nor final recertification reads the diagnostic copy.

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
   and post-refit feasibility guards, and then the previous objective gate.
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
   eligible for commit. When the ordinary component list is nonempty, run its
   global audit/salvage only if the previous audit objective is available;
   otherwise reject the remaining selections as exhausted. With previous
   available, attempt cooperative rescue even if there were no ordinary
   components. Build maximal eligible boundary components from the
   safe accepted state and the best finite objective-rejected proposal retained
   for each rejected cluster. Components containing at least one such proposal
   receive the same endpoint, active-column joint correction, and common-factor
   checks. Cooperative rescue permits a member to deteriorate only within
   `1e-8 + 1e-3 * abs(previous member objective)`; the component and tentative
   assembled global objectives must both improve strictly, and the global
   historical-best tolerance still applies. A failed correction retains a valid
   rescue endpoint, and any failed rescue atomically leaves the safe accepted
   state unchanged.
9. Only when previous is available and cooperative reconciliation accepts a
   result, call global previous/best audit/salvage on the complete assembled
   state. Both call sites follow the [conditional selection audit](#conditional-selection-audit).
   If tolerance-level deterioration from independent reconciliation units causes
   aggregate rejection, remove non-improving units from worst to best until the
   first passing subset is found. If only strictly improving units remain but the
   historical best gate still fails, roll back the attempt without weakening the
   gate.
10. Stage Active/Frozen quarantine for the next iteration without changing
    the audited state. Publish quarantine and trust-radius updates together. An accepted state's adaptive rebuild can queue a
    new partition for the next attempt; pending changes block convergence.
    For an accepted result, update the global best and audit patience using
    same-background scores, reusing `final_audit_objective` or recomputing it
    when absent. This post-commit scoring is not an acceptance gate; all-rejected
    attempts return before it.

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

Shape and offset endpoint availability are tracked separately for each selected
atom. A non-hard-failure, valid joint-offset model can supply an available soft
endpoint without solver qualification. Missing shape evidence makes both shape
residual coordinates infinite; missing offset evidence makes its coordinate
infinite. `operator_complete` requires every nominal atom's availability masks;
transformed finiteness and solver qualification are separate checks. The p99
predicate applies independently to log peak, log width, and offset coordinates.

Production requires every nominal shape solve to report `SUCCESS` and every
nominal offset solve to report `Converged`, using the status of that endpoint's
actual solve. Ordinary accepted active coordinates retain their qualification
checks as well. Both require a full undamped, non-fallback endpoint.
A usable soft endpoint may continue through candidate selection without being
solver qualified. Historical cluster rollups and active proposal residuals are
not evaluated by the current runtime and do not define production.

The logarithmic coordinates keep amplitude and width positive when a candidate
is decoded. A candidate is invalid when its amplitude or width is not finite
and positive, its offset is not finite, or its transformed coordinates cannot
be decoded to a valid Gaussian model.

Joint-offset and joint-polish conditioning normalize design columns and use
LDLT `min(D)/max(D)` of the normalized Gram matrix as a conditioning proxy before
ridge. A ratio at or below `1e-8` requires a ridge multiplier floor of `10`.
Empty/invalid columns, failed factorization, and nonpositive/nonfinite pivots
return the zero sentinel and require the guard. Other safeguards may increase
the multiplier further. This proxy is neither a singular-value condition number
nor an effective-rank certificate.

## Cluster objective

Every selected raw sample belongs to exactly one owner cluster: the cluster
containing the sample's selected target atom. Unselected contributors never own
objective rows. `src/core/detail/gaussian_fit/FittingRanges.hpp` defines three internal
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

Candidate scoring uses only production objective references. Extra observation
copies existing evidence; it performs no historical re-evaluation or snapshot
retention. The production global `best_audit_state` remains part of acceptance,
patience and final-state selection. See the [audit contract](second-stage-audit.md).

At each background refresh, re-evaluate both the previous selected state and
the retained global best under the new cache before comparing or choosing
them. An unavailable retained-best objective discards that best entry.
Background-only refresh does not rebuild the sampling domain or its fixed
robust scales. Refresh itself is not an improvement: audit patience uses the
candidate's strict improvement over the recomputed historical best, alongside
the existing domain/quarantine/radius reset conditions.

All second-stage audit tolerances use:

```text
tolerance(reference) = absolute tolerance + relative tolerance * abs(reference)
```

Progress and deterioration comparisons use
`1e-8 + 1e-3 * abs(reference)`. Strict best, tie, and polish-improvement
comparisons use `1e-10 + 1e-8 * abs(reference)`. Global gate comparisons against previous and best retain their existing tolerances.
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
order. Trust-inadmissible trials do not run guard or objective evaluation.
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
keys never grow: accepted steps otherwise keep their radius. The controller stores
an unsigned shrink level, `0..4`, rather than a floating-point radius. `GetRadius`
returns the exact table value `[1.0, 0.5, 0.25, 0.125, 0.0625]`. Shrink increments
the level until 4, after which it reports saturation; minimum-radius recovery sets
level 4 directly. Reconciliation preserves surviving keys, drops removed keys and
starts new keys at level 0. Update ordering, repeated keys, exhausted-key exclusion,
missing-key exceptions and changed/saturated reporting are unchanged.
Cooperative rescue retains its existing acceptance and lifecycle policy.
The production controller uses neither actual-reduction growth nor rho.

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

After the accepted-only pass and its conditional audit/salvage, cooperative
rescue runs only with an available previous audit objective, including when
there were no ordinary components. Every rejected cluster that produced a finite
objective proposal retains its lowest-objective trust-region patch in memory.
The rejected proposals and safe accepted clusters form maximal eligible
boundary components. Suspicious and hard-failure atoms are not removed from a
component; their inactive parameter blocks remain fixed while active neighbors
can supply correction columns. A cooperative endpoint may contain a member
whose objective is slightly worse than its previous value, but only within the
normal progress tolerance. The endpoint or its joint/backtracked replacement
must strictly improve the component audit, and the tentatively assembled state
must then strictly improve the previous global audit without violating the
global historical-best tolerance. Successful rescue
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

When the [conditional selection audit](#conditional-selection-audit) is triggered,
the complete assembled state must pass the unchanged global previous/best gates.
On aggregate failure, independent components
and singleton clusters are scored by their exact global objective delta. Only
non-improving units are removed, from worst delta to best, with the full audit
recomputed after each removal. The first passing subset is committed atomically.
If all remaining units strictly improve the previous objective but cannot satisfy
the historical best gate, the complete remaining attempt is marked exhausted;
objective tolerances are never relaxed.

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
audit. If it fails, each currently accepted component is tried as a single
removal using a full-state audit. Choose the available removal with the lowest
objective, requiring improvement over the assembled objective when available;
ties retain the first traversal position. Repeat until the base-improvement
gate passes or no improving removal remains. Reuse the chosen removal's already
evaluated objective. Unless the surviving state strictly improves
the pre-polish global objective, the complete polish is discarded and the base
state is written unchanged. When the selected base state stopped by production
convergence, an objective-accepted polish remains provisional until the strict
fixed-point operator is evaluated again on the polished state. The polish is
applied only when every solver is qualified, the nominal-DOF operator is
complete, and every operator-residual p99 is below `1e-4`. An unavailable,
failed, or above-threshold certificate discards the polish and writes the
already converged base state unchanged.

For all non-convergence stops, the chosen base state is persisted directly,
without polish, but with its own nominal operator evaluation before persistence.
That certificate describes the chosen parameters under the last frozen background;
it never substitutes the final attempted iteration's evidence. Maximum residual remains diagnostic. Applied
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
- Active targets freeze after five consecutive observations of the same
  quarantine failure reason: solver hard failure, invalid candidate, or guard
  infeasibility. Targets remain shape blocks, offset blocks or hard-failure clusters.
  A changed reason or absent observation resets active tracking; a background
  refresh alone does not clear the consecutive-failure count.
  `objective-exhausted` remains a search terminal diagnostic but is not quarantine
  evidence. On its own it breaks an Active target's failure streak; other guard,
  invalid, or solver hard failures in the same attempt still count.
- Frozen targets record the last attempted Frozen-recovery revision, independent
  of the objective-domain revision used by the phase observer. Recovery advances
  on an applied partition/domain change or changed background response, not on an unchanged rebuild or merely queued topology. Each new
  revision permits one retry; retries are transient, not a third lifecycle.
  Background updates can consequently permit a retry every iteration.
  Objective reevaluation itself does not schedule recovery; the independent
  revision ownership and triggers are specified below.
- Retry keys use minimum radius `0.0625` and `10x` ridge. Non-retrying Frozen
  masks retain priority; overlapping retry targets do not unlock coordinates
  still frozen by another target. No cooldown, attempt limit or Exhausted state
  remains. Frozen targets awaiting a recovery revision change do not hold audit
  patience open as a scheduled recovery. Actual transitions retain their patience
  reset; active failure tracking alone does not.
- Recovery needs no affecting failure and all required target coordinates active.
  It requires a finally accepted material change, or a complete, guard-safe,
  solver-qualified nonmaterial unrestricted endpoint. Missing evidence cannot
  release a target. A failed retry records the revision and remains Frozen.
  Objective exhaustion does not block recovery, but does not itself release a
  Frozen target or bypass these domain-retry and endpoint requirements.
- Newly Frozen targets affect only the next proposal. Lifecycle updates never
  rewrite the audited model/provenance, so quarantine fallback re-audit is removed.
  Immediate block isolation and solver fallback remain. Final activity still
  fixes unresolved targets; transition convergence blockers remain.
- Rejected-cluster debug output distinguishes failures before objective
  evaluation from objective rejection. Pre-objective failures report their
  proposal reason, radius, available step norm, and `objective =
  not-evaluated`; `objective-unavailable` is reserved for an objective that was
  actually attempted but could not be calculated.

### Objective and Frozen-recovery revisions

`objective_domain_revision` identifies the objective context, including for
phase observation. `frozen_recovery_revision` independently schedules Frozen
retries; it is never copied or derived from the objective revision. Each target's
`last_recovery_revision` records freezing, retry initiation, and failed retry.
Only the recovery revision is passed to quarantine retry and publication.

| Event | Objective revision | Frozen-recovery revision |
| --- | --- | --- |
| Initialization | 1 | 1 |
| Queued partition applied at the next attempt | +1 | +1 |
| Same partition, changed `response_by_atom` | +1 | +1 |
| Partition applied together with background change | +1 once | +1 once |
| Same partition and exactly equal background response | Unchanged | Unchanged |
| Topology merely queued | Unchanged | Unchanged |

Objective reset owns the objective increment. Recovery advances explicitly after
a successful partition reset, or in the background-change branch. Objective
reevaluation alone is not a recovery event. Existing `domain-retry`/probation
diagnostic labels refer to recovery scheduling; their names do not imply shared
counter ownership.

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

Rescue remains enabled and shares the normal component pipeline and
result-application entry. Local outcomes stay provisional per key through
reconciliation and any triggered global salvage, after which accepted/rejected
lists are materialized once.

## Production member-best ownership

`IterationState::member_best` retains one committed parameter patch per cluster.
Local search and ordinary boundary members check both their previous objective and
this historical reference. Cooperative rescue keeps its existing member rules.
The historical comparison replaces only that member's parameters in the candidate
overlay, preserving every other candidate parameter and the same background/domain.
No saved scalar crosses backgrounds. Updates occur after successful transactions;
new cluster keys initialize from the committed state and obsolete keys are removed.

## Global audit and stopping

The global audit uses the fixed per-cluster fit/tail scales and retains the
earliest state that improves the best objective beyond the strict tolerance.

Pre-commit acceptance follows the [conditional selection audit](#conditional-selection-audit).
After commit, an all-rejected attempt has already published quarantine and
radius updates and retained the previous model; the runner returns before candidate scoring or best/patience updates and schedules controlled recovery for the next outer attempt. For an accepted
result, it advances accepted progress and checks adaptive topology, then reuses
`selection.final_audit_objective` when available. Otherwise it builds a snapshot
of the committed assembled state and evaluates its audit objective in the
current domain and frozen background.

An available score is used to update retained best and to assess strict
improvement over the same-background historical best for patience. If the
score remains unavailable, the improvement flag remains false and the existing
patience reset/increment rules still apply. This scoring does not reject or
roll back the committed candidate and does not substitute for a pre-commit
acceptance gate.

One internal `ConvergenceCertificate` is the sole source of convergence truth.
`ProductionConverged()` requires solver qualification, accepted active-DOF p99
below `1e-4`, a complete nominal-DOF operator, nominal fixed-point residual p99
below `1e-4`, and clear orthogonal blockers. Fixed and
quarantined coordinates are excluded only from the accepted population; they
remain in the nominal operator population. An empty accepted population passes
its percentile check vacuously, but an all-fixed state still needs qualified,
complete, sufficiently small nominal operator evidence. Every nominal shape and
offset endpoint carries the status of the solve that produced it, including the
separate unrestricted shape refit when quarantine changes the proposal. Offset
qualification also checks each active atom against its owning cluster's solve status. An
unavailable endpoint makes the operator incomplete instead of substituting
the previous state as a zero residual. `StrictOperatorPassed()` reuses the same certificate for converged
final-polish certification without the accepted-movement or orthogonal-blocker
terms.

Orthogonal blockers cover objective-domain changes, quarantine transitions,
suspicious block fallback, and rejected clusters. `suspicious_block_fallback`
includes shape and hard-failure evidence as well as offsets. Maximum values are
not certificate predicates; their uses in nonmaterial search/backtracking and
topology drift remain independent of convergence.

The stage stops on the first applicable condition:

- no valid initial seed is available for every selected atom;
- accepted active-DOF p99 and complete nominal-DOF fixed-point residual p99 are
  both below `1e-4`, every nominal endpoint and active coordinate is solver-qualified, all clusters
  are accepted, and no orthogonal blocker is present;
- controlled fixed-point recovery fails after patience or an all-rejected attempt;
- final persisted-state certification fails after provisional convergence;
- `kLocalFittingMaximumIterations` outer attempts are reached.

Three accepted iterations without strict historical-best improvement schedule recovery,
as does an all-rejected ordinary attempt. Actual domain/quarantine transitions and
finite rejected-radius shrink actions retain their reset behavior; merely having
an active quarantine tracker does not indefinitely reset patience. Recovery is a
persistent production mode, bounded by the same 100 outer attempts. Its eight-step
search, qualified nominal residual requirement and historical objective envelope
are specified in [Production fitting](production-fitting.md).

Convergence writes the current accepted state. Recovery failure, audit-patience, all-rejected,
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
on the existing fast path. More precisely, the [conditional selection audit](#conditional-selection-audit)
is skipped when the ordinary component list is empty and rescue accepts no
result. A committed local selection can still require the runner's post-commit
objective evaluation for best/patience tracking.

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

Objective-domain startup reports weights and scalar sample/cluster counts.
Warnings report cumulative quarantine entries, releases, failed retries and
unresolved targets. Progress, necessary warnings, final summary and basic timings
remain available without extra audit payloads.

Enabled, non-quiet Debug sessions emit `Second-stage audit: schema=2, payload={...}`.
There is one start record, one summary per attempt, and one terminal record.
Only actual production evidence is recorded; no solver/operator replay, historical
model scoring, complete snapshots, atom dump, or alternate coupling threshold scan
remains. Each attempt and finalization retain at most five abnormal details at
collection time, with complete category totals and omission counts.
See [Second-stage audit](second-stage-audit.md) for fields, references, failure
isolation, parser usage and verification tiers.

The one-time atom-cutoff text remains `Local-fitting atom cutoff: atoms=N,
limit=100, clusters=C, max-atoms=M, cutoff-edges=E.` and counts selected atoms only.
Formal graph thresholds, hysteresis, binary fallback and cutoff are unchanged.

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
