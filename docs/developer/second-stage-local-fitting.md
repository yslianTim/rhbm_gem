# Second-stage local fitting

## Scope

`RunSecondStageLocalFitting(ModelObject&, const FitOptions&)` refines the
Gaussian estimates of the selected atoms after the first-stage local and group
fits. It accounts for overlapping responses from neighboring selected atoms
while updating each atom's amplitude, width, and offset.

The stage keeps candidate states in memory and writes one validated terminal
state to `ModelObject`. Individual outer iterations do not partially update the
stored atom estimates.

## Model context and initialization

The fitting context contains, for each selected atom:

- its local-potential samples and trained `alpha_r`;
- a width prior, preferring the group prior over the same-group median and the
  atom's current valid width;
- the selected neighboring atoms that contribute to each sample.

Neighbor candidates are searched within `kNeighborAtomSearchRange`. A neighbor
contributes to a sample only when its distance from that sample does not exceed
`kNeighborContributionDistanceMax`.

The current local MDPDE estimate is used as the initial Gaussian seed when it is
valid. An invalid seed is repaired from the first valid source in this order:

1. group posterior;
2. group prior;
3. local OLS estimate;
4. same-group parameter median;
5. global parameter median.

If a valid seed cannot be obtained for every selected atom, the stage exits
without changing the stored estimates.

The valid initial state is used to build a weighted coupling topology. The
topology records atoms that jointly affect objective samples. Each iteration
partitions the active portion of this topology into deterministic cluster keys,
so independent clusters can be evaluated and accepted separately.

## Iteration flow

Each outer attempt performs the following sequence:

1. Build the active atom list by excluding terminal-fallback atoms.
2. Partition the active coupling topology and reconcile the per-cluster
   objective and trust-region states.
3. Jointly estimate one shared offset per represented group within each cluster
   using robust IRLS and the fixed `kJointOffsetRidgeRatio`.
4. Build component-wise group-median models from the post-solve snapshot. For
   each active atom, subtract the group-median responses of its selected
   neighbors from the observed sample responses.
5. Refit the atom's local Gaussian with its trained `alpha_r`, using its
   group-median model as the fixed offset model. These refits form the raw
   fixed-point state.
6. Limit each cluster's raw proposal to its trust region and score the resulting
   base candidate.
7. For a stationarity-eligible cluster without suspicious atoms, attempt one
   joint amplitude/width/offset polish. Keep the polish only when it strictly
   improves the base candidate on the same objective scale.
8. If clusters share boundary samples, validate the assembled candidates with
   the combined-objective guard.
9. Update trust radii, persistent-failure state, the global audit state, and the
   stopping conditions.

The neighbor-adjusted response for atom `i` and one of its samples is:

```text
adjusted response = observed response
                  - sum(neighbor fitted responses at the sample position)
```

This offset solve and neighbor-adjusted local refit constitute one raw
fixed-point update.

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

Each cluster tracks:

- a scale-reference tracker;
- the previous accepted candidate's objective samples;
- the best objective samples and the corresponding maximum transformed change
  used to break objective ties.

The objective combines:

- Cauchy robust residual loss;
- a width-prior penalty weighted by
  `kLocalFittingWidthPriorPenaltyWeight`;
- an offset-plausibility penalty weighted by
  `kLocalFittingOffsetPlausibilityPenaltyWeight`.

The width penalty measures the log-width displacement from the atom's prior.
The offset penalty applies when the offset response is too large relative to
the fitted peak and residual scale.

Candidate scoring uses a provisional copy of the cluster objective state. A
rejected base candidate or a combined-objective rejection does not advance the
scale warmup, previous candidate, or best candidate. Candidate, previous, and
best objectives are evaluated on the same provisional scale. Joint polish uses
the committed scale of its accepted base candidate.

## Trust region

Each base proposal requests the full raw step. Its effective damping is limited
by the transformed step norm and the cluster radius:

```text
effective damping = min(1.0, trust radius / transformed step norm)
```

The polish step is limited by the radius remaining after the base movement. A
rejected cluster shrinks its own radius. An accepted cluster grows its radius
only when the objective improves and the accepted step is close to the current
boundary. Trust-region updates are isolated by cluster.

## Numerical defenses and terminal isolation

- Joint-offset conditioning and column-collinearity guards can increase the
  ridge multiplier for affected atoms.
- A suspicious joint-offset or post-refit model is rolled back to the previous
  validated atom state. On the next attempt, affected active atoms receive
  `kSuspiciousJointOffsetRidgeMultiplier`.
- A failed or invalid local refit falls back to the previous Gaussian
  parameters with the newly estimated offset when that model remains valid.
  The owning cluster is then stationarity-ineligible for that attempt. If the
  fallback is also invalid, the affected post-refit cluster is rolled back.
- Repeated stable hard solver failures or suspicious rollbacks become terminal
  cluster fallbacks. Terminal atoms retain their previous validated states and
  no longer participate in later solves.
- Terminal isolation removes only the affected cluster, allowing independent
  active clusters to continue fitting.

When terminal isolation changes the active domain, the global audit removes the
terminal samples and penalty atoms and reconciles its fallback state with the
validated non-terminal progress.

## Global audit and stopping

The global audit uses a fixed objective scale and retains the earliest state
that improves the best objective beyond the tie tolerance.

The stage stops on the first applicable condition:

- no valid initial seed is available for every selected atom;
- every selected atom has become terminal;
- the accepted and raw fixed-point transformed changes both converge, all
  active objective references are locked, and no active cluster is rejected,
  suspicious, or unhealthy;
- `kLocalFittingAuditPatience` accepted iterations produce no strict global
  audit improvement;
- all clusters reject at their minimum trust radius;
- `kLocalFittingMaximumIterations` outer attempts are reached.

Trust-radius shrink retries do not consume audit patience while a rejected
cluster's radius is still changing.

Convergence writes the current accepted state. Audit-patience, minimum-radius
all-reject, and iteration-limit stops write the best validated audit state when
one is available, otherwise the current validated fallback. Terminal
reconciliation preserves validated progress from non-terminal clusters.

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
| `Polish E/A/R/S` | Eligible / accepted / rejected / skipped polish clusters after the combined-objective guard; `E = A + R + S` |
| `Suspicious` | Atoms rolled back by the suspicious-offset checks in this attempt |
| `dMax A/R` | Maximum transformed change in the accepted/raw state; accepted is `-` on an all-rejected attempt |

Debug rejection diagnostics finish the active progress line before printing
their details. Terminal, convergence, and summary messages do the same before
normal line output.

Non-quiet runs end with this summary format:

```text
Second-stage local fitting summary: accepted_iterations=<N>, best_iteration=<initial|N|unavailable>, stop_reason=<reason>, best_audit_objective=<value|unavailable>, final_uses_polish=<yes|no|unavailable>.
```

`final_uses_polish` describes the state actually written to `ModelObject`. It is
`yes` when at least one atom's most recent transformed-parameter update in that
state came from an accepted polish, `no` when none did, and `unavailable` when
the stage exits before a valid state can be formed.

`quiet_mode` suppresses the second-stage informational logging.
