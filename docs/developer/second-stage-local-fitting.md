# Second-stage local fitting

## Scope

`RunSecondStageLocalFitting(ModelObject&, const FitOptions&)` refines the
selected atoms after first-stage local fitting. Its public signature and
workflow position are stable. The implementation keeps a separate in-memory
state until a validated terminal state is chosen, so an iteration never
partially mutates `ModelObject`.

The second stage is intentionally a single fixed-point path. It does not use
Anderson acceleration, damping ladders, freeze/thaw tracking, stagnation
regimes, or adaptive cluster/global solver ridge.

## Safety invariants

The simplification does not remove the numerical defenses that protect model
validity or cluster isolation:

- Invalid initial Gaussian seeds are repaired from the established fallback
  sources before iteration. If every selected atom cannot receive a valid
  seed, the stage is skipped without changing stored estimates.
- The weighted coupling graph partitions active atoms into deterministic
  cluster keys. Boundary samples are represented in every cluster that owns
  them, and a combined-objective guard protects cross-cluster interactions.
- Joint offsets are solved with robust IRLS and a fixed
  `kJointOffsetRidgeRatio`. Conditioning and column-collinearity guards may
  raise an atom-local ridge multiplier.
- A suspicious joint-offset or post-refit model is rolled back to the previous
  validated atom state. On the next iteration only the affected atoms receive
  `kSuspiciousJointOffsetRidgeMultiplier`.
- Persistent hard solver failures and persistent suspicious rollbacks become
  terminal cluster fallbacks. Terminal atoms no longer participate in later
  solves; healthy remote clusters remain active.
- Amplitude and width must remain finite and positive. Offset and transformed
  coordinates must remain finite. Candidate construction failure rejects only
  the owning cluster.
- The scientific objective, terminal isolation, and global best-audit fallback
  remain independent of trust-radius bookkeeping.

## Iteration flow

Each outer iteration performs the following sequence:

1. Build the active atom list from every non-terminal selected atom.
2. Rebuild the weighted coupling partition and reconcile the per-cluster
   objective and trust-radius maps.
3. Solve all active joint offsets once using the fixed ridge ratio. Apply
   suspicious rollback and run the local Gaussian refits to produce one raw
   fixed-point state.
4. For every cluster, request the full raw step (`1.0`). The trust region may
   truncate that step once; there is no secondary damping search.
5. Score the base candidate on the cluster's fixed objective scale. A reject
   leaves both the assembled state and objective tracker unchanged.
6. If the base candidate passed and the cluster is stationarity-eligible,
   attempt one joint amplitude/width/offset polish. The polish requests a full
   step, is limited by the remaining trust radius, and replaces the base only
   when it strictly improves the same fixed-scale objective.
7. If boundary samples exist, apply the combined-objective guard to the
   assembled cluster candidates. Cluster objective state is committed only
   after this guard passes.
8. Grow a cluster radius only for objective improvement near its boundary.
   Shrink only rejected cluster radii.
9. Update persistent terminal-failure state, the fixed global audit, and the
   stopping conditions.

This is the complete candidate selection path. There is no candidate-kind
dispatch, forced fixed-point mode, history suppression, discrete damping
ladder, or adaptive objective ridge.

## Cluster objective state

The production state is a direct map owned by the implementation file. Each
cluster stores only:

- a scale-reference tracker;
- the previous transformed-change statistics and objective samples;
- the best tracked candidate's transformed-change statistics and objective
  samples, when a valid objective has been observed.

Candidate scoring works on a copy of this map. A rejected candidate, including
a combined-objective rejection, cannot advance the scale warmup, previous
candidate, or best candidate. Candidate, previous, and best objectives are
recomputed on the same provisional objective scale for each attempt rather than
being stored in multiple state layers.

The objective contains robust residual loss plus width-prior and offset
plausibility penalties. Joint polish uses the same committed objective scale as
its base candidate and therefore cannot win by changing normalization.

## Trust region

Every base and polish proposal requests damping `1.0`. The effective step is
computed directly from the transformed parameter norm and the cluster radius:

```text
effective damping = min(1.0, trust radius / transformed step norm)
```

Polish is limited by the radius remaining after the base movement. Objective
rejection shrinks only the rejected cluster. Improvement grows a radius only
when the accepted step is close to the boundary. At the minimum radius there
is no solver-ridge escalation or global retry regime.

## Audit patience and stopping

The global audit uses one fixed scale and retains the earliest state that
improves beyond the tie tolerance. After terminal isolation changes the audit
domain, terminal samples and penalty atoms are removed and the current
validated non-terminal state becomes the new fallback baseline.

The stage stops on the first applicable condition:

- transformed accepted and raw fixed-point changes both converge, with all
  active objective references locked and no suspicious or unhealthy cluster;
- three accepted iterations produce no strict global audit improvement;
- all clusters reject at their minimum trust radius;
- every cluster has become terminal;
- the maximum outer-iteration count is reached.

Trust-radius shrink attempts do not consume audit patience while they are still
changing a rejected cluster's radius. This prevents a partially accepted
cluster from terminating recovery of an independent rejected cluster.

Convergence applies the current accepted state. Audit patience and iteration
limits apply the best validated audit state. A minimum-radius all-reject stop
also uses the best audit fallback. Terminal reconciliation preserves validated
non-terminal cluster progress rather than restoring an obsolete whole-model
snapshot.

## Logging

After a valid seed is available, non-quiet runs print one iteration-table
header and update one progress row per outer attempt with
`Logger::ProgressLine`. An outer attempt includes both an accepted iteration
and an all-rejected trust-radius retry, so `Try` can advance without `Acc`.

| Column | Meaning after the current outer attempt |
|---|---|
| `Try/Acc` | One-based outer attempt / cumulative accepted iterations |
| `Atom A/T` | Remaining active / cumulative terminal atoms |
| `Cmp/Max/R` | Partition components / largest component atoms / largest-component-to-active-atom ratio |
| `Cand C;A A/R` | Accepted/rejected candidate clusters; accepted/rejected candidate atoms after the objective gate |
| `Pol` | Clusters whose strict-improvement polish replaced the base candidate |
| `TR G/S/M` | Trust-radius grows / successful shrinks / rejected clusters already at minimum radius |
| `Guard S/U/C` | Suspicious atoms / stationarity-ineligible clusters / combined-objective rejection (`0` or `1`) |
| `dMax A/R` | Maximum transformed change in the accepted/raw state; accepted is `-` on an all-rejected attempt |
| `Audit B/P` | Best audit iteration (`I`, accepted-iteration number, or `-`) / current audit patience |

Rows have a precomputed fixed width and overwrite the previous row. Debug
rejection diagnostics finish the active progress line before printing their
details, then the current attempt row is displayed. Terminal, convergence, and
summary messages likewise finish the progress line before normal line output.

Non-quiet runs end with one stable summary line:

```text
Second-stage local fitting summary: accepted_iterations=<N>, best_iteration=<initial|N|unavailable>, stop_reason=<reason>, best_audit_objective=<value|unavailable>.
```

The fold-168 benchmark parses this line and requires at most ten accepted
iterations. `quiet_mode` continues to suppress second-stage informational
logging.

## Verification

The focused estimator suite covers seed repair, non-finite fallback,
conditioning/collinearity ridge, suspicious rollback, remote-cluster
isolation, joint polish, trust-region limits, and intensity-scale invariance.
Shared transformed-change, robust-loss, scale-tracking, and
weighted-ridge utilities retain focused algorithm tests because the production
path still depends on them.

The external fold-168 regression checks all 168 serial IDs, finite and valid
parameters, truth-based amplitude/width/offset RMSE, maximum absolute offset,
and the accepted-iteration limit. Wall time is reported but is not a blocking
gate.
