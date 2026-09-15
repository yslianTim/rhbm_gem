# Production fitting and controlled recovery

## Inputs and statistical contract

Production fitting consumes sampled map responses, selected atomic geometry,
explicit fitting options and previously estimated state. It never loads the
simulation manifest, generated charges, truth widths or scoring output. Atomic
number is not an amplitude truth initializer. `potential_analysis --simulation`
continues to skip normalization; `-r` supplies model metadata only. SQLite and
the schema 7 offline truth definition are unchanged.

This first repair preserves the MDPDE equations, joint-offset IRLS equations,
inner iteration budgets, transformed p99 threshold `1e-4` and 100 outer attempts.
Parameter accuracy does not establish convergence. Inner numerical failures are
reported and reproduced rather than reclassified as success.

## Evidence and ownership

Each nominal shape endpoint stores its own MDPDE status, iterations, terminal
squared beta change, relative variance change and variance. Each offset endpoint
stores its joint solve status, iterations, normalized change and robust scale.
Quarantine may require a separate unrestricted shape solve: its status cannot be
borrowed from the constrained proposal. All nominal endpoints must be qualified,
including coordinates excluded from accepted-movement statistics.

`EvaluateNominalOperator` uses the same proposal implementation with unrestricted
activity and scratch solver workspaces. It changes neither model parameters,
quarantine, trust radii nor historical bests. Production invokes it for recovery
and final certification regardless of audit enablement. Every persisted state,
including a restored best, receives its own operator assessment under the last
frozen background also used by persistent peeling. Failed provisional final
certification reports `final-certificate-failed`; non-convergence stop reasons
are not upgraded just because a final residual happens to be small.

`member_best` is production-owned, one parameter patch per cluster. Ordinary
local/boundary acceptance reevaluates that patch in the candidate overlay,
replacing only the target member. Historical scalars are never reused across
backgrounds. Successful commits update history; new partition keys initialize
from committed state and removed keys are discarded. Cooperative rescue retains
its previous member policy.

## Recovery

Three accepted iterations without strict historical-best improvement, or an
all-rejected ordinary attempt, enter recovery on the next outer attempt. Actual
domain/quarantine transitions and finite rejected-radius shrink actions retain
their patience resets. The existence of active quarantine tracking is not progress.
Recovery remains active until convergence or failure; ordinary updates cannot
immediately undo its progress.

Recovery freezes the background, domain, partition and ridge settings during
each search. Its direction is the unrestricted nominal endpoint, with the existing
Gaussian interpolation and guards. It tries at most `1, 1/2, ..., 1/128`, checking
each member's current trust radius. Both current and trial operators must be
complete and qualified. Inactive coordinates or an unqualified operator block
recovery explicitly. An already fixed current state may be certified without a
search step, provided it satisfies the objective bound.

For all selected atoms and three transformed coordinates, let `M` be the mean
squared nominal residual divided by `1e-4` squared. A trial at factor `lambda`
must satisfy both:

```text
M_trial <= (1 - 1e-3 * lambda) * M_current
J_trial <= J_best + 1e-8 + 1e-3 * abs(J_best)
```

The best objective is reevaluated in this same environment. The tolerance is
anchored to historical best and cannot accumulate after each permitted increase.
This recovery policy is separate from ordinary member-best gates. An exhausted
or blocked recovery stops with `recovery-failed` and preserves the best validated
state. All recovery attempts count toward the outer limit; basic final-state
output reports extra recovery and certificate operator evaluations separately.

## Verification and failure replay

Use matched compiler, feature flags, worker count, verbosity and quiet setting
when comparing audit OFF/ON. Debug verbosity has an existing scheduling effect;
a `-v 4` diagnosis must not be presented as the identical `-v 3` benchmark run.
The numerical probe compares work, commits, terminal state, parameters and peeling.
The external runner compares 168 atoms and retains the 25-iteration budget;
`convergence_acceptance` additionally requires a qualified persisted-state
certificate and a `converged` stop within 25 outer attempts. Quality remains
`uncalibrated`, so its overall exit status remains 1.

Testing builds can save the first failure of each solver/status pair:

```sh
RHBM_TEST_SOLVER_CAPTURE_DIR=build/production-fitting/failures \
  build/production-fitting/on/bin/RHBM-GEM potential_analysis [normal arguments]
RHBM_TEST_REPLAY_DIR=build/production-fitting/failures \
  build/production-fitting/on/bin/RHBM-GEM-TEST \
  --gtest_filter=ProductionFittingTest.ReplaysCapturedSolverFailureWithoutSimulationInputs
```

These environment variables exist only in `BUILD_TESTING` binaries. Capture is
independent of the observer, never changes a result, and the validation harness
must verify that requested artifacts were produced. No map or manifest is needed
to replay a captured subproblem. Retain executable/library hashes and build/source
provenance beside the fixtures; replay uses the same compiled solver constants.

Fixture version 1 is whitespace-delimited with round-trip double precision:

- `shape 1`: alpha, threads, iteration limit, tolerance, data-weight floor;
  dense X and y; expected status/variance; expected OLS and MDPDE beta vectors.
- `offset 1`: sparse X dimensions/nonzero count and row/column/value entries
  (including explicit zeros); y, anchor and ridge vectors; expected status and
  offset vector. The compiled IRLS constants remain unchanged.
- Dense matrices/vectors contain row and column counts followed by row-major
  values. Replay requires exact parameters and matching terminal status.

The first-repair report distinguishes implementation verification, convergence
acceptance and quality calibration. An early failure with a reproducible solver
case controls the cycle but does not meet convergence acceptance.
