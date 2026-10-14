# Joint stage integration

This guide owns the current integration between Joint Component fitting and
estimator-neutral stage results. Numerical fitting rules are owned by the
[runtime guide](joint-component-runtime.md); SQLite and neutral-document
semantics are owned by the
[data I/O architecture](architecture/dataobject-io-architecture.md).

## Stage estimates and provenance

A local stage record separates its optional point, source, convergence, reason,
and estimator-neutral uncertainty. Its source identifies estimator method,
atom, component, run, and target role. An unavailable estimate is not encoded as
a zero point, and Joint points are not stored in OLS or local MDPDE diagnostics.

The Joint adapter maps returned component states through contributor identities.
It preserves an available point even when runtime convergence fails, with the
failure status alongside it. A missing component state yields an unavailable
estimate with its stop reason. For recorded selection domains, target and halo
roles come from the saved domain; inputs without selection metadata retain
not-recorded. Observable-contribution-only halos do not receive a point.

Joint points occupy the neutral Second-stage estimate. The stage summary,
group fitting, displays, and feature/export consumers include only Joint target
points. Available unconverged target points remain visible and are marked
unconverged. Charge C is descriptive; group inference uses the eligible target
parameter evidence and covariance, and does not infer charge uncertainty.
Joint point estimates do not carry a full-component covariance. Group
ParameterEvidence stores the separate local conditional covariance used by
Group Fitting.

## Workflow ownership

The map-aware Joint workflow builds one immutable problem and initializes all
eligible contributors using their model identities. First-stage initialization
supplies widths; only those widths seed Joint fitting. Initialization sampling
may read outside the target observation rows, while the Joint objective remains
on the saved fixed domain. The second-stage refinement option is for the
two-stage estimator; the model-only workflow does not accept Joint requests.

After fitting, the workflow publishes the Joint stage estimates, builds
target-only post-fit peeling, computes local conditional covariance at each
eligible fixed Joint endpoint, and builds GroupParameterEvidence in temporary
results. Peeling and evidence are validated and published together as one
batch before the Joint snapshot is stored and group fitting runs. The batch
checks target and run identity, raw/peeled sample counts, finite responses and
evidence identity before changing the model; it collects affected groups once
and invalidates each once. The ordinary single-atom setters keep their own
invalidation behavior.

Peeling responses paired with raw sample geometry supply each atom's local
observations; neighboring Joint parameters and the atom's Joint charge are
fixed. Peeling reads the raw sample list through a const view while the model
is unchanged, so it does not copy the full list per target. The local estimator
only computes the A/B curvature and never updates the Joint point. No
full-component Joint covariance or cross-atom covariance is computed.

For each eligible target, Group Fitting receives the encoded Joint point and a
2-by-2 covariance in its group parameter coordinates. That covariance is a
conditional working covariance under an i.i.d. Gaussian residual model; it does
not include uncertainty in neighboring parameters. Group inference uses a
by-atom block-diagonal working covariance.

## Post-fit peeling

Peeling subtracts the fitted model on the Joint observation rows, including
contributions from halo atoms. It preserves each raw sample and stores an
optional peeled response with a reason when coverage or a contributor state is
missing. Ratios require complete paired coverage in their interval and retain
signed values without clamping. Peeling does not extend the observation domain.

Replacing a Joint endpoint invalidates derived peeling, local group evidence,
posterior results, and the saved Joint snapshot for the affected run. A
second-stage endpoint change clears its estimator-neutral uncertainty and
group evidence; group summaries affected by the changed endpoint are
invalidated. Replacing raw samples or peeling responses clears local group
evidence and affected posterior results while preserving Joint point estimates.
The editor enforces these source and invalidation relationships.

## Persistence boundary

The in-memory model owns the neutral stage record and the Joint result snapshot.
Storage validates and round-trips these values but does not recompute numerical
evidence. Current database versions, migration, transaction, and JSON payload
contracts are defined in the
[data I/O architecture](architecture/dataobject-io-architecture.md) and the
[potential-analysis command guide](commands/potential_analysis.md).

## Permanent coverage

Permanent tests cover Joint-to-stage mapping and point preservation, target and
halo roles, local covariance against a finite-difference reference, unavailable
local designs, target-only summaries and group inference, peeling coverage,
persistence, and invalidation after endpoint or sample changes. The current test lane is
documented in [tests/README.md](../../tests/README.md). Historical integration
results and their boundaries are indexed in the
[canonical evidence document](joint-component-evidence.md).
