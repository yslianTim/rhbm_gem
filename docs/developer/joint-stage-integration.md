# Joint stage integration

This guide owns the current integration between Joint Component fitting and
estimator-neutral stage results. Numerical fitting rules are owned by the
[runtime guide](joint-component-runtime.md); SQLite and neutral-document
semantics are owned by the
[data I/O architecture](architecture/dataobject-io-architecture.md).

## Stage estimates and provenance

A local stage record separates its optional point, source, convergence, reason,
and uncertainty. Its source identifies estimator method, atom, component, run,
and target role. An unavailable estimate is not encoded as a zero point, and
Joint points are not stored in OLS or local MDPDE diagnostics.

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

## Workflow ownership

The map-aware Joint workflow builds one immutable problem and initializes all
eligible contributors using their model identities. First-stage initialization
supplies widths; only those widths seed Joint fitting. Initialization sampling
may read outside the target observation rows, while the Joint objective remains
on the saved fixed domain. The second-stage refinement option is for the
two-stage estimator; the model-only workflow does not accept Joint requests.

After fitting, the workflow publishes the Joint stage estimates, builds
target-only post-fit peeling from the fixed observation rows, computes eligible
target uncertainty and group evidence, stores the Joint diagnostic snapshot,
then runs group fitting. These steps consume the returned Joint state; no
separate local fit substitutes for a missing component.

## Post-fit peeling

Peeling subtracts the fitted model on the Joint observation rows, including
contributions from halo atoms. It preserves each raw sample and stores an
optional peeled response with a reason when coverage or a contributor state is
missing. Ratios require complete paired coverage in their interval and retain
signed values without clamping. Peeling does not extend the observation domain.

Replacing a Joint endpoint invalidates derived peeling, parameter evidence,
posterior results, and the saved Joint snapshot for the affected run. A
second-stage endpoint change clears its uncertainty and group evidence; group
summaries affected by the changed endpoint are invalidated. Replacing raw
samples clears peeling while preserving point estimates. The editor enforces
these source and invalidation relationships.

## Persistence boundary

The in-memory model owns the neutral stage record and the Joint result snapshot.
Storage validates and round-trips these values but does not recompute numerical
evidence. Current database versions, migration, transaction, and JSON payload
contracts are defined in the
[data I/O architecture](architecture/dataobject-io-architecture.md) and the
[potential-analysis command guide](commands/potential_analysis.md).

## Permanent coverage

Permanent tests cover Joint-to-stage mapping, target and halo roles, unavailable
states, target-only summaries and group inference, peeling coverage, and
invalidation after endpoint or sample changes. The current test lane is
documented in [tests/README.md](../../tests/README.md). Historical integration
results and their boundaries are indexed in the
[canonical evidence document](joint-component-evidence.md).
