# Component-local certification and runtime parity

This work starts at `6f30510c`. The exact-component v1 and certification v2
contracts remain the numerical baseline. Local certificates supplement, and
never replace, monolithic/global and globally restricted certificates.

## Component-local evidence

`joint-abc-components-local-audit DATASET RUN OUTPUT [CASE]` reads frozen
component fits and their immutable search context. It does not read an assembled
fit or a final global context. `joint-abc-rerun-local-component DATASET CASE ID
CONTEXT OUTPUT` requires only `search_context` in its context bundle.

The Python component runner exposes these through `audit --local-only` and
`rerun-component --local-only`. Ordinary full audits additionally produce
`local-components/`, preserving the existing `components/` restricted evidence.

Each local scope records the stable component ID, parent snapshot identity and
exact eta/beta of the last trusted state used by assembly. Assessment uses those
coefficients, with a separately reprofiled consistency control. Normalization
and active-set response norm still come from all parent observations. Rank
sizes remain the existing independent-component sizes.

Inherited directions are discarded. The local endpoint supplies normalized
all-ones, alternating and weakest projected-width directions, with direction
norms and Jacobian response norms saved. Missing trusted states or width spectra
are explicit limitations, not unit-axis substitutes. Search completion, usable
state, derivative evidence and regular qualification remain separate.

The two-step derivative checks and registered Richardson/50-100 digit audits
retain their existing thresholds. Truth never selects directions, intervals or
qualification. Local certificates additionally require a usable actual state,
weak-direction evidence and agreement with the profiled control. Missing
components cannot be promoted by another component's certificate.

Independence means identical component input and immutable parent context:
sibling execution order and sibling endpoint availability cannot affect local
evidence. It does not assert invariance when parent observations/scale change.

## Implementation and validation records

The implementation is delivered in three stages: local certification, typed
runtime extraction, and a fresh Map/Model entry point. Validation results are
recorded under `figures/joint-component-runtime/`; historical archives are
read-only reference inputs. A successful local scope does not change any
historical global qualification.
