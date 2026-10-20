# Joint Peeling Chemical-label Audit

The baseline section records the production dependency before the Halo policy
change. The reference tree was `develop` at
`76d52c94096f38fa8eeaba7e8245ddc92c62087f`.

## Call graph and label use

```text
RunJointComponentWorkflow
  ├─ BuildPureHaloSupportExclusionProblem
  │    └─ selects target and halo roles; keeps the existing non-hydrogen policy
  ├─ FitJointComponentsInitializedWithProgress
  │    └─ qualifies target FullABC endpoints
  ├─ BuildHaloMedianReconstruction
  │    └─ groups qualified target donors and original-domain Halos by GroupKey
  └─ JointPostprocessingContext
       ├─ validates saved Halo GroupKey against ModelObject
       ├─ BuildPostFitPeelingSamples
       │    ├─ uses original-domain support and saved Halo A/B/C
       │    └─ reconstructs off-domain nodes from the full model atom list
       ├─ BuildPostJointLocalSamples
       │    └─ pairs peeled samples and subtracts the target's own Joint charge
       └─ RunJointGroupPotentialFitting
            └─ groups local results by GroupKey
```

| Location | Read | Effect | Classification |
| --- | --- | --- | --- |
| `HaloMedianReconstruction.cpp` | `GetComponentKey()` and `GetAtomKey()` via `GetGroupKey()` | Chooses target donors for each original-domain Halo's A/B/C; changes the predicted field | Direct numerical dependency to remove |
| `PostFitPeeling.cpp` snapshot validation | `GetGroupKey()` | Rejects a snapshot if the Halo's current chemical grouping differs | Chemical-label validation; remove with the old provenance contract |
| `PostFitPeeling.cpp` contributor universe | `GetElement()` | Omits Hydrogen from reconstruction-only contributors | Preserve the upstream non-hydrogen policy, but carry its result into Peeling as immutable eligibility |
| `PostFitPeeling.cpp` extra contributors | `GetGroupKey()` | Assigns the same-group median or global median to atoms outside the Joint input domain | Direct numerical dependency to remove |
| `PostFitPeeling.cpp` identity lookup | Atom serial ID | Associates a contributor with its Joint index and excludes the target's own signal | Allowed identity use; do not use the ID to recover a chemical class |
| `BuildPostJointLocalSamples` | Target position, Joint B/C, raw geometry | Builds the target's own charge profile and Gaussian-only samples | Numeric endpoint/geometry use; no chemical classification |
| `RunJointGroupPotentialFitting` | `GetGroupKey()` | Defines the requested downstream chemical groups | Preserve |
| UMAP and result export | Residue, element, spot, atom ID | Labels rows and plots | Preserve as display metadata; do not feed labels into Tail Ratio |

`AtomClassifier::GetGroupKey()` packs `ComponentKey` and `AtomKey`. Residue,
Element and Spot are not used by the median helper itself, but can affect which
values those keys contain. The element filter is a separate contributor
eligibility decision and must stay upstream of the label-free Peeling API.

## Invariance contract

For a Peeling-level comparison, freeze map values, positions, atom identities,
target/halo roles, the eligible contributor universe, qualified target A/B/C,
sampling geometry, support radius and tricubic interpolation. Change only
Residue, ComponentKey, AtomKey, Spot, or Element while retaining the already
fixed eligibility universe. Halo parameters and reconstructed predictions,
Peeling samples, Tail Ratio, own-charge correction and Post-Joint Local Gaussian
samples must then be equal. Atom IDs may identify atoms and break exact-distance
ties, but may not be looked up to recover a chemical label.

Changing a non-hydrogen atom to Hydrogen before constructing the upstream
contributor universe changes eligibility by policy; that is not a Peeling
invariance comparison. Likewise, rerunning name-dependent upstream target
selection is outside this fixed-universe contract.

## Current baseline experiments

The deterministic helper fixture has three qualified donors with points
`(A, B, C) = (1,1,2), (3,2,5), (5,3,8)`. When all donors share a GroupKey, the
current global and same-group median is `(3,2,5)`. Relabeling only the third
donor so it no longer shares the Halo's GroupKey leaves the global median at
`(3,2,5)`, but changes that Halo's same-group median to `(2,1.5,3.5)`. The
regression is `CurrentGroupKeyPolicyChangesHaloEstimateWhenOnlyTargetLabelsChange`.

The reconstruction-only fixture freezes both targets, all four atom positions,
the map, and roles while changing only the two off-domain contributors'
ComponentKey and AtomKey. The regression verifies that current extra-domain
Peeling responses differ, demonstrating the separate GroupKey dependency in
the node-reconstruction path.

Existing postprocessing metrics already record qualified donor count,
same-group/fallback Halo counts, fallback support rows, reconstruction-only
node evaluations, cache peak, contributor lookup count and estimated
workspace bytes. Tail coverage and unavailable reasons are tracked per run.
Halo contribution magnitude and per-target coverage are not currently
persisted as a separate production audit. The Phase 2 research runner records
synthetic per-Halo signal errors, donor coverage, availability, runtime and
workspace estimates in its CSV output.

## Current production boundary

`BuildHaloMedianReconstruction` now builds its donor pool only from qualified
Joint Target point estimates, then uses a KD-tree over target positions to
select the four nearest donors for each Original-Domain Halo. Componentwise
medians produce A/B/C; equal-distance ordering uses Atom ID. The exact search
policy is recorded in the Joint JSON schema-18 reconstruction snapshot.

`BuildJointPeelingContributorUniverse` is run before postprocessing and carries
every model atom, its position and identity, any Original-Domain index and
Target/Halo role, and the already decided non-hydrogen eligibility. The
universe includes eligible contributors outside the Joint domain. Peeling
uses the same qualified Target donor index for those contributors and the
saved Halo profiles, so the two paths share one policy. Missing qualified
donors leave reconstruction unavailable; no parameter or sample is filled
with zero.

The production Joint Halo and Peeling files no longer include or call
`AtomClassifier`; `GetElement()` is read only while the upstream contributor
universe freezes the existing hydrogen-exclusion decision. Atom ID is used
for identity, attribution, self-contribution exclusion, and deterministic
ties. `RunJointGroupPotentialFitting` still uses GroupKey after local fitting.
UMAP continues to export and display residue/element/spot labels, while its
Tail Peeling Ratio uses the already paired, label-free peeling values.

The Phase 4 invariance tests change ComponentKey, AtomKey, ComponentID,
Element (between non-hydrogen elements), and other contributor labels while
holding map, coordinates, identity, roles, eligible universe, endpoint and
sampling fixed. They compare Halo A/B/C, peeling response, Tail Ratio,
charge-corrected Gaussian response, and Post-Joint local samples. Legacy
GroupKey snapshots decode and round-trip as schema 17; the new Peeling path
rejects them instead of interpreting them under the geometry contract.
