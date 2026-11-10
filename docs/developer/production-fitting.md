# Production fitting

Potential Analysis runs the Joint Component Estimator. It consumes map responses
and atomic geometry. Simulation metadata supplies the model resolution, and a
simulation request skips map normalization. The fit uses neither simulation
manifests nor truth parameters.

1. Prepare the selected non-hydrogen Targets, check the model and map
   fingerprints, and build the Joint observation domains.
2. Sample FirstStage responses on deterministic spheres. The native local
   OLS/MDPDE fit supplies initial widths. Global calibration uses only qualified
   Target donors and supplies the width bound.
3. Partition the Joint problem by structural component. Gauss–Seidel block
   coordinate descent uses constrained A/C variable projection and bounded
   active-set LM for width updates.
4. Verify global practical convergence, endpoint feasibility, and independent
   replay. Record Local Conditioning and Target Point Quality.
5. Publish qualified Target points and reconstruct Halo profiles. Pair
   Original-grid raw and peeled responses, correct the target charge, and fit
   Post-Joint Local MDPDE and eligible Group results.
6. Save the Joint snapshot and stage records. An unqualified run saves its
   diagnostics with the corresponding command outcome.

The [command contract](commands/potential_analysis.md) defines input, selection,
and outcome semantics. [Runtime](joint-component-runtime.md) and
[Numerics](joint-component-numerics.md) specify the solver and qualification.
[Stage Integration](joint-stage-integration.md) specifies FirstStage and
downstream data ownership.
