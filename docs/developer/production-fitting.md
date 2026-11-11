# Production fitting

Potential Analysis runs the Joint Component Estimator. It consumes stored voxel
Map Values and atomic geometry. Simulation metadata supplies the model resolution, and a
simulation request skips map normalization. The fit uses neither simulation
manifests nor truth parameters.

1. Prepare the selected non-hydrogen Targets, check the model and map
   fingerprints, and build the Joint observation domains.
2. **InitialStage** obtains Sampling Values by Tricubic interpolation of Map Values
   on deterministic spheres and fits native OLS/MDPDE estimates that supply initial widths.
3. **JointStage** uses Original-grid voxel observations read directly from
   `MapObject::GetMapValue()`, without interpolation. It calibrates the global
   width bound using qualified Target donors and partitions the problem by
   structural component. Gauss–Seidel block
   coordinate descent uses constrained A/C variable projection and bounded
   active-set LM for width updates. The stage verifies global practical convergence,
   endpoint feasibility, and independent replay, and records Local Conditioning
   and Target Point Quality. After qualification it reconstructs Halo profiles
   and publishes Target points, peels the Original grid, corrects target charge,
   and fits Post-Joint Local MDPDE.
4. **GroupingStage** uses eligible Post-Joint Local results for group alpha
   training, group MDPDE, and WEB posterior inference.
5. Save the Joint snapshot and stage records. An unqualified run saves its
   diagnostics with the corresponding command outcome.

The CLI and Python use the same public Command API. Installed C++ consumers
can run Potential Analysis through `CommandSystem.hpp` and `CommandTypes.hpp`.
The workflow and standalone Gaussian fitting routines live in the private
`src/core/estimator/` module. Map sampling and shared interpolation live in
`src/core/sampling/`.

`JointAnalysisResult`, the analysis views and editor, and the Gaussian and Joint
result types are part of the installed data API. Atom-local views expose
published points through `GetInitialStageModel()` and `GetJointStageModel()`.
Native InitialStage diagnostics are available through
`GetInitialStageGaussianResult()`.

The [command contract](commands/potential_analysis.md) defines input, selection,
and outcome semantics. [Runtime](joint-component-runtime.md) and
[Numerics](joint-component-numerics.md) specify the solver and qualification.
[Stage Integration](joint-stage-integration.md) specifies InitialStage and
downstream data ownership.

Map Value always refers to the current stored voxel value, including a normalized
map. Sampling Value always refers to an interpolation output. Observations and
responses are neutral solver/container terms; peeled responses, model predictions,
charge contributions and residuals retain their specific derived-value names.
