#pragma once
#include <rhbm_gem/data/object/JointAnalysisResult.hpp>
#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace rhbm_gem {
class MapObject;
class ModelObject;
namespace core {
namespace joint_component {struct ProblemData;}

struct JointSupport {std::size_t row{}; double squared_distance{};};
struct JointProblemInput
{
    std::vector<double> observations;
    std::vector<std::vector<JointSupport>> support;
    std::vector<std::string> atom_ids,row_ids;
    std::optional<JointSelectionDomain> selection_domain;
};
// Owns a validated immutable snapshot. Construction never reads historical fits.
// BuildJointProblem is the historical low-level fixed-halo research builder;
// initial target widths alone do not create its required FirstStage halo snapshot.
class JointProblem
{
public:
    explicit JointProblem(JointProblemInput);
    const JointProblemInput & Input() const;
    const JointProblemInput & OriginalInput() const;
    const std::vector<std::size_t> & FittingToOriginalContributors() const;
    const std::vector<std::size_t> & FittingToOriginalRows() const;
    const std::vector<JointFixedHaloParameter> & FixedHalos() const;
    // Attaches exactly one qualified FirstStage (A,B,0) value for every
    // supported non-target contributor. The returned problem keeps that
    // snapshot immutable during target-only fitting.
    JointProblem WithFixedHalos(std::vector<JointFixedHaloParameter>) const;
    JointProblem WithReconstructionDomain(JointProblemInput original,
        std::vector<std::size_t> fitting_to_original_contributors,
        std::vector<std::size_t> fitting_to_original_rows) const;
    double ObservationScale() const;
    const JointParameterLayout & ParameterLayout() const;
private:
    explicit JointProblem(std::shared_ptr<const joint_component::ProblemData>);
    std::shared_ptr<const joint_component::ProblemData> m_data;
    friend struct JointProblemAccess;
};
using ::rhbm_gem::JointCheckStatus;
using ::rhbm_gem::JointEvidenceScope;
using ::rhbm_gem::JointCheck;
using ::rhbm_gem::JointInitializationAtom;
using ::rhbm_gem::JointInitialization;
using ::rhbm_gem::JointState;
using ::rhbm_gem::JointFitCosts;
struct JointComponentResult : JointComponentData
{
    JointCheckStatus RuntimeConvergence() const;
    JointCheckStatus TargetRuntimeConvergence() const;
};
struct JointFitResult
{
    std::optional<JointProblem> problem;
    JointInitialization initialization;
    JointFitCosts costs;
    std::vector<JointComponentResult> components;
    std::optional<JointState> assembled_state;
    std::optional<std::vector<double>> prediction;
    // Includes constant rows; unavailable when any component state is missing.
    std::optional<double> objective;
    std::vector<bool> available_row_mask;
    std::vector<JointCheck> evidence;
    bool search_completed{};
    double observation_scale{1};
    JointCheckStatus regular_certificate{JointCheckStatus::NotRun};
    JointSolverProvenance solver_provenance;
    // Actual-state numerical evidence only; independent of search termination and offline audits.
    JointCheckStatus RuntimeConvergence() const;
    JointCheckStatus TargetRuntimeConvergence() const;
    std::optional<JointTargetEvidence> target_evidence;
};
// Captures results without retaining the problem snapshot or running numerical checks.
JointAnalysisResult CaptureJointAnalysisResult(const JointFitResult &, JointAnalysisMetadata = {});
// Historical fixed-halo problem builder. New research call sites should use
// the explicit BuildFixedHaloResearchProblem name.
JointProblem BuildJointProblem(const MapObject &,const ModelObject &);
JointProblem BuildFixedHaloResearchProblem(const MapObject &,const ModelObject &);
// Production target-only problem. Halo support is used only to exclude fitting
// rows; reconstruction keeps the original-domain input in OriginalInput().
JointProblem BuildPureHaloSupportExclusionProblem(const MapObject &,const ModelObject &);
// Fits only selected targets. For a fixed-halo problem, initial_b has one
// width per selected target and the problem must already own its complete
// qualified fixed-halo snapshot.
JointFitResult FitJointComponents(const JointProblem &,const std::vector<double> & initial_b);
// Validates FirstStage initialization and snapshots its qualified fixed halos
// before fitting selected targets. A raw BuildJointProblem result without
// qualified initialization cannot proceed to target optimization.
JointFitResult FitJointComponentsInitialized(const JointProblem &,const JointInitialization &);
// Builds the Production Pure Exclusion problem, performs FirstStage
// initialization, and then fits targets only. Inspect initialization validity
// separately from target convergence.
// Updates successful targets' first-stage analysis only. Joint estimates are returned, not written
// into the existing second-stage or group-fitting result fields. Replacing raw samples
// invalidates sample-derived peeling, while fixed Second and group results are retained.
JointFitResult EstimateJointComponents(MapObject &,ModelObject &);
// Explicit historical Fixed Halo research/compatibility entry point.
JointFitResult EstimateJointComponentsFixedHalo(MapObject &,ModelObject &);
}
}
