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
class JointProblem
{
public:
    explicit JointProblem(JointProblemInput);
    const JointProblemInput & Input() const;
    double ObservationScale() const;
private:
    std::shared_ptr<const joint_component::ProblemData> m_data;
    friend struct JointProblemAccess;
};
using ::rhbm_gem::JointCheckStatus;
using ::rhbm_gem::JointEvidenceScope;
using ::rhbm_gem::JointCheck;
using ::rhbm_gem::JointRankEvidence;
using ::rhbm_gem::JointInitializationAtom;
using ::rhbm_gem::JointInitialization;
using ::rhbm_gem::JointState;
using ::rhbm_gem::JointFitCosts;
struct JointComponentResult : JointComponentData
{
    JointCheckStatus RuntimeConvergence() const;
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
    std::vector<JointRankEvidence> ranks;
    bool search_completed{};
    double observation_scale{1};
    JointCheckStatus regular_certificate{JointCheckStatus::NotRun};
    // Actual-state numerical evidence only; independent of search termination and offline audits.
    JointCheckStatus RuntimeConvergence() const;
};
// Captures results without retaining the problem snapshot or running numerical checks.
JointAnalysisResult CaptureJointAnalysisResult(const JointFitResult &, JointAnalysisMetadata = {});
// Selected non-hydrogen atoms fix the rows; all non-hydrogen contributors on
// those rows are fitted, including unselected halo atoms.
JointProblem BuildJointProblem(const MapObject &,const ModelObject &);
JointFitResult FitJointComponents(const JointProblem &,const std::vector<double> & initial_b);
// Updates successful targets' first-stage analysis only. Joint estimates are returned, not written
// into the existing second-stage or group-fitting result fields. Replacing raw samples
// invalidates sample-derived peeling, while fixed Second and group results are retained.
JointFitResult EstimateJointComponents(MapObject &,ModelObject &);
}
}
