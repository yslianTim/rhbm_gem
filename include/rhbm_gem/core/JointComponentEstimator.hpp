#pragma once
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
enum class JointCheckStatus {Passed,Failed,Unavailable,NotRun};
enum class JointEvidenceScope {ComponentLocal,AssembledGlobal};
struct JointCheck
{
    std::string name;
    JointCheckStatus status{JointCheckStatus::NotRun};
    JointEvidenceScope scope{JointEvidenceScope::ComponentLocal};
    std::optional<double> value,threshold;
    std::string reason;
};
struct JointRankEvidence
{
    std::string name;
    JointEvidenceScope scope{JointEvidenceScope::ComponentLocal};
    std::size_t rank{};
    double threshold{};
    std::vector<double> singular_values;
};
struct JointInitializationAtom
{
    std::string id;
    std::array<double,3> ols{},mdpde{};
    double alpha{};
    std::size_t sample_count{};
    std::optional<int> native_status;
};
struct JointInitialization
{
    bool valid{};
    std::string reason;
    std::vector<double> b;
    std::vector<JointInitializationAtom> atoms;
};
struct JointState
{
    std::vector<double> ac,b,log_b,width_gradient;
    // 0.5 * squared residual norm / parent ObservationScale() squared.
    double objective{};
};
struct JointComponentResult
{
    std::string id,stop_reason;
    std::vector<std::size_t> atoms,rows;
    bool search_completed{};
    int profile_evaluations{},reference_evaluations{},accepted_updates{},native_status{};
    std::optional<JointState> state;
    std::vector<JointCheck> evidence;
    std::vector<JointRankEvidence> ranks;
    JointCheckStatus regular_certificate{JointCheckStatus::NotRun};
    // Actual-state numerical evidence only; independent of search termination and offline audits.
    JointCheckStatus RuntimeConvergence() const;
};
struct JointFitCosts
{
    double initialization_seconds{},search_seconds{},search_reference_seconds{},assessment_seconds{},assembly_seconds{};
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
// V1 includes every non-hydrogen atom. A partial non-hydrogen selection is rejected.
JointProblem BuildJointProblem(const MapObject &,const ModelObject &);
JointFitResult FitJointComponents(const JointProblem &,const std::vector<double> & initial_b);
// Updates first-stage analysis only. Joint estimates are returned, not written
// into the existing second-stage or group-fitting result fields.
JointFitResult EstimateJointComponents(MapObject &,ModelObject &);
}
}
