#pragma once

#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/data/object/JointAnalysisResult.hpp>
#include <cstddef>
#include <array>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace rhbm_gem {
class ModelObject;
namespace core { class JointProblem; }
namespace core::detail {

struct HaloMedianCoverageRecord
{
    std::string halo_atom_id;
    std::size_t target_donor_count{};
    std::size_t affected_support_rows{};
    std::array<double,3> median{};
    std::array<double,3> dispersion{};
};

struct HaloMedianCoverageAudit
{
    std::size_t valid_target_count{};
    std::size_t geometry_local_halos{};
    std::array<double,3> global_median{};
    std::vector<HaloMedianCoverageRecord> halos;
};

struct HaloMedianTargetDonor
{
    std::size_t contributor_index{};
    std::string atom_id;
    std::array<double,3> position{};
    std::array<double,3> point{};
};

struct HaloMedianBuildResult
{
    JointHaloMedianSnapshot snapshot;
    HaloMedianCoverageAudit audit;
    std::vector<HaloMedianTargetDonor> target_donors;
};

enum class HaloFullABCRole {FixedGlobalTargetMedian,FreeFullABC};
enum class HaloFullABCParameterSource
{
    Unavailable,PreFixedGlobalMedian,FullABC,ObservationACFixedWidthRecovery,PostFailureGlobalMedianRecovery
};
enum class HaloFullABCRecoveryMethod
{
    None,FixedWidthACPerHalo,FixedWidthACPerComponent,PostFailureGlobalMedian,IsolateFailedComponent
};
struct HaloFullABCClassificationRecord
{
    std::string atom_id;
    std::optional<std::size_t> original_contributor_index;
    std::size_t distinct_original_observations{};
    HaloFullABCRole role{HaloFullABCRole::FixedGlobalTargetMedian};
    std::vector<std::size_t> original_support_rows;
    std::optional<std::array<double,3>> parameter_profile;
    bool operator==(const HaloFullABCClassificationRecord &) const = default;
};
struct HaloFullABCPreclassification
{
    bool available{};
    std::string unavailable_reason;
    std::optional<std::array<double,3>> global_target_median;
    std::size_t fixed_halo_count{},free_halo_count{};
    std::vector<HaloFullABCClassificationRecord> halos;
    bool operator==(const HaloFullABCPreclassification &) const = default;
};
struct HaloFullABCReconstructionHalo
{
    HaloFullABCClassificationRecord classification;
    std::optional<std::array<double,3>> reconstructed_parameters;
    HaloFullABCParameterSource parameter_source{HaloFullABCParameterSource::Unavailable};
    bool qualified{};
    std::string qualification_reason;
};
struct HaloFullABCBlockFailureDiagnostic
{
    std::size_t component_index{},sweep_index{},block_index{},affected_rows{};
    std::string component_id,status,reason,local_search_stop_reason,initial_profile_reason,
        solver_reason,lm_subproblem_status,first_invalid_atom_id,first_invalid_parameter;
    std::vector<std::string> core_atom_ids;
    std::size_t local_design_rows{},local_design_columns{};
    std::optional<std::size_t> local_design_rank,first_invalid_coefficient;
    std::vector<std::size_t> active_amplitude_indices;
    std::vector<double> local_column_norms,local_design_singular_values,initial_widths,
        width_derivative_column_norms,projected_width_singular_values,projected_width_column_norms,proposed_eta_step;
    std::optional<double> local_design_rank_threshold,condition_indicator,local_response_norm,
        observation_scale,primal_feasibility,scaled_projected_kkt,width_upper_bound,width_step_inf,
        width_gradient_inf,lm_actual_decrease,lm_predicted_decrease,lm_ratio,weak_width_step_fraction,
        projected_width_rank_threshold,local_objective_before,local_objective_after,global_objective_before,
        global_objective_after,global_replay_delta,local_global_delta_error,local_final_ac_kkt,
        global_ac_kkt,global_width_gradient_inf_norm,replay_prediction_discrepancy,
        replay_gradient_discrepancy,replay_kkt_discrepancy;
    bool initial_profile_valid{},initial_profile_trusted{},solver_valid{},width_bound_active{};
};
struct HaloFullABCBlockSummary
{
    std::size_t sweep{},block{},affected_rows{};
    std::string status,reason,local_search_stop_reason;
    int profile_evaluations{},accepted_updates{};
    std::optional<double> local_objective_before,local_objective_after,global_objective_before,
        global_objective_after,global_replay_delta,local_global_delta_error,local_final_ac_kkt;
};
struct HaloFullABCCandidateProfile
{
    std::string atom_id;
    std::array<double,3> parameters{};
};
struct HaloFullABCComponentSearchDiagnostic
{
    std::string component_id,reason,endpoint_failure,trust_failure;
    bool search_converged{};
    std::size_t sweep_count{};
    std::optional<double> objective,global_ac_kkt,global_width_gradient_inf_norm,eta_change_inf,
        objective_replay_error,cache_replay_error;
    bool endpoint_verified{},endpoint_qualified{};
    std::vector<HaloFullABCBlockSummary> blocks;
    std::vector<JointBlockConditioning> local_conditioning;
    std::vector<HaloFullABCCandidateProfile> candidate_profiles;
};
struct HaloFullABCRecoveryDiagnostic
{
    std::string component_id,method,status,reason;
    std::size_t attempted_halos{},qualified_halos{};
    std::optional<double> primary_ac_kkt,independent_ac_kkt,objective_replay_error,ac_condition,fixed_width;
    std::optional<std::size_t> ac_rank,ac_columns;
    bool full_ac_rank{},independent_replay_passed{},fullabc_refit_attempted{},
        fullabc_refit_qualified{},qualified{};
    std::vector<std::string> atom_ids;
};
struct HaloFullABCJointDiagnostics
{
    bool solver_invoked{},solver_qualified{},reconstruction_qualified{},assembly_replay_qualified{};
    JointCheckStatus solver_status{JointCheckStatus::NotRun};
    std::string unavailable_reason;
    std::size_t original_observation_rows{},free_support_entries{},zero_sensitivity_rows{},
        qualified_first_stage_width_donors{};
    std::optional<double> first_stage_width_p95,width_upper_bound,objective,observation_scale;
    std::optional<double> assembled_prediction_replay_error,assembled_objective_replay_error,
        original_domain_objective,original_domain_normalization_scale;
    std::optional<std::size_t> estimated_peak_memory_bytes;
    std::optional<HaloFullABCBlockFailureDiagnostic> first_failing_block;
    std::vector<HaloFullABCComponentSearchDiagnostic> component_searches;
    std::string recovery_method;
    std::vector<HaloFullABCRecoveryDiagnostic> recovery;
    JointSolverProvenance solver_provenance;
    std::vector<JointCheck> evidence;
};
struct HaloFullABCReconstructionResult
{
    bool available{};
    std::string unavailable_reason,joint_run_id;
    std::optional<std::array<double,3>> global_target_median;
    std::size_t fixed_halo_count{},free_halo_count{};
    std::vector<HaloFullABCReconstructionHalo> halos;
    HaloFullABCJointDiagnostics diagnostics;
};
std::vector<std::size_t> DistinctOriginalSupportRows(const std::vector<JointSupport> &);

class GeometryLocalTargetDonorIndex
{
public:
    explicit GeometryLocalTargetDonorIndex(const std::vector<HaloMedianTargetDonor> &);
    ~GeometryLocalTargetDonorIndex();
    GeometryLocalTargetDonorIndex(const GeometryLocalTargetDonorIndex &) = delete;
    GeometryLocalTargetDonorIndex & operator=(const GeometryLocalTargetDonorIndex &) = delete;
    std::vector<const HaloMedianTargetDonor *> Query(const std::array<double,3> &) const;

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};
std::array<double,3> GeometryLocalTargetMedian(
    const std::vector<const HaloMedianTargetDonor *> &);
std::vector<HaloMedianTargetDonor> QualifiedJointTargetDonors(
    const ModelObject &,const JointProblem &,const JointAnalysisResult &);
HaloFullABCPreclassification BuildHaloFullABCPreclassification(
    const ModelObject &,const JointProblem &,const JointAnalysisResult &);
HaloFullABCReconstructionResult EstimateHaloFullABCJoint(
    const ModelObject &,const JointProblem &,const JointAnalysisResult &,std::string_view joint_run_id,
    bool collect_solver_diagnostics=false);
HaloFullABCReconstructionResult EstimateHaloFullABCJointWithRecovery(
    const ModelObject &,const JointProblem &,const JointAnalysisResult &,std::string_view joint_run_id,
    HaloFullABCRecoveryMethod);

// Builds the immutable post-Joint Halo reconstruction snapshot. Only target
// point evidence from this JointAnalysisResult may enter the donor pool.
HaloMedianBuildResult BuildHaloMedianReconstruction(
    const ModelObject &, const JointProblem &, const JointAnalysisResult &,
    std::string_view joint_run_id);

} // namespace core::detail
} // namespace rhbm_gem
