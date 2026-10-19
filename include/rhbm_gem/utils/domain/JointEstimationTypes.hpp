#pragma once
#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <vector>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>

namespace rhbm_gem {
inline constexpr std::string_view FullABCHaloContributorPolicy="fullabc-eligible-halo-pruning-v1";
inline constexpr std::string_view FixedHaloContributorPolicy="first-stage-fixed-halo-v1";
inline constexpr std::string_view PureHaloSupportExclusionContributorPolicy="pure-halo-support-exclusion-v1";
inline constexpr std::string_view JointPureHaloSupportExclusionParameterizationContract=
    "target-only-pure-halo-support-exclusion-v1";
inline constexpr std::string_view HaloMedianReconstructionContractVersion=
    "group-key-target-median-v1";
struct ExcludedHaloRecord
{
    std::string atom_id;
    std::size_t support_rows{},distinct_radial_distances{};
    std::string reason;
    bool operator==(const ExcludedHaloRecord &) const = default;
};
struct JointHaloEligibilityProvenance
{
    std::size_t minimum_observations{3},minimum_distinct_radial_information{3};
    std::size_t candidate_halos{},retained_halos{};
    std::vector<ExcludedHaloRecord> excluded_halos;
};
struct JointFixedHaloSelectionProvenance
{
    std::size_t minimum_observations{1};
    std::size_t candidate_halos{},retained_halos{},zero_support_halos{};
    bool operator==(const JointFixedHaloSelectionProvenance &) const = default;
};
struct JointPureHaloSupportExclusionProvenance
{
    std::size_t original_rows{},fit_rows{},excluded_rows{};
    std::size_t candidate_halos{},retained_halos{},zero_support_halos{};
    bool operator==(const JointPureHaloSupportExclusionProvenance &) const = default;
};
struct JointSelectionDomain
{
    std::string contract{"fixed-selected-voxel-closure-v1"};
    std::vector<std::size_t> target_indices;
    double observation_radius{2.5},support_radius{2.5};
    std::string contributor_policy{"all-non-hydrogen"};
    std::optional<JointHaloEligibilityProvenance> halo_eligibility;
    std::optional<JointFixedHaloSelectionProvenance> fixed_halo_selection;
    std::optional<JointPureHaloSupportExclusionProvenance> pure_halo_support_exclusion;
    bool IsValid(std::size_t atoms) const
    {
        if(contract!="fixed-selected-voxel-closure-v1" ||
            (contributor_policy!="all-non-hydrogen" && contributor_policy!=FullABCHaloContributorPolicy &&
                contributor_policy!=FixedHaloContributorPolicy &&
                contributor_policy!=PureHaloSupportExclusionContributorPolicy) ||
            observation_radius!=2.5 || support_radius!=2.5 || target_indices.empty()) return false;
        for(std::size_t i=0;i<target_indices.size();++i)
            if(target_indices[i]>=atoms || (i && target_indices[i-1]>=target_indices[i])) return false;
        if(contributor_policy=="all-non-hydrogen") return !halo_eligibility && !fixed_halo_selection;
        if(contributor_policy==FixedHaloContributorPolicy)
        {
            if(halo_eligibility || pure_halo_support_exclusion || !fixed_halo_selection) return false;
            const auto & policy=*fixed_halo_selection;
            return policy.minimum_observations==1 &&
                policy.retained_halos==atoms-target_indices.size() &&
                policy.candidate_halos>=policy.retained_halos &&
                policy.zero_support_halos==policy.candidate_halos-policy.retained_halos;
        }
        if(contributor_policy==PureHaloSupportExclusionContributorPolicy)
        {
            if(halo_eligibility || fixed_halo_selection || !pure_halo_support_exclusion) return false;
            const auto & policy=*pure_halo_support_exclusion;
            return policy.original_rows>=policy.fit_rows && policy.excluded_rows==policy.original_rows-policy.fit_rows &&
                policy.retained_halos==atoms-target_indices.size() &&
                policy.candidate_halos>=policy.retained_halos &&
                policy.zero_support_halos==policy.candidate_halos-policy.retained_halos &&
                policy.fit_rows>0;
        }
        if(!halo_eligibility || fixed_halo_selection || pure_halo_support_exclusion) return false;
        const auto & policy=*halo_eligibility;
        if(policy.minimum_observations<3 || policy.minimum_distinct_radial_information!=3 ||
            policy.retained_halos!=atoms-target_indices.size() ||
            policy.candidate_halos<policy.retained_halos ||
            policy.excluded_halos.size()!=policy.candidate_halos-policy.retained_halos) return false;
        std::string_view previous_id;
        for(const auto & excluded:policy.excluded_halos)
        {
            if(excluded.atom_id.empty() || (!previous_id.empty() && previous_id>=excluded.atom_id)) return false;
            previous_id=excluded.atom_id;
            if(excluded.reason=="insufficient-observations")
            {if(excluded.support_rows>=policy.minimum_observations) return false;}
            else if(excluded.reason=="insufficient-radial-information")
            {if(excluded.support_rows<policy.minimum_observations ||
                excluded.distinct_radial_distances>=policy.minimum_distinct_radial_information) return false;}
            else return false;
        }
        return true;
    }
};
enum class JointCheckStatus {Passed,Failed,Unavailable,NotRun};
enum class JointEvidenceScope {ComponentLocal,AssembledGlobal};
struct JointNuisanceGroup
{
    std::size_t row{};
    std::vector<std::size_t> atoms;
};
// All indices refer to the original immutable problem, including component layouts.
// Newly constructed problems use every atom and observation row; groups remain
// solely to decode and re-export historical singleton-profile results.
struct JointParameterLayout
{
    std::vector<std::size_t> full_atoms,informative_rows;
    std::vector<JointNuisanceGroup> groups;
};
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
    std::string context_identity;
};
inline constexpr std::string_view JointTargetPointValidationContractVersion="target-point-validation-v1";
struct JointAtomIdentifiability
{
    std::size_t atom{};
    JointCheckStatus status{JointCheckStatus::Unavailable};
    std::optional<double> null_space_leakage;
    std::string reason;
};
struct JointAtomPointValidation
{
    std::size_t atom{};
    JointCheckStatus component_convergence{JointCheckStatus::Unavailable},parameter_validity{JointCheckStatus::Unavailable};
    std::string local_conditioning{"NotAssessed"},parameter_stability{"NotAssessed"},point_quality{"Unavailable"};
    std::vector<std::string> local_conditioning_blocks,parameter_stability_blocks;
    std::string reason;
};
struct JointTargetEvidence
{
    std::string contract{JointTargetPointValidationContractVersion};
    std::vector<JointCheck> checks;
    std::vector<JointRankEvidence> ranks;
    // Historical target-estimability evidence retains its original meaning.
    std::vector<JointAtomIdentifiability> atoms;
    std::vector<JointAtomPointValidation> point_atoms;
    std::vector<double> column_scales;
    std::size_t original_rows{};
};
struct JointInitializationAtom
{
    std::string id;
    std::array<double,3> ols{},mdpde{};
    double alpha{};
    std::size_t sample_count{};
    std::optional<int> native_status;
    bool target{},first_stage_fit_qualified{};
    std::string reason;
    std::optional<double> original_b,used_b;
    std::string seed_source{"not-recorded"};
    std::size_t donor_count{};
    std::optional<RHBMSolveQualification> first_stage_qualification{};
    std::optional<RHBMFailedOnlyRecoveryDiagnostics> first_stage_recovery{};
};
inline constexpr std::string_view FixedHaloParameterSourceContractVersion="first-stage-mdpde-fixed-halo-v1";
struct JointFixedHaloParameter
{
    std::size_t contributor_index{};
    std::string atom_id;
    double amplitude{},width{},charge{};
    std::size_t first_stage_sample_count{};
    std::optional<int> first_stage_status;
    std::string source_contract{FixedHaloParameterSourceContractVersion};
    bool first_stage_fit_qualified{};
    std::optional<RHBMSolveQualification> first_stage_qualification{};
    std::optional<RHBMFailedOnlyRecoveryDiagnostics> first_stage_recovery{};
    bool operator==(const JointFixedHaloParameter & other) const;
};
inline RHBMSolveQualification FirstStageQualification(
    std::optional<int> native_status, bool fit_qualified,
    const std::optional<RHBMSolveQualification> & qualification,
    const std::optional<RHBMFailedOnlyRecoveryDiagnostics> & recovery)
{
    if (qualification) return *qualification;
    if (recovery) return recovery->accepted ? RHBMSolveQualification::RefinedSuccess :
        RHBMSolveQualification::Unqualified;
    return fit_qualified && native_status && *native_status == static_cast<int>(RHBMEstimationStatus::SUCCESS) ?
        RHBMSolveQualification::NativeSuccess : RHBMSolveQualification::Unqualified;
}
inline bool IsQualifiedFirstStageEstimate(
    std::optional<int> native_status, bool fit_qualified,
    const std::optional<RHBMSolveQualification> & qualification,
    const std::optional<RHBMFailedOnlyRecoveryDiagnostics> & recovery)
{
    if (!fit_qualified || !native_status) return false;
    const auto effective = FirstStageQualification(native_status, fit_qualified, qualification, recovery);
    if (effective == RHBMSolveQualification::NativeSuccess)
        return *native_status == static_cast<int>(RHBMEstimationStatus::SUCCESS);
    return effective == RHBMSolveQualification::RefinedSuccess &&
        *native_status != static_cast<int>(RHBMEstimationStatus::SUCCESS) &&
        recovery && recovery->accepted;
}
inline bool JointFixedHaloParameter::operator==(const JointFixedHaloParameter & other) const
{
    return contributor_index==other.contributor_index && atom_id==other.atom_id &&
        amplitude==other.amplitude && width==other.width && charge==other.charge &&
        first_stage_sample_count==other.first_stage_sample_count &&
        first_stage_status==other.first_stage_status && source_contract==other.source_contract &&
        first_stage_fit_qualified==other.first_stage_fit_qualified &&
        FirstStageQualification(first_stage_status,first_stage_fit_qualified,
            first_stage_qualification,first_stage_recovery)==
            FirstStageQualification(other.first_stage_status,other.first_stage_fit_qualified,
                other.first_stage_qualification,other.first_stage_recovery) &&
        first_stage_recovery==other.first_stage_recovery;
}
inline bool IsQualifiedFirstStageFixedHalo(const JointFixedHaloParameter & fixed)
{
    return IsQualifiedFirstStageEstimate(fixed.first_stage_status,
        fixed.first_stage_fit_qualified, fixed.first_stage_qualification, fixed.first_stage_recovery);
}
enum class HaloMedianSource { GroupKeyMedian, GlobalTargetMedian };
struct JointHaloMedianReconstruction
{
    std::size_t halo_contributor_index{};
    std::string atom_id;
    GroupKey group_key{};
    double amplitude{},width{},charge{};
    HaloMedianSource source{HaloMedianSource::GlobalTargetMedian};
    std::size_t donor_count{};
    std::string joint_run_id;
    std::string contract{HaloMedianReconstructionContractVersion};
    bool operator==(const JointHaloMedianReconstruction &) const = default;
};
struct JointHaloMedianSnapshot
{
    std::string contract{HaloMedianReconstructionContractVersion};
    std::string joint_run_id;
    std::array<double,3> global_median{};
    std::size_t valid_target_count{};
    std::vector<JointHaloMedianReconstruction> halos;
    bool operator==(const JointHaloMedianSnapshot &) const = default;
};
struct JointInitialization
{
    bool valid{};
    std::string reason;
    std::vector<double> b;
    std::vector<JointInitializationAtom> atoms;
    std::string data_scope{"caller-provided-widths"};
};
struct JointState
{
    std::vector<double> ac,b,log_b,width_gradient;
    // 0.5 * squared residual norm / parent ObservationScale() squared.
    double objective{};
    // Historical singleton-profile amplitudes; current FullABC states leave this empty.
    std::vector<double> nuisance_amplitudes{};
};
// A parameter-only snapshot of the latest accepted search candidate. This is
// diagnostic data; the component's `state` remains the trusted scientific state.
struct JointCandidateDiagnostic
{
    bool candidate_available{};
    std::vector<double> candidate_eta,candidate_ac;
    std::optional<double> candidate_objective,candidate_ac_kkt,candidate_width_gradient_inf_norm;
    std::string convergence_criteria_id;
    std::optional<double> global_ac_kkt_threshold,global_width_gradient_threshold,
        confirmed_log_width_change_threshold;
    std::optional<double> candidate_log_width_change_inf_norm;
    std::optional<bool> candidate_feasible;
    JointCheckStatus certification_status{JointCheckStatus::NotRun};
    JointCheckStatus endpoint_verification_status{JointCheckStatus::NotRun};
    JointCheckStatus numerical_qualification_status{JointCheckStatus::NotRun};
    std::string verification_reason;
    std::string resource_limit_reason;
};
struct JointBlockConditioning
{
    std::string block_id,status{"not-assessed"},reason;
    std::string conditioning_classification{"NotAssessed"};
    std::vector<std::size_t> core_atoms;
    std::size_t core_atom_count{},affected_row_count{},free_ac_columns{},active_amplitude_count{};
    std::optional<std::size_t> free_ac_rank,projected_width_rank,profile_rank,active_width_bound_targets,
        workspace_peak_memory_estimate_bytes;
    std::optional<double> free_ac_sigma_min,free_ac_sigma_max,free_ac_condition,
        free_ac_rank_threshold,free_ac_rank_threshold_margin_relative,
        local_ac_projected_kkt,local_projected_width_kkt,
        projected_width_sigma_min,projected_width_sigma_max,projected_width_condition,
        projected_width_rank_threshold,projected_width_rank_threshold_margin_relative,
        profile_sigma_min,profile_sigma_max,profile_condition,
        profile_rank_threshold,profile_rank_threshold_margin_relative;
    std::optional<bool> free_ac_threshold_sensitive,projected_width_threshold_sensitive,profile_threshold_sensitive;
    std::string free_ac_status{"not-assessed"},free_ac_reason,
        projected_width_status{"not-assessed"},projected_width_reason,
        profile_status{"not-assessed"},profile_reason;
    std::string free_ac_conditioning{"NotAssessed"},projected_width_conditioning{"NotAssessed"},
        profile_conditioning{"NotAssessed"};
};
struct JointStabilityProbe
{
    std::string direction,source,status;
    std::string minus_seed_status{"NotAssessed"},plus_seed_status{"NotAssessed"},
        minus_seed_reason,plus_seed_reason;
    std::vector<double> unit_direction;
    std::size_t evaluations{};
    std::optional<double> parameter_change,prediction_change,objective_change;
    std::optional<double> minus_requested_displacement_inf,minus_effective_displacement_inf,
        plus_requested_displacement_inf,plus_effective_displacement_inf;
};
struct JointBlockParameterStability
{
    std::string block_id,status{"NotAssessed"},reason;
    bool probe_selected{};
    std::vector<std::size_t> core_atoms;
    std::vector<JointStabilityProbe> probes;
    std::string prediction_status{"NotAssessed"},objective_consistency_status{"NotAssessed"};
    std::optional<double> log_width_parameter_change_inf,scaled_ac_parameter_change_inf,
        local_prediction_change_relative,objective_difference;
    std::optional<double> baseline_local_ac_kkt,baseline_projected_width_kkt,baseline_width_bound_slack,
        baseline_local_objective,baseline_local_residual_rmse_normalized,conditional_reprofile_local_ac_kkt,
        conditional_reprofile_projected_width_kkt,conditional_reprofile_objective_difference,
        conditional_reprofile_scaled_ac_parameter_difference,conditional_reprofile_prediction_change_relative;
    std::optional<std::size_t> baseline_active_ac_constraints,baseline_active_width_bound_targets;
    std::string baseline_numerical_validity_status{"NotAssessed"},
        baseline_conditional_optimality_status{"NotAssessed"},baseline_eligibility_status{"NotAssessed"},
        baseline_eligibility_reason;
    std::string baseline_independent_replay_status{"NotRun"},baseline_independent_replay_reason,
        conditional_reprofile_status{"NotRun"},conditional_reprofile_independent_replay_status{"NotRun"},
        conditional_reprofile_independent_replay_reason;
    double cpu_seconds{};
    std::size_t perturbation_directions{},perturbation_seeds{},perturbation_solves{},perturbation_evaluations{};
};
struct JointParameterStability
{
    std::size_t recent_sweep_count{};
    std::optional<double> maximum_eta_change,maximum_scaled_ac_change;
    std::string trajectory_stability_status{"insufficient-history"},perturbation_status{"NotAssessed"};
    std::vector<JointBlockParameterStability> blocks;
    double seconds{},cpu_seconds{};
    std::size_t perturbation_directions{},perturbation_seeds{},perturbation_solves{},perturbation_evaluations{};
    std::optional<std::size_t> workspace_peak_memory_estimate_bytes;
};
struct JointComponentData
{
    std::string id,stop_reason;
    std::vector<std::size_t> atoms,rows;
    bool search_completed{};
    std::size_t search_sweeps{},first_order_stationarity_sweep{},confirmed_stationarity_sweep{},
        endpoint_certification_sweep{};
    int profile_evaluations{},reference_evaluations{},accepted_updates{},native_status{};
    std::optional<JointState> state;
    std::optional<JointCandidateDiagnostic> candidate_diagnostic;
    std::vector<JointCheck> evidence;
    std::vector<JointBlockConditioning> local_conditioning;
    double local_conditioning_seconds{};
    std::optional<std::size_t> local_conditioning_peak_memory_estimate_bytes;
    std::optional<JointParameterStability> parameter_stability;
    JointCheckStatus regular_certificate{JointCheckStatus::NotRun};
    // Absent means the historical identity layout (all component atoms are FullABC).
    std::optional<JointParameterLayout> layout;
    std::optional<JointTargetEvidence> target_evidence;
};
struct JointFitCosts
{
    double initialization_seconds{},search_seconds{},search_reference_seconds{},assessment_seconds{},assembly_seconds{};
    double construction_seconds{};
};
}
