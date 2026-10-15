#pragma once
#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace rhbm_gem {
inline constexpr std::string_view FullABCHaloContributorPolicy="fullabc-eligible-halo-pruning-v1";
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
struct JointSelectionDomain
{
    std::string contract{"fixed-selected-voxel-closure-v1"};
    std::vector<std::size_t> target_indices;
    double observation_radius{2.5},support_radius{2.5};
    std::string contributor_policy{"all-non-hydrogen"};
    std::optional<JointHaloEligibilityProvenance> halo_eligibility;
    bool IsValid(std::size_t atoms) const
    {
        if(contract!="fixed-selected-voxel-closure-v1" ||
            (contributor_policy!="all-non-hydrogen" && contributor_policy!=FullABCHaloContributorPolicy) ||
            observation_radius!=2.5 || support_radius!=2.5 || target_indices.empty()) return false;
        for(std::size_t i=0;i<target_indices.size();++i)
            if(target_indices[i]>=atoms || (i && target_indices[i-1]>=target_indices[i])) return false;
        if(contributor_policy=="all-non-hydrogen") return !halo_eligibility;
        if(!halo_eligibility) return false;
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
struct JointAtomIdentifiability
{
    std::size_t atom{};
    JointCheckStatus status{JointCheckStatus::Unavailable};
    std::optional<double> null_space_leakage;
    std::string reason;
};
struct JointTargetEvidence
{
    std::string contract{"target-estimability-v4"};
    std::vector<JointCheck> checks;
    std::vector<JointRankEvidence> ranks;
    std::vector<JointAtomIdentifiability> atoms;
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
    std::string reason;
    std::optional<double> original_b,used_b;
    std::string seed_source{"not-recorded"};
    std::size_t donor_count{};
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
    std::optional<bool> candidate_feasible;
    JointCheckStatus reference_status{JointCheckStatus::NotRun};
    JointCheckStatus rank_assessment_status{JointCheckStatus::NotRun};
    JointCheckStatus certification_status{JointCheckStatus::NotRun};
    JointCheckStatus endpoint_verification_status{JointCheckStatus::NotRun};
    JointCheckStatus numerical_qualification_status{JointCheckStatus::NotRun};
    std::string verification_reason;
    std::string resource_limit_reason;
};
struct JointBlockConditioning
{
    std::string block_id,status{"not-assessed"},reason;
    std::size_t core_atom_count{},affected_row_count{},free_ac_columns{},active_amplitude_count{};
    std::optional<std::size_t> free_ac_rank,projected_width_rank,profile_rank,workspace_peak_memory_estimate_bytes;
    std::optional<double> free_ac_sigma_min,free_ac_sigma_max,free_ac_condition,
        projected_width_condition,profile_condition;
    std::string free_ac_status{"not-assessed"},free_ac_reason,
        projected_width_status{"not-assessed"},projected_width_reason,
        profile_status{"not-assessed"},profile_reason;
};
struct JointBlockParameterStability
{
    std::string block_id,status{"NotAssessed"},reason;
    std::string prediction_status{"NotAssessed"},objective_consistency_status{"NotAssessed"};
    std::optional<double> log_width_parameter_change_inf,scaled_ac_parameter_change_inf,
        local_prediction_change_relative,objective_difference;
    double cpu_seconds{};
    std::size_t perturbation_solves{},perturbation_evaluations{};
};
struct JointParameterStability
{
    std::size_t recent_sweep_count{};
    std::optional<double> maximum_eta_change,maximum_scaled_ac_change;
    std::string trajectory_stability_status{"insufficient-history"},perturbation_status{"NotAssessed"};
    std::vector<JointBlockParameterStability> blocks;
    double seconds{},cpu_seconds{};
    std::size_t perturbation_solves{},perturbation_evaluations{};
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
    std::vector<JointRankEvidence> ranks;
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
