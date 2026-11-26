#pragma once
#include <algorithm>
#include <array>
#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>

namespace rhbm_gem {
enum class FittingRole { NotRecorded, OutputTarget, BufferTarget, ExternalHalo };
inline constexpr bool IsFittingTarget(FittingRole role)
{
    return role==FittingRole::OutputTarget || role==FittingRole::BufferTarget;
}

inline constexpr std::string_view BufferExpansionContributorPolicy="buffer-expansion-external-halo-exclusion-v1";
inline constexpr std::string_view JointBufferExpansionParameterizationContract=
    "fixed-output-domain-buffer-expansion-v1";
inline constexpr std::string_view ContributorReconstructionContractVersion=
    "quality-qualified-contributor-reconstruction-v1";
inline constexpr std::size_t GeometryLocalFittingMedianNeighborCount=4;
inline constexpr std::string_view ContributorReconstructionDonorPolicy=
    "interior-adequate-joint-fitting-points-knn4-unbounded-euclidean-angstrom-distance-atom-id-tie-v1";
struct JointBufferExpansionProvenance
{
    std::size_t original_rows{},fit_rows{},excluded_rows{};
    std::size_t non_output_candidates{},external_halos{},zero_support_halos{};
    std::size_t recovered_rows{};
    bool operator==(const JointBufferExpansionProvenance &) const = default;
};
struct JointSelectionDomain
{
    std::string contract{"fixed-output-voxel-closure-v2"};
    // Both role sets use the contributor indices of the snapshot that owns them.
    std::vector<std::size_t> output_indices,buffer_indices;
    double observation_radius{2.5},support_radius{2.5};
    std::string contributor_policy{"all-non-hydrogen"};
    std::optional<JointBufferExpansionProvenance> buffer_expansion;
    bool operator==(const JointSelectionDomain &) const = default;

    FittingRole RoleOf(std::size_t atom) const
    {
        if(std::binary_search(output_indices.begin(),output_indices.end(),atom))
            return FittingRole::OutputTarget;
        if(!buffer_indices.empty() && buffer_indices.front()==atom) return FittingRole::BufferTarget;
        return FittingRole::ExternalHalo;
    }

    std::size_t FittingCount() const { return output_indices.size()+buffer_indices.size(); }

    bool MatchesFittingContributors(std::span<const std::size_t> mapping) const
    {
        if(mapping.size()!=FittingCount()) return false;
        for(std::size_t i=0;i<mapping.size();++i)
            if((i && mapping[i-1]>=mapping[i]) || !IsFittingTarget(RoleOf(mapping[i]))) return false;
        return true;
    }

    bool IsValid(std::size_t atoms) const
    {
        if(contract!="fixed-output-voxel-closure-v2" ||
            (contributor_policy!="all-non-hydrogen" &&
                contributor_policy!=BufferExpansionContributorPolicy) ||
            observation_radius!=2.5 || support_radius!=2.5 || output_indices.empty()) return false;
        for(std::size_t i=0;i<output_indices.size();++i)
            if(output_indices[i]>=atoms || (i && output_indices[i-1]>=output_indices[i])) return false;
        if(buffer_indices.size()>1 || (!buffer_indices.empty() &&
            (buffer_indices.front()>=atoms || RoleOf(buffer_indices.front())!=FittingRole::BufferTarget)))
            return false;
        if(contributor_policy=="all-non-hydrogen") return !buffer_expansion;
        if(contributor_policy==BufferExpansionContributorPolicy)
        {
            if(!buffer_expansion) return false;
            const auto & policy=*buffer_expansion;
            return policy.original_rows>=policy.fit_rows && policy.excluded_rows==policy.original_rows-policy.fit_rows &&
                policy.external_halos==atoms-FittingCount() &&
                policy.non_output_candidates>=atoms-output_indices.size() &&
                policy.zero_support_halos==policy.non_output_candidates-(atoms-output_indices.size()) &&
                policy.recovered_rows<=policy.fit_rows &&
                (policy.recovered_rows>0)==!buffer_indices.empty() &&
                policy.fit_rows>0;
        }
        return false;
    }
};
enum class JointCheckStatus {Passed,Failed,Unavailable,NotRun};
enum class JointEvidenceScope {ComponentLocal,AssembledGlobal};
// All indices refer to the original immutable problem, including component layouts.
struct JointParameterLayout
{
    std::vector<std::size_t> full_atoms,informative_rows;
};
struct JointCheck
{
    std::string name;
    JointCheckStatus status{JointCheckStatus::NotRun};
    JointEvidenceScope scope{JointEvidenceScope::ComponentLocal};
    std::optional<double> value,threshold;
    std::string reason;
};
inline constexpr std::string_view JointTargetPointValidationContractVersion="target-point-validation-v2";
struct JointAtomPointValidation
{
    std::size_t atom{};
    JointCheckStatus component_convergence{JointCheckStatus::Unavailable},parameter_validity{JointCheckStatus::Unavailable};
    std::string local_conditioning{"NotAssessed"},point_quality{"Unavailable"};
    std::vector<std::string> local_conditioning_blocks;
    std::string reason;
};
struct JointTargetEvidence
{
    std::string contract{JointTargetPointValidationContractVersion};
    std::vector<JointCheck> checks;
    std::vector<JointAtomPointValidation> point_atoms;
    std::size_t original_rows{};
};
struct JointInitializationAtom
{
    std::string id;
    std::array<double,3> ols{},mdpde{};
    double alpha{};
    std::size_t sample_count{};
    std::optional<int> native_status;
    FittingRole role{FittingRole::NotRecorded};
    bool initial_stage_fit_qualified{};
    std::string reason;
    std::optional<double> original_tau,used_tau;
    std::string seed_source{"not-recorded"};
    std::size_t donor_count{};
    std::optional<RHBMSolveQualification> initial_stage_qualification{};
    std::optional<RHBMFailedOnlyRecoveryDiagnostics> initial_stage_recovery{};
};
inline RHBMSolveQualification InitialStageQualification(
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
inline bool IsQualifiedInitialStageEstimate(
    std::optional<int> native_status, bool fit_qualified,
    const std::optional<RHBMSolveQualification> & qualification,
    const std::optional<RHBMFailedOnlyRecoveryDiagnostics> & recovery)
{
    if (!fit_qualified || !native_status) return false;
    const auto effective = InitialStageQualification(native_status, fit_qualified, qualification, recovery);
    if (effective == RHBMSolveQualification::NativeSuccess)
        return *native_status == static_cast<int>(RHBMEstimationStatus::SUCCESS);
    return effective == RHBMSolveQualification::RefinedSuccess &&
        *native_status != static_cast<int>(RHBMEstimationStatus::SUCCESS) &&
        recovery && recovery->accepted;
}
enum class ReconstructionSource { JointPoint, GeometryLocalFittingMedian, Unavailable };
enum ReconstructionRejectionReason : unsigned
{
    ExternalHaloReconstruction = 1,
    WidthLowerBoundReconstruction = 2,
    WidthUpperBoundReconstruction = 4,
    ConditioningNotAdequateReconstruction = 8
};
struct JointContributorReplacement
{
    std::size_t contributor_index{};
    std::string atom_id;
    FittingRole role{FittingRole::NotRecorded};
    unsigned rejection_reasons{};
    double amplitude{},width{},charge{};
    ReconstructionSource source{ReconstructionSource::GeometryLocalFittingMedian};
    std::size_t donor_count{};
    std::string joint_run_id;
    std::string contract{ContributorReconstructionContractVersion};
    std::vector<std::string> donor_atom_ids;
    bool operator==(const JointContributorReplacement &) const = default;
};
struct JointContributorReconstructionSnapshot
{
    std::string contract{ContributorReconstructionContractVersion};
    std::string joint_run_id;
    std::size_t qualified_donor_count{};
    std::vector<JointContributorReplacement> replacements;
    std::size_t neighbor_count{GeometryLocalFittingMedianNeighborCount};
    std::string donor_policy{ContributorReconstructionDonorPolicy};
    bool operator==(const JointContributorReconstructionSnapshot &) const = default;
};
struct JointInitialization
{
    bool valid{};
    std::string reason;
    std::vector<double> tau;
    std::vector<JointInitializationAtom> atoms;
    std::string data_scope{"caller-provided-widths"};
};
struct JointState
{
    // tau = exp(eta); log_tau = eta; eta_gradient = dF/deta (not dF/dtau).
    std::vector<double> ac,tau,log_tau,eta_gradient;
    // 0.5 * squared residual norm / parent ObservationScale() squared.
    double objective{};
};
// A parameter-only snapshot of the latest accepted search candidate. This is
// diagnostic data; the component's `state` remains the trusted scientific state.
struct JointCandidateDiagnostic
{
    bool candidate_available{};
    std::vector<double> candidate_eta,candidate_ac;
    std::optional<double> candidate_objective,candidate_ac_kkt,candidate_eta_gradient_inf_norm;
    std::string convergence_criteria_id;
    std::optional<double> global_ac_kkt_threshold,global_eta_gradient_threshold;
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
    std::vector<double> local_tau,projected_width_column_norms;
    std::string free_ac_status{"not-assessed"},free_ac_reason,
        projected_width_status{"not-assessed"},projected_width_reason,
        profile_status{"not-assessed"},profile_reason;
    std::string free_ac_conditioning{"NotAssessed"},projected_width_conditioning{"NotAssessed"},
        profile_conditioning{"NotAssessed"};
};
struct JointComponentData
{
    std::string id,stop_reason;
    std::vector<std::size_t> atoms,rows;
    bool search_completed{};
    std::size_t search_sweeps{},first_order_stationarity_sweep{},confirmed_stationarity_sweep{},
        endpoint_certification_sweep{};
    int profile_evaluations{},accepted_updates{},native_status{};
    std::optional<JointState> state;
    std::optional<JointCandidateDiagnostic> candidate_diagnostic;
    std::vector<JointCheck> evidence;
    std::vector<JointBlockConditioning> local_conditioning;
    double local_conditioning_seconds{};
    std::optional<std::size_t> local_conditioning_peak_memory_estimate_bytes;
    // Absent means component atom indices map directly to the state.
    std::optional<JointParameterLayout> layout;
    std::optional<JointTargetEvidence> target_evidence;
};
struct JointFitCosts
{
    double initialization_seconds{},search_seconds{},assessment_seconds{},assembly_seconds{};
    double construction_seconds{};
};
}
