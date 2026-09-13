#pragma once

#include "core/detail/second_stage/CandidateSelection.hpp"

#include <memory>
#include <string>
#include <string_view>

namespace rhbm_gem::core::detail {

struct ClusterHistoryDiagnostic;

struct JointCandidateObjectiveDiagnostic
{
    std::size_t history_observation{ 0 };
    std::size_t locally_deteriorated_member_count{ 0 };
    double maximum_local_deterioration{ 0.0 };
    std::string_view source{};
    std::size_t round{ 0 };
    std::size_t candidate_number{ 0 };
    std::optional<double> factor{};
    ClusterKey member_key{};
    std::optional<ObjectiveBreakdown> previous{};
    std::optional<ObjectiveBreakdown> best{};
    std::optional<ObjectiveBreakdown> candidate{};
    std::optional<ObjectiveBreakdown> stored_best{};
    bool best_checked{ false };
    std::string_view outcome{ "accepted" };
    std::string best_source_id{};
    std::vector<std::string> best_comparison_lines{};
};

// Numerical evidence is copied here for output; production never reads this copy.
struct ObjectiveAttemptDiagnostic : CandidateDecisionEvidence
{
    std::optional<double> pre_objective_attempted_step_norm{};
    std::optional<ObjectiveScale> scale{};
    std::size_t fit_sample_count{ 0 };
    std::size_t tail_sample_count{ 0 };
    std::shared_ptr<const ClusterHistoryDiagnostic> history{};
    double trust_region_radius{ 0.0 };
    double trust_region_step_norm{ 0.0 };
    std::size_t trial_count{ 0 };
    std::size_t trust_skipped_trial_count{ 0 };
};

struct ClusterCandidateDiagnostic
{
    ClusterKey key{};
    ObjectiveAttemptDiagnostic attempt{};

    bool boundary_rescued{ false };
};

struct BoundaryComponentReconciliationDiagnostic
{
    std::vector<ClusterKey> key_list{};
    std::size_t atom_count{ 0 };
    std::size_t boundary_sample_count{ 0 };
    std::size_t trial_count{ 1 };
    std::optional<double> accepted_factor{};
    BoundaryComponentAcceptedSource accepted_source{ BoundaryComponentAcceptedSource::None };
    std::optional<BoundaryJointCorrectionStatus> joint_correction_status{};
    std::size_t interface_atom_count{ 0 };
    std::size_t shape_active_atom_count{ 0 };
    std::size_t offset_active_atom_count{ 0 };
    std::size_t suspicious_candidate_atom_count{ 0 };
    std::size_t joint_parameter_count{ 0 };
    std::optional<double> joint_damping{};
    std::optional<double> maximum_normalized_trust_step{};
    std::optional<double> previous_component_objective{};
    std::optional<double> endpoint_component_objective{};
    std::optional<double> joint_reference_component_objective{};
    std::optional<double> joint_candidate_component_objective{};
    std::optional<double> candidate_component_objective{};
    std::size_t accepted_cluster_count{ 0 };
    std::size_t rescue_candidate_cluster_count{ 0 };
    std::size_t rescued_cluster_count{ 0 };
    std::size_t locally_deteriorated_member_count{ 0 };
    double maximum_local_deterioration{ 0.0 };
    std::optional<double> component_improvement{};
    std::optional<double> global_improvement{};
    std::vector<JointCandidateObjectiveDiagnostic> objective_diagnostic_list{};
    bool is_rescue_attempt{ false };
    bool exhausted{ false };
};

struct FinalDependencyPolishDiagnostic
{
    std::size_t component_count{ 0 };
    std::size_t attempted_component_count{ 0 };
    std::size_t accepted_component_count{ 0 };
    std::size_t atom_count{ 0 };
    std::size_t parameter_count{ 0 };
    std::size_t round_count{ 0 };
    std::size_t suspicious_candidate_atom_count{ 0 };
    std::optional<double> objective_before{};
    std::optional<double> objective_after{};
    double elapsed_milliseconds{ 0.0 };
    struct Component
    {
        std::vector<ClusterKey> key_list{};
        std::size_t atom_count{ 0 };
        std::size_t parameter_count{ 0 };
        std::size_t round_count{ 0 };
        std::size_t suspicious_candidate_atom_count{ 0 };
        std::size_t symbolic_analysis_count{ 0 };
        std::optional<double> objective_before{};
        std::optional<double> objective_after{};
        double elapsed_milliseconds{ 0.0 };
        std::vector<JointCandidateObjectiveDiagnostic> objective_diagnostic_list{};
        bool accepted{ false };
    };
    std::vector<Component> component_list{};
};

struct IterationDiagnostics
{
    std::optional<double> accepted_maximum_transformed_change{};
    double proposal_maximum_transformed_change{ 0.0 };
};

struct IterationObservation
{
    std::map<ClusterKey, ClusterCandidateDiagnostic> candidate_by_key{};
    std::vector<ClusterCandidateDiagnostic> accepted_cluster_diagnostic_list{};
    std::vector<ClusterCandidateDiagnostic> rejected_cluster_diagnostic_list{};
    std::vector<BoundaryComponentReconciliationDiagnostic> boundary_reconciliation_diagnostic_list{};
    IterationDiagnostics diagnostics{};
};

} // namespace rhbm_gem::core::detail
