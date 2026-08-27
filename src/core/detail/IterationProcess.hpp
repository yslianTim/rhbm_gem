#pragma once

#include "core/detail/CandidateSelection.hpp"

#include <array>
#include <cstddef>
#include <map>
#include <optional>
#include <span>
#include <variant>
#include <vector>

namespace rhbm_gem::core::detail {

enum class SecondStageSeedSource
{
    GroupPosterior,
    GroupPrior,
    GroupMedian,
    GlobalMedian
};

struct SecondStageSeedCandidates
{
    std::optional<GaussianModel3DWithUncertainty> group_posterior{};
    std::optional<GaussianModel3DWithUncertainty> group_prior{};
    std::optional<GaussianModel3DWithUncertainty> group_median{};
    std::optional<GaussianModel3DWithUncertainty> global_median{};
};

struct SecondStageSeedSelection
{
    SecondStageSeedSource source{ SecondStageSeedSource::GlobalMedian };
    GaussianModel3DWithUncertainty model{};
};

std::optional<SecondStageSeedSelection> SelectSecondStageSeed(const SecondStageSeedCandidates & candidates);

constexpr double kAdaptiveTopologyRebuildDriftThreshold{ 0.10 };
constexpr std::size_t kAdaptiveTopologyRebuildAcceptedIterationInterval{ 3 };

enum class AdaptiveTopologyRebuildTrigger
{
    None,
    Drift,
    Interval
};

struct AdaptiveTopologyRebuildDecision
{
    AdaptiveTopologyRebuildTrigger trigger{ AdaptiveTopologyRebuildTrigger::None };
    double maximum_transformed_drift{ 0.0 };
};

AdaptiveTopologyRebuildDecision EvaluateAdaptiveTopologyRebuildTrigger(
    const FitState & accepted_state,
    const FitState & topology_reference_state,
    const std::vector<std::size_t> & active_index_list,
    std::size_t accepted_iterations_since_rebuild);

struct ConvergenceChangePredicates
{
    bool percentile_converged{ false };
    bool maximum_converged{ false };

    bool Converged() const
    {
        return percentile_converged && maximum_converged;
    }
};

struct ConvergenceSafeguardPredicates
{
    bool stationarity_eligible{ false };
    ConvergenceChangePredicates accepted{};
    ConvergenceChangePredicates raw{};

    bool Converged() const
    {
        return stationarity_eligible && accepted.Converged() && raw.Converged();
    }
};

ConvergenceSafeguardPredicates EvaluateConvergenceSafeguardPredicates(
    bool stationarity_eligible,
    const TransformedChangeSummary & accepted_change,
    const TransformedChangeSummary & raw_change);

TransformedChangeIndexListByParameter BuildActiveBlockChangeIndexLists(
    const std::vector<std::size_t> & atom_index_list,
    const SuspiciousBlockActivity & block_activity);

struct ActiveCoordinateAuditPopulation
{
    TransformedChangeIndexListByParameter member_index_list_by_parameter{};
    std::vector<ClusterKey> active_offset_group_atom_index_list{};
    std::vector<std::size_t> active_offset_group_size_list{};
    std::vector<char> mixed_offset_group_mask{};
    std::size_t total_offset_group_count{ 0 };
    std::size_t fixed_offset_group_count{ 0 };
    std::size_t quarantined_offset_group_count{ 0 };
    std::size_t mixed_offset_group_count{ 0 };
};

ActiveCoordinateAuditPopulation BuildActiveCoordinateAuditPopulation(
    const std::vector<std::size_t> & atom_index_list,
    const std::vector<ClusterKey> & cluster_key_list,
    const std::vector<std::size_t> & group_id_by_atom_index,
    const SuspiciousBlockActivity & block_activity,
    const SuspiciousBlockActivity & quarantine_activity);

struct ActiveCoordinateChangeAudit
{
    TransformedChangeSummary member{};
    TransformedChangeSummary shared_dof{};
};

ActiveCoordinateChangeAudit EvaluateActiveCoordinateChangeAudit(
    const std::vector<algorithm::ParameterChange> & change_list,
    const ActiveCoordinateAuditPopulation & population);

ActiveCoordinateChangeAudit EvaluateActiveCoordinateChangeAudit(
    const FitState & current_state,
    const FitState & previous_state,
    const ActiveCoordinateAuditPopulation & population);

struct ConvergenceStationarityAudit
{
    bool active_block_eligible{ false };
    bool full_cluster_eligible{ false };
    std::size_t active_block_ineligible_cluster_count{ 0 };
    std::size_t refit_ineligible_cluster_count{ 0 };
    std::size_t soft_joint_nonconverged_cluster_count{ 0 };
    std::size_t hard_joint_failure_cluster_count{ 0 };
    std::array<std::size_t, 7> joint_offset_status_count{};
};

ConvergenceStationarityAudit EvaluateConvergenceStationarityAudit(
    const ClusterHealthMap & health_by_key);

struct StrictConvergenceStationarityAudit
{
    bool current_eligible{ false };
    bool strict_eligible{ true };
    bool restricted_active_set{ false };
    bool all_fixed{ false };
    std::size_t active_shape_count{ 0 };
    std::size_t qualified_shape_count{ 0 };
    std::size_t soft_nonstationary_shape_count{ 0 };
    std::size_t hard_failure_shape_count{ 0 };
    std::size_t fixed_shape_count{ 0 };
    std::size_t quarantined_shape_count{ 0 };
    std::size_t active_offset_group_count{ 0 };
    std::size_t qualified_offset_group_count{ 0 };
    std::size_t soft_nonstationary_offset_group_count{ 0 };
    std::size_t hard_failure_offset_group_count{ 0 };
    std::size_t fixed_offset_group_count{ 0 };
    std::size_t quarantined_offset_group_count{ 0 };
    std::size_t mixed_offset_group_count{ 0 };
};

StrictConvergenceStationarityAudit EvaluateStrictConvergenceStationarityAudit(
    const std::vector<std::size_t> & atom_index_list,
    const std::vector<ClusterKey> & cluster_key_list,
    const std::vector<std::size_t> & group_id_by_atom_index,
    const SuspiciousBlockActivity & block_activity,
    const SuspiciousBlockActivity & quarantine_activity,
    const SuspiciousUpdateMask & shape_stationarity_eligible_atom_mask,
    const SuspiciousUpdateMask & offset_stationarity_eligible_atom_mask,
    std::span<const std::optional<RHBMEstimationStatus>> local_refit_status_by_atom,
    const ClusterHealthMap & health_by_key);

constexpr std::size_t kPersistentQuarantineFailureIterationLimit{ 5 };
constexpr std::size_t kQuarantineProbationCooldown{ 2 };
constexpr std::size_t kQuarantineMaximumProbationCount{ 3 };

enum class QuarantineTargetKind
{
    ShapeAtom,
    OffsetGroup,
    HardFailureCluster
};

struct QuarantineTarget
{
    QuarantineTargetKind kind{ QuarantineTargetKind::ShapeAtom };
    std::vector<std::size_t> atom_index_list{};

    friend auto operator<=>(const QuarantineTarget &, const QuarantineTarget &) = default;
};

using QuarantineFailureReason =
    std::variant<SuspiciousGaussianReason, JointOffsetSolveStatus>;

struct QuarantineFailureObservation
{
    QuarantineTarget target{};
    QuarantineFailureReason reason{};
};

struct QuarantineFailureState
{
    QuarantineFailureReason reason{};
    std::size_t stable_iteration_count{ 0 };
    std::size_t probation_count{ 0 };
    std::size_t next_probation_iteration{ 0 };
    bool quarantined{ false };
    bool probation_active{ false };
    bool probation_exhausted{ false };
};

using QuarantineFailureStateMap = std::map<QuarantineTarget, QuarantineFailureState>;

struct QuarantineStateTransition
{
    std::vector<QuarantineTarget> entered_target_list{};
    std::vector<QuarantineTarget> released_target_list{};
    std::vector<QuarantineTarget> failed_probation_target_list{};
};

QuarantineStateTransition UpdateQuarantineFailureState(
    const std::vector<QuarantineFailureObservation> & observation_list,
    const std::vector<QuarantineTarget> & successful_probation_target_list,
    std::size_t accepted_iteration_count,
    QuarantineFailureStateMap & state_by_target);

} // namespace rhbm_gem::core::detail
