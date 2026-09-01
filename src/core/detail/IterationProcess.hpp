#pragma once

#include "core/detail/CandidateSelection.hpp"
#include "core/detail/GaussianModelOperations.hpp"

#include <array>
#include <cstddef>
#include <map>
#include <optional>
#include <span>
#include <string_view>
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

enum class FixedPointResidualInterpretation
{
    Restricted,
    UnqualifiedSmall,
    StepLimited,
    PostprocessedMovement,
    BulkFixedPointWithTail,
    FixedPointConverged,
    Progressing
};

FixedPointResidualInterpretation EvaluateFixedPointResidualInterpretation(
    bool operator_complete,
    bool qualification_passed,
    const TransformedChangeSummary & accepted_change,
    const TransformedChangeSummary & fixed_point_residual);

std::string_view GetFixedPointResidualInterpretationText(
    FixedPointResidualInterpretation interpretation);

struct ActiveCoordinatePopulation
{
    TransformedChangeIndexListByParameter active_atom_index_list_by_parameter{};
    std::vector<ClusterKey> active_offset_group_atom_index_list{};
    std::vector<char> mixed_offset_group_mask{};
    std::size_t total_offset_group_count{ 0 };
    std::size_t fixed_offset_group_count{ 0 };
    std::size_t quarantined_offset_group_count{ 0 };
    std::size_t mixed_offset_group_count{ 0 };
};

ActiveCoordinatePopulation BuildActiveCoordinatePopulation(
    const std::vector<std::size_t> & atom_index_list,
    const std::vector<ClusterKey> & cluster_key_list,
    const std::vector<std::size_t> & group_id_by_atom_index,
    const SuspiciousBlockActivity & block_activity,
    const SuspiciousBlockActivity & quarantine_activity);

TransformedChangeSummary SummarizeActiveDofChanges(
    const std::vector<TransformedChange> & change_list,
    const ActiveCoordinatePopulation & population);

TransformedChangeSummary SummarizeActiveDofChanges(
    const FitState & current_state,
    const FitState & previous_state,
    const ActiveCoordinatePopulation & population);

SuspiciousUpdateMask BuildSuspiciousFailureAtomMask(
    const SuspiciousBlockActivity & block_activity,
    std::span<const SuspiciousGaussianAssessment> assessment_by_atom);

struct SolverQualificationAudit
{
    bool production_qualified{ false };
    bool solver_qualified{ true };
    bool restricted_active_set{ false };
    bool all_fixed{ false };
    std::size_t active_shape_count{ 0 };
    std::size_t qualified_shape_count{ 0 };
    std::size_t soft_unqualified_shape_count{ 0 };
    std::size_t hard_failure_shape_count{ 0 };
    std::size_t fixed_shape_count{ 0 };
    std::size_t quarantined_shape_count{ 0 };
    std::size_t active_offset_group_count{ 0 };
    std::size_t qualified_offset_group_count{ 0 };
    std::size_t soft_unqualified_offset_group_count{ 0 };
    std::size_t hard_failure_offset_group_count{ 0 };
    std::size_t fixed_offset_group_count{ 0 };
    std::size_t quarantined_offset_group_count{ 0 };
    std::size_t mixed_offset_group_count{ 0 };
    std::array<std::size_t, 7> joint_offset_status_count{};
};

struct ConvergenceOrthogonalBlockers
{
    bool objective_domain_changed{ false };
    bool quarantine_transition{ false };
    bool suspicious_offset_fallback{ false };
    bool rejected_cluster{ false };

    bool Clear() const;
};

struct ConvergenceCertificate
{
    ActiveCoordinatePopulation accepted_active_population{};
    TransformedChangeSummary accepted_active_movement{};
    ActiveCoordinatePopulation operator_nominal_population{};
    TransformedChangeSummary operator_nominal_residual{};
    SolverQualificationAudit solver_qualification{};
    std::array<
        std::size_t,
        GaussianModel3D::TransformedCoordinateSize()> operator_unavailable_count{};
    std::array<std::size_t, 3> operator_unavailable_reason_count{};
    std::array<
        std::size_t,
        GaussianModel3D::TransformedCoordinateSize()> operator_tail_count{};
    bool operator_shadow_shape_refit_performed{ false };
    ConvergenceOrthogonalBlockers blockers{};

    bool AcceptedPercentilePassed() const;
    bool OperatorPercentilePassed() const;
    bool OperatorComplete() const;
    bool InvariantsClear() const;
    bool MixedSharedGroup() const;
    bool StrictOperatorPassed() const;
    bool ProductionConverged() const;
};

SolverQualificationAudit EvaluateSolverQualificationAudit(
    const std::vector<std::size_t> & atom_index_list,
    const std::vector<ClusterKey> & cluster_key_list,
    const std::vector<std::size_t> & group_id_by_atom_index,
    const SuspiciousBlockActivity & block_activity,
    const SuspiciousBlockActivity & quarantine_activity,
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

struct StabilizationTerminalFailure
{
    StabilizationTerminalReason category{ StabilizationTerminalReason::None };
    std::optional<SuspiciousGaussianReason> guard_reason{};

    friend auto operator<=>(
        const StabilizationTerminalFailure &,
        const StabilizationTerminalFailure &) = default;
};

using QuarantineFailureReason =
    std::variant<JointOffsetSolveStatus, StabilizationTerminalFailure>;

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

bool HasPendingQuarantineLifecycle(
    const QuarantineFailureStateMap & state_by_target);

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
