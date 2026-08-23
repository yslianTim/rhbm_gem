#pragma once

#include "core/detail/FittingModel.hpp"
#include "core/detail/CouplingGraph.hpp"
#include "core/detail/JointFitting.hpp"
#include "core/detail/CandidateSelection.hpp"

#include <array>
#include <cstddef>
#include <map>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>
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

std::optional<SecondStageSeedSelection> SelectValidSecondStageSeedCandidate(
    SecondStageSeedSource source,
    const std::optional<GaussianModel3DWithUncertainty> & candidate);

std::optional<SecondStageSeedSelection> SelectSecondStageSeed(const SecondStageSeedCandidates & candidates);

constexpr double kNeighborContributionDistanceMax{ 2.5 };
constexpr double kNeighborAtomSearchRange{ 2.0 * kNeighborContributionDistanceMax };

struct SecondStageSeedSelectionRecord
{
    detail::SecondStageSeedSource source{ detail::SecondStageSeedSource::GlobalMedian };
    GaussianModel3D original_model{};
    GaussianModel3D selected_model{};
};

struct UnselectedSecondStageSeedSelectionRecord
{
    int atom_serial_id{ 0 };
    detail::SecondStageSeedSource source{ detail::SecondStageSeedSource::GlobalMedian };
    GaussianModel3D selected_model{};
};

struct SecondStageInitialStateBuildResult
{
    FitState state{};
    std::vector<SecondStageSeedSelectionRecord> selection_record_list{};
    std::vector<UnselectedSecondStageSeedSelectionRecord> unselected_selection_record_list{};
    enum class Failure
    {
        None,
        SelectedSeedUnavailable,
        UnselectedSeedUnavailable
    } failure{ Failure::None };
};

SecondStageContext BuildSecondStageContext(
    const ModelObject & model_object,
    const FitOptions & options);

void StoreSecondStageNeighborCounts(
    ModelObject & model_object,
    const SecondStageContext & context);


std::optional<GaussianModel3DWithUncertainty> BuildValidGaussianParameterMedian(
    const std::vector<GaussianModel3D> & model_list);

SecondStageInitialStateBuildResult BuildInitialFitState(SecondStageContext & context);

const char * GetSecondStageSeedSourceText(SecondStageSeedSource source);

void LogSecondStageSeedSelections(
    const std::vector<SecondStageSeedSelectionRecord> & selection_record_list,
    bool quiet_mode);

void LogUnselectedSecondStageSeedSelections(
    const std::vector<UnselectedSecondStageSeedSelectionRecord> & selection_record_list,
    bool quiet_mode);


constexpr std::size_t kPersistentTerminalFailureIterationLimit{ 5 };

using PersistentSuspiciousRollbackReason = std::vector<std::size_t>;
using PersistentTerminalFailureReason = std::variant<PersistentSuspiciousRollbackReason, JointOffsetSolveStatus>;

struct PersistentTerminalFailureState
{
    PersistentTerminalFailureReason reason{};
    std::size_t stable_iteration_count{ 0 };
};

using PersistentTerminalFailureStateMap = std::map<ClusterKey, PersistentTerminalFailureState>;
using TerminalPersistentFailureMap = std::map<ClusterKey, PersistentTerminalFailureReason>;

struct TerminalSummary
{
    std::size_t suspicious_cluster_count{ 0 };
    std::size_t suspicious_atom_count{ 0 };
    std::size_t joint_offset_failure_cluster_count{ 0 };
    std::size_t joint_offset_failure_atom_count{ 0 };
    std::map<JointOffsetSolveStatus, std::size_t> joint_offset_failure_status_count{};

    std::size_t AtomCount() const
    {
        return suspicious_atom_count + joint_offset_failure_atom_count;
    }

    bool HasFailures() const
    {
        return AtomCount() > 0;
    }
};

void AppendTerminalSummary(std::ostream & stream, const TerminalSummary & summary);

std::vector<ClusterKey> AccumulateTerminalFailureSummary(
    const TerminalPersistentFailureMap & terminal_failure_by_key,
    TerminalSummary & terminal_summary);

TerminalPersistentFailureMap UpdatePersistentTerminalFailureState(
    const std::vector<ClusterKey> & accepted_key_list,
    const SuspiciousUpdateMask & suspicious_atom_mask,
    const ClusterHealthMap & health_by_key,
    const FitState & assembled_state,
    const FitState & previous_state,
    PersistentTerminalFailureStateMap & state_by_key);

void ApplyTerminalFallbackClusters(
    const std::vector<ClusterKey> & terminal_key_list,
    const FitState & previous_state,
    const PolishProvenance & previous_polish_provenance,
    std::vector<char> & terminal_atom_mask,
    FitState & assembled_state,
    PolishProvenance & assembled_polish_provenance);

struct TerminalFailureState
{
    PersistentTerminalFailureStateMap persistent_state_by_key{};
    std::vector<char> terminal_atom_mask{};
    TerminalSummary terminal_summary{};
    TerminalFailureState() = default;

    explicit TerminalFailureState(std::size_t atom_count) : terminal_atom_mask(atom_count, 0)
    {
    }

    bool HasFailures() const { return terminal_summary.HasFailures(); }
    std::size_t AtomCount() const { return terminal_summary.AtomCount(); }

    std::vector<std::size_t> BuildEligibleActiveIndexList() const;

    bool IsolatePersistentFailures(
        const std::vector<ClusterKey> & accepted_key_list,
        SuspiciousUpdateMask & suspicious_atom_mask,
        const ClusterHealthMap & health_by_key,
        FitState & assembled_state,
        const FitState & previous_state,
        const PolishProvenance & previous_polish_provenance,
        PolishProvenance & assembled_polish_provenance);
};


constexpr std::size_t kMaximumIterations{ 100 };
constexpr std::size_t kAuditPatience{ 3 };

struct IterationDiagnostics
{
    std::vector<ClusterCandidateDiagnostic> accepted_cluster_diagnostic_list{};
    std::vector<ClusterCandidateDiagnostic> rejected_cluster_diagnostic_list{};
    std::size_t combined_backtracking_trial_count{ 0 };
    std::optional<double> combined_backtracking_factor{};
    bool combined_backtracking_exhausted{ false };
    TrustRegionIterationUpdate trust_region_update{};
};

struct IterationState
{
    FitState previous_state{};
    PolishProvenance previous_polish_provenance{};
    SuspiciousUpdateMask rollback_atom_mask{};
    std::vector<std::size_t> active_index_list{};
    CouplingGraphPartition graph_partition{};
    ClusterSolverWorkspaceMap solver_workspace_by_key{};
    ObjectiveDomain objective_domain{};
    BestAuditState best_audit_state{};
    TerminalFailureState terminal_failure_state{};
    ClusterObjectiveStateMap cluster_objective_state{};
    TrustRegionStateSet trust_region_state{};
    std::vector<ClusterKey> unchanged_state_exhausted_key_list{};
    std::size_t accepted_iteration_count{ 0 };
    std::size_t audit_patience_count{ 0 };
};

struct IterationProgress
{
    std::size_t attempt_number{ 0 };
    std::size_t accepted_iteration_count{ 0 };
    std::size_t active_atom_count{ 0 };
    std::size_t terminal_atom_count{ 0 };
    std::size_t accepted_cluster_count{ 0 };
    std::size_t rejected_cluster_count{ 0 };
    PolishProgress polish_progress{};
    std::size_t suspicious_atom_count{ 0 };
    std::optional<double> accepted_maximum_transformed_change{};
    double raw_maximum_transformed_change{ 0.0 };
};

using ProgressColumnWidths = std::array<std::size_t, 6>;

struct IterationResult
{
    IterationDiagnostics diagnostics{};
    IterationProgress progress{};
    std::optional<AllRejectedResolution> all_rejected_resolution{};
    bool objective_domain_changed{ false };
    bool converged{ false };
    bool audit_patience_exhausted{ false };
    algorithm::ParameterChangeStats transformed_change_stats{};
};

void AppendObjectiveBreakdown(
    std::ostringstream & stream,
    const std::optional<ObjectiveBreakdown> & breakdown);

std::string_view GetPreObjectiveFailureReasonText(PreObjectiveFailureReason reason);

void LogRejectedClusterDiagnostics(
    bool quiet_mode,
    const std::vector<ClusterCandidateDiagnostic> & diagnostic_list);

std::string_view GetAllRejectedResolutionText(AllRejectedResolution resolution);

void LogAllRejectedResolution(
    bool quiet_mode,
    const TrustRegionIterationUpdate & trust_region_update,
    AllRejectedResolution resolution);

void LogAcceptedBacktrackingDiagnostics(
    bool quiet_mode,
    const IterationDiagnostics & diagnostics);

std::string FormatProgressMaximum(double value);

ProgressColumnWidths BuildProgressColumnWidths(std::size_t atom_size);

void LogProgressHeader(bool quiet_mode, const ProgressColumnWidths & column_widths);

void LogIterationProgress(
    bool quiet_mode,
    const ProgressColumnWidths & column_widths,
    const IterationProgress & progress);

struct RawIterationResult
{
    FitState state{};
    SuspiciousUpdateMask rollback_atom_mask{};
    ClusterHealthMap health_by_key{};
};

struct LocalAtomRefitResult
{
    LocalGaussianResult result{};
    bool is_stationarity_eligible{ false };
};

std::optional<LocalAtomRefitResult> FitAtomWithJointOffsetFallback(
    const AtomContext & atom_context,
    const LocalGaussianResult & previous_result,
    const GaussianModel3D & offset_model,
    const std::vector<double> & adjusted_response_list,
    const FitOptions & options);

RawIterationResult RunRawIteration(
    const SecondStageContext & context,
    const std::vector<ClusterKey> & cluster_key_list,
    const FitState & previous_state,
    const FitOptions & options,
    const std::vector<double> & ridge_multiplier_list,
    ClusterSolverWorkspaceMap & solver_workspace_by_key);

IterationState BuildIterationState(
    const SecondStageContext & context,
    const GraphTopology & graph_topology,
    FitState initial_state,
    const FitOptions & options);

IterationResult RunIteration(
    const SecondStageContext & context,
    const GraphTopology & graph_topology,
    const FitOptions & options,
    std::size_t attempt_number,
    IterationState & iteration_state,
    PerformanceCounters & performance_counters);

bool RunSecondStageIterationProcess(
    ModelObject & model_object,
    const FitOptions & options);

} // namespace rhbm_gem::core::detail
