#pragma once

#include "core/detail/CandidateSelection.hpp"
#include "core/detail/GaussianModelOperations.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace rhbm_gem {
class ModelObject;
}

namespace rhbm_gem::core {
struct FitOptions;
}

namespace rhbm_gem::core::detail {

enum class SecondStageSeedSource
{
    LocalMdpde,
    GlobalMedian
};

struct SecondStageSeedSelection
{
    SecondStageSeedSource source{ SecondStageSeedSource::GlobalMedian };
    GaussianModel3DWithUncertainty model{};
};

struct SecondStageSeedSelectionRecord
{
    SecondStageSeedSource source{ SecondStageSeedSource::GlobalMedian };
    GaussianModel3D original_model{};
    GaussianModel3D selected_model{};
};

std::optional<SecondStageSeedSelection> SelectSecondStageSeed(
    const GaussianModel3DWithUncertainty & local_mdpde,
    const std::optional<GaussianModel3D> & global_median);

enum class SecondStageStopReason
{
    None,
    Quarantine,
    Converged,
    AuditPatience,
    AllRejectedBacktrackingExhausted,
    AllRejectedAtMaximumIterations,
    MaximumIterations
};

constexpr std::size_t kMaximumIterations{ 100 };

struct IterationDiagnostics
{
    std::optional<double> accepted_maximum_transformed_change{};
    double proposal_maximum_transformed_change{ 0.0 };
};

struct IterationResult
{
    std::vector<ClusterCandidateDiagnostic> accepted_cluster_diagnostic_list{};
    std::vector<ClusterCandidateDiagnostic> rejected_cluster_diagnostic_list{};
    std::vector<BoundaryComponentReconciliationDiagnostic>
        boundary_reconciliation_diagnostic_list{};
    TrustRegionRadiusUpdate trust_region_update{};
    std::size_t attempt_number{ 0 };
    std::size_t accepted_iteration_count{ 0 };
    std::size_t active_atom_count{ 0 };
    std::size_t quarantine_atom_count{ 0 };
    PolishProgress polish_progress{};
    std::size_t suspicious_atom_count{ 0 };
    IterationDiagnostics diagnostics{};
    SecondStageStopReason stop_reason{ SecondStageStopReason::None };
    bool objective_domain_changed{ false };
    TransformedChange transformed_change_percentile{};
};

constexpr double kAdaptiveTopologyRebuildDriftThreshold{ 0.10 };

double CalculateAdaptiveTopologyDrift(
    const FitState & accepted_state,
    const FittedGaussianSnapshot & topology_reference_state,
    const std::vector<std::size_t> & active_index_list);

struct ActiveCoordinatePopulation
{
    std::vector<std::size_t> active_shape_atom_index_list{};
    std::vector<std::size_t> active_offset_atom_index_list{};
};

ActiveCoordinatePopulation BuildActiveCoordinatePopulation(
    const std::vector<std::size_t> & atom_index_list,
    const SuspiciousBlockActivity & block_activity);

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

struct ConvergenceCertificate
{
    TransformedChange accepted_active_p99{};
    TransformedChange operator_nominal_p99{};
    bool solver_qualified{ true };
    bool operator_complete{ true };
    bool objective_domain_changed{ false };
    bool quarantine_transition{ false };
    bool suspicious_offset_fallback{ false };
    bool rejected_cluster{ false };

    bool StrictOperatorPassed() const;
    bool ProductionConverged() const;
};

struct ConvergenceDiagnostics
{
    TransformedChangeSummary accepted_active_movement{};
    TransformedChangeSummary operator_nominal_residual{};
};

struct ConvergenceAssessment
{
    ConvergenceCertificate certificate{};
    ConvergenceDiagnostics diagnostics{};
};

enum class FinalPolishCertificationPolicy
{
    RequireResidualNonRegression,
    RequireStrictFixedPoint
};

enum class FinalPolishResidualSafetyStatus
{
    NotEvaluated,
    AbsolutePassed,
    RelativePassed,
    Failed,
    Error
};

bool AreActiveCoordinatesSolverQualified(
    const std::vector<std::size_t> & atom_index_list,
    const std::vector<ClusterKey> & cluster_key_list,
    const SuspiciousBlockActivity & block_activity,
    std::span<const std::optional<RHBMEstimationStatus>> local_refit_status_by_atom,
    const ClusterHealthMap & health_by_key);

void RunSecondStageIterations(ModelObject & model_object, const FitOptions & options);

} // namespace rhbm_gem::core::detail
