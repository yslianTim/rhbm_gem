#pragma once

#include "core/detail/second_stage/JointFitting.hpp"
#include "core/detail/second_stage/SuspiciousUpdate.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace rhbm_gem::core::detail {

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

struct ConvergenceCertificate
{
    TransformedChange accepted_active_p99{};
    TransformedChange operator_nominal_p99{};
    bool solver_qualified{ true };
    bool operator_complete{ true };
    bool objective_domain_changed{ false };
    bool quarantine_transition{ false };
    bool suspicious_block_fallback{ false };
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

enum class FinalPolishResidualSafetyStatus
{
    NotEvaluated,
    AbsolutePassed,
    Failed,
    Error
};

bool AreActiveCoordinatesSolverQualified(
    const std::vector<std::size_t> & atom_index_list,
    const std::vector<ClusterKey> & cluster_key_list,
    const SuspiciousBlockActivity & block_activity,
    std::span<const std::optional<RHBMEstimationStatus>> local_refit_status_by_atom,
    const ClusterHealthMap & health_by_key);

} // namespace rhbm_gem::core::detail
