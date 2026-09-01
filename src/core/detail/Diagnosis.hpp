#pragma once

#include "core/detail/IterationProcess.hpp"

#include <array>
#include <cstddef>
#include <span>
#include <string_view>
#include <vector>

namespace rhbm_gem::core::detail {

struct SecondStageSeedSelectionRecord
{
    SecondStageSeedSource source{ SecondStageSeedSource::GlobalMedian };
    GaussianModel3D original_model{};
    GaussianModel3D selected_model{};
};

struct UnselectedSecondStageSeedSelectionRecord
{
    int atom_serial_id{ 0 };
    SecondStageSeedSource source{ SecondStageSeedSource::GlobalMedian };
    GaussianModel3D selected_model{};
};

using ProgressColumnWidths = std::array<std::size_t, 6>;

bool IsDebugLogLevelEnabled();
void FinishProgressLine(bool quiet_mode);

void LogSecondStageStart(bool quiet_mode);
void LogSecondStageInitializationFailure(
    bool quiet_mode,
    SecondStageInitializationFailure failure);
void LogNoSelectedAtoms(bool quiet_mode);

void LogSecondStageSeedSelections(
    const std::vector<SecondStageSeedSelectionRecord> & selection_record_list,
    bool quiet_mode);
void LogUnselectedSecondStageSeedSelections(
    const std::vector<UnselectedSecondStageSeedSelectionRecord> & selection_record_list,
    bool quiet_mode);

#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
void LogTrustModelShadowDiagnostics(
    bool quiet_mode,
    const IterationResult & iteration_result);
#endif

void LogRejectedClusterDiagnostics(
    bool quiet_mode,
    const std::vector<ClusterCandidateDiagnostic> & diagnostic_list);
void LogAllRejectedResolution(
    bool quiet_mode,
    const IterationResult & iteration_result);
void LogAcceptedCandidateSearchDiagnostics(
    bool quiet_mode,
    const IterationResult & iteration_result);

ProgressColumnWidths BuildProgressColumnWidths(std::size_t atom_count);
void LogProgressHeader(
    bool quiet_mode,
    const ProgressColumnWidths & column_widths);
void LogIterationProgress(
    bool quiet_mode,
    const ProgressColumnWidths & column_widths,
    const IterationResult & iteration_result);

void LogUnrestrictedOperatorAssessments(
    bool quiet_mode,
    std::span<const SuspiciousGaussianAssessment> assessment_by_atom,
    const SuspiciousBlockActivity & block_activity);
void LogConvergenceSafeguardAudit(
    bool quiet_mode,
    const IterationResult & iteration_result,
    const ConvergenceCertificate & certificate,
    const ActiveCoordinatePopulation & active_population,
    const SuspiciousBlockActivity & block_activity,
    std::span<const std::optional<RHBMEstimationStatus>> local_refit_status_by_atom,
    const FitState & assembled_state,
    const FitState & operator_proposal_state,
    const std::vector<std::size_t> & active_index_list,
    bool assembled_uses_polish);

void LogAdaptiveTopologyRebuild(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    const AdaptiveTopologyRebuildDecision & decision,
    const GraphTopology & previous_topology,
    const GraphTopology & rebuilt_topology,
    const CouplingGraphPartition & previous_partition,
    const CouplingGraphPartition & rebuilt_partition,
    bool partition_changed);

void LogFinalDependencyPolish(
    bool quiet_mode,
    const FinalDependencyPolishResult & polish_result,
    FinalPolishCertificationPolicy certification_policy,
    FinalPolishResidualSafetyStatus safety_status,
    bool applied,
    const ConvergenceCertificate * base_certificate = nullptr,
    const ConvergenceCertificate * candidate_certificate = nullptr);

void LogSecondStageAuditTerminal(
    bool quiet_mode,
    const SecondStageContext & context,
    SecondStageStopReason reason,
    std::size_t attempt_number,
    std::size_t accepted_iteration_count,
    const FitState & finalized_state,
    const ObjectiveDomain & comparison_objective_domain);
void LogQuarantineFallback(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    std::size_t entered_target_count,
    std::size_t released_target_count,
    std::size_t failed_probation_count,
    std::size_t unresolved_target_count,
    const FitState & finalized_state);
void LogConverged(
    bool quiet_mode,
    const IterationResult & iteration_result,
    const FitState & finalized_state);
void LogMaximumIterations(
    bool quiet_mode,
    std::size_t entered_target_count,
    std::size_t released_target_count,
    std::size_t failed_probation_count,
    std::size_t unresolved_target_count,
    const BestAuditState & best_audit_state,
    const FitState & latest_state);
void LogSecondStageSummary(
    bool quiet_mode,
    std::size_t accepted_iteration_count,
    const BestAuditState & best_audit_state,
    const PolishProvenance & latest_polish_provenance,
    SecondStageStopReason stop_reason,
    bool final_uses_best_audit);

std::string_view GetFixedPointResidualInterpretationText(
    bool operator_complete,
    bool qualification_passed,
    const TransformedChangeSummary & accepted_change,
    const TransformedChangeSummary & fixed_point_residual);

} // namespace rhbm_gem::core::detail
