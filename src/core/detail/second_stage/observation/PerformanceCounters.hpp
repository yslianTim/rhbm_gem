#pragma once

#include "core/detail/second_stage/JointFitting.hpp"

#include <atomic>
#include <chrono>
#include <cstddef>

namespace rhbm_gem::core::detail {

class PerformanceCounters
{
    friend void LogSecondStagePerformance(const PerformanceCounters &, std::size_t, double);

    const bool m_quiet_mode;
    const ClusterSolverWorkspaceMap & m_solver_workspace_by_key;
    const BoundaryJointCorrectionWorkspaceMap & m_boundary_joint_correction_workspace_by_key;
    const std::chrono::steady_clock::time_point m_start_time;
    const std::size_t m_cached_sample_count;
    std::atomic<std::size_t> m_full_state_materialization_count{ 0 };
    std::atomic<std::size_t> m_gaussian_cache_hit_count{ 0 };
    std::atomic<std::size_t> m_gaussian_cache_miss_count{ 0 };
    std::atomic<std::size_t> m_objective_recomputed_sample_count{ 0 };
    std::atomic<std::size_t> m_objective_reused_sample_count{ 0 };
    std::size_t m_retired_solver_symbolic_analysis_count{ 0 };
    std::size_t m_topology_rebuild_attempt_count{ 0 };
    std::size_t m_topology_partition_change_count{ 0 };
    std::size_t m_boundary_reconciliation_attempt_count{ 0 };
    std::size_t m_boundary_reconciliation_backtracked_count{ 0 };
    std::size_t m_boundary_reconciliation_rejected_count{ 0 };
    std::size_t m_boundary_joint_correction_attempt_count{ 0 };
    std::size_t m_boundary_joint_correction_accepted_count{ 0 };
    std::size_t m_boundary_joint_correction_fallback_count{ 0 };
    std::size_t m_boundary_rescue_attempt_count{ 0 };
    std::size_t m_boundary_rescue_accepted_count{ 0 };
    std::size_t m_boundary_rescue_fallback_count{ 0 };
    std::size_t m_boundary_rescue_rejected_count{ 0 };
    std::size_t m_boundary_rescue_hard_failure_exclusion_count{ 0 };
    std::size_t m_boundary_rescue_invalid_proposal_exclusion_count{ 0 };
    std::size_t m_boundary_rescue_objective_unavailable_exclusion_count{ 0 };
    std::size_t m_dependency_polish_component_count{ 0 };
    std::size_t m_dependency_polish_attempt_count{ 0 };
    std::size_t m_dependency_polish_accepted_count{ 0 };
    std::size_t m_dependency_polish_fallback_count{ 0 };
    std::size_t m_dependency_polish_atom_count{ 0 };
    std::size_t m_dependency_polish_parameter_count{ 0 };
    std::size_t m_dependency_polish_round_count{ 0 };
    double m_iteration_phase_milliseconds{ 0.0 };
    double m_candidate_phase_milliseconds{ 0.0 };
    double m_topology_rebuild_milliseconds{ 0.0 };
    double m_boundary_reconciliation_milliseconds{ 0.0 };
    double m_boundary_joint_correction_milliseconds{ 0.0 };
    double m_dependency_polish_milliseconds{ 0.0 };

public:
    PerformanceCounters(
        bool quiet_mode,
        const SecondStageContext & context,
        const ClusterSolverWorkspaceMap & solver_workspace_by_key,
        const BoundaryJointCorrectionWorkspaceMap & boundary_joint_correction_workspace_by_key);

    ~PerformanceCounters();

    void RecordFullStateMaterialization();
    void RecordGaussianCacheMisses();
    void RecordGaussianCacheHits();
    void RecordObjectiveSampleEvaluation(std::size_t recomputed_sample_count, std::size_t total_sample_count);
    void FinishIterationPhase(std::chrono::steady_clock::time_point start_time);
    void FinishCandidatePhase(std::chrono::steady_clock::time_point start_time);
    void RecordSolverWorkspaceReset();
    void RecordTopologyRebuild(double elapsed_milliseconds, bool partition_changed);
    void RecordBoundaryReconciliation(
        std::size_t attempt_count,
        std::size_t backtracked_count,
        std::size_t rejected_count,
        double elapsed_milliseconds);
    void RecordBoundaryJointCorrection(bool accepted, double elapsed_milliseconds);
    void RecordBoundaryRescue(bool accepted, bool used_fallback);
    void RecordBoundaryRescueExclusions(
        std::size_t hard_failure_count,
        std::size_t invalid_proposal_count,
        std::size_t objective_unavailable_count);
    void RecordDependencyPolish(
        std::size_t component_count,
        std::size_t attempt_count,
        std::size_t accepted_count,
        std::size_t fallback_count,
        std::size_t atom_count,
        std::size_t parameter_count,
        std::size_t round_count,
        double elapsed_milliseconds);

private:
    static double CalculateElapsedMilliseconds(std::chrono::steady_clock::time_point start_time);
    static std::size_t CountRawSamplingEntries(const SecondStageContext & context);
    std::size_t CountCurrentSolverSymbolicAnalyses() const;
};

} // namespace rhbm_gem::core::detail
