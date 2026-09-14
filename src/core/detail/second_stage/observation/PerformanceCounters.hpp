#pragma once

#include "core/detail/second_stage/JointFitting.hpp"

#include <array>
#include <atomic>
#include <chrono>
#include <cstddef>

namespace rhbm_gem::core::detail {

class SecondStageObservationSession;

class PerformanceCounters
{
    friend void LogSecondStagePerformance(const PerformanceCounters &, std::size_t, double);

    const bool m_quiet_mode;
    const SecondStageObservationSession * m_observation;
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
        const BoundaryJointCorrectionWorkspaceMap & boundary_joint_correction_workspace_by_key,
        const SecondStageObservationSession * observation = nullptr);

    ~PerformanceCounters();
    bool AuditEnabled() const noexcept;
    std::array<std::size_t, 6> AuditCounts() const;

    void RecordFullStateMaterialization();
    void RecordGaussianCacheMisses();
    void RecordGaussianCacheHits();
    void RecordObjectiveSampleEvaluation(std::size_t recomputed_sample_count, std::size_t total_sample_count);
    void FinishIterationPhase(std::chrono::steady_clock::time_point start_time);
    void FinishCandidatePhase(std::chrono::steady_clock::time_point start_time);
    void RecordSolverWorkspaceReset();
    void RecordTopologyRebuild(double elapsed_milliseconds);
    void RecordBoundaryReconciliation(
        double elapsed_milliseconds);
    void RecordBoundaryJointCorrection( double elapsed_milliseconds);
    void RecordDependencyPolish(
        double elapsed_milliseconds);

private:
    static double CalculateElapsedMilliseconds(std::chrono::steady_clock::time_point start_time);
    static std::size_t CountRawSamplingEntries(const SecondStageContext & context);
    std::size_t CountCurrentSolverSymbolicAnalyses() const;
};

} // namespace rhbm_gem::core::detail
