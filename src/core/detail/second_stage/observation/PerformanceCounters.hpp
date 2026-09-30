#pragma once

#include <chrono>

namespace rhbm_gem::core::detail {

class PerformanceCounters
{
    friend void LogSecondStagePerformance(const PerformanceCounters &, double);

    const bool m_quiet_mode;
    const std::chrono::steady_clock::time_point m_start_time;
    double m_iteration_phase_milliseconds{ 0.0 };
    double m_candidate_phase_milliseconds{ 0.0 };
    double m_topology_rebuild_milliseconds{ 0.0 };
    double m_boundary_reconciliation_milliseconds{ 0.0 };
    double m_boundary_joint_correction_milliseconds{ 0.0 };
    double m_dependency_polish_milliseconds{ 0.0 };

public:
    explicit PerformanceCounters(bool quiet_mode);

    ~PerformanceCounters();
    void FinishIterationPhase(std::chrono::steady_clock::time_point start_time);
    void FinishCandidatePhase(std::chrono::steady_clock::time_point start_time);
    void RecordTopologyRebuild(double elapsed_milliseconds);
    void RecordBoundaryReconciliation(double elapsed_milliseconds);
    void RecordBoundaryJointCorrection(double elapsed_milliseconds);
    void RecordDependencyPolish(double elapsed_milliseconds);

private:
    static double CalculateElapsedMilliseconds(std::chrono::steady_clock::time_point start_time);
};

} // namespace rhbm_gem::core::detail
