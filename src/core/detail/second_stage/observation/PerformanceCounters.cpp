#include "core/detail/second_stage/observation/PerformanceCounters.hpp"
#include "core/detail/second_stage/observation/SecondStageLogging.hpp"
#include "core/detail/second_stage/observation/SecondStageObservation.hpp"

namespace rhbm_gem::core::detail {

PerformanceCounters::PerformanceCounters(
    bool quiet_mode,
    const SecondStageContext & context,
    const ClusterSolverWorkspaceMap & solver_workspace_by_key,
    const BoundaryJointCorrectionWorkspaceMap &
        boundary_joint_correction_workspace_by_key,
    const SecondStageObservationSession * observation)
    : m_quiet_mode{ quiet_mode },
      m_observation{ observation },
      m_solver_workspace_by_key{ solver_workspace_by_key },
      m_boundary_joint_correction_workspace_by_key{ boundary_joint_correction_workspace_by_key },
      m_start_time{ std::chrono::steady_clock::now() },
      m_cached_sample_count{ AuditEnabled() ? CountRawSamplingEntries(context) : 0 }
{
}

bool PerformanceCounters::AuditEnabled() const noexcept
{
    return m_observation && m_observation->Enabled();
}

PerformanceCounters::~PerformanceCounters()
{
    if (m_quiet_mode) return;

    const auto total_milliseconds{
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - m_start_time).count()
    };
    LogSecondStagePerformance(*this, 0, total_milliseconds);
}

void PerformanceCounters::RecordFullStateMaterialization()
{
    if (!AuditEnabled()) return;
    m_full_state_materialization_count.fetch_add(1, std::memory_order_relaxed);
}

void PerformanceCounters::RecordGaussianCacheMisses()
{
    if (!AuditEnabled()) return;
    m_gaussian_cache_miss_count.fetch_add(m_cached_sample_count, std::memory_order_relaxed);
}

void PerformanceCounters::RecordGaussianCacheHits()
{
    if (!AuditEnabled()) return;
    m_gaussian_cache_hit_count.fetch_add(m_cached_sample_count, std::memory_order_relaxed);
}

void PerformanceCounters::RecordObjectiveSampleEvaluation(
    std::size_t recomputed_sample_count,
    std::size_t total_sample_count)
{
    if (!AuditEnabled()) return;
    m_objective_recomputed_sample_count.fetch_add(
        recomputed_sample_count,
        std::memory_order_relaxed);
    m_objective_reused_sample_count.fetch_add(
        total_sample_count > recomputed_sample_count ?
            total_sample_count - recomputed_sample_count : 0,
        std::memory_order_relaxed);
}

void PerformanceCounters::FinishIterationPhase(std::chrono::steady_clock::time_point start_time)
{
    m_iteration_phase_milliseconds += CalculateElapsedMilliseconds(start_time);
}

void PerformanceCounters::FinishCandidatePhase(std::chrono::steady_clock::time_point start_time)
{
    m_candidate_phase_milliseconds += CalculateElapsedMilliseconds(start_time);
}

void PerformanceCounters::RecordSolverWorkspaceReset()
{
    if (!AuditEnabled()) return;
    m_retired_solver_symbolic_analysis_count += CountCurrentSolverSymbolicAnalyses();
}

void PerformanceCounters::RecordTopologyRebuild(double elapsed_milliseconds)
{
    m_topology_rebuild_milliseconds += elapsed_milliseconds;
}

void PerformanceCounters::RecordBoundaryReconciliation(double elapsed_milliseconds)
{
    m_boundary_reconciliation_milliseconds += elapsed_milliseconds;
}

void PerformanceCounters::RecordBoundaryJointCorrection(double elapsed_milliseconds)
{
    m_boundary_joint_correction_milliseconds += elapsed_milliseconds;
}





void PerformanceCounters::RecordDependencyPolish(double elapsed_milliseconds)
{
    m_dependency_polish_milliseconds += elapsed_milliseconds;
}

double PerformanceCounters::CalculateElapsedMilliseconds(std::chrono::steady_clock::time_point start_time)
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start_time).count();
}

std::size_t PerformanceCounters::CountRawSamplingEntries(const SecondStageContext & context)
{
    std::size_t count{ 0 };
    for (const auto & atom_context : context.atom_list)
    {
        count += atom_context.raw_sampling_entries.size();
    }
    return count;
}

std::size_t PerformanceCounters::CountCurrentSolverSymbolicAnalyses() const
{
    std::size_t count{ 0 };
    for (const auto & [key, workspace] : m_solver_workspace_by_key)
    {
        static_cast<void>(key);
        count += workspace.joint_offset.GetSymbolicAnalysisCount();
        count += workspace.joint_polish.GetSymbolicAnalysisCount();
    }
    for (const auto & [key, workspace] :
        m_boundary_joint_correction_workspace_by_key)
    {
        static_cast<void>(key);
        count += workspace.GetSymbolicAnalysisCount();
    }
    return count;
}

std::array<std::size_t, 6> PerformanceCounters::AuditCounts() const
{
    if (!AuditEnabled()) return {};
    return {m_full_state_materialization_count.load(), m_gaussian_cache_hit_count.load(),
        m_gaussian_cache_miss_count.load(), m_objective_recomputed_sample_count.load(),
        m_objective_reused_sample_count.load(), m_retired_solver_symbolic_analysis_count + CountCurrentSolverSymbolicAnalyses()};
}
} // namespace rhbm_gem::core::detail
