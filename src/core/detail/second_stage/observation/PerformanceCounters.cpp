#include "core/detail/second_stage/observation/PerformanceCounters.hpp"
#include "core/detail/second_stage/observation/SecondStageLogging.hpp"

namespace rhbm_gem::core::detail {

PerformanceCounters::PerformanceCounters(
    bool quiet_mode,
    const SecondStageContext & context,
    const ClusterSolverWorkspaceMap & solver_workspace_by_key,
    const BoundaryJointCorrectionWorkspaceMap &
        boundary_joint_correction_workspace_by_key)
    : m_quiet_mode{ quiet_mode },
      m_solver_workspace_by_key{ solver_workspace_by_key },
      m_boundary_joint_correction_workspace_by_key{ boundary_joint_correction_workspace_by_key },
      m_start_time{ std::chrono::steady_clock::now() },
      m_cached_sample_count{ CountRawSamplingEntries(context) }
{
}

PerformanceCounters::~PerformanceCounters()
{
    if (m_quiet_mode) return;

    const auto symbolic_analysis_count{
        m_retired_solver_symbolic_analysis_count + CountCurrentSolverSymbolicAnalyses()
    };
    const auto total_milliseconds{
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - m_start_time).count()
    };
    LogSecondStagePerformance(*this, symbolic_analysis_count, total_milliseconds);
}

void PerformanceCounters::RecordFullStateMaterialization()
{
    m_full_state_materialization_count.fetch_add(1, std::memory_order_relaxed);
}

void PerformanceCounters::RecordGaussianCacheMisses()
{
    m_gaussian_cache_miss_count.fetch_add(m_cached_sample_count, std::memory_order_relaxed);
}

void PerformanceCounters::RecordGaussianCacheHits()
{
    m_gaussian_cache_hit_count.fetch_add(m_cached_sample_count, std::memory_order_relaxed);
}

void PerformanceCounters::RecordObjectiveSampleEvaluation(
    std::size_t recomputed_sample_count,
    std::size_t total_sample_count)
{
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
    m_retired_solver_symbolic_analysis_count += CountCurrentSolverSymbolicAnalyses();
}

void PerformanceCounters::RecordTopologyRebuild(double elapsed_milliseconds, bool partition_changed)
{
    m_topology_rebuild_attempt_count++;
    if (partition_changed) m_topology_partition_change_count++;
    m_topology_rebuild_milliseconds += elapsed_milliseconds;
}

void PerformanceCounters::RecordBoundaryReconciliation(
    std::size_t attempt_count,
    std::size_t backtracked_count,
    std::size_t rejected_count,
    double elapsed_milliseconds)
{
    m_boundary_reconciliation_attempt_count += attempt_count;
    m_boundary_reconciliation_backtracked_count += backtracked_count;
    m_boundary_reconciliation_rejected_count += rejected_count;
    m_boundary_reconciliation_milliseconds += elapsed_milliseconds;
}

void PerformanceCounters::RecordBoundaryJointCorrection(bool accepted, double elapsed_milliseconds)
{
    m_boundary_joint_correction_attempt_count++;
    if (accepted)
    {
        m_boundary_joint_correction_accepted_count++;
    }
    else
    {
        m_boundary_joint_correction_fallback_count++;
    }
    m_boundary_joint_correction_milliseconds += elapsed_milliseconds;
}

void PerformanceCounters::RecordBoundaryRescue(bool accepted, bool used_fallback)
{
    m_boundary_rescue_attempt_count++;
    if (accepted)
    {
        m_boundary_rescue_accepted_count++;
    }
    else
    {
        m_boundary_rescue_rejected_count++;
    }
    if (used_fallback) m_boundary_rescue_fallback_count++;
}

void PerformanceCounters::RecordBoundaryRescueExclusions(
    std::size_t hard_failure_count,
    std::size_t invalid_proposal_count,
    std::size_t objective_unavailable_count)
{
    m_boundary_rescue_hard_failure_exclusion_count += hard_failure_count;
    m_boundary_rescue_invalid_proposal_exclusion_count += invalid_proposal_count;
    m_boundary_rescue_objective_unavailable_exclusion_count += objective_unavailable_count;
}

void PerformanceCounters::RecordDependencyPolish(
    std::size_t component_count,
    std::size_t attempt_count,
    std::size_t accepted_count,
    std::size_t fallback_count,
    std::size_t atom_count,
    std::size_t parameter_count,
    std::size_t round_count,
    double elapsed_milliseconds)
{
    m_dependency_polish_component_count += component_count;
    m_dependency_polish_attempt_count += attempt_count;
    m_dependency_polish_accepted_count += accepted_count;
    m_dependency_polish_fallback_count += fallback_count;
    m_dependency_polish_atom_count += atom_count;
    m_dependency_polish_parameter_count += parameter_count;
    m_dependency_polish_round_count += round_count;
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

} // namespace rhbm_gem::core::detail
