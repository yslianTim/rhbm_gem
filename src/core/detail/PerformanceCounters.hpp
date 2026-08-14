#pragma once

#include <atomic>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <sstream>

#include <rhbm_gem/utils/domain/Logger.hpp>

#include "core/detail/ReusableWeightedRidgeSolver.hpp"
#include "core/detail/SecondStageContext.hpp"

namespace rhbm_gem::core::detail {

class PerformanceCounters
{
    const FitOptions & m_options;
    const ClusterSolverWorkspaceMap & m_solver_workspace_by_key;
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
    
public:
    PerformanceCounters(
        const FitOptions & options,
        const SecondStageContext & context,
        const ClusterSolverWorkspaceMap & solver_workspace_by_key)
        : m_options{ options },
          m_solver_workspace_by_key{ solver_workspace_by_key },
          m_start_time{ std::chrono::steady_clock::now() },
          m_cached_sample_count{ CountRawSamplingEntries(context) }
    {
    }

    ~PerformanceCounters()
    {
        if (m_options.quiet_mode) return;

        const auto symbolic_analysis_count{
            m_retired_solver_symbolic_analysis_count + CountCurrentSolverSymbolicAnalyses()
        };
        const auto total_milliseconds{
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - m_start_time).count()
        };
        std::ostringstream message;
        message
            << "Second-stage local fitting performance: full_state_materializations="
            << m_full_state_materialization_count.load()
            << ", gaussian_cache_hit/miss="
            << m_gaussian_cache_hit_count.load() << "/"
            << m_gaussian_cache_miss_count.load()
            << ", objective_recomputed/reused_samples="
            << m_objective_recomputed_sample_count.load() << "/"
            << m_objective_reused_sample_count.load()
            << ", solver_symbolic_analyses=" << symbolic_analysis_count
            << ", iteration/candidate/total_ms="
            << std::fixed << std::setprecision(3)
            << m_iteration_phase_milliseconds << "/"
            << m_candidate_phase_milliseconds << "/"
            << total_milliseconds << ".";
        Logger::Log(LogLevel::Info, message.str());
    }

    void RecordFullStateMaterialization()
    {
        m_full_state_materialization_count.fetch_add(1, std::memory_order_relaxed);
    }

    void RecordGaussianCacheMisses()
    {
        m_gaussian_cache_miss_count.fetch_add(m_cached_sample_count, std::memory_order_relaxed);
    }

    void RecordGaussianCacheHits()
    {
        m_gaussian_cache_hit_count.fetch_add(m_cached_sample_count, std::memory_order_relaxed);
    }

    void RecordObjectiveSampleEvaluation(std::size_t recomputed_sample_count, std::size_t total_sample_count)
    {
        m_objective_recomputed_sample_count.fetch_add(recomputed_sample_count, std::memory_order_relaxed);
        m_objective_reused_sample_count.fetch_add(
            total_sample_count > recomputed_sample_count ?
                total_sample_count - recomputed_sample_count : 0,
            std::memory_order_relaxed);
    }

    [[nodiscard]] std::chrono::steady_clock::time_point StartIterationPhase() const
    {
        return std::chrono::steady_clock::now();
    }

    void FinishIterationPhase(std::chrono::steady_clock::time_point start_time)
    {
        m_iteration_phase_milliseconds += CalculateElapsedMilliseconds(start_time);
    }

    [[nodiscard]] std::chrono::steady_clock::time_point StartCandidatePhase() const
    {
        return std::chrono::steady_clock::now();
    }

    void FinishCandidatePhase(std::chrono::steady_clock::time_point start_time)
    {
        m_candidate_phase_milliseconds += CalculateElapsedMilliseconds(start_time);
    }

    void RecordSolverWorkspaceReset()
    {
        m_retired_solver_symbolic_analysis_count += CountCurrentSolverSymbolicAnalyses();
    }

private:
    static double CalculateElapsedMilliseconds(std::chrono::steady_clock::time_point start_time)
    {
        return std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - start_time).count();
    }

    static std::size_t CountRawSamplingEntries(const SecondStageContext & context)
    {
        std::size_t count{ 0 };
        for (const auto & atom_context : context)
        {
            count += atom_context.raw_sampling_entries.size();
        }
        return count;
    }

    std::size_t CountCurrentSolverSymbolicAnalyses() const
    {
        std::size_t count{ 0 };
        for (const auto & [key, workspace] : m_solver_workspace_by_key)
        {
            static_cast<void>(key);
            count += workspace.joint_offset.GetSymbolicAnalysisCount();
            count += workspace.joint_polish.GetSymbolicAnalysisCount();
        }
        return count;
    }
};

} // namespace rhbm_gem::core::detail
