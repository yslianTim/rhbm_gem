#pragma once

#include <atomic>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <sstream>

#include <rhbm_gem/utils/domain/Logger.hpp>

#include "core/detail/ReusableWeightedRidgeSolver.hpp"

namespace rhbm_gem::core::detail {

struct LocalFittingPerformanceCounters
{
    const FitOptions & options;
    const LocalFittingClusterSolverWorkspaceMap & solver_workspace_by_key;
    std::chrono::steady_clock::time_point start_time{
        std::chrono::steady_clock::now()
    };
    std::atomic<std::size_t> full_state_materialization_count{ 0 };
    std::atomic<std::size_t> gaussian_cache_hit_count{ 0 };
    std::atomic<std::size_t> gaussian_cache_miss_count{ 0 };
    std::atomic<std::size_t> objective_recomputed_sample_count{ 0 };
    std::atomic<std::size_t> objective_reused_sample_count{ 0 };
    double iteration_phase_milliseconds{ 0.0 };
    double candidate_phase_milliseconds{ 0.0 };

    ~LocalFittingPerformanceCounters()
    {
        if (options.quiet_mode) return;
        std::size_t symbolic_analysis_count{ 0 };
        for (const auto & [key, workspace] : solver_workspace_by_key)
        {
            static_cast<void>(key);
            symbolic_analysis_count +=
                workspace.joint_offset.GetSymbolicAnalysisCount();
            symbolic_analysis_count +=
                workspace.joint_polish.GetSymbolicAnalysisCount();
        }
        const auto total_milliseconds{
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - start_time).count()
        };
        std::ostringstream message;
        message
            << "Second-stage local fitting performance: full_state_materializations="
            << full_state_materialization_count.load()
            << ", gaussian_cache_hit/miss="
            << gaussian_cache_hit_count.load() << "/"
            << gaussian_cache_miss_count.load()
            << ", objective_recomputed/reused_samples="
            << objective_recomputed_sample_count.load() << "/"
            << objective_reused_sample_count.load()
            << ", solver_symbolic_analyses=" << symbolic_analysis_count
            << ", iteration/candidate/total_ms="
            << std::fixed << std::setprecision(3)
            << iteration_phase_milliseconds << "/"
            << candidate_phase_milliseconds << "/"
            << total_milliseconds << ".";
        Logger::Log(LogLevel::Info, message.str());
    }
};

} // namespace rhbm_gem::core::detail
