#pragma once

#include <cstddef>
#include <exception>
#include <optional>
#include <utility>
#include <vector>

#include <rhbm_gem/utils/domain/Logger.hpp>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "core/detail/JointOffset.hpp"
#include "core/detail/ClusterHealth.hpp"
#include "core/detail/ResidualEvaluation.hpp"
#include "core/detail/FitStateView.hpp"
#include "core/detail/SuspiciousUpdate.hpp"
#include "core/detail/ReusableWeightedRidgeSolver.hpp"
#include "core/detail/ScopedEigenThreadCount.hpp"
#include "core/detail/SecondStageContext.hpp"

namespace rhbm_gem::core::detail {

struct LocalFittingIterationResult
{
    FitState state{};
    SuspiciousUpdateMask rollback_atom_mask{};
    ClusterHealthMap health_by_key{};
};

struct LocalAtomRefitResult
{
    LocalGaussianResult result{};
    bool is_stationarity_eligible{ false };
};

inline std::optional<LocalAtomRefitResult> FitAtomWithJointOffsetFallback(
    const SecondStageContext & context,
    std::size_t atom_index,
    const LocalGaussianResult & previous_result,
    const GaussianModel3D & offset_model,
    const std::vector<double> & adjusted_response_list,
    const FitOptions & options)
{
    auto adjusted_sampling_entries{
        BuildSecondStageAdjustedSamples(context, atom_index, adjusted_response_list)
    };
    const auto & previous_model{ previous_result.mdpde.GetModel() };
    const auto previous_baseline{
        local_fitting_suspicious_internal::BuildPreviousSuspiciousProfileBaseline(
            adjusted_sampling_entries,
            previous_model,
            options)
    };
    const auto is_post_refit_candidate_acceptable = [&](const GaussianModel3D & model)
    {
        const auto reason{ local_fitting_suspicious_internal::EvaluateSuspiciousGaussianUpdate(
                adjusted_sampling_entries,
                previous_model,
                model,
                options,
                previous_baseline,
                true) };
        return reason == SuspiciousGaussianReason::None;
    };
    const auto is_offset_only_fallback_acceptable =
        [&](const GaussianModel3D & model)
    {
        const auto reason{ local_fitting_suspicious_internal::EvaluateSuspiciousGaussianUpdate(
                adjusted_sampling_entries,
                previous_model,
                model,
                options,
                previous_baseline,
                false) };
        return reason == SuspiciousGaussianReason::None;
    };
    try
    {
        auto candidate_result{
            EstimateLocalGaussianPrepared(
                context.at(atom_index).refit_design_template,
                adjusted_response_list,
                context.at(atom_index).alpha_r,
                options,
                offset_model)
        };
        const auto candidate_model{ candidate_result.mdpde.GetModel() };
        if (is_post_refit_candidate_acceptable(candidate_model))
        {
            const auto is_stationarity_eligible{
                candidate_result.fit_result.has_value() &&
                IsLocalRefitStatusStationarityEligible(candidate_result.fit_result->status)
            };
            return LocalAtomRefitResult{
                std::move(candidate_result),
                is_stationarity_eligible
            };
        }
    }
    catch (const std::exception &)
    {
    }

    auto result{ previous_result };
    result.ols = GaussianModel3DWithUncertainty{
        result.ols.GetModel().WithOffset(offset_model.GetOffset()),
        result.ols.GetStandardDeviationModel()
    };
    result.mdpde = GaussianModel3DWithUncertainty{
        result.mdpde.GetModel().WithOffset(offset_model.GetOffset()),
        result.mdpde.GetStandardDeviationModel()
    };
    if (!is_offset_only_fallback_acceptable(result.mdpde.GetModel()))
    {
        return std::nullopt;
    }
    return LocalAtomRefitResult{ std::move(result), false };
}

inline LocalFittingIterationResult RunLocalFittingIteration(
    const SecondStageContext & context,
    const std::vector<ClusterKey> & cluster_key_list,
    const FitState & previous_state,
    const FitOptions & options,
    const std::vector<double> & ridge_multiplier_list,
    ClusterSolverWorkspaceMap & solver_workspace_by_key)
{
    const auto selected_atom_size{ context.size() };
    auto current_model_snapshot{
        BuildSecondStageModelSnapshot(
            context,
            BuildFittedGaussianSnapshot(previous_state))
    };
    const auto is_debug_logging_enabled{
        Logger::GetLogLevel() >= LogLevel::Debug
    };
    const auto log_debug_diagnostics{
        !options.quiet_mode && is_debug_logging_enabled
    };
    std::vector<JointOffsetSolveResult> joint_offset_result_list(cluster_key_list.size());
    std::vector<std::exception_ptr> joint_offset_exception_list(cluster_key_list.size());
    const auto solve_joint_offset = [&](std::size_t cluster_position)
    {
        try
        {
            joint_offset_result_list.at(cluster_position) =
                EstimateJointOffsets(
                context,
                cluster_key_list.at(cluster_position),
                current_model_snapshot,
                ridge_multiplier_list,
                solver_workspace_by_key.at(
                    cluster_key_list.at(cluster_position)).joint_offset,
                log_debug_diagnostics);
        }
        catch (...)
        {
            joint_offset_exception_list.at(cluster_position) = std::current_exception();
        }
    };
#ifdef USE_OPENMP
    const bool parallel_joint_offsets{
        !is_debug_logging_enabled &&
        options.thread_size > 1 &&
        cluster_key_list.size() > 1
    };
    if (parallel_joint_offsets)
    {
        ScopedEigenThreadCount eigen_thread_guard{ 1 };
#pragma omp parallel for schedule(dynamic) num_threads(options.thread_size)
        for (std::size_t cluster_position = 0; cluster_position < cluster_key_list.size(); cluster_position++)
        {
            solve_joint_offset(cluster_position);
        }
    }
    else
#endif
    {
        for (std::size_t cluster_position = 0; cluster_position < cluster_key_list.size(); cluster_position++)
        {
            solve_joint_offset(cluster_position);
        }
    }
    for (const auto & exception : joint_offset_exception_list)
    {
        if (exception) std::rethrow_exception(exception);
    }

    ClusterHealthMap health_by_key;
    for (std::size_t cluster_position = 0; cluster_position < cluster_key_list.size(); cluster_position++)
    {
        const auto & key{ cluster_key_list.at(cluster_position) };
        const auto & result{ joint_offset_result_list.at(cluster_position) };
        for (std::size_t i = 0; i < key.size(); i++)
        {
            const auto atom_index{ key.at(i) };
            current_model_snapshot.selected.at(atom_index) =
                current_model_snapshot.selected.at(atom_index)
                    .WithOffset(result.offset(static_cast<Eigen::Index>(i)));
        }
        health_by_key.emplace(key, ClusterHealth{ result.status });
    }

    auto iteration_state{ previous_state };
    SuspiciousUpdateMask rollback_atom_mask(selected_atom_size, 0);
    std::vector<GroupKey> group_key_by_atom_index;
    group_key_by_atom_index.reserve(context.size());
    for (const auto & atom_context : context)
    {
        group_key_by_atom_index.emplace_back(atom_context.group_key);
    }
    for (std::size_t cluster_position = 0; cluster_position < cluster_key_list.size(); cluster_position++)
    {
        const auto & key{ cluster_key_list.at(cluster_position) };
        std::vector<GroupKey> group_key_by_position;
        SuspiciousUpdateMask suspicious_seed_mask(key.size(), 0);
        group_key_by_position.reserve(key.size());
        for (std::size_t position = 0; position < key.size(); position++)
        {
            const auto atom_index{ key.at(position) };
            group_key_by_position.emplace_back(group_key_by_atom_index.at(atom_index));
            if (EvaluateSuspiciousOffsetUpdate(
                    context.at(atom_index).raw_sampling_entries,
                    previous_state.at(atom_index).mdpde.GetModel(),
                    current_model_snapshot.selected.at(atom_index),
                    options) != SuspiciousGaussianReason::None)
            {
                suspicious_seed_mask.at(position) = 1;
            }
        }
        const auto cluster_rollback_mask{
            ExpandSuspiciousSharedOffsetGroups(group_key_by_position, suspicious_seed_mask)
        };
        for (std::size_t position = 0; position < key.size(); position++)
        {
            if (cluster_rollback_mask.at(position) == 0) continue;
            rollback_atom_mask.at(key.at(position)) = 1;
        }
    }
    for (std::size_t atom_index = 0; atom_index < rollback_atom_mask.size(); atom_index++)
    {
        if (rollback_atom_mask.at(atom_index) == 0) continue;
        current_model_snapshot.selected.at(atom_index) = previous_state.at(atom_index).mdpde.GetModel();
    }

    auto refit_model_snapshot{
        BuildLocalFittingGroupMedianModelList(
            group_key_by_atom_index,
            current_model_snapshot.selected)
    };
    const auto refit_model_bundle{
        BuildSecondStageModelSnapshot(context, std::move(refit_model_snapshot))
    };
    const auto refit_response_cache{
        BuildSecondStageAdjustedResponseCache(context, refit_model_bundle)
    };
    std::vector<std::size_t> refit_atom_index_list;
    for (const auto & key : cluster_key_list)
    {
        for (const auto atom_index : key)
        {
            if (rollback_atom_mask.at(atom_index) == 0)
            {
                refit_atom_index_list.emplace_back(atom_index);
            }
        }
    }
    std::vector<std::optional<LocalAtomRefitResult>> refit_result_list(refit_atom_index_list.size());
    std::vector<std::exception_ptr> refit_exception_list(refit_atom_index_list.size());
#ifdef USE_OPENMP
    const bool parallel_refits{
        !is_debug_logging_enabled &&
        options.thread_size > 1 &&
        refit_atom_index_list.size() > 1
    };
#else
    const bool parallel_refits{ false };
#endif
    FitOptions refit_options{ options };
    if (parallel_refits)
    {
        refit_options.thread_size = 1;
    }
    const auto run_refit = [&](std::size_t refit_position)
    {
        const auto atom_index{ refit_atom_index_list.at(refit_position) };
        try
        {
            refit_result_list.at(refit_position) =
                FitAtomWithJointOffsetFallback(
                    context,
                    atom_index,
                    previous_state.at(atom_index),
                    refit_model_bundle.selected.at(atom_index),
                    refit_response_cache.at(atom_index),
                    refit_options);
        }
        catch (...)
        {
            refit_exception_list.at(refit_position) = std::current_exception();
        }
    };
#ifdef USE_OPENMP
    if (parallel_refits)
    {
        ScopedEigenThreadCount eigen_thread_guard{ 1 };
#pragma omp parallel for schedule(dynamic) num_threads(options.thread_size)
        for (std::size_t refit_position = 0; refit_position < refit_atom_index_list.size(); refit_position++)
        {
            run_refit(refit_position);
        }
    }
    else
#endif
    {
        for (std::size_t refit_position = 0; refit_position < refit_atom_index_list.size(); refit_position++)
        {
            run_refit(refit_position);
        }
    }
    for (const auto & exception : refit_exception_list)
    {
        if (exception) std::rethrow_exception(exception);
    }

    std::vector<std::size_t> post_refit_suspicious_seed_atom_index_list;
    std::size_t refit_position{ 0 };
    for (const auto & key : cluster_key_list)
    {
        auto & health{ health_by_key.at(key) };
        for (const auto atom_index : key)
        {
            if (rollback_atom_mask.at(atom_index) != 0) continue;

            auto refit_result{ std::move(refit_result_list.at(refit_position++)) };
            if (!refit_result.has_value())
            {
                health.is_refit_stationarity_eligible = false;
                post_refit_suspicious_seed_atom_index_list.emplace_back(atom_index);
                continue;
            }
            if (!refit_result->is_stationarity_eligible)
            {
                health.is_refit_stationarity_eligible = false;
            }
            iteration_state.at(atom_index) = std::move(refit_result->result);
        }
    }
    ExpandPostRefitSuspiciousClusters(
        cluster_key_list,
        post_refit_suspicious_seed_atom_index_list,
        rollback_atom_mask);
    for (std::size_t atom_index = 0; atom_index < rollback_atom_mask.size(); atom_index++)
    {
        if (rollback_atom_mask.at(atom_index) == 0) continue;
        iteration_state.at(atom_index) = previous_state.at(atom_index);
    }

    LocalFittingIterationResult iteration_result;
    iteration_result.state = std::move(iteration_state);
    iteration_result.rollback_atom_mask = std::move(rollback_atom_mask);
    iteration_result.health_by_key = std::move(health_by_key);
    return iteration_result;
}

} // namespace rhbm_gem::core::detail
