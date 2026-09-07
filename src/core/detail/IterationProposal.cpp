#include "core/detail/IterationProposal.hpp"

#include "core/detail/Diagnosis.hpp"

#include <exception>
#include <limits>
#include <utility>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/utils/math/EigenHelper.hpp>

namespace rhbm_gem::core::detail {

namespace {

struct LocalAtomRefitResult
{
    LocalGaussianResult result{};
    std::optional<GaussianModel3D> unrestricted_model{};
    SuspiciousGaussianAssessment assessment{};
    std::optional<RHBMEstimationStatus> attempted_refit_status{};
};

static std::optional<LocalAtomRefitResult> FitAtomWithJointOffsetFallback(
    const AtomContext & atom_context,
    const LocalGaussianResult & previous_result,
    const GaussianModel3D & offset_model,
    const std::vector<double> & adjusted_response_list,
    const FitOptions & options)
{
    auto adjusted_sampling_entries{
        BuildSecondStageAdjustedSamples(atom_context, adjusted_response_list)
    };
    const auto & previous_model{ previous_result.mdpde.GetModel() };
    const auto previous_baseline{
        BuildPreviousSuspiciousProfileBaseline(adjusted_sampling_entries, previous_model, options)
    };
    SuspiciousGaussianAssessment failed_shape_assessment;
    std::optional<RHBMEstimationStatus> attempted_refit_status;
    std::optional<GaussianModel3D> unrestricted_model;
    try
    {
        auto candidate_result{
            atom_context.refit_design.Estimate(
                adjusted_response_list,
                atom_context.alpha_r,
                options.thread_size,
                offset_model)
        };
        if (candidate_result.fit_result.has_value())
        {
            attempted_refit_status = candidate_result.fit_result->status;
        }
        if (IsValidSecondStageGaussianModel(candidate_result.mdpde.GetModel()))
        {
            unrestricted_model = candidate_result.mdpde.GetModel();
        }
        const auto assessment{
            AssessSuspiciousGaussianUpdate(
                adjusted_sampling_entries,
                candidate_result.mdpde.GetModel(),
                options,
                previous_baseline,
                SuspiciousUpdateMode::PostRefit)
        };
        if (unrestricted_model.has_value())
        {
            return LocalAtomRefitResult{
                std::move(candidate_result),
                std::move(unrestricted_model),
                assessment,
                attempted_refit_status
            };
        }
        failed_shape_assessment = assessment;
    }
    catch (const std::exception &)
    {
        failed_shape_assessment = SuspiciousGaussianAssessment{
            .reason = SuspiciousGaussianReason::InvalidModel,
            .normalized_margin = std::numeric_limits<double>::infinity()
        };
    }

    auto result{ previous_result };
    SetLocalResultOffset(result, offset_model.GetOffset());
    auto fallback_assessment{
        AssessSuspiciousGaussianUpdate(
            adjusted_sampling_entries,
            result.mdpde.GetModel(),
            options,
            previous_baseline,
            SuspiciousUpdateMode::OffsetOnly)
    };
    if (fallback_assessment.IsSuspicious())
    {
        return std::nullopt;
    }
    return LocalAtomRefitResult{
        std::move(result),
        std::move(unrestricted_model),
        failed_shape_assessment,
        attempted_refit_status
    };
}

static std::vector<std::optional<GaussianModel3D>>
RunUnrestrictedShapeRefits(
    const SecondStageContext & context,
    const FittedGaussianSnapshot & operator_offset_state,
    const FitOptions & options)
{
    const auto operator_model_bundle{
        BuildSecondStageModelSnapshot(context, operator_offset_state)
    };
    const auto adjusted_response_cache{
        BuildSecondStageAdjustedResponseCache(
            context,
            operator_model_bundle)
    };
    std::vector<std::optional<GaussianModel3D>> result(
        context.atom_list.size());
    int refit_thread_size{ options.thread_size };
#ifdef USE_OPENMP
    const bool parallel_refits{
        !IsDebugLogLevelEnabled() &&
        options.thread_size > 1 && context.atom_list.size() > 1
    };
    if (parallel_refits) refit_thread_size = 1;
#pragma omp parallel for schedule(dynamic) if(parallel_refits) num_threads(options.thread_size)
#endif
    for (std::size_t atom_index = 0; atom_index < context.atom_list.size(); atom_index++)
    {
        try
        {
            const auto candidate{
                context.atom_list.at(atom_index).refit_design.Estimate(
                    adjusted_response_cache.at(atom_index),
                    context.atom_list.at(atom_index).alpha_r,
                    refit_thread_size,
                    GetFitModel(operator_model_bundle.node, atom_index))
            };
            if (IsValidSecondStageGaussianModel(candidate.mdpde.GetModel()))
            {
                result.at(atom_index) = candidate.mdpde.GetModel();
            }
        }
        catch (const std::exception &)
        {
        }
    }
    return result;
}

} // namespace

IterationProposalResult BuildIterationProposal(
    const SecondStageContext & context,
    const std::vector<ClusterKey> & cluster_key_list,
    const FitState & previous_state,
    const FitOptions & options,
    const std::vector<double> & ridge_multiplier_list,
    const SuspiciousBlockActivity & quarantine_activity,
    ClusterSolverWorkspaceMap & solver_workspace_by_key,
    std::string_view diagnostic_phase)
{
    auto current_model_snapshot{
        BuildSecondStageModelSnapshot(context, previous_state)
    };
    const auto is_debug_logging_enabled{ IsDebugLogLevelEnabled() };
    const auto log_debug_diagnostics{ !options.quiet_mode && is_debug_logging_enabled };
    std::vector<JointOffsetSolveResult> joint_offset_result_list(cluster_key_list.size());
    std::vector<std::exception_ptr> joint_offset_exception_list(cluster_key_list.size());
    const auto solve_joint_offset = [&](std::size_t cluster_position)
    {
        try
        {
            joint_offset_result_list.at(cluster_position) = EstimateJointOffsets(
                context,
                cluster_key_list.at(cluster_position),
                current_model_snapshot,
                ridge_multiplier_list,
                solver_workspace_by_key.at(
                    cluster_key_list.at(cluster_position)).joint_offset,
                log_debug_diagnostics,
                diagnostic_phase);
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
        eigen_helper::ScopedEigenThreadCount eigen_thread_guard{ 1 };
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

    FixedPointOperatorEvidence fixed_point_operator;
    fixed_point_operator.state = current_model_snapshot.node;
    fixed_point_operator.shape_available_atom_mask.assign(context.atom_list.size(), 0);
    fixed_point_operator.offset_available_atom_mask.assign(context.atom_list.size(), 0);
    for (std::size_t cluster_position = 0;
        cluster_position < cluster_key_list.size(); cluster_position++)
    {
        const auto & key{ cluster_key_list.at(cluster_position) };
        const auto & offset_result{
            joint_offset_result_list.at(cluster_position)
        };
        for (std::size_t position = 0; position < key.size(); position++)
        {
            const auto atom_index{ key.at(position) };
            if (IsJointOffsetSolveHardFailure(offset_result.status))
            {
                continue;
            }
            const auto proposed_offset{
                offset_result.offset(static_cast<Eigen::Index>(position))
            };
            const auto operator_model{
                previous_state.at(atom_index).mdpde.GetModel().WithOffset(proposed_offset)
            };
            if (!IsValidSecondStageGaussianModel(operator_model))
            {
                continue;
            }
            fixed_point_operator.state.at(atom_index) = operator_model;
            fixed_point_operator.offset_available_atom_mask.at(atom_index) = 1;
        }
    }

    auto proposal_state{ previous_state };
    SuspiciousBlockActivity block_activity{
        SuspiciousUpdateMask(context.atom_list.size(), 0),
        SuspiciousUpdateMask(context.atom_list.size(), 0),
        SuspiciousUpdateMask(context.atom_list.size(), 0)
    };
    std::vector<std::optional<RHBMEstimationStatus>>
        local_refit_status_by_atom(context.atom_list.size());
    ClusterHealthMap health_by_key;
    for (std::size_t cluster_position = 0;
        cluster_position < cluster_key_list.size();
        cluster_position++)
    {
        const auto & key{ cluster_key_list.at(cluster_position) };
        const auto & offset_result{ joint_offset_result_list.at(cluster_position) };
        auto [health_iter, inserted]{
            health_by_key.emplace(key, ClusterHealth{ offset_result.status })
        };
        static_cast<void>(inserted);
        auto & health{ health_iter->second };
        if (IsJointOffsetSolveHardFailure(offset_result.status))
        {
            health.all_local_refits_solver_qualified = false;
            for (const auto atom_index : key)
            {
                block_activity.offset_fixed_atom_mask.at(atom_index) = 1;
                block_activity.hard_failure_atom_mask.at(atom_index) = 1;
            }
            continue;
        }

        for (const auto atom_index : key)
        {
            if (!quarantine_activity.HasActiveOffset(atom_index) ||
                fixed_point_operator.offset_available_atom_mask.at(atom_index) == 0)
            {
                health.all_local_refits_solver_qualified = false;
                block_activity.offset_fixed_atom_mask.at(atom_index) = 1;
                continue;
            }
            current_model_snapshot.node.at(atom_index) =
                fixed_point_operator.state.at(atom_index);
        }
    }

    bool operator_offsets_complete{ true };
    bool operator_offsets_match_proposal{ true };
    for (const auto & key : cluster_key_list)
    {
        for (const auto atom_index : key)
        {
            if (fixed_point_operator.offset_available_atom_mask.at(atom_index) == 0)
            {
                operator_offsets_complete = false;
                operator_offsets_match_proposal = false;
                continue;
            }
            operator_offsets_match_proposal =
                operator_offsets_match_proposal &&
                fixed_point_operator.state.at(atom_index).GetOffset() ==
                    current_model_snapshot.node.at(atom_index).GetOffset();
        }
    }

    const auto refit_response_cache{
        BuildSecondStageAdjustedResponseCache(context, current_model_snapshot)
    };
    std::vector<std::optional<LocalAtomRefitResult>> refit_result_list(context.atom_list.size());
    std::vector<std::exception_ptr> refit_exception_list(context.atom_list.size());
#ifdef USE_OPENMP
    const bool parallel_refits{
        !is_debug_logging_enabled &&
        options.thread_size > 1 &&
        context.atom_list.size() > 1
    };
#else
    const bool parallel_refits{ false };
#endif
    FitOptions refit_options{ options };
    if (parallel_refits)
    {
        refit_options.thread_size = 1;
    }
    const auto run_refit = [&](std::size_t atom_index)
    {
        try
        {
            refit_result_list.at(atom_index) =
                FitAtomWithJointOffsetFallback(
                    context.atom_list.at(atom_index),
                    previous_state.at(atom_index),
                    GetFitModel(current_model_snapshot.node, atom_index),
                    refit_response_cache.at(atom_index),
                    refit_options);
        }
        catch (...)
        {
            refit_exception_list.at(atom_index) = std::current_exception();
        }
    };
#ifdef USE_OPENMP
    if (parallel_refits)
    {
        eigen_helper::ScopedEigenThreadCount eigen_thread_guard{ 1 };
#pragma omp parallel for schedule(dynamic) num_threads(options.thread_size)
        for (std::size_t atom_index = 0; atom_index < context.atom_list.size(); atom_index++)
        {
            run_refit(atom_index);
        }
    }
    else
#endif
    {
        for (std::size_t atom_index = 0; atom_index < context.atom_list.size(); atom_index++)
        {
            run_refit(atom_index);
        }
    }
    for (const auto & exception : refit_exception_list)
    {
        if (exception) std::rethrow_exception(exception);
    }

    std::vector<SuspiciousGaussianAssessment> assessment_by_atom(context.atom_list.size());
    for (const auto & key : cluster_key_list)
    {
        auto & health{ health_by_key.at(key) };
        for (const auto atom_index : key)
        {
            auto refit_result{ std::move(refit_result_list.at(atom_index)) };
            if (!refit_result.has_value())
            {
                health.all_local_refits_solver_qualified = false;
                block_activity.shape_fixed_atom_mask.at(atom_index) = 1;
                block_activity.offset_fixed_atom_mask.at(atom_index) = 1;
                assessment_by_atom.at(atom_index) = SuspiciousGaussianAssessment{
                    .reason = SuspiciousGaussianReason::InvalidModel,
                    .normalized_margin = std::numeric_limits<double>::infinity()
                };
                continue;
            }
            if (operator_offsets_complete && operator_offsets_match_proposal)
            {
                if (refit_result->unrestricted_model.has_value())
                {
                    fixed_point_operator.state.at(atom_index) = *refit_result->unrestricted_model;
                    fixed_point_operator.shape_available_atom_mask.at(atom_index) = 1;
                }
            }
            local_refit_status_by_atom.at(atom_index) = refit_result->attempted_refit_status;
            const auto shape_solver_qualified{
                refit_result->attempted_refit_status.has_value() &&
                IsLocalRefitStatusSolverQualified(*refit_result->attempted_refit_status)
            };
            if (!shape_solver_qualified)
            {
                health.all_local_refits_solver_qualified = false;
            }
            if (!refit_result->unrestricted_model.has_value())
            {
                block_activity.shape_fixed_atom_mask.at(atom_index) = 1;
            }
            assessment_by_atom.at(atom_index) = refit_result->assessment;
            proposal_state.at(atom_index) = std::move(refit_result->result);
        }
    }
    if (operator_offsets_complete && !operator_offsets_match_proposal)
    {
        const auto unrestricted_shape_list{
            RunUnrestrictedShapeRefits(
                context,
                fixed_point_operator.state,
                options)
        };
        for (std::size_t atom_index = 0; atom_index < context.atom_list.size(); atom_index++)
        {
            if (unrestricted_shape_list.at(atom_index).has_value())
            {
                fixed_point_operator.state.at(atom_index) = *unrestricted_shape_list.at(atom_index);
                fixed_point_operator.shape_available_atom_mask.at(atom_index) = 1;
            }
        }
    }
    for (std::size_t atom_index = 0; atom_index < context.atom_list.size(); atom_index++)
    {
        if (block_activity.hard_failure_atom_mask.at(atom_index) != 0)
        {
            proposal_state.at(atom_index) = previous_state.at(atom_index);
            continue;
        }
        if (block_activity.offset_fixed_atom_mask.at(atom_index) == 0)
        {
            continue;
        }
        const auto previous_offset{
            previous_state.at(atom_index).mdpde.GetModel().GetOffset()
        };
        SetLocalResultOffset(proposal_state.at(atom_index), previous_offset);
    }

    for (std::size_t atom_index = 0; atom_index < context.atom_list.size(); atom_index++)
    {
        if (!quarantine_activity.HasActiveShape(atom_index))
        {
            const auto accepted_offset{
                proposal_state.at(atom_index).mdpde.GetModel().GetOffset()
            };
            proposal_state.at(atom_index).ols = WithPreservedUncertaintyOffset(
                previous_state.at(atom_index).ols,
                accepted_offset);
            proposal_state.at(atom_index).mdpde = WithPreservedUncertaintyOffset(
                previous_state.at(atom_index).mdpde,
                accepted_offset);
            block_activity.shape_fixed_atom_mask.at(atom_index) = 1;
        }
        if (!quarantine_activity.HasActiveOffset(atom_index))
        {
            const auto previous_offset{
                previous_state.at(atom_index).mdpde.GetModel().GetOffset()
            };
            SetLocalResultOffset(proposal_state.at(atom_index), previous_offset);
            block_activity.offset_fixed_atom_mask.at(atom_index) = 1;
        }
        if (quarantine_activity.hard_failure_atom_mask.at(atom_index) != 0)
        {
            proposal_state.at(atom_index) = previous_state.at(atom_index);
            block_activity.hard_failure_atom_mask.at(atom_index) = 1;
        }
    }
    return IterationProposalResult{
        std::move(proposal_state),
        std::move(fixed_point_operator),
        std::move(block_activity),
        std::move(assessment_by_atom),
        std::move(local_refit_status_by_atom),
        std::move(health_by_key)
    };
}

} // namespace rhbm_gem::core::detail
