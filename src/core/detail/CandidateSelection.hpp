#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "core/detail/CouplingGraph.hpp"
#include "core/detail/JointPolish.hpp"
#include "core/detail/BacktrackingWorkspace.hpp"
#include "core/detail/LocalFittingCandidateEvaluationOverlay.hpp"
#include "core/detail/LocalFittingGroupMedian.hpp"
#include "core/detail/LocalFittingObjective.hpp"
#include "core/detail/LocalFittingPerformanceCounters.hpp"
#include "core/detail/FitStateView.hpp"
#include "core/detail/TransformedChange.hpp"
#include "core/detail/TrustRegion.hpp"
#include "core/detail/ReusableWeightedRidgeSolver.hpp"
#include "core/detail/ScopedEigenThreadCount.hpp"
#include "core/detail/SecondStageContext.hpp"

namespace rhbm_gem::core::detail {

struct RejectedClusterDiagnostic
{
    ClusterKey key{};
    LocalFittingObjectiveAttemptDiagnostic attempt{};
};

struct PolishProgress
{
    std::size_t eligible_count{ 0 };
    std::size_t accepted_count{ 0 };
    std::size_t rejected_count{ 0 };
    std::size_t skipped_count{ 0 };
};

struct CandidateSelection
{
    FitState assembled_state{};
    PolishProvenance assembled_polish_provenance{};
    std::vector<ClusterKey> accepted_key_list{};
    std::vector<ClusterKey> rejected_key_list{};
    std::vector<RejectedClusterDiagnostic> accepted_cluster_diagnostic_list{};
    std::vector<RejectedClusterDiagnostic> rejected_cluster_diagnostic_list{};
    std::vector<ClusterKey> grow_trust_region_key_list{};
    std::vector<ClusterKey> backtracking_exhausted_key_list{};
    std::size_t combined_backtracking_trial_count{ 0 };
    std::optional<double> combined_backtracking_factor{};
    std::optional<LocalFittingObjectiveBreakdown> combined_backtracking_objective{};
    bool combined_backtracking_exhausted{ false };
    PolishProgress polish_progress{};
};


namespace candidate_internal {

struct BaseProposal
{
    FitStatePatch patch{};
    double effective_damping{ 0.0 };
    double step_norm{ 0.0 };
};

struct BaseProposalBuildResult
{
    std::optional<BaseProposal> proposal{};
    LocalFittingPreObjectiveFailureReason failure_reason{ LocalFittingPreObjectiveFailureReason::None };
    std::optional<double> attempted_step_norm{};
};

inline BaseProposalBuildResult BuildSharedOffsetBaseProposal(
    const SecondStageContext & context,
    const FitState & outer_previous_state,
    const FitState & raw_state,
    const ClusterKey & key,
    double trust_region_radius)
{
    if (key.empty())
    {
        return BaseProposalBuildResult{
            std::nullopt,
            LocalFittingPreObjectiveFailureReason::InvalidModel,
            std::nullopt
        };
    }

    std::vector<GroupKey> group_key_by_atom_position;
    std::vector<GaussianModel3D> previous_model_list;
    std::vector<GaussianModel3D> raw_model_list;
    group_key_by_atom_position.reserve(key.size());
    previous_model_list.reserve(key.size());
    raw_model_list.reserve(key.size());
    for (const auto atom_index : key)
    {
        group_key_by_atom_position.emplace_back(context.at(atom_index).group_key);
        previous_model_list.emplace_back(outer_previous_state.at(atom_index).mdpde.GetModel());
        raw_model_list.emplace_back(raw_state.at(atom_index).mdpde.GetModel());
    }
    const auto previous_shared_offset_model_list{
        BuildLocalFittingGroupMedianModelList(group_key_by_atom_position, previous_model_list)
    };
    const auto raw_shared_offset_model_list{
        BuildLocalFittingGroupMedianModelList(group_key_by_atom_position, raw_model_list)
    };

    const auto seed_model_list{
        BuildLocalFittingSharedOffsetDampedModelList(
            previous_model_list,
            raw_model_list,
            previous_shared_offset_model_list,
            raw_shared_offset_model_list,
            0.0)
    };
    if (!seed_model_list.has_value())
    {
        return BaseProposalBuildResult{
            std::nullopt,
            LocalFittingPreObjectiveFailureReason::InvalidModel,
            std::nullopt
        };
    }
    const auto seed_step_norm{
        CalculateClusterModelTrustRegionStepNorm(outer_previous_state, key, *seed_model_list)
    };
    if (!seed_step_norm.has_value())
    {
        return BaseProposalBuildResult{
            std::nullopt,
            LocalFittingPreObjectiveFailureReason::InvalidModel,
            std::nullopt
        };
    }
    if (!IsTrustRegionStepWithinRadius(*seed_step_norm, trust_region_radius))
    {
        return BaseProposalBuildResult{
            std::nullopt,
            LocalFittingPreObjectiveFailureReason::PreviousSharedOffsetProjectionOutsideTrustRegion,
            *seed_step_norm
        };
    }

    double damping{ 1.0 };
    std::optional<double> attempted_step_norm;
    while (damping >= std::numeric_limits<double>::epsilon())
    {
        auto candidate_model_list{
            BuildLocalFittingSharedOffsetDampedModelList(
                previous_model_list,
                raw_model_list,
                previous_shared_offset_model_list,
                raw_shared_offset_model_list,
                damping)
        };
        if (candidate_model_list.has_value())
        {
            const auto step_norm{
                CalculateClusterModelTrustRegionStepNorm(
                    outer_previous_state,
                    key,
                    *candidate_model_list)
            };
            if (step_norm.has_value() &&
                IsTrustRegionStepWithinRadius(*step_norm, trust_region_radius))
            {
                BaseProposal proposal;
                proposal.patch.atom_index_list = key;
                proposal.patch.mdpde_list.reserve(key.size());
                proposal.effective_damping = damping;
                proposal.step_norm = *step_norm;
                for (std::size_t atom_position = 0; atom_position < key.size(); atom_position++)
                {
                    const auto atom_index{ key.at(atom_position) };
                    proposal.patch.mdpde_list.emplace_back(
                        GaussianModel3DWithUncertainty{
                            candidate_model_list->at(atom_position),
                            raw_state.at(atom_index).mdpde
                                .GetStandardDeviationModel()
                        });
                }
                return BaseProposalBuildResult{
                    std::move(proposal),
                    LocalFittingPreObjectiveFailureReason::None,
                    *step_norm
                };
            }
            if (step_norm.has_value()) attempted_step_norm = *step_norm;
        }
        damping *= 0.5;
    }
    return BaseProposalBuildResult{
        std::nullopt,
        LocalFittingPreObjectiveFailureReason::NoCandidateWithinTrustRegion,
        attempted_step_norm
    };
}


struct ClusterCandidateResult
{
    ClusterKey key{};
    std::optional<FitStatePatch> accepted_patch{};
    std::vector<char> polish_provenance{};
    LocalFittingClusterObjectiveState objective_state{};
    LocalFittingObjectiveAttemptDiagnostic diagnostic{};
    PolishProgress polish_progress{};
    bool accepted{ false };
    bool grow_trust_region{ false };
};

} // namespace candidate_internal


inline void RejectCombinedCandidate(
    const FitState & previous_state,
    const PolishProvenance & previous_polish_provenance,
    const std::vector<ClusterKey> & cluster_key_list,
    CandidateSelection & selection)
{
    selection.assembled_state = previous_state;
    selection.assembled_polish_provenance = previous_polish_provenance;
    selection.accepted_key_list.clear();
    selection.rejected_key_list = cluster_key_list;
    selection.grow_trust_region_key_list.clear();
    selection.combined_backtracking_objective.reset();
    selection.polish_progress.rejected_count += selection.polish_progress.accepted_count;
    selection.polish_progress.accepted_count = 0;
}

inline std::vector<Eigen::Vector3d> InterpolateTransformedEstimations(
    const std::vector<Eigen::Vector3d> & previous_estimation_list,
    const std::vector<Eigen::Vector3d> & candidate_estimation_list,
    double damping)
{
    if (previous_estimation_list.size() != candidate_estimation_list.size() ||
        !std::isfinite(damping) || damping < 0.0 || damping > 1.0)
    {
        throw std::invalid_argument("Local fitting transformed interpolation inputs are invalid.");
    }
    std::vector<Eigen::Vector3d> interpolated_list;
    interpolated_list.reserve(previous_estimation_list.size());
    for (std::size_t i = 0; i < previous_estimation_list.size(); i++)
    {
        interpolated_list.emplace_back(
            (previous_estimation_list.at(i) + damping * (candidate_estimation_list.at(i) - previous_estimation_list.at(i))).eval());
    }
    return interpolated_list;
}

inline std::optional<FitStatePatch> BuildCandidatePatch(
    const FitState & previous_state,
    const std::vector<Eigen::Vector3d> & previous_transformed_estimation_list,
    const std::vector<Eigen::Vector3d> & candidate_transformed_estimation_list,
    const FitState & uncertainty_state,
    const std::vector<std::size_t> & active_index_list)
{
    if (previous_transformed_estimation_list.size() != active_index_list.size() ||
        candidate_transformed_estimation_list.size() != active_index_list.size())
    {
        throw std::invalid_argument(
            "Local fitting candidate transformed coordinate count is inconsistent.");
    }
    FitStatePatch candidate_patch;
    candidate_patch.atom_index_list = active_index_list;
    candidate_patch.mdpde_list.reserve(active_index_list.size());
    for (std::size_t local_position = 0; local_position < active_index_list.size(); local_position++)
    {
        const auto active_index{ active_index_list.at(local_position) };
        const auto & previous_transformed_estimation{
            previous_transformed_estimation_list.at(local_position)
        };
        const auto & candidate_transformed_estimation{
            candidate_transformed_estimation_list.at(local_position)
        };
        if (!previous_transformed_estimation.allFinite() || !candidate_transformed_estimation.allFinite())
        {
            return std::nullopt;
        }
        if ((candidate_transformed_estimation.array() == previous_transformed_estimation.array()).all())
        {
            candidate_patch.mdpde_list.emplace_back(
                GaussianModel3DWithUncertainty{
                previous_state.at(active_index).mdpde.GetModel(),
                uncertainty_state.at(active_index).mdpde.GetStandardDeviationModel()
                });
            continue;
        }
        const auto candidate_model{
            DecodeTransformedCoordinates(candidate_transformed_estimation)
        };
        if (!candidate_model.has_value())
        {
            return std::nullopt;
        }
        if (!IsValidSecondStageGaussianModel(*candidate_model))
        {
            return std::nullopt;
        }
        candidate_patch.mdpde_list.emplace_back(
            GaussianModel3DWithUncertainty{
            *candidate_model,
            uncertainty_state.at(active_index).mdpde.GetStandardDeviationModel()
            });
    }
    return candidate_patch;
}

inline bool ShouldGrowTrustRegion(const LocalFittingObjectiveAttemptDiagnostic & diagnostic)
{
    return diagnostic.candidate_objective.has_value() &&
        diagnostic.previous_objective.has_value() &&
        IsTrustRegionGrowthEligible(
            diagnostic.trust_region_step_norm,
            diagnostic.trust_region_radius,
            IsBetterLocalFittingAuditObjective(
                diagnostic.candidate_objective->total_objective,
                diagnostic.previous_objective->total_objective,
                kLocalFittingObjectiveStrictTolerance));
}

inline candidate_internal::ClusterCandidateResult SelectClusterCandidate(
    const SecondStageContext & context,
    const SecondStageModelSnapshot & previous_model_snapshot,
    const LocalFittingResidualBaseline & residual_baseline,
    const ClusterKey & key,
    const std::vector<LocalFittingObjectiveSampleRef> & objective_sample_ref_list,
    bool is_polish_eligible,
    bool is_unchanged_state_exhausted,
    const FitState & previous_state,
    const PolishProvenance & previous_polish_provenance,
    const FitState & raw_state,
    const std::vector<Eigen::Vector3d> & previous_transformed_estimation_list,
    const std::vector<Eigen::Vector3d> & raw_transformed_estimation_list,
    const std::vector<char> & rollback_atom_mask,
    const std::vector<double> & ridge_multiplier_list,
    const LocalFittingObjectiveDomain & objective_domain,
    const LocalFittingObjectiveByKey & previous_objective_by_key,
    const LocalFittingClusterObjectiveState & previous_objective_state,
    double trust_region_radius,
    LocalFittingClusterSolverWorkspace & solver_workspace,
    LocalFittingPerformanceCounters & performance_counters)
{
    using candidate_internal::BuildSharedOffsetBaseProposal;
    using candidate_internal::ClusterCandidateResult;
    using candidate_internal::BaseProposal;
    ClusterCandidateResult result;
    result.key = key;
    result.objective_state = previous_objective_state;
    result.polish_provenance.reserve(key.size());
    for (const auto atom_index : key)
    {
        result.polish_provenance.emplace_back(
            previous_polish_provenance.at(atom_index));
    }
    if (is_polish_eligible) result.polish_progress.eligible_count = 1;
    result.diagnostic.trust_region_radius = trust_region_radius;
    if (is_unchanged_state_exhausted)
    {
        if (is_polish_eligible) result.polish_progress.skipped_count = 1;
        result.diagnostic.backtracking_exhausted = true;
        return result;
    }

    LocalFittingClusterObjectiveStateMap objective_state;
    objective_state.emplace(key, result.objective_state);
    const FitStateView previous_state_view{ previous_state };
    const auto contains_suspicious_atom{
        std::any_of(
            key.begin(),
            key.end(),
            [&](std::size_t atom_index)
            {
                return rollback_atom_mask.at(atom_index) != 0;
            })
    };
    std::optional<BaseProposal> base_proposal;
    if (!contains_suspicious_atom)
    {
        auto proposal_result{
            BuildSharedOffsetBaseProposal(
                context,
                previous_state,
                raw_state,
                key,
                trust_region_radius)
        };
        result.diagnostic.pre_objective_failure_reason =
            proposal_result.failure_reason;
        if (proposal_result.attempted_step_norm.has_value())
        {
            result.diagnostic.pre_objective_attempted_step_norm =
                proposal_result.attempted_step_norm;
            result.diagnostic.trust_region_step_norm =
                *proposal_result.attempted_step_norm;
        }
        base_proposal = std::move(proposal_result.proposal);
    }
    else
    {
        std::vector<Eigen::Vector3d> previous_cluster_estimation_list;
        std::vector<Eigen::Vector3d> raw_cluster_estimation_list;
        previous_cluster_estimation_list.reserve(key.size());
        raw_cluster_estimation_list.reserve(key.size());
        for (const auto atom_index : key)
        {
            previous_cluster_estimation_list.emplace_back(
                previous_transformed_estimation_list.at(atom_index));
            raw_cluster_estimation_list.emplace_back(
                raw_transformed_estimation_list.at(atom_index));
        }
        const auto trust_region_damping{
            LimitTrustRegionDamping(
                previous_cluster_estimation_list,
                raw_cluster_estimation_list,
                kTrustRegionParameterScale,
                1.0,
                trust_region_radius)
        };
        const auto base_cluster_estimation_list{
            InterpolateTransformedEstimations(
                previous_cluster_estimation_list,
                raw_cluster_estimation_list,
                trust_region_damping.effective_damping)
        };
        auto base_patch{
            BuildCandidatePatch(
                previous_state,
                previous_cluster_estimation_list,
                base_cluster_estimation_list,
                raw_state,
                key)
        };
        if (base_patch.has_value())
        {
            base_proposal = BaseProposal{
                std::move(*base_patch),
                trust_region_damping.effective_damping,
                trust_region_damping.step_norm
            };
        }
        else
        {
            result.diagnostic.is_invalid_model = true;
            result.diagnostic.pre_objective_failure_reason =
                LocalFittingPreObjectiveFailureReason::InvalidModel;
        }
    }
    if (!base_proposal.has_value())
    {
        if (is_polish_eligible) result.polish_progress.skipped_count = 1;
        return result;
    }

    result.diagnostic.effective_damping = base_proposal->effective_damping;
    result.diagnostic.trust_region_step_norm = base_proposal->step_norm;
    result.diagnostic.backtracking_trial_count = 1;
    result.diagnostic.accepted_backtracking_factor = 1.0;
    auto & base_patch{ base_proposal->patch };
    FitStateView base_state_view{ previous_state, base_patch };
    BacktrackingWorkspace backtracking_workspace{
        context,
        previous_state,
        base_state_view,
        key,
        kTransformedChangeTolerance
    };
    LocalFittingCandidateEvaluationOverlay base_overlay{
        context,
        previous_model_snapshot,
        residual_baseline,
        base_state_view,
        base_patch
    };
    auto accepted_base_candidate{
        TryCommitLocalFittingClusterCandidate(
            context,
            base_state_view,
            base_overlay,
            previous_state_view,
            key,
            objective_sample_ref_list,
            previous_objective_by_key.at(key),
            false,
            objective_domain,
            objective_state,
            result.diagnostic,
            &performance_counters)
    };
    auto accepted_by_backtracking{ false };
    const PolishProvenance non_polished_endpoint_provenance(key.size(), 0);
    if (!accepted_base_candidate)
    {
        result.diagnostic.accepted_backtracking_factor.reset();
        const auto step{ backtracking_workspace.FindAcceptedCandidate(
            [&](const BacktrackingStep & step)
            {
                const auto factor{ step.factor };
                const auto * backtracked_patch{ step.candidate_patch };
                LocalFittingObjectiveAttemptDiagnostic trial_diagnostic;
                trial_diagnostic.effective_damping = base_proposal->effective_damping * factor;
                trial_diagnostic.trust_region_radius = trust_region_radius;
                trial_diagnostic.trust_region_step_norm = base_proposal->step_norm * factor;
                trial_diagnostic.backtracking_trial_count = step.trial_number;
                const FitStateView backtracked_state_view{
                    previous_state,
                    *backtracked_patch
                };
                const LocalFittingCandidateEvaluationOverlay backtracked_overlay{
                    context,
                    previous_model_snapshot,
                    residual_baseline,
                    backtracked_state_view,
                    *backtracked_patch
                };
                if (TryCommitLocalFittingClusterCandidate(
                        context,
                        backtracked_state_view,
                        backtracked_overlay,
                        previous_state_view,
                        key,
                        objective_sample_ref_list,
                        previous_objective_by_key.at(key),
                        false,
                        objective_domain,
                        objective_state,
                        trial_diagnostic,
                        &performance_counters))
                {
                    trial_diagnostic.accepted_backtracking_factor = factor;
                    result.diagnostic = std::move(trial_diagnostic);
                    return true;
                }
                result.diagnostic = std::move(trial_diagnostic);
                return false;
            }) };
        if (step.status == BacktrackingStepStatus::Exhausted)
        {
            result.diagnostic.backtracking_exhausted = true;
        }
        else if (step.status == BacktrackingStepStatus::InvalidCandidate)
        {
            result.diagnostic.is_invalid_model = true;
        }
        else
        {
            result.polish_provenance =
                backtracking_workspace.BuildActiveCandidatePolishProvenance(
                    result.polish_provenance,
                    non_polished_endpoint_provenance);
            base_patch = backtracking_workspace.TakeCandidatePatch();
            accepted_base_candidate = true;
            accepted_by_backtracking = true;
        }
    }
    else
    {
        result.polish_provenance =
            backtracking_workspace.BuildActiveCandidatePolishProvenance(
                result.polish_provenance,
                non_polished_endpoint_provenance);
    }
    if (!accepted_base_candidate)
    {
        if (is_polish_eligible) result.polish_progress.skipped_count = 1;
        result.objective_state = objective_state.at(key);
        return result;
    }

    result.accepted = true;
    result.grow_trust_region = !accepted_by_backtracking && ShouldGrowTrustRegion(result.diagnostic);
    result.accepted_patch = base_patch;
    if (is_polish_eligible)
    {
        auto polished_candidate{
            BuildJointPolishProposal(
                context,
                previous_state,
                base_state_view,
                key,
                objective_sample_ref_list,
                ridge_multiplier_list,
                solver_workspace.joint_polish,
                trust_region_radius)
        };
        if (!polished_candidate.has_value())
        {
            result.polish_progress.skipped_count = 1;
        }
        else
        {
            LocalFittingObjectiveAttemptDiagnostic polish_diagnostic;
            polish_diagnostic.effective_damping = polished_candidate->effective_damping;
            polish_diagnostic.trust_region_radius = trust_region_radius;
            polish_diagnostic.trust_region_step_norm = polished_candidate->step_norm;
            const FitStateView polished_state_view{
                previous_state,
                polished_candidate->patch
            };
            const LocalFittingCandidateEvaluationOverlay polished_overlay{
                context,
                previous_model_snapshot,
                residual_baseline,
                polished_state_view,
                polished_candidate->patch
            };
            if (!TryCommitLocalFittingClusterCandidate(
                    context,
                    polished_state_view,
                    polished_overlay,
                    base_state_view,
                    key,
                    objective_sample_ref_list,
                    result.diagnostic.candidate_objective,
                    true,
                    objective_domain,
                    objective_state,
                    polish_diagnostic,
                    &performance_counters))
            {
                result.polish_progress.rejected_count = 1;
            }
            else
            {
                result.polish_progress.accepted_count = 1;
                result.accepted_patch = std::move(polished_candidate->patch);
                for (const auto atom_index : polished_candidate->changed_atom_index_list)
                {
                    const auto position{
                        static_cast<std::size_t>(std::distance(
                            key.begin(),
                            std::find(key.begin(), key.end(), atom_index)))
                    };
                    result.polish_provenance.at(position) = 1;
                }
                if (!accepted_by_backtracking && ShouldGrowTrustRegion(polish_diagnostic))
                {
                    result.grow_trust_region = true;
                }
            }
        }
    }
    result.objective_state = objective_state.at(key);
    return result;
}

inline CandidateSelection SelectClusterCandidates(
    const SecondStageContext & context,
    const SecondStageModelSnapshot & previous_model_snapshot,
    const LocalFittingResidualBaseline & residual_baseline,
    const CouplingGraphPartition & partition,
    const std::vector<ClusterKey> & polish_eligible_key_list,
    const FitState & previous_state,
    const PolishProvenance & previous_polish_provenance,
    const FitState & raw_state,
    const std::vector<Eigen::Vector3d> & previous_transformed_estimation_list,
    const std::vector<Eigen::Vector3d> & raw_transformed_estimation_list,
    const std::vector<char> & rollback_atom_mask,
    const std::vector<double> & ridge_multiplier_list,
    const std::vector<ClusterKey> & unchanged_state_exhausted_key_list,
    const LocalFittingObjectiveDomain & objective_domain,
    const LocalFittingObjectiveByKey & previous_objective_by_key,
    LocalFittingClusterObjectiveStateMap & cluster_objective_state,
    const TrustRegionStateSet & trust_region_state,
    LocalFittingClusterSolverWorkspaceMap & solver_workspace_by_key,
    int thread_size,
    LocalFittingPerformanceCounters & performance_counters)
{
    std::vector<const ClusterKey *> key_list;
    std::vector<const std::vector<LocalFittingObjectiveSampleRef> *> sample_ref_list_by_position;
    std::vector<LocalFittingClusterSolverWorkspace *> solver_workspace_list;
    std::vector<const LocalFittingClusterObjectiveState *> objective_state_list;
    std::vector<char> polish_eligible_list;
    std::vector<char> exhausted_list;
    key_list.reserve(partition.sample_id_list_by_key.size());
    sample_ref_list_by_position.reserve(partition.sample_id_list_by_key.size());
    solver_workspace_list.reserve(partition.sample_id_list_by_key.size());
    objective_state_list.reserve(partition.sample_id_list_by_key.size());
    for (const auto & [key, sample_ref_list] : partition.sample_id_list_by_key)
    {
        key_list.emplace_back(&key);
        sample_ref_list_by_position.emplace_back(&sample_ref_list);
        solver_workspace_list.emplace_back(&solver_workspace_by_key.at(key));
        objective_state_list.emplace_back(&cluster_objective_state.at(key));
        polish_eligible_list.emplace_back(
            std::find(
                polish_eligible_key_list.begin(),
                polish_eligible_key_list.end(),
                key) != polish_eligible_key_list.end());
        exhausted_list.emplace_back(
            std::find(
                unchanged_state_exhausted_key_list.begin(),
                unchanged_state_exhausted_key_list.end(),
                key) != unchanged_state_exhausted_key_list.end());
    }

    using candidate_internal::ClusterCandidateResult;
    std::vector<ClusterCandidateResult> result_list(key_list.size());
    std::vector<std::exception_ptr> exception_list(key_list.size());
    const auto select_candidate = [&](std::size_t position)
    {
        try
        {
            result_list.at(position) = SelectClusterCandidate(
                context,
                previous_model_snapshot,
                residual_baseline,
                *key_list.at(position),
                *sample_ref_list_by_position.at(position),
                polish_eligible_list.at(position) != 0,
                exhausted_list.at(position) != 0,
                previous_state,
                previous_polish_provenance,
                raw_state,
                previous_transformed_estimation_list,
                raw_transformed_estimation_list,
                rollback_atom_mask,
                ridge_multiplier_list,
                objective_domain,
                previous_objective_by_key,
                *objective_state_list.at(position),
                trust_region_state.GetRadius(*key_list.at(position)),
                *solver_workspace_list.at(position),
                performance_counters);
        }
        catch (...)
        {
            exception_list.at(position) = std::current_exception();
        }
    };
#ifdef USE_OPENMP
    if (thread_size > 1 && key_list.size() > 1)
    {
        ScopedEigenThreadCount eigen_thread_guard{ 1 };
#pragma omp parallel for schedule(dynamic) num_threads(thread_size)
        for (std::size_t position = 0; position < key_list.size(); position++)
        {
            select_candidate(position);
        }
    }
    else
#endif
    {
        for (std::size_t position = 0; position < key_list.size(); position++)
        {
            select_candidate(position);
        }
    }
    for (const auto & exception : exception_list)
    {
        if (exception) std::rethrow_exception(exception);
    }

    CandidateSelection selection;
    selection.assembled_state = previous_state;
    selection.assembled_polish_provenance = previous_polish_provenance;
    for (auto & result : result_list)
    {
        cluster_objective_state.at(result.key) = std::move(result.objective_state);
        selection.polish_progress.eligible_count += result.polish_progress.eligible_count;
        selection.polish_progress.accepted_count += result.polish_progress.accepted_count;
        selection.polish_progress.rejected_count += result.polish_progress.rejected_count;
        selection.polish_progress.skipped_count += result.polish_progress.skipped_count;
        if (!result.accepted)
        {
            selection.rejected_key_list.emplace_back(result.key);
            if (result.diagnostic.backtracking_exhausted)
            {
                selection.backtracking_exhausted_key_list.emplace_back(result.key);
            }
            selection.rejected_cluster_diagnostic_list.emplace_back(
                RejectedClusterDiagnostic{
                    result.key,
                    std::move(result.diagnostic)
                });
            continue;
        }
        selection.accepted_key_list.emplace_back(result.key);
        selection.accepted_cluster_diagnostic_list.emplace_back(
            RejectedClusterDiagnostic{
                result.key,
                result.diagnostic
            });
        if (result.grow_trust_region)
        {
            selection.grow_trust_region_key_list.emplace_back(result.key);
        }
        result.accepted_patch->ApplyTo(selection.assembled_state);
        for (std::size_t position = 0; position < result.key.size(); position++)
        {
            selection.assembled_polish_provenance.at(result.key.at(position)) = result.polish_provenance.at(position);
        }
    }
    return selection;
}

inline bool TryBacktrackCombinedCandidate(
    const SecondStageContext & context,
    const SecondStageModelSnapshot & previous_model_snapshot,
    const LocalFittingResidualBaseline & residual_baseline,
    const CouplingGraphPartition & partition,
    const FitState & previous_state,
    const PolishProvenance & previous_polish_provenance,
    const LocalFittingObjectiveDomain & objective_domain,
    const LocalFittingObjectiveByKey & previous_objective_by_key,
    const std::optional<LocalFittingObjectiveBreakdown> & previous_audit_objective,
    const LocalFittingBestAuditState & best_audit_state,
    const LocalFittingClusterObjectiveStateMap & committed_objective_state,
    LocalFittingClusterObjectiveStateMap & working_objective_state,
    CandidateSelection & selection,
    LocalFittingPerformanceCounters & performance_counters)
{
    const auto endpoint_state{ selection.assembled_state };
    const auto endpoint_provenance{ selection.assembled_polish_provenance };
    const auto accepted_key_list{ selection.accepted_key_list };
    std::vector<std::size_t> changed_atom_index_list;
    for (const auto & key : accepted_key_list)
    {
        changed_atom_index_list.insert(
            changed_atom_index_list.end(),
            key.begin(),
            key.end());
    }
    const auto affected_sample_ref_list{
        BuildGraphAffectedSampleUnion(partition, accepted_key_list)
    };
    BacktrackingWorkspace backtracking_workspace{
        context,
        previous_state,
        endpoint_state,
        changed_atom_index_list,
        kTransformedChangeTolerance
    };
    const FitStateView previous_state_view{ previous_state };

    selection.combined_backtracking_trial_count = 1;
    LocalFittingClusterObjectiveStateMap accepted_trial_objective_state;
    const auto step{ backtracking_workspace.FindAcceptedCandidate(
        [&](const BacktrackingStep & step)
        {
            const auto factor{ step.factor };
            const auto * candidate_patch{ step.candidate_patch };
            const FitStateView candidate_state_view{
                previous_state,
                *candidate_patch
            };
            const LocalFittingCandidateEvaluationOverlay candidate_overlay{
                context,
                previous_model_snapshot,
                residual_baseline,
                candidate_state_view,
                *candidate_patch
            };
            selection.combined_backtracking_trial_count = step.trial_number;
            auto trial_objective_state{ committed_objective_state };
            auto local_criteria_accepted{ true };
            for (const auto & key : accepted_key_list)
            {
                const auto sample_iter{
                    partition.sample_id_list_by_key.find(key)
                };
                if (sample_iter == partition.sample_id_list_by_key.end())
                {
                    local_criteria_accepted = false;
                    break;
                }
                LocalFittingObjectiveAttemptDiagnostic diagnostic;
                diagnostic.backtracking_trial_count = selection.combined_backtracking_trial_count;
                diagnostic.accepted_backtracking_factor = factor;
                if (!TryCommitLocalFittingClusterCandidate(
                        context,
                        candidate_state_view,
                        candidate_overlay,
                        previous_state_view,
                        key,
                        sample_iter->second,
                        previous_objective_by_key.at(key),
                        false,
                        objective_domain,
                        trial_objective_state,
                        diagnostic,
                        &performance_counters))
                {
                    local_criteria_accepted = false;
                    break;
                }
            }
            const auto combined_check{
                local_criteria_accepted ?
                    EvaluateLocalFittingCombinedObjective(
                        context,
                        previous_state,
                        candidate_state_view,
                        residual_baseline,
                        candidate_overlay,
                        candidate_patch->atom_index_list,
                        affected_sample_ref_list,
                        objective_domain,
                        best_audit_state,
                        previous_audit_objective,
                        &performance_counters) :
                    LocalFittingCombinedObjectiveCheck{}
            };
            if (combined_check.accepted)
            {
                selection.combined_backtracking_factor = factor;
                selection.combined_backtracking_objective = combined_check.candidate_objective;
                accepted_trial_objective_state = std::move(trial_objective_state);
                return true;
            }
            return false;
        }) };
    if (step.status == BacktrackingStepStatus::Exhausted)
    {
        selection.combined_backtracking_exhausted = true;
        return false;
    }
    if (step.status == BacktrackingStepStatus::InvalidCandidate)
    {
        return false;
    }

    auto candidate_state{ backtracking_workspace.MaterializeCandidateState() };
    performance_counters.full_state_materialization_count++;
    selection.assembled_state = std::move(candidate_state);
    selection.assembled_polish_provenance =
        backtracking_workspace.BuildCandidatePolishProvenance(
            previous_polish_provenance,
            endpoint_provenance);
    selection.grow_trust_region_key_list.clear();
    working_objective_state = std::move(accepted_trial_objective_state);
    return true;
}



} // namespace rhbm_gem::core::detail
