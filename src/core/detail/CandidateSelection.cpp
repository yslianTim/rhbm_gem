#include "core/detail/ClusterHistoryObserver.hpp"
#include "utils/hrl/EstimationAudit.hpp"
#include "core/detail/PhaseAudit.hpp"
#include "core/detail/TrustModelAudit.hpp"
#include "core/detail/CandidateSelection.hpp"
#include "core/detail/CandidateTransaction.hpp"
#include "core/detail/CandidateEvaluation.hpp"
#include "core/detail/Diagnosis.hpp"
#include "core/detail/GaussianModelOperations.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
#include <utility>

#include <Eigen/Dense>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/utils/math/EigenHelper.hpp>

namespace rhbm_gem::core::detail {
namespace {
constexpr double kTrustRegionInitialRadius{ 1.0 };
constexpr double kTrustRegionMinimumRadius{ 0.0625 };
constexpr double kTrustRegionShrinkFactor{ 0.5 };

struct ClusterCandidateResult
{
    std::optional<FitStatePatch> accepted_patch{};
    std::optional<FitStatePatch> rescue_patch{};
    PolishProvenance polish_provenance{};
    ObjectiveAttemptDiagnostic diagnostic{};
    PolishProgress polish_progress{};
    bool shrink_trust_region{ false };

};

} // namespace

void TrustRegionStateSet::Reconcile(
    const std::vector<ClusterKey> & key_list)
{
    std::map<ClusterKey, double> next_radius_by_key;
    for (const auto & key : key_list)
    {
        const auto iter{ m_radius_by_key.find(key) };
        next_radius_by_key.emplace(
            key,
            iter == m_radius_by_key.end() ?
                kTrustRegionInitialRadius : iter->second);
    }
    m_radius_by_key = std::move(next_radius_by_key);
}

double TrustRegionStateSet::GetRadius(const ClusterKey & key) const
{
    const auto iter{ m_radius_by_key.find(key) };
    if (iter == m_radius_by_key.end())
    {
        throw std::invalid_argument(
            "Local fitting trust-region state is missing.");
    }
    return iter->second;
}

void TrustRegionStateSet::ResetToMinimum(const std::vector<ClusterKey> & key_list)
{
    for (const auto & key : key_list)
    {
        const auto iter{ m_radius_by_key.find(key) };
        if (iter == m_radius_by_key.end())
        {
            throw std::invalid_argument("Local fitting trust-region state is missing.");
        }
        iter->second = kTrustRegionMinimumRadius;
    }
}

TrustRegionRadiusUpdate TrustRegionStateSet::ApplyRadiusUpdates(
    const std::vector<ClusterKey> & accepted_shrink_key_list,
    const std::vector<ClusterKey> & rejected_key_list,
    const std::vector<ClusterKey> & exhausted_key_list)
{
    TrustRegionRadiusUpdate update;
    const auto shrink = [&](const std::vector<ClusterKey> & key_list)
    {
        for (const auto & key : key_list)
        {
            auto iter{ m_radius_by_key.find(key) };
            if (iter == m_radius_by_key.end())
            {
                throw std::invalid_argument(
                    "Local fitting trust-region state is missing.");
            }
            if (iter->second <= kTrustRegionMinimumRadius)
            {
                update.saturated_key_list.emplace_back(key);
                continue;
            }
            iter->second = std::max(
                kTrustRegionMinimumRadius,
                iter->second * kTrustRegionShrinkFactor);
            update.changed_key_list.emplace_back(key);
        }
    };
    shrink(accepted_shrink_key_list);
    std::vector<ClusterKey> retryable_key_list;
    for (const auto & key : rejected_key_list)
    {
        if (std::ranges::find(exhausted_key_list, key) == exhausted_key_list.end())
        {
            retryable_key_list.emplace_back(key);
        }
    }
    shrink(retryable_key_list);
    return update;
}

BacktrackingStep BacktrackingWorkspace::BuildNextCandidate()
{
    const auto factor{ m_next_factor };
    if (!std::isfinite(factor) || factor < std::numeric_limits<double>::epsilon())
    {
        return BacktrackingStep{
            BacktrackingStepStatus::Exhausted,
            factor,
            m_trial_number
        };
    }
    m_next_factor *= 0.5;

    if (!BuildCandidate(factor))
    {
        return BacktrackingStep{
            BacktrackingStepStatus::InvalidCandidate,
            factor,
            m_trial_number
        };
    }
    const auto maximum_transformed_change{ GetMaximumTransformedChange() };
    if (maximum_transformed_change < m_minimum_transformed_change)
    {
        return BacktrackingStep{
            BacktrackingStepStatus::Exhausted,
            factor,
            m_trial_number
        };
    }
    m_trial_number++;
    return BacktrackingStep{
        BacktrackingStepStatus::CandidateReady,
        factor,
        m_trial_number
    };
}

PolishProvenance BacktrackingWorkspace::BuildCandidatePolishProvenance(
    const PolishProvenance & previous_provenance,
    const PolishProvenance & endpoint_provenance) const
{
    if (previous_provenance.size() != m_previous_state_size ||
        endpoint_provenance.size() != m_previous_state_size)
    {
        throw std::invalid_argument(
            "Local fitting backtracking provenance sizes are inconsistent.");
    }
    auto provenance{ previous_provenance };
    for (std::size_t i = 0; i < m_candidate_patch.atom_index_list.size(); i++)
    {
        if (IsTransformedChangeMaterial(
                CalculateTransformedChange(
                    m_candidate_patch.mdpde_list.at(i).GetModel(),
                    m_previous_model_list.at(i)),
                m_minimum_transformed_change))
        {
            provenance.at(m_candidate_patch.atom_index_list.at(i)) =
                endpoint_provenance.at(m_candidate_patch.atom_index_list.at(i));
        }
    }
    return provenance;
}

bool BacktrackingWorkspace::BuildCandidate(double factor)
{
    const auto candidate_model_list{
        BuildDampedModelList(
            m_previous_model_list,
            m_endpoint_model_list,
            factor)
    };
    if (!candidate_model_list.has_value())
    {
        return false;
    }

    for (std::size_t i = 0; i < m_candidate_patch.atom_index_list.size(); i++)
    {
        const auto endpoint_uncertainty{
            m_candidate_patch.mdpde_list.at(i).GetStandardDeviationModel()
        };
        m_candidate_patch.mdpde_list.at(i) =
            GaussianModel3DWithUncertainty{
                candidate_model_list->at(i),
                endpoint_uncertainty
            };
    }
    return true;
}

double BacktrackingWorkspace::GetMaximumTransformedChange() const
{
    double maximum_change{ 0.0 };
    for (std::size_t i = 0; i < m_candidate_patch.atom_index_list.size(); i++)
    {
        maximum_change = std::max(
            maximum_change,
            std::ranges::max(
                CalculateTransformedChange(
                    m_candidate_patch.mdpde_list.at(i).GetModel(),
                    m_previous_model_list.at(i))));
    }
    return maximum_change;
}

BacktrackingWorkspace::BacktrackingWorkspace(
    const FitState & previous_state,
    const FitStatePatch & endpoint_patch,
    double minimum_transformed_change)
    : m_previous_state_size{ previous_state.size() },
      m_minimum_transformed_change{ minimum_transformed_change }
{
    if (!std::isfinite(m_minimum_transformed_change) ||
        m_minimum_transformed_change < 0.0)
    {
        throw std::invalid_argument(
            "Local fitting backtracking minimum transformed change is invalid.");
    }
    m_candidate_patch = endpoint_patch;
    m_previous_model_list.reserve(m_candidate_patch.atom_index_list.size());
    m_endpoint_model_list.reserve(m_candidate_patch.atom_index_list.size());
    for (std::size_t atom_position = 0;
        atom_position < m_candidate_patch.atom_index_list.size();
        atom_position++)
    {
        const auto atom_index{ m_candidate_patch.atom_index_list.at(atom_position) };
        m_previous_model_list.emplace_back(previous_state.at(atom_index).mdpde.GetModel());
        const auto & endpoint_mdpde{
            m_candidate_patch.mdpde_list.at(atom_position)
        };
        m_endpoint_model_list.emplace_back(endpoint_mdpde.GetModel());
    }
}

static std::optional<FitStateProposal> BuildAtomProposal(
    const FitState & outer_previous_state,
    const FitStatePatch & endpoint_patch,
    const SuspiciousBlockActivity & block_activity,
    double factor)
{
    const auto & key{ endpoint_patch.atom_index_list };
    std::vector<GaussianModel3D> previous_model_list;
    std::vector<GaussianModel3D> raw_model_list;
    previous_model_list.reserve(key.size());
    raw_model_list.reserve(key.size());
    for (std::size_t atom_position = 0; atom_position < key.size(); atom_position++)
    {
        const auto atom_index{ key.at(atom_position) };
        const auto & previous_model{ outer_previous_state.at(atom_index).mdpde.GetModel() };
        auto raw_model{ endpoint_patch.mdpde_list.at(atom_position).GetModel() };
        if (!block_activity.HasActiveShape(atom_index))
        {
            raw_model = previous_model.WithOffset(raw_model.GetOffset());
        }
        if (!block_activity.HasActiveOffset(atom_index))
        {
            raw_model = raw_model.WithOffset(previous_model.GetOffset());
        }
        previous_model_list.emplace_back(previous_model);
        raw_model_list.emplace_back(raw_model);
    }
    const auto candidate_model_list{
        BuildDampedModelList(previous_model_list, raw_model_list, factor)
    };
    if (!candidate_model_list.has_value())
    {
        return std::nullopt;
    }
    FitStateProposal proposal{
        .patch{ .atom_index_list = key },
        .effective_damping = factor
    };
    proposal.patch.mdpde_list.reserve(key.size());
    for (std::size_t atom_position = 0; atom_position < key.size(); atom_position++)
    {
        proposal.patch.mdpde_list.emplace_back(
            GaussianModel3DWithUncertainty{
                candidate_model_list->at(atom_position),
                endpoint_patch.mdpde_list.at(atom_position).GetStandardDeviationModel()
            });
    }
    const auto step_norm{ CalculateModelTrustRegionStepNorm(previous_model_list, *candidate_model_list) };
    if (!step_norm.has_value()) return std::nullopt;
    proposal.step_norm = *step_norm;
    return proposal;
}

static bool ContainsClusterKey(const std::vector<ClusterKey> & key_list, const ClusterKey & key)
{
    return std::ranges::find(key_list, key) != key_list.end();
}

bool ShouldShrinkAcceptedTrustRegionRadius(
    std::optional<double> first_objective_evaluated_factor,
    const ObjectiveAttemptDiagnostic & diagnostic)
{
    return first_objective_evaluated_factor.has_value() &&
        diagnostic.accepted_factor.has_value() &&
        *diagnostic.accepted_factor < *first_objective_evaluated_factor;
}

static ClusterCandidateResult SelectClusterCandidate(
    const CandidateSelectionInputs & inputs,
    const ClusterKey & key,
    const std::vector<SampleRef> & objective_sample_ref_list,
    ClusterSolverWorkspace & solver_workspace)
{
    const auto & context{ inputs.context };
    const auto & residual_baseline{ inputs.residual_baseline };
    const auto & previous_state{ inputs.previous_state };
    const auto & previous_polish_provenance{ inputs.previous_polish_provenance };
    const auto & ridge_multiplier_list{ inputs.ridge_multiplier_list };
    const auto & objective_domain{ inputs.objective_domain };
    const auto & previous_objective_entry{ inputs.previous_objective_by_key.at(key) };
    const auto * previous_objective{
        previous_objective_entry.has_value() ? &*previous_objective_entry : nullptr
    };
    const auto trust_region_radius{ inputs.trust_region_state.GetRadius(key) };
    auto & performance_counters{ inputs.performance_counters };
    ClusterCandidateResult result;
    result.polish_provenance.reserve(key.size());
    for (const auto atom_index : key)
    {
        result.polish_provenance.emplace_back(previous_polish_provenance.at(atom_index));
    }
    auto search_endpoint_patch{ FitStatePatch::FromState(inputs.proposal_state, key) };
    auto search_block_activity{ inputs.block_activity };
    std::vector<StabilizationTerminalDiagnostic> terminal_diagnostic_list;
    double best_rescue_objective{ std::numeric_limits<double>::infinity() };
    std::optional<double> first_objective_evaluated_factor;
    bool is_polish_eligible{ false };
    TrustModelTrialObserver observer(inputs, key, objective_sample_ref_list);
    for (;;)
    {
        observer.BeginSearch();
        if (context.cluster_history) context.cluster_history->BeginSearch(key);
        result.diagnostic = ObjectiveAttemptDiagnostic{};
        result.diagnostic.trust_region_radius = trust_region_radius;
        result.polish_progress = PolishProgress{};
        is_polish_eligible =
            inputs.health_by_key.at(key).IsSolverQualified() &&
            std::ranges::all_of(
                key,
                [&](const auto atom_index)
                {
                    return search_block_activity.HasActiveShape(atom_index) &&
                        search_block_activity.HasActiveOffset(atom_index);
                });
        if (is_polish_eligible) result.polish_progress.eligible_count = 1;
        const auto has_active_parameter{
            std::ranges::any_of(
                key,
                [&](const auto atom_index)
                {
                    return search_block_activity.HasActiveShape(atom_index) ||
                        (search_block_activity.HasActiveOffset(atom_index));
                })
        };
        if (!has_active_parameter)
        {
            result.accepted_patch = FitStatePatch::FromState(previous_state, key);
            result.diagnostic.previous_objective = previous_objective_entry;
            result.diagnostic.candidate_objective = previous_objective_entry;
            result.diagnostic.trust_region_step_norm = 0.0;
            result.diagnostic.accepted_factor = 0.0;
            if (is_polish_eligible) result.polish_progress.skipped_count = 1;
            break;
        }

        first_objective_evaluated_factor.reset();
        std::optional<StabilizationTerminalDiagnostic> last_guard_failure;
        for (double factor{ 1.0 };
            factor >= std::numeric_limits<double>::epsilon(); factor *= 0.5)
        {
            result.diagnostic.trial_count++;
            observer.Generated();
            auto proposal_result{
                BuildAtomProposal(
                    previous_state,
                    search_endpoint_patch,
                    search_block_activity,
                    factor)
            };
            if (!proposal_result.has_value())
            {
                result.diagnostic.invalid_trial_count++;
                observer.Invalid();
                result.diagnostic.pre_objective_failure_reason =
                    PreObjectiveFailureReason::InvalidModel;
                continue;
            }
            auto proposal{ std::move(*proposal_result) };
            result.diagnostic.pre_objective_attempted_step_norm = proposal.step_norm;
            result.diagnostic.trust_region_step_norm = proposal.step_norm;
            const CandidateEvaluationOverlay candidate_overlay{
                context,
                residual_baseline,
                previous_state,
                proposal.patch
            };
            const auto maximum_change{
                std::ranges::max(
                    SummarizeTransformedChanges(
                        candidate_overlay.GetState(),
                        residual_baseline.model_snapshot.node,
                        key).maximum_list)
            };
            if (maximum_change < kTransformedChangeTolerance)
            {
                observer.Nonmaterial();
                if (factor == 1.0 && is_polish_eligible)
                {
                    result.diagnostic.accepted_factor = 1.0;
                    result.diagnostic.previous_objective = previous_objective_entry;
                    result.diagnostic.candidate_objective = previous_objective_entry;
                    result.accepted_patch = std::move(proposal.patch);
                }
                break;
            }
            const auto preflight{ EvaluateCandidate(candidate_overlay,
                CandidatePreflightReference{key, search_block_activity, proposal.step_norm, trust_region_radius}) };
            if (preflight.failure_stage == CandidateFailureStage::Trust)
            {
                result.diagnostic.trust_skipped_trial_count++;
                observer.TrustSkipped();
                result.diagnostic.pre_objective_failure_reason =
                    PreObjectiveFailureReason::NoCandidateWithinTrustRegion;
                continue;
            }
            if (preflight.guard_failure)
            {
                result.diagnostic.guard_rejected_trial_count++;
                observer.GuardRejected();
                last_guard_failure = preflight.guard_failure;
                continue;
            }
            if (!first_objective_evaluated_factor.has_value())
            {
                first_objective_evaluated_factor = factor;
            }

            ObjectiveAttemptDiagnostic trial_diagnostic{
                .accepted_factor = factor,
                .trust_region_radius = trust_region_radius,
                .trust_region_step_norm = proposal.step_norm,
                .trial_count = result.diagnostic.trial_count,
                .invalid_trial_count = result.diagnostic.invalid_trial_count,
                .trust_skipped_trial_count =
                    result.diagnostic.trust_skipped_trial_count,
                .guard_rejected_trial_count =
                    result.diagnostic.guard_rejected_trial_count,
                .objective_rejected_trial_count =
                    result.diagnostic.objective_rejected_trial_count
            };
            const auto evaluation{ EvaluateCandidate(candidate_overlay, CandidateScope::LocalSearch,
                LocalCandidateReference{key, objective_sample_ref_list, previous_objective,
                    objective_domain, trial_diagnostic, performance_counters}) };
            trial_diagnostic = evaluation.diagnostic;
            const auto committed{ evaluation.accepted };
            if (context.cluster_history)
                trial_diagnostic.history = context.cluster_history->Local(candidate_overlay, key,
                    objective_sample_ref_list, objective_domain, "local-candidate", committed, trial_diagnostic);
            observer.Trial(proposal.patch, trial_diagnostic, false, factor, committed);
            if (committed)
            {
                result.diagnostic = std::move(trial_diagnostic);
                result.accepted_patch = std::move(proposal.patch);
                break;
            }
            trial_diagnostic.objective_rejected_trial_count++;
            if (trial_diagnostic.candidate_objective.has_value())
            {
                const auto rescue_objective{
                    trial_diagnostic.candidate_objective->GetTotalObjective()
                };
                if (std::isfinite(rescue_objective) &&
                    rescue_objective < best_rescue_objective)
                {
                    result.rescue_patch = proposal.patch;
                    best_rescue_objective = rescue_objective;
                }
            }
            result.diagnostic = std::move(trial_diagnostic);
        }
        if (result.accepted_patch.has_value()) break;

        if (result.diagnostic.objective_rejected_trial_count == 0 &&
            result.diagnostic.guard_rejected_trial_count != 0 &&
            last_guard_failure.has_value())
        {
            if (!last_guard_failure->guard_atom_index.has_value() ||
                !last_guard_failure->guard_mode.has_value())
            {
                terminal_diagnostic_list.emplace_back(*last_guard_failure);
                result.diagnostic.terminal_diagnostic_list = std::move(terminal_diagnostic_list);
                ObservePhaseMissing(context, "local-search", key, "guard-failed");
                return result;
            }
            const auto atom_index{ *last_guard_failure->guard_atom_index };
            const auto mode{ *last_guard_failure->guard_mode };
            auto & endpoint_model{ search_endpoint_patch.mdpde_list.at(
                static_cast<std::size_t>(std::ranges::lower_bound(key, atom_index) - key.begin())) };
            if (mode == SuspiciousUpdateMode::OffsetOnly)
            {
                const auto previous_offset{
                    previous_state.at(atom_index).mdpde.GetModel().GetOffset()
                };
                endpoint_model = WithPreservedUncertaintyOffset(
                    endpoint_model, previous_offset);
                search_block_activity.offset_fixed_atom_mask.at(atom_index) = 1;
            }
            else
            {
                const auto retained_offset{
                    endpoint_model.GetModel().GetOffset()
                };
                endpoint_model =
                    WithPreservedUncertaintyOffset(
                        previous_state.at(atom_index).mdpde,
                        retained_offset);
                search_block_activity.shape_fixed_atom_mask.at(atom_index) = 1;
            }
            terminal_diagnostic_list.emplace_back(*last_guard_failure);
            continue;
        }
        if (result.diagnostic.objective_rejected_trial_count != 0)
        {
            terminal_diagnostic_list.emplace_back(
                StabilizationTerminalDiagnostic{
                    StabilizationTerminalReason::ObjectiveExhausted });
        }
        else if (result.diagnostic.invalid_trial_count != 0)
        {
            terminal_diagnostic_list.emplace_back(
                StabilizationTerminalDiagnostic{
                    StabilizationTerminalReason::InvalidCandidate });
        }
        result.diagnostic.terminal_diagnostic_list =
            std::move(terminal_diagnostic_list);
        if (is_polish_eligible) result.polish_progress.skipped_count = 1;
        ObservePhaseMissing(context, "local-search", key, "search-exhausted");
        return result;
    }
    result.diagnostic.terminal_diagnostic_list =
        std::move(terminal_diagnostic_list);

    const FitStateView base_state_view{ previous_state, *result.accepted_patch };
    ObservePhaseCandidate(context, "local-search", key, base_state_view,
        nullptr, result.diagnostic.accepted_factor.value_or(1.0), "accepted");
    for (std::size_t position = 0; position < key.size(); position++)
    {
        if (IsTransformedChangeMaterial(
                CalculateTransformedChange(
                    base_state_view.GetModel(key.at(position)),
                    previous_state.at(key.at(position)).mdpde.GetModel()),
                kTransformedChangeTolerance))
        {
            result.polish_provenance.at(position) = 0;
        }
    }
    result.shrink_trust_region = ShouldShrinkAcceptedTrustRegionRadius(
        first_objective_evaluated_factor,
        result.diagnostic);
    if (is_polish_eligible)
    {
        auto polished_candidate{
            BuildJointPolishProposal(
                context,
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
            ObservePhaseMissing(context, "local-polish", key, "no-polish-proposal");
        }
        else
        {
            ObjectiveAttemptDiagnostic polish_diagnostic;
            polish_diagnostic.accepted_factor = polished_candidate->effective_damping;
            polish_diagnostic.trust_region_radius = trust_region_radius;
            polish_diagnostic.trust_region_step_norm = polished_candidate->step_norm;
            const CandidateEvaluationOverlay polished_overlay{
                context,
                residual_baseline,
                previous_state,
                polished_candidate->patch
            };
            const auto evaluation{ EvaluateCandidate(polished_overlay, CandidateScope::LocalPolish,
                LocalCandidateReference{key, objective_sample_ref_list,
                    result.diagnostic.candidate_objective ? &*result.diagnostic.candidate_objective : nullptr,
                    objective_domain, polish_diagnostic, performance_counters}) };
            polish_diagnostic = evaluation.diagnostic;
            const auto polish_committed{ evaluation.accepted };
            if (context.cluster_history)
                polish_diagnostic.history = context.cluster_history->Local(polished_overlay, key,
                    objective_sample_ref_list, objective_domain, "local-polish", polish_committed, polish_diagnostic);
            ObservePhaseLocalPolish(context, key, polished_overlay.GetState(), base_state_view,
                polished_candidate->effective_damping, polish_committed, polish_diagnostic);
            observer.Trial(polished_candidate->patch, polish_diagnostic, true,
                polished_candidate->effective_damping, polish_committed);
            if (!polish_committed)
            {
                result.polish_progress.rejected_count = 1;
            }
            else
            {
                result.polish_progress.accepted_count = 1;
                for (std::size_t position = 0; position < key.size(); position++)
                {
                    const auto base_coordinates{
                        base_state_view.GetModel(key.at(position)).ToTransformedCoordinates()
                    };
                    const auto candidate_coordinates{
                        polished_candidate->patch.mdpde_list.at(position)
                            .GetModel().ToTransformedCoordinates()
                    };
                    if ((base_coordinates->array() != candidate_coordinates->array()).any())
                    {
                        result.polish_provenance.at(position) = 1;
                    }
                }
                result.accepted_patch = std::move(polished_candidate->patch);

            }
        }
    }
    observer.Finish(result.shrink_trust_region, first_objective_evaluated_factor, result.diagnostic);
    return result;
}

void CandidateTransactionBuilder::Select(const CandidateSelectionInputs & inputs)
{
    const auto candidate_phase_start{ std::chrono::steady_clock::now() };
    const auto & partition{ inputs.partition };
    auto & solver_workspace_by_key{ inputs.solver_workspace_by_key };
    const auto cluster_key_list{ BuildGraphClusterKeyList(partition) };
    // Validate workspace coverage before any parallel solver mutates state.
    for (const auto & key : cluster_key_list)
    {
        static_cast<void>(solver_workspace_by_key.at(key));
    }

    std::vector<ClusterCandidateResult> result_list(cluster_key_list.size());
    std::vector<std::exception_ptr> exception_list(cluster_key_list.size());
    const auto audit_context{ estimation_audit::current };
    const auto select_candidate = [&](std::size_t position)
    {
        try
        {
            estimation_audit::Scope audit_scope(audit_context.attempt, audit_context.source);
            const auto & key{ cluster_key_list.at(position) };
            result_list.at(position) = SelectClusterCandidate(
                inputs,
                key,
                partition.sample_id_list_by_key.at(key),
                solver_workspace_by_key.at(key));
        }
        catch (...)
        {
            exception_list.at(position) = std::current_exception();
        }
    };
#ifdef USE_OPENMP
    if (inputs.options.thread_size > 1 && cluster_key_list.size() > 1)
    {
        eigen_helper::ScopedEigenThreadCount eigen_thread_guard{ 1 };
#pragma omp parallel for schedule(dynamic) num_threads(inputs.options.thread_size)
        for (std::size_t position = 0;
            position < cluster_key_list.size(); position++)
        {
            select_candidate(position);
        }
    }
    else
#endif
    {
        for (std::size_t position = 0;
            position < cluster_key_list.size(); position++)
        {
            select_candidate(position);
        }
    }
    for (const auto & exception : exception_list)
    {
        if (exception) std::rethrow_exception(exception);
    }

    m_selection = CandidateSelection{
        .block_activity = inputs.block_activity,
        .assembled_state = inputs.previous_state,
        .assembled_polish_provenance = inputs.previous_polish_provenance
    };
    auto & selection{ m_selection };
    m_candidate_by_key.clear();
    m_next_rejection_order = 0;
    std::vector<ClusterKey> locally_polished_key_list;
    locally_polished_key_list.reserve(result_list.size());
    for (std::size_t position = 0; position < result_list.size(); position++)
    {
        auto & result{ result_list.at(position) };
        const auto & key{ cluster_key_list.at(position) };
        for (const auto & terminal_diagnostic :
            result.diagnostic.terminal_diagnostic_list)
        {
            if (terminal_diagnostic.reason !=
                    StabilizationTerminalReason::GuardInfeasible ||
                !terminal_diagnostic.guard_atom_index.has_value() ||
                !terminal_diagnostic.guard_mode.has_value())
            {
                continue;
            }
            const auto atom_index{ *terminal_diagnostic.guard_atom_index };
            const auto mode{ *terminal_diagnostic.guard_mode };
            if (mode == SuspiciousUpdateMode::OffsetOnly)
            {
                selection.block_activity.offset_fixed_atom_mask.at(atom_index) = 1;
            }
            else
            {
                selection.block_activity.shape_fixed_atom_mask.at(atom_index) = 1;
            }
        }
        selection.polish_progress.eligible_count +=
            result.polish_progress.eligible_count;
        selection.polish_progress.accepted_count +=
            result.polish_progress.accepted_count;
        selection.polish_progress.rejected_count +=
            result.polish_progress.rejected_count;
        selection.polish_progress.skipped_count +=
            result.polish_progress.skipped_count;
        if (result.polish_progress.accepted_count != 0)
        {
            locally_polished_key_list.emplace_back(key);
        }

        const auto is_accepted{ result.accepted_patch.has_value() };
        auto & pending{ m_candidate_by_key[key] };
        pending.selected = is_accepted;
        if (!is_accepted)
        {
            if (result.rescue_patch.has_value())
            {
                pending.cooperative_patch = std::move(result.rescue_patch);
            }
            else if (IsJointOffsetSolveHardFailure(
                inputs.health_by_key.at(key).joint_offset_status))
            {
                inputs.performance_counters.RecordBoundaryRescueExclusions(
                    1,
                    0,
                    0);
            }
            else if (result.diagnostic.pre_objective_failure_reason !=
                PreObjectiveFailureReason::None)
            {
                inputs.performance_counters.RecordBoundaryRescueExclusions(
                    0,
                    1,
                    0);
            }
            else
            {
                inputs.performance_counters.RecordBoundaryRescueExclusions(
                    0,
                    0,
                    1);
            }
        }

        pending.diagnostic = ClusterCandidateDiagnostic{ key, std::move(result.diagnostic) };
        if (!is_accepted)
        {
            pending.rejection_order = m_next_rejection_order++;
            continue;
        }
        pending.shrink_trust_region = result.shrink_trust_region;
        result.accepted_patch->ApplyTo(selection.assembled_state);
        for (std::size_t key_position = 0;
            key_position < key.size(); key_position++)
        {
            selection.assembled_polish_provenance.at(key.at(key_position)) =
                result.polish_provenance.at(key_position);
        }
    }
    inputs.performance_counters.FinishCandidatePhase(candidate_phase_start);
    inputs.performance_counters.RecordFullStateMaterialization();

    ObservePhaseSearchAssembly(inputs.context, selection.assembled_state);
    ReconcileSelectedBoundaries(inputs);
    ObservePhaseState(inputs.context, "boundary-final", selection.assembled_state);
    for (const auto & key : locally_polished_key_list)
    {
        if (ContainsClusterKey(selection.accepted_key_list, key)) continue;
        selection.polish_progress.accepted_count--;
        selection.polish_progress.rejected_count++;
    }
    FinalizeTrustModelAudit(inputs.context, selection);

}

} // namespace rhbm_gem::core::detail
