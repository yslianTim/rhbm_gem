#include "core/detail/CandidateSelection.hpp"

#include "core/detail/BoundaryReconciliation.hpp"
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
#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>
#include <rhbm_gem/utils/math/EigenHelper.hpp>

namespace rhbm_gem::core::detail {

namespace {

constexpr double kTrustRegionInitialRadius{ 1.0 };
constexpr double kTrustRegionMinimumRadius{ 0.0625 };
constexpr double kTrustRegionMaximumRadius{ 4.0 };
constexpr double kTrustRegionShrinkFactor{ 0.5 };
constexpr double kTrustRegionGrowthFactor{ 2.0 };
constexpr double kTrustRegionGrowthBoundaryRatio{ 0.8 };

bool IsTrustRegionStepAtGrowthBoundary(double step_norm, double radius)
{
    return std::isfinite(step_norm) &&
        std::isfinite(radius) &&
        radius > 0.0 &&
        step_norm >= kTrustRegionGrowthBoundaryRatio * radius;
}

struct ClusterCandidateResult
{
    std::optional<FitStatePatch> accepted_patch{};
    std::optional<FitStatePatch> rescue_patch{};
    PolishProvenance polish_provenance{};
    ClusterObjectiveState objective_state{};
    ObjectiveAttemptDiagnostic diagnostic{};
    PolishProgress polish_progress{};
    TrustRegionRadiusAction radius_action{ TrustRegionRadiusAction::Keep };
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
    std::vector<TrustModelShadowDiagnostic> trust_model_shadow_trial_list{};
    TrustModelCandidateFunnel trust_model_candidate_funnel{};
#endif
};

#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
std::optional<double> EvaluateTrustModelResponseDirection(
    const GaussianModel3D & previous_model,
    const GaussianModel3D & candidate_model,
    double distance)
{
    const auto previous_coordinates{ previous_model.ToTransformedCoordinates() };
    const auto candidate_coordinates{ candidate_model.ToTransformedCoordinates() };
    if (!previous_coordinates.has_value() || !candidate_coordinates.has_value())
    {
        return std::nullopt;
    }
    const auto direction{ *candidate_coordinates - *previous_coordinates };
    if (!direction.allFinite()) return std::nullopt;
    if (direction.isZero()) return 0.0;
    const auto invariants{ BuildTransformedModelInvariants(previous_model) };
    if (!invariants.has_value()) return std::nullopt;
    const auto jacobian{ EvaluateTransformedJacobian(*invariants, distance) };
    if (!jacobian.has_value()) return std::nullopt;
    const auto response_direction{ jacobian->dot(direction) };
    return std::isfinite(response_direction) ?
        std::optional<double>{ response_direction } : std::nullopt;
}
#endif

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
    const std::vector<ClusterKey> & grow_key_list,
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
    for (const auto & key : grow_key_list)
    {
        auto iter{ m_radius_by_key.find(key) };
        if (iter == m_radius_by_key.end())
        {
            throw std::invalid_argument(
                "Local fitting trust-region state is missing.");
        }
        iter->second = std::min(
            kTrustRegionMaximumRadius,
            iter->second * kTrustRegionGrowthFactor);
    }
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

static bool ShouldGrowTrustRegion(
    const ObjectiveAttemptDiagnostic & diagnostic)
{
    return diagnostic.candidate_objective.has_value() &&
        diagnostic.previous_objective.has_value() &&
        IsTrustRegionStepAtGrowthBoundary(
            diagnostic.trust_region_step_norm,
            diagnostic.trust_region_radius) &&
        IsBetterAuditObjective(
            diagnostic.candidate_objective->GetTotalObjective(),
            diagnostic.previous_objective->GetTotalObjective(),
            kObjectiveProgressTolerance);
}

TrustRegionRadiusAction DetermineAcceptedTrustRegionRadiusAction(
    std::optional<double> first_objective_evaluated_factor,
    const ObjectiveAttemptDiagnostic & diagnostic)
{
    if (first_objective_evaluated_factor.has_value() &&
        diagnostic.accepted_factor.has_value() &&
        *diagnostic.accepted_factor < *first_objective_evaluated_factor)
    {
        return TrustRegionRadiusAction::Shrink;
    }
    return ShouldGrowTrustRegion(diagnostic) ?
        TrustRegionRadiusAction::Grow : TrustRegionRadiusAction::Keep;
}

#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
TrustRegionRadiusAction DetermineTrustModelShadowAction(
    const TrustModelShadowDiagnostic & diagnostic)
{
    if (diagnostic.objective_backtracked)
    {
        return TrustRegionRadiusAction::Shrink;
    }
    if (diagnostic.status != TrustModelPredictionStatus::Available ||
        !diagnostic.rho.has_value())
    {
        return diagnostic.current_action;
    }
    if (*diagnostic.rho < 0.25)
    {
        return TrustRegionRadiusAction::Shrink;
    }
    if (*diagnostic.rho > 0.75 &&
        std::isfinite(diagnostic.boundary_utilization) &&
        diagnostic.boundary_utilization >= 0.8)
    {
        return TrustRegionRadiusAction::Grow;
    }
    return TrustRegionRadiusAction::Keep;
}

TrustModelShadowDiagnostic EvaluateTrustModelShadow(
    const SecondStageContext & context,
    const ResidualBaseline & residual_baseline,
    const FitState & previous_state,
    const FitStatePatch & candidate_patch,
    const ClusterKey & key,
    const std::vector<SampleRef> & objective_sample_ref_list,
    const ObjectiveDomain & objective_domain,
    const std::optional<ObjectiveBreakdown> & previous_objective,
    const std::optional<ObjectiveBreakdown> & candidate_objective,
    double trust_region_radius,
    TrustRegionRadiusAction current_action,
    TrustModelCandidateSource candidate_source,
    bool objective_backtracked)
{
    TrustModelShadowDiagnostic result{
        .candidate_source = candidate_source,
        .current_action = current_action,
        .objective_backtracked = objective_backtracked
    };
    if (!previous_objective.has_value() || !candidate_objective.has_value())
    {
        result.status = TrustModelPredictionStatus::ObjectiveUnavailable;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }
    const auto previous_objective_value{ previous_objective->GetTotalObjective() };
    const auto candidate_objective_value{ candidate_objective->GetTotalObjective() };
    if (!std::isfinite(previous_objective_value) ||
        !std::isfinite(candidate_objective_value))
    {
        result.status = TrustModelPredictionStatus::Nonfinite;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }
    result.actual_reduction = previous_objective_value - candidate_objective_value;

    const FitStateView candidate_state{ previous_state, candidate_patch };
    std::vector<GaussianModel3D> previous_model_list;
    std::vector<GaussianModel3D> candidate_model_list;
    previous_model_list.reserve(key.size());
    candidate_model_list.reserve(key.size());
    for (const auto atom_index : key)
    {
        previous_model_list.emplace_back(previous_state.at(atom_index).mdpde.GetModel());
        candidate_model_list.emplace_back(candidate_state.GetModel(atom_index));
    }
    const auto step_norm{
        CalculateModelTrustRegionStepNorm(previous_model_list, candidate_model_list)
    };
    if (!step_norm.has_value())
    {
        result.status = TrustModelPredictionStatus::ModelUnavailable;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }
    result.step_norm = *step_norm;
    result.boundary_utilization =
        std::isfinite(trust_region_radius) && trust_region_radius > 0.0 ?
            *step_norm / trust_region_radius : 0.0;
    if (*step_norm < kTransformedChangeTolerance)
    {
        result.status = TrustModelPredictionStatus::NonmaterialStep;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }
    if (objective_domain.active_atom_count == 0)
    {
        result.status = TrustModelPredictionStatus::ResidualUnavailable;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }

    SecondStageModelSnapshot candidate_snapshot;
    try
    {
        candidate_snapshot = BuildSecondStageModelSnapshot(context, candidate_state);
    }
    catch (const std::exception &)
    {
        result.status = TrustModelPredictionStatus::ModelUnavailable;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }

    double predicted_residual_reduction{ 0.0 };
    for (const auto & sample_ref : objective_sample_ref_list)
    {
        const auto & owner_key{
            objective_domain.owner_key_by_atom_index.at(sample_ref.atom_index)
        };
        if (owner_key.empty()) continue;
        const auto in_fit{
            objective_domain.fit_sample_mask_by_atom.at(sample_ref.atom_index).at(sample_ref.sample_index) != 0
        };
        const auto in_tail{
            objective_domain.tail_sample_mask_by_atom.at(sample_ref.atom_index).at(sample_ref.sample_index) != 0
        };
        if (!in_fit && !in_tail) continue;
        const auto owner_iter{ objective_domain.cluster_by_key.find(owner_key) };
        const auto previous_residual{ residual_baseline(sample_ref) };
        if (owner_iter == objective_domain.cluster_by_key.end() ||
            !owner_iter->second.scale.has_value() ||
            !previous_residual.has_value())
        {
            result.status = TrustModelPredictionStatus::ResidualUnavailable;
            result.shadow_action = DetermineTrustModelShadowAction(result);
            return result;
        }
        const auto & atom_context{ context.atom_list.at(sample_ref.atom_index) };
        const auto target_direction{
            EvaluateTrustModelResponseDirection(
                residual_baseline.model_snapshot.node.at(sample_ref.atom_index),
                candidate_snapshot.node.at(sample_ref.atom_index),
                atom_context.raw_sampling_entries.at(sample_ref.sample_index)
                    .point.distance)
        };
        if (!target_direction.has_value())
        {
            result.status = TrustModelPredictionStatus::ModelUnavailable;
            result.shadow_action = DetermineTrustModelShadowAction(result);
            return result;
        }
        double residual_direction{ -*target_direction };
        for (const auto & neighbor : atom_context.Neighbors(sample_ref.sample_index))
        {
            const auto & previous_neighbor{
                GetFitModel(residual_baseline.model_snapshot.node, neighbor.atom_index)
            };
            const auto & candidate_neighbor{
                GetFitModel(candidate_snapshot.node, neighbor.atom_index)
            };
            const auto neighbor_direction{
                EvaluateTrustModelResponseDirection(
                    previous_neighbor,
                    candidate_neighbor,
                    neighbor.distance)
            };
            if (!neighbor_direction.has_value())
            {
                result.status = TrustModelPredictionStatus::ModelUnavailable;
                result.shadow_action = DetermineTrustModelShadowAction(result);
                return result;
            }
            residual_direction -= *neighbor_direction;
        }
        const auto linearized_residual{
            previous_residual->residual + residual_direction
        };
        for (const bool is_fit_range : { true, false })
        {
            if (!(is_fit_range ? in_fit : in_tail)) continue;
            const auto sample_count{ is_fit_range ?
                owner_iter->second.fit_sample_ref_list.size() :
                owner_iter->second.tail_sample_ref_list.size()
            };
            const auto scale{ is_fit_range ?
                owner_iter->second.scale->fit : owner_iter->second.scale->tail
            };
            if (sample_count == 0 || !std::isfinite(scale) || scale <= 0.0)
            {
                result.status = TrustModelPredictionStatus::ResidualUnavailable;
                result.shadow_action = DetermineTrustModelShadowAction(result);
                return result;
            }

            const auto weight{
                algorithm::CalculateCauchyWeight(
                    previous_residual->residual,
                    scale,
                    kObjectiveRobustLossCutoffMultiplier)
            };
            const auto coefficient{
                CalculateClusterAtomWeight(
                    owner_iter->second.selected_atom_count,
                    objective_domain.active_atom_count) /
                static_cast<double>(sample_count)
            };
            const auto range_weight{ is_fit_range ? kFitRangeWeight : kTailValidationWeight };
            const auto previous_normalized{ previous_residual->residual / scale };
            const auto linearized_normalized{ linearized_residual / scale };
            const auto contribution{
                0.5 * range_weight * coefficient * weight *
                (previous_normalized * previous_normalized -
                    linearized_normalized * linearized_normalized)
            };
            if (!std::isfinite(contribution))
            {
                result.status = TrustModelPredictionStatus::Nonfinite;
                result.shadow_action = DetermineTrustModelShadowAction(result);
                return result;
            }
            predicted_residual_reduction += contribution;
        }
    }
    result.predicted_residual_reduction = predicted_residual_reduction;
    result.predicted_penalty_reduction =
        previous_objective->offset_plausibility_penalty -
        candidate_objective->offset_plausibility_penalty;
    const auto predicted_reduction{
        *result.predicted_residual_reduction +
        *result.predicted_penalty_reduction
    };
    if (!std::isfinite(predicted_reduction) ||
        !result.actual_reduction.has_value() ||
        !std::isfinite(*result.actual_reduction))
    {
        result.status = TrustModelPredictionStatus::Nonfinite;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }
    result.predicted_reduction = predicted_reduction;
    if (predicted_reduction <= 0.0)
    {
        result.status = TrustModelPredictionStatus::NonpositivePrediction;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }
    if (predicted_reduction <= CalculateObjectiveTolerance(
            previous_objective_value,
            kObjectiveProgressTolerance))
    {
        result.status = TrustModelPredictionStatus::NonmaterialPrediction;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }
    const auto rho{ *result.actual_reduction / predicted_reduction };
    if (!std::isfinite(rho))
    {
        result.status = TrustModelPredictionStatus::Nonfinite;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }
    result.status = TrustModelPredictionStatus::Available;
    result.rho = rho;
    result.shadow_action = DetermineTrustModelShadowAction(result);
    return result;
}
#endif

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
    const auto & previous_objective_state{ inputs.cluster_objective_state.at(key) };
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
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
    std::optional<std::size_t> final_trust_model_trial_index;
    std::size_t trust_model_search_pass{ 0 };
    const auto record_trust_model_trial = [&]
        (const FitStatePatch & patch,
         const ObjectiveAttemptDiagnostic & diagnostic,
         TrustModelCandidateSource source,
         std::size_t search_pass,
         std::size_t trial_number,
         double factor,
         bool accepted,
         bool rejected_by_strict_polish)
    {
        const auto start{ std::chrono::steady_clock::now() };
        auto shadow{ EvaluateTrustModelShadow(
            context,
            residual_baseline,
            previous_state,
            patch,
            key,
            objective_sample_ref_list,
            objective_domain,
            previous_objective_entry,
            diagnostic.candidate_objective,
            trust_region_radius,
            TrustRegionRadiusAction::Keep,
            source,
            false) };
        shadow.search_pass = search_pass;
        shadow.trial_number = trial_number;
        shadow.factor = factor;
        shadow.trial_disposition = accepted ?
            TrustModelTrialDisposition::Accepted :
            TrustModelTrialDisposition::ObjectiveRejected;
        shadow.rejected_by_previous = diagnostic.rejected_by_previous;
        shadow.rejected_by_best = diagnostic.rejected_by_best;
        shadow.rejected_by_strict_polish = rejected_by_strict_polish;
        if (source == TrustModelCandidateSource::Polish &&
            diagnostic.previous_objective.has_value() &&
            diagnostic.candidate_objective.has_value())
        {
            shadow.polish_reduction =
                diagnostic.previous_objective->GetTotalObjective() -
                diagnostic.candidate_objective->GetTotalObjective();
        }
        shadow.shadow_action.reset();
        shadow.elapsed_milliseconds =
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - start).count();
        result.trust_model_shadow_trial_list.emplace_back(std::move(shadow));
        return result.trust_model_shadow_trial_list.size() - 1;
    };
#endif
    for (;;)
    {
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
        trust_model_search_pass++;
#endif
        result.objective_state = previous_objective_state;
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
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
            result.trust_model_candidate_funnel.generated_count++;
#endif
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
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
                result.trust_model_candidate_funnel.invalid_count++;
#endif
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
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
                result.trust_model_candidate_funnel.nonmaterial_count++;
#endif
                if (factor == 1.0 && is_polish_eligible)
                {
                    result.diagnostic.accepted_factor = 1.0;
                    result.diagnostic.previous_objective = previous_objective_entry;
                    result.diagnostic.candidate_objective = previous_objective_entry;
                    result.accepted_patch = std::move(proposal.patch);
                }
                break;
            }
            if (!IsTrustRegionStepWithinRadius(
                    proposal.step_norm, trust_region_radius))
            {
                result.diagnostic.trust_skipped_trial_count++;
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
                result.trust_model_candidate_funnel.trust_skipped_count++;
#endif
                result.diagnostic.pre_objective_failure_reason =
                    PreObjectiveFailureReason::NoCandidateWithinTrustRegion;
                continue;
            }
            const auto guard_failure{
                EvaluateClusterCandidateGuard(
                    context,
                    residual_baseline.model_snapshot,
                    key,
                    candidate_overlay.GetState(),
                    search_block_activity)
            };
            if (guard_failure.has_value())
            {
                result.diagnostic.guard_rejected_trial_count++;
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
                result.trust_model_candidate_funnel.guard_rejected_count++;
#endif
                last_guard_failure = guard_failure;
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
            const auto committed{ TryCommitClusterCandidate(
                    candidate_overlay,
                    key,
                    objective_sample_ref_list,
                    previous_objective,
                    false,
                    objective_domain,
                    result.objective_state,
                    trial_diagnostic,
                    performance_counters) };
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
            result.trust_model_candidate_funnel.objective_evaluated_count++;
            const auto trust_model_trial_index{ record_trust_model_trial(
                proposal.patch,
                trial_diagnostic,
                TrustModelCandidateSource::Base,
                trust_model_search_pass,
                trial_diagnostic.trial_count,
                factor,
                committed,
                false) };
#endif
            if (committed)
            {
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
                final_trust_model_trial_index = trust_model_trial_index;
#endif
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
        result.radius_action = TrustRegionRadiusAction::Keep;
        if (is_polish_eligible) result.polish_progress.skipped_count = 1;
        return result;
    }
    result.diagnostic.terminal_diagnostic_list =
        std::move(terminal_diagnostic_list);

    const FitStateView base_state_view{ previous_state, *result.accepted_patch };
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
    result.radius_action = DetermineAcceptedTrustRegionRadiusAction(
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
            const auto polish_committed{ TryCommitClusterCandidate(
                    polished_overlay,
                    key,
                    objective_sample_ref_list,
                    result.diagnostic.candidate_objective.has_value() ?
                        &*result.diagnostic.candidate_objective : nullptr,
                    true,
                    objective_domain,
                    result.objective_state,
                    polish_diagnostic,
                    performance_counters) };
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
            result.trust_model_candidate_funnel.polish_objective_evaluated_count++;
            const auto rejected_by_strict_polish{
                !polish_committed &&
                polish_diagnostic.candidate_objective.has_value() &&
                polish_diagnostic.previous_objective.has_value() &&
                !polish_diagnostic.rejected_by_previous &&
                !polish_diagnostic.rejected_by_best
            };
            const auto trust_model_trial_index{ record_trust_model_trial(
                polished_candidate->patch,
                polish_diagnostic,
                TrustModelCandidateSource::Polish,
                trust_model_search_pass,
                1,
                polished_candidate->effective_damping,
                polish_committed,
                rejected_by_strict_polish) };
#endif
            if (!polish_committed)
            {
                result.polish_progress.rejected_count = 1;
            }
            else
            {
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
                final_trust_model_trial_index = trust_model_trial_index;
#endif
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
                if (result.radius_action != TrustRegionRadiusAction::Shrink &&
                    ShouldGrowTrustRegion(polish_diagnostic))
                {
                    result.radius_action = TrustRegionRadiusAction::Grow;
                }
            }
        }
    }
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
    if (final_trust_model_trial_index.has_value())
    {
        auto & final_shadow{
            result.trust_model_shadow_trial_list.at(*final_trust_model_trial_index)
        };
        final_shadow.final_local_candidate = true;
        final_shadow.readiness_eligible = true;
        final_shadow.current_action = result.radius_action;
        final_shadow.objective_backtracked =
            first_objective_evaluated_factor.has_value() &&
            result.diagnostic.accepted_factor.has_value() &&
            *result.diagnostic.accepted_factor < *first_objective_evaluated_factor;
        final_shadow.shadow_action = DetermineTrustModelShadowAction(final_shadow);
    }
#endif
    return result;
}

#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
static void FinalizeTrustModelShadowDisposition(CandidateSelection & selection)
{
    const auto update = [&](ClusterCandidateDiagnostic & diagnostic, bool accepted)
    {
        diagnostic.boundary_touched = std::ranges::any_of(
            selection.boundary_reconciliation_diagnostic_list,
            [&](const auto & boundary_diagnostic)
            {
                return ContainsClusterKey(
                    boundary_diagnostic.key_list,
                    diagnostic.key);
            });
        for (auto & shadow : diagnostic.trust_model_shadow_trial_list)
        {
            shadow.readiness_eligible =
                shadow.final_local_candidate &&
                accepted &&
                !diagnostic.boundary_touched &&
                !diagnostic.boundary_rescued;
            if (!shadow.readiness_eligible) shadow.shadow_action.reset();
        }
    };
    for (auto & diagnostic : selection.accepted_cluster_diagnostic_list)
    {
        update(diagnostic, true);
    }
    for (auto & diagnostic : selection.rejected_cluster_diagnostic_list)
    {
        update(diagnostic, false);
    }
}
#endif

CandidateSelection SelectClusterCandidates(const CandidateSelectionInputs & inputs)
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
    const auto select_candidate = [&](std::size_t position)
    {
        try
        {
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

    CandidateSelection selection{
        .block_activity = inputs.block_activity,
        .cluster_objective_state = inputs.cluster_objective_state,
        .assembled_state = inputs.previous_state,
        .assembled_polish_provenance = inputs.previous_polish_provenance
    };
    std::map<ClusterKey, FitStatePatch> rescue_patch_by_key;
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
        selection.cluster_objective_state.at(key) = std::move(result.objective_state);

        const auto is_accepted{ result.accepted_patch.has_value() };
        if (!is_accepted)
        {
            if (result.rescue_patch.has_value())
            {
                rescue_patch_by_key.emplace(
                    key,
                    std::move(*result.rescue_patch));
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

        ClusterCandidateDiagnostic cluster_diagnostic{
            key,
            std::move(result.diagnostic)
        };
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
        cluster_diagnostic.trust_model_shadow_trial_list =
            std::move(result.trust_model_shadow_trial_list);
        cluster_diagnostic.trust_model_candidate_funnel =
            result.trust_model_candidate_funnel;
#endif
        if (!is_accepted)
        {
            selection.rejected_key_list.emplace_back(key);
            selection.rejected_cluster_diagnostic_list.emplace_back(
                std::move(cluster_diagnostic));
            continue;
        }

        selection.accepted_key_list.emplace_back(key);
        if (result.radius_action == TrustRegionRadiusAction::Grow)
        {
            selection.grow_trust_region_key_list.emplace_back(key);
        }
        else if (result.radius_action == TrustRegionRadiusAction::Shrink)
        {
            selection.shrink_trust_region_key_list.emplace_back(key);
        }
        selection.accepted_cluster_diagnostic_list.emplace_back(
            std::move(cluster_diagnostic));
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

    ReconcileSelectedBoundaries(inputs, rescue_patch_by_key, selection);
    for (const auto & key : locally_polished_key_list)
    {
        if (ContainsClusterKey(selection.accepted_key_list, key)) continue;
        selection.polish_progress.accepted_count--;
        selection.polish_progress.rejected_count++;
    }
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
    FinalizeTrustModelShadowDisposition(selection);
#endif
    return selection;
}

} // namespace rhbm_gem::core::detail
