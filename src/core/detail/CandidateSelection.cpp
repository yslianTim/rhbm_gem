#include "core/detail/CandidateSelection.hpp"

#include "core/detail/GaussianModelOperations.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iomanip>
#include <limits>
#include <ranges>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <Eigen/Dense>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/EigenHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

namespace rhbm_gem::core::detail {

namespace {

constexpr double kSuspiciousProfileInnermostSignFlipRatio{ 0.25 };
constexpr double kSuspiciousProfileNoiseScaleMultiplier{ 3.0 };
constexpr double kSuspiciousProfileScaleMin{ 1.0e-12 };
constexpr std::size_t kSuspiciousProfileMinimumRadiusCount{ 3 };
constexpr double kSuspiciousProfileDistanceTolerance{ 1.0e-6 };
constexpr double kSuspiciousProfileReboundCenterRatio{ 1.5 };
constexpr double kSuspiciousProfileReboundReferenceRatio{ 0.25 };
constexpr double kSuspiciousProfileUpwardExcursionReferenceRatio{ 0.20 };
constexpr int kSuspiciousProfileMaximumUpwardExcursions{ 1 };
constexpr double kSuspiciousWidthGrowthLimit{ 1.5 };
constexpr double kSuspiciousWidthRangeLimitRatio{ 1.5 };
constexpr double kSuspiciousCompensationResponseRatio{ 2.0 };

constexpr double kObjectiveResidualScaleFloorRatio{ 1.0e-6 };
constexpr double kObjectiveResidualScaleMin{ 1.0e-12 };
constexpr double kFitRangeWeight{ 1.0 };
constexpr double kOffsetPlausibilityPenaltyWeight{ 1.0e-2 };
constexpr double kOffsetPeakRatioMax{ 1.0 };
constexpr double kObjectiveRobustLossCutoffMultiplier{ 1.345 };
constexpr ObjectiveTolerance kObjectiveStrictTolerance{ 1.0e-10, 1.0e-8 };
constexpr ObjectiveTolerance kObjectiveProgressTolerance{ 1.0e-8, 1.0e-3 };
constexpr double kTrustRegionInitialRadius{ 1.0 };
constexpr double kTrustRegionMinimumRadius{ 0.0625 };
constexpr double kTrustRegionMaximumRadius{ 4.0 };
constexpr double kTrustRegionShrinkFactor{ 0.5 };
constexpr double kTrustRegionGrowthFactor{ 2.0 };
constexpr double kTrustRegionGrowthBoundaryRatio{ 0.8 };

using PartitionEntry = std::pair<const ClusterKey, std::vector<SampleRef>>;
using CandidateWork = std::pair<const PartitionEntry *, ClusterSolverWorkspace *>;

enum class SuspiciousProfileAnalysisMode
{
    Candidate,
    PreviousBaseline
};

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

template<typename ResidualEvaluator>
std::optional<ObjectiveBreakdown> EvaluateResidualObjectiveContribution(
    const std::vector<SampleRef> & sample_ref_list,
    const ObjectiveDomain & domain,
    const ResidualEvaluator & residual_evaluator)
{
    if (domain.active_atom_count == 0) return std::nullopt;
    ObjectiveBreakdown contribution;
    for (const auto & sample_ref : sample_ref_list)
    {
        const auto & owner_key{
            domain.owner_key_by_atom_index.at(sample_ref.atom_index)
        };
        if (owner_key.empty()) continue;
        const auto owner_iter{ domain.cluster_by_key.find(owner_key) };
        if (owner_iter == domain.cluster_by_key.end() ||
            !owner_iter->second.scale.has_value())
        {
            return std::nullopt;
        }
        const auto residual_sample{ residual_evaluator(sample_ref) };
        if (!residual_sample.has_value()) return std::nullopt;
        const auto is_fit_range{
            domain.fit_sample_mask_by_atom.at(sample_ref.atom_index).at(sample_ref.sample_index) != 0
        };
        const auto sample_count{
            is_fit_range ? owner_iter->second.fit_sample_ref_list.size() :
                owner_iter->second.tail_sample_ref_list.size()
        };
        if (sample_count == 0) return std::nullopt;
        const auto scale{
            is_fit_range ? owner_iter->second.scale->fit : owner_iter->second.scale->tail
        };
        const auto loss{
            algorithm::CalculateCauchyLoss(
                residual_sample->residual / scale,
                kObjectiveRobustLossCutoffMultiplier)
        };
        const auto coefficient{
            CalculateClusterAtomWeight(owner_key.size(), domain.active_atom_count) /
            static_cast<double>(sample_count)
        };
        if (is_fit_range)
        {
            contribution.fit_range_residual_objective += kFitRangeWeight * coefficient * loss;
        }
        else
        {
            contribution.tail_validation_loss += coefficient * loss;
        }
    }
    if (!std::isfinite(contribution.fit_range_residual_objective) ||
        !std::isfinite(contribution.tail_validation_loss))
    {
        return std::nullopt;
    }
    return contribution;
}

template<typename State>
std::optional<double> EvaluateOffsetPlausibilityPenalty(
    const State & state,
    const ClusterKey & changed_key,
    const ObjectiveDomain & domain)
{
    if (domain.active_atom_count == 0) return std::nullopt;
    double penalty{ 0.0 };
    for (const auto atom_index : changed_key)
    {
        const auto owner_iter{
            domain.cluster_by_key.find(domain.owner_key_by_atom_index.at(atom_index))
        };
        if (owner_iter == domain.cluster_by_key.end() || !owner_iter->second.scale.has_value())
        {
            return std::nullopt;
        }
        const auto & model{ GetFitModel(state, atom_index) };
        if (!IsValidSecondStageGaussianModel(model)) return std::nullopt;
        const auto peak_signal{ model.SignalAtDistance(0.0) };
        const auto offset_peak{ model.GetOffset() * model.OffsetBasisAtDistance(0.0) };
        if (!std::isfinite(peak_signal) || !std::isfinite(offset_peak))
        {
            return std::nullopt;
        }
        const auto offset_ratio{
            std::abs(offset_peak) /
            std::max({
                std::abs(peak_signal),
                owner_iter->second.scale->fit,
                kObjectiveResidualScaleMin
            })
        };
        const auto offset_excess{ std::max(0.0, offset_ratio - kOffsetPeakRatioMax) };
        penalty +=
            kOffsetPlausibilityPenaltyWeight * offset_excess * offset_excess /
            static_cast<double>(domain.active_atom_count);
    }
    return std::isfinite(penalty) ? std::optional<double>{ penalty } : std::nullopt;
}

template<typename Evaluator>
std::optional<ObjectiveBreakdown> EvaluateObjectiveContribution(
    const Evaluator & evaluator,
    const ClusterKey & changed_key,
    const std::vector<SampleRef> & sample_ref_list,
    const ObjectiveDomain & domain)
{
    const auto residual_contribution{
        EvaluateResidualObjectiveContribution(sample_ref_list, domain, evaluator)
    };
    if (!residual_contribution.has_value()) return std::nullopt;
    const auto offset_penalty{
        EvaluateOffsetPlausibilityPenalty(evaluator.GetState(), changed_key, domain)
    };
    if (!offset_penalty.has_value()) return std::nullopt;
    return BuildObjectiveBreakdown(
        residual_contribution->fit_range_residual_objective,
        residual_contribution->tail_validation_loss,
        *offset_penalty);
}

template<typename Evaluator>
std::optional<ObjectiveBreakdown> EvaluateAuditObjectiveImpl(
    const ObjectiveDomain & domain,
    const Evaluator & evaluator)
{
    double fit_range_residual_objective{ 0.0 };
    double tail_validation_loss{ 0.0 };
    double offset_plausibility_penalty{ 0.0 };
    for (const auto & [key, cluster_domain] : domain.cluster_by_key)
    {
        const auto fit_contribution{
            EvaluateResidualObjectiveContribution(
                cluster_domain.fit_sample_ref_list,
                domain,
                evaluator)
        };
        if (!fit_contribution.has_value()) return std::nullopt;
        const auto tail_contribution{
            EvaluateResidualObjectiveContribution(
                cluster_domain.tail_sample_ref_list,
                domain,
                evaluator)
        };
        if (!tail_contribution.has_value()) return std::nullopt;
        const auto offset_contribution{
            EvaluateOffsetPlausibilityPenalty(
                evaluator.GetState(),
                key,
                domain)
        };
        if (!offset_contribution.has_value()) return std::nullopt;
        fit_range_residual_objective += fit_contribution->fit_range_residual_objective;
        fit_range_residual_objective += tail_contribution->fit_range_residual_objective;
        tail_validation_loss += fit_contribution->tail_validation_loss;
        tail_validation_loss += tail_contribution->tail_validation_loss;
        offset_plausibility_penalty += *offset_contribution;
    }
    return BuildObjectiveBreakdown(
        fit_range_residual_objective,
        tail_validation_loss,
        offset_plausibility_penalty);
}

template<typename Evaluator>
ObjectiveByKey BuildObjectiveByKeyImpl(
    const CouplingGraphPartition & partition,
    const ObjectiveDomain & domain,
    const Evaluator & evaluator)
{
    ObjectiveByKey objective_by_key;
    for (const auto & [key, sample_ref_list] :
        partition.sample_id_list_by_key)
    {
        objective_by_key.emplace(
            key,
            EvaluateObjectiveContribution(
                evaluator,
                key,
                sample_ref_list,
                domain));
    }
    return objective_by_key;
}

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
    std::ostringstream message_info;
    message_info
        << " Second-Stage Local Fitting Performance :\n"
        << " - boundary_reconciliation_ms = "
        << std::fixed << std::setprecision(3) << m_boundary_reconciliation_milliseconds << "\n"
        << " - boundary_joint_correction_ms = " << m_boundary_joint_correction_milliseconds << "\n"
        << " - dependency_polish_ms = " << m_dependency_polish_milliseconds << "\n"
        << " - iteration/candidate/topology/total_ms = "
        << std::fixed << std::setprecision(3)
        << m_iteration_phase_milliseconds << "/"
        << m_candidate_phase_milliseconds << "/"
        << m_topology_rebuild_milliseconds << "/"
        << total_milliseconds << "\n";

    std::ostringstream message_debug;
    message_debug
        << " - full_state_materializations = " << m_full_state_materialization_count.load() << "\n"
        << " - gaussian_cache_hit/miss = "
        << m_gaussian_cache_hit_count.load() << "/" << m_gaussian_cache_miss_count.load() << "\n"
        << " - objective_recomputed/reused_samples = "
        << m_objective_recomputed_sample_count.load() << "/"
        << m_objective_reused_sample_count.load() <<"\n"
        << " - solver_symbolic_analyses = " << symbolic_analysis_count << "\n"
        << " - topology_rebuilds/partition_changes = "
        << m_topology_rebuild_attempt_count << "/" << m_topology_partition_change_count<< "\n"
        << " - boundary_reconciliations/backtracked/rejected = "
        << m_boundary_reconciliation_attempt_count << "/"
        << m_boundary_reconciliation_backtracked_count << "/"
        << m_boundary_reconciliation_rejected_count << "\n"
        << " - boundary_joint_correction_attempts/accepted/fallback = "
        << m_boundary_joint_correction_attempt_count << "/"
        << m_boundary_joint_correction_accepted_count << "/"
        << m_boundary_joint_correction_fallback_count << "\n"
        << " - boundary_rescues/accepted/fallback/rejected = "
        << m_boundary_rescue_attempt_count << "/"
        << m_boundary_rescue_accepted_count << "/"
        << m_boundary_rescue_fallback_count << "/"
        << m_boundary_rescue_rejected_count << "\n"
        << " - boundary_rescue_exclusions_hard/invalid/no-objective = "
        << m_boundary_rescue_hard_failure_exclusion_count << "/"
        << m_boundary_rescue_invalid_proposal_exclusion_count << "/"
        << m_boundary_rescue_objective_unavailable_exclusion_count << "\n"
        << " - dependency_polish_components/attempted/accepted/fallback = "
        << m_dependency_polish_component_count << "/"
        << m_dependency_polish_attempt_count << "/"
        << m_dependency_polish_accepted_count << "/"
        << m_dependency_polish_fallback_count << "\n"
        << " - dependency_polish_atoms/parameters/rounds = "
        << m_dependency_polish_atom_count << "/"
        << m_dependency_polish_parameter_count << "/"
        << m_dependency_polish_round_count << "\n";
        
    Logger::Log(LogLevel::Info, message_info.str());
    Logger::Log(LogLevel::Debug, message_debug.str());
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
    std::ostringstream audit_performance;
    audit_performance << std::scientific << std::setprecision(17)
        << "Trust-model performance: schema=1"
        << ", candidate-ms=" << m_candidate_phase_milliseconds
        << ", total-ms=" << total_milliseconds;
    Logger::Log(LogLevel::Debug, audit_performance.str());
#endif
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
    for (const auto & atom_context : context)
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

CandidateEvaluationOverlay::CandidateEvaluationOverlay(
    const SecondStageContext & context,
    const ResidualBaseline & baseline,
    const FitStateView & candidate_state)
    : m_context{ context },
      m_baseline{ baseline },
      m_candidate_state{ candidate_state },
      m_changed_group_mask(context.selected_atom_index_list_by_group.size(), 0),
      m_changed_group_median(context.selected_atom_index_list_by_group.size())
{
    for (const auto atom_index : m_candidate_state.GetOverrideAtomIndexList())
    {
        m_changed_group_mask.at(m_context.at(atom_index).group_id) = 1;
    }
    std::vector<GaussianModel3D> model_list;
    for (std::size_t group_id = 0; group_id < m_changed_group_mask.size(); group_id++)
    {
        if (m_changed_group_mask.at(group_id) == 0) continue;
        const auto & atom_index_list{
            m_context.selected_atom_index_list_by_group.at(group_id)
        };
        model_list.clear();
        model_list.reserve(atom_index_list.size());
        for (const auto atom_index : atom_index_list)
        {
            model_list.emplace_back(m_candidate_state.GetModel(atom_index));
        }
        m_changed_group_median.at(group_id) = BuildGaussianParameterMedian(model_list);
    }
}

std::optional<ResidualSample> CandidateEvaluationOverlay::operator()(const SampleRef & sample_ref) const
{
    const auto & baseline{
        m_baseline.sample_list.at(sample_ref.atom_index).at(sample_ref.sample_index)
    };
    if (!baseline.has_value()) return std::nullopt;
    const auto & atom_context{ m_context.at(sample_ref.atom_index) };
    const auto & sample{
        atom_context.raw_sampling_entries.at(sample_ref.sample_index)
    };
    auto adjusted_response{ baseline->adjusted_response };
    for (const auto & neighbor_atom_sample :
        atom_context.Neighbors(sample_ref.sample_index))
    {
        const GaussianModel3D * candidate_model{ nullptr };
        const GaussianModel3D * baseline_model{ nullptr };
        if (neighbor_atom_sample.is_selected)
        {
            if (m_candidate_state.FindOverride(neighbor_atom_sample.atom_index) == nullptr)
            {
                continue;
            }
            baseline_model = &GetFitModel(
                m_baseline.model_snapshot.selected,
                neighbor_atom_sample.atom_index);
            candidate_model = &m_candidate_state.GetModel(neighbor_atom_sample.atom_index);
        }
        else
        {
            const auto & unselected_atom_contributor{
                m_context.unselected_atom_list.at(neighbor_atom_sample.atom_index)
            };
            if (!unselected_atom_contributor.selected_group_id.has_value() ||
                m_changed_group_mask.at(*unselected_atom_contributor.selected_group_id) == 0)
            {
                continue;
            }
            baseline_model = &GetFitModel(
                m_baseline.model_snapshot.unselected,
                neighbor_atom_sample.atom_index);
            const auto & median{
                m_changed_group_median.at(*unselected_atom_contributor.selected_group_id)
            };
            if (median.has_value())
            {
                candidate_model = &*median;
            }
            else
            {
                candidate_model = &unselected_atom_contributor.initial_seed;
            }
        }
        adjusted_response +=
            baseline_model->ResponseAtDistance(neighbor_atom_sample.distance) -
            candidate_model->ResponseAtDistance(neighbor_atom_sample.distance);
    }
    const auto expected_response{
        m_candidate_state.FindOverride(sample_ref.atom_index) == nullptr ?
            baseline->adjusted_response - baseline->residual :
            m_candidate_state.GetModel(sample_ref.atom_index).ResponseAtDistance(
                sample.point.distance)
    };
    const auto residual{ adjusted_response - expected_response };
    if (!std::isfinite(adjusted_response) || !std::isfinite(residual))
    {
        return std::nullopt;
    }
    return ResidualSample{ adjusted_response, residual };
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
        BuildSharedOffsetDampedModelList(
            m_previous_model_list,
            m_endpoint_model_list,
            m_previous_shared_offset_list,
            m_endpoint_shared_offset_list,
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

SuspiciousUpdateMask SuspiciousBlockActivity::BuildCombinedFixedAtomMask() const
{
    if (shape_fixed_atom_mask.size() != offset_fixed_atom_mask.size() ||
        shape_fixed_atom_mask.size() != hard_failure_atom_mask.size())
    {
        throw std::invalid_argument("Suspicious block activity mask sizes are inconsistent.");
    }
    SuspiciousUpdateMask mask(shape_fixed_atom_mask.size(), 0);
    for (std::size_t atom_index = 0; atom_index < mask.size(); atom_index++)
    {
        mask.at(atom_index) = shape_fixed_atom_mask.at(atom_index) != 0 ||
                offset_fixed_atom_mask.at(atom_index) != 0 ||
                hard_failure_atom_mask.at(atom_index) != 0 ?
            1 : 0;
    }
    return mask;
}

bool SuspiciousBlockActivity::HasActiveShape(std::size_t atom_index) const
{
    return shape_fixed_atom_mask.at(atom_index) == 0;
}

bool SuspiciousBlockActivity::HasActiveOffset(std::size_t atom_index) const
{
    return offset_fixed_atom_mask.at(atom_index) == 0 &&
        hard_failure_atom_mask.at(atom_index) == 0;
}

static bool IsSameSuspiciousProfileRadius(double lhs, double rhs)
{
    const auto scale{ std::max({ std::abs(lhs), std::abs(rhs), 1.0 }) };
    return std::abs(lhs - rhs) <= kSuspiciousProfileDistanceTolerance * scale;
}

static SuspiciousProfileAnalysis BuildSuspiciousProfileAnalysis(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & model,
    const FitOptions & options,
    SuspiciousProfileAnalysisMode mode)
{
    SuspiciousProfileAnalysis analysis;
    std::vector<std::pair<double, double>> profile_samples;
    std::vector<double> residual_list;
    const auto calculate_residual_scale{
        mode == SuspiciousProfileAnalysisMode::PreviousBaseline
    };
    profile_samples.reserve(sample_entries.size());
    if (calculate_residual_scale) residual_list.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        const auto distance{ sample.point.distance };
        const auto response{
            sample.response -
            (model.ResponseAtDistance(distance) - model.SignalAtDistance(distance))
        };
        if (!std::isfinite(response))
        {
            analysis.all_responses_finite = false;
            continue;
        }
        if (distance < options.distance_min || distance > options.distance_max) continue;
        profile_samples.emplace_back(distance, response);
        if (calculate_residual_scale)
        {
            const auto residual{ response - model.SignalAtDistance(distance) };
            if (std::isfinite(residual)) residual_list.emplace_back(residual);
        }
    }

    if (profile_samples.empty()) return analysis;
    std::ranges::sort(
        profile_samples,
        {},
        &std::pair<double, double>::first);

    analysis.distance_range = profile_samples.back().first - profile_samples.front().first;
    for (std::size_t i = 0; i < profile_samples.size();)
    {
        const auto radius{ profile_samples.at(i).first };
        std::vector<double> response_list;
        while (i < profile_samples.size() && IsSameSuspiciousProfileRadius(profile_samples.at(i).first, radius))
        {
            const auto response{ profile_samples.at(i).second };
            analysis.max_abs_response = std::max(
                analysis.max_abs_response,
                std::abs(response));
            response_list.emplace_back(response);
            i++;
        }
        analysis.radius_response_median_list.emplace_back(
            array_helper::ComputeMedian(response_list));
    }
    if (calculate_residual_scale)
    {
        analysis.robust_residual_scale =
            array_helper::ComputeMedianAbsoluteDeviationScale(residual_list);
    }
    return analysis;
}

static bool HasUsableSuspiciousProfileBaseline(
    const GaussianModel3D & previous_model,
    const SuspiciousProfileAnalysis & previous_profile)
{
    if (previous_profile.radius_response_median_list.size() < kSuspiciousProfileMinimumRadiusCount)
    {
        return false;
    }
    if (!std::isfinite(previous_model.GetAmplitude()) ||
        !std::isfinite(previous_model.GetWidth()) ||
        !std::isfinite(previous_model.GetOffset()) ||
        previous_model.GetWidth() <= 0.0 ||
        previous_profile.max_abs_response <= kSuspiciousProfileScaleMin ||
        !std::isfinite(previous_profile.robust_residual_scale))
    {
        return false;
    }
    const auto innermost_scale{
        std::max(
            std::abs(previous_profile.radius_response_median_list.front()),
            kSuspiciousProfileScaleMin)
    };
    for (std::size_t i = 1; i < previous_profile.radius_response_median_list.size(); i++)
    {
        const auto current_scale{
            std::abs(previous_profile.radius_response_median_list.at(i))
        };
        if (current_scale > kSuspiciousProfileReboundCenterRatio * innermost_scale)
        {
            return false;
        }
    }
    return true;
}

static double CalculateNormalizedRatioMargin(double observed, double limit)
{
    if (std::isfinite(observed) && std::isfinite(limit) && limit > 0.0)
    {
        return observed / limit - 1.0;
    }
    return observed > limit ?
        std::numeric_limits<double>::infinity() :
        -std::numeric_limits<double>::infinity();
}

static double CalculateOffsetMagnitudeMargin(
    const GaussianModel3D & previous_model,
    const GaussianModel3D & candidate_model,
    double previous_profile_max_abs_response)
{
    const auto previous_offset_response{
        previous_model.GetOffset() * previous_model.OffsetBasisAtDistance(0.0)
    };
    const auto candidate_offset_response{
        candidate_model.GetOffset() * candidate_model.OffsetBasisAtDistance(0.0)
    };
    if (!std::isfinite(previous_offset_response) ||
        !std::isfinite(candidate_offset_response))
    {
        return std::numeric_limits<double>::infinity();
    }
    const auto reference_scale{
        std::max({
            std::abs(previous_model.SignalAtDistance(0.0)),
            std::abs(previous_offset_response),
            previous_profile_max_abs_response,
            kSuspiciousProfileScaleMin
        })
    };
    return CalculateNormalizedRatioMargin(
        std::abs(candidate_offset_response),
        kSuspiciousCompensationResponseRatio * reference_scale);
}

static double CalculateCenterSignFlipMargin(
    double previous_innermost_response,
    double candidate_innermost_response,
    double previous_residual_scale)
{
    if (!std::isfinite(previous_innermost_response) ||
        !std::isfinite(candidate_innermost_response) ||
        !std::isfinite(previous_residual_scale) || previous_residual_scale < 0.0)
    {
        return -std::numeric_limits<double>::infinity();
    }
    const auto noise_threshold{
        std::max(
            kSuspiciousProfileNoiseScaleMultiplier * previous_residual_scale,
            kSuspiciousProfileScaleMin)
    };
    const auto negative_threshold{
        std::max(
            kSuspiciousProfileInnermostSignFlipRatio * previous_innermost_response,
            noise_threshold)
    };
    return std::min(
        CalculateNormalizedRatioMargin(previous_innermost_response, noise_threshold),
        CalculateNormalizedRatioMargin(-candidate_innermost_response, negative_threshold));
}

static double CalculateRadialReboundMargin(
    const SuspiciousProfileAnalysis & previous_profile,
    const SuspiciousProfileAnalysis & candidate_profile)
{
    if (candidate_profile.radius_response_median_list.size() < kSuspiciousProfileMinimumRadiusCount)
    {
        return -std::numeric_limits<double>::infinity();
    }
    const auto noise_threshold{
        std::max(
            kSuspiciousProfileNoiseScaleMultiplier * previous_profile.robust_residual_scale,
            kSuspiciousProfileScaleMin)
    };
    const auto rebound_threshold{
        std::max(
            kSuspiciousProfileReboundReferenceRatio *
                std::abs(previous_profile.radius_response_median_list.front()),
            noise_threshold)
    };
    const auto excursion_threshold{
        std::max(
            kSuspiciousProfileUpwardExcursionReferenceRatio *
                std::abs(previous_profile.radius_response_median_list.front()),
            noise_threshold)
    };
    const auto candidate_center_scale{
        std::max(
            std::abs(candidate_profile.radius_response_median_list.front()),
            kSuspiciousProfileScaleMin)
    };
    double margin{ -std::numeric_limits<double>::infinity() };
    int excursion_count{ 0 };
    auto previous_abs_response{ std::abs(candidate_profile.radius_response_median_list.front()) };
    for (std::size_t i = 1; i < candidate_profile.radius_response_median_list.size(); i++)
    {
        const auto current_abs_response{
            std::abs(candidate_profile.radius_response_median_list.at(i))
        };
        margin = std::max(
            margin,
            std::min(
                CalculateNormalizedRatioMargin(
                    current_abs_response,
                    kSuspiciousProfileReboundCenterRatio * candidate_center_scale),
                CalculateNormalizedRatioMargin(current_abs_response, rebound_threshold)));
        if (current_abs_response > previous_abs_response + excursion_threshold)
        {
            excursion_count++;
        }
        previous_abs_response = current_abs_response;
    }
    margin = std::max(
        margin,
        static_cast<double>(excursion_count) /
                static_cast<double>(kSuspiciousProfileMaximumUpwardExcursions) -
            1.0);
    return margin;
}

static double CalculateWidthGrowthMargin(
    const GaussianModel3D & previous_model,
    const GaussianModel3D & candidate_model,
    const SuspiciousProfileAnalysis & previous_profile)
{
    const auto width{ candidate_model.GetWidth() };
    if (!std::isfinite(width) || width <= 0.0)
    {
        return std::numeric_limits<double>::infinity();
    }
    auto margin{ CalculateNormalizedRatioMargin(
        width,
        kSuspiciousWidthGrowthLimit * previous_model.GetWidth()) };
    if (!previous_profile.radius_response_median_list.empty() &&
        previous_profile.distance_range > 0.0)
    {
        margin = std::max(
            margin,
            CalculateNormalizedRatioMargin(
                width,
                kSuspiciousWidthRangeLimitRatio * previous_profile.distance_range));
    }
    return margin;
}

static double CalculateAmplitudeOffsetCompensationMargin(
    const GaussianModel3D & previous_model,
    const GaussianModel3D & candidate_model,
    const SuspiciousProfileAnalysis & previous_profile)
{
    const auto signal_delta{
        candidate_model.SignalAtDistance(0.0) - previous_model.SignalAtDistance(0.0)
    };
    const auto offset_delta_response{
        candidate_model.GetOffset() * candidate_model.OffsetBasisAtDistance(0.0) -
            previous_model.GetOffset() * previous_model.OffsetBasisAtDistance(0.0)
    };
    if (!(signal_delta * offset_delta_response < 0.0)) return -1.0;
    const auto reference_scale{
        std::max({
            previous_profile.radius_response_median_list.empty() ?
                0.0 :
                std::abs(previous_profile.radius_response_median_list.front()),
            std::abs(previous_model.SignalAtDistance(0.0)),
            kSuspiciousProfileScaleMin
        })
    };
    const auto limit{ kSuspiciousCompensationResponseRatio * reference_scale };
    return std::min(
        CalculateNormalizedRatioMargin(std::abs(signal_delta), limit),
        CalculateNormalizedRatioMargin(std::abs(offset_delta_response), limit));
}

SuspiciousGaussianAssessment AssessSuspiciousGaussianUpdate(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & candidate_model,
    const FitOptions & options,
    const SuspiciousUpdateBaseline & previous_baseline,
    SuspiciousUpdateMode mode)
{
    const auto & previous_model{ previous_baseline.previous_model };
    const auto & previous_analysis{ previous_baseline.previous_analysis };
    SuspiciousGaussianAssessment assessment;
    if (mode == SuspiciousUpdateMode::PostRefit && !IsValidSecondStageGaussianModel(candidate_model))
    {
        assessment.reason = SuspiciousGaussianReason::InvalidModel;
        assessment.normalized_margin = std::numeric_limits<double>::infinity();
        return assessment;
    }
    const auto candidate_analysis{
        BuildSuspiciousProfileAnalysis(
            sample_entries,
            candidate_model,
            options,
            SuspiciousProfileAnalysisMode::Candidate)
    };
    if (!candidate_analysis.all_responses_finite)
    {
        assessment.reason = SuspiciousGaussianReason::NonFiniteResponse;
        assessment.normalized_margin = std::numeric_limits<double>::infinity();
        return assessment;
    }
    const auto offset_margin{ CalculateOffsetMagnitudeMargin(
            previous_model,
            candidate_model,
            previous_analysis.max_abs_response) };
    assessment.normalized_margin = offset_margin;
    if (offset_margin > 0.0)
    {
        assessment.reason = SuspiciousGaussianReason::OffsetMagnitude;
        return assessment;
    }
    const auto has_usable_radial_baseline{
        previous_analysis.all_responses_finite &&
        !previous_analysis.radius_response_median_list.empty() &&
        HasUsableSuspiciousProfileBaseline(previous_model, previous_analysis)
    };
    if (has_usable_radial_baseline)
    {
        if (candidate_analysis.radius_response_median_list.empty())
        {
            assessment.reason = SuspiciousGaussianReason::NonFiniteResponse;
            assessment.normalized_margin = std::numeric_limits<double>::infinity();
            return assessment;
        }
        const auto sign_flip_margin{
            CalculateCenterSignFlipMargin(
                previous_analysis.radius_response_median_list.front(),
                candidate_analysis.radius_response_median_list.front(),
                previous_analysis.robust_residual_scale)
        };
        assessment.normalized_margin = std::max(
            assessment.normalized_margin,
            sign_flip_margin);
        if (sign_flip_margin > 0.0)
        {
            assessment.reason = SuspiciousGaussianReason::CenterSignFlip;
            assessment.normalized_margin = sign_flip_margin;
            return assessment;
        }
        const auto rebound_margin{
            CalculateRadialReboundMargin(
                previous_analysis,
                candidate_analysis)
        };
        assessment.normalized_margin = std::max(
            assessment.normalized_margin,
            rebound_margin);
        if (rebound_margin > 0.0)
        {
            assessment.reason = SuspiciousGaussianReason::RadialRebound;
            assessment.normalized_margin = rebound_margin;
            return assessment;
        }
    }
    if (mode == SuspiciousUpdateMode::OffsetOnly)
    {
        return assessment;
    }
    const auto width_margin{
        CalculateWidthGrowthMargin(previous_model, candidate_model, previous_analysis)
    };
    assessment.normalized_margin = std::max(assessment.normalized_margin, width_margin);
    if (width_margin > 0.0)
    {
        assessment.reason = SuspiciousGaussianReason::WidthGrowth;
        assessment.normalized_margin = width_margin;
        return assessment;
    }
    const auto compensation_margin{
        CalculateAmplitudeOffsetCompensationMargin(
            previous_model,
            candidate_model,
            previous_analysis)
    };
    assessment.normalized_margin = std::max(
        assessment.normalized_margin,
        compensation_margin);
    if (compensation_margin > 0.0)
    {
        assessment.reason = SuspiciousGaussianReason::AmplitudeOffsetCompensation;
        assessment.normalized_margin = compensation_margin;
        return assessment;
    }
    return assessment;
}

SuspiciousUpdateBaseline BuildPreviousSuspiciousProfileBaseline(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & previous_model,
    const FitOptions & options)
{
    return SuspiciousUpdateBaseline{
        previous_model,
        BuildSuspiciousProfileAnalysis(
            sample_entries,
            previous_model,
            options,
            SuspiciousProfileAnalysisMode::PreviousBaseline)
    };
}

double CalculateClusterAtomWeight(std::size_t cluster_atom_count, std::size_t active_atom_count)
{
    if (cluster_atom_count == 0 || active_atom_count == 0 || cluster_atom_count > active_atom_count)
    {
        throw std::invalid_argument("Local fitting cluster atom counts are invalid.");
    }
    return static_cast<double>(cluster_atom_count) / static_cast<double>(active_atom_count);
}

static void ValidateObjectiveTolerance(ObjectiveTolerance tolerance)
{
    if (!std::isfinite(tolerance.absolute_tolerance) || tolerance.absolute_tolerance < 0.0 ||
        !std::isfinite(tolerance.relative_tolerance) || tolerance.relative_tolerance < 0.0)
    {
        throw std::invalid_argument(
            "Local fitting audit objective tolerances must be finite and non-negative.");
    }
}

static double CalculateObjectiveTolerance(
    double reference,
    ObjectiveTolerance tolerance)
{
    return tolerance.absolute_tolerance + tolerance.relative_tolerance * std::abs(reference);
}

static bool IsObjectiveDeteriorated(
    double candidate,
    double reference,
    ObjectiveTolerance tolerance)
{
    if (!std::isfinite(reference)) return true;
    return candidate > reference + CalculateObjectiveTolerance(reference, tolerance);
}

std::optional<ObjectiveBreakdown> BuildObjectiveBreakdown(
    double fit_range_residual_objective,
    double tail_validation_loss,
    double offset_plausibility_penalty)
{
    if (!std::isfinite(fit_range_residual_objective) ||
        !std::isfinite(tail_validation_loss) ||
        !std::isfinite(offset_plausibility_penalty))
    {
        return std::nullopt;
    }

    const ObjectiveBreakdown breakdown{
        fit_range_residual_objective,
        tail_validation_loss,
        offset_plausibility_penalty
    };
    if (!std::isfinite(breakdown.GetTotalObjective())) return std::nullopt;
    return breakdown;
}

bool IsBetterAuditObjective(double candidate, double best, ObjectiveTolerance tolerance)
{
    ValidateObjectiveTolerance(tolerance);
    if (!std::isfinite(candidate)) return false;
    if (!std::isfinite(best)) return true;
    return candidate < best - CalculateObjectiveTolerance(best, tolerance);
}

bool IsAuditObjectiveAcceptableForProgress(
    double candidate,
    double previous,
    const ObjectiveBreakdown * best,
    ObjectiveTolerance tolerance)
{
    ValidateObjectiveTolerance(tolerance);
    if (!std::isfinite(candidate) || !std::isfinite(previous)) return false;
    return !IsObjectiveDeteriorated(candidate, previous, tolerance) &&
        (best == nullptr ||
            !IsObjectiveDeteriorated(candidate, best->GetTotalObjective(), tolerance));
}

static void AppendObjectiveScaleSummary(std::ostringstream & message, const std::vector<double> & scale_list)
{
    if (scale_list.empty())
    {
        message << "unavailable";
        return;
    }
    message
        << array_helper::ComputePercentile(scale_list, 0.5) << "/"
        << array_helper::ComputePercentile(scale_list, 0.99) << "/"
        << std::ranges::max(scale_list);
}

void LogObjectiveDomain(
    const ObjectiveDomain & domain,
    bool quiet_mode,
    bool is_terminal_reset)
{
    if (quiet_mode) return;
    std::vector<double> fit_scale_list;
    std::vector<double> tail_scale_list;
    for (const auto & entry : domain.cluster_by_key)
    {
        const auto & cluster_domain{ entry.second };
        if (!cluster_domain.scale.has_value()) continue;
        fit_scale_list.emplace_back(cluster_domain.scale->fit);
        if (!cluster_domain.tail_sample_ref_list.empty())
        {
            tail_scale_list.emplace_back(cluster_domain.scale->tail);
        }
    }
    std::ostringstream message;
    message
        << (is_terminal_reset ?
            "Reset second-stage objective domain" : "Initialize second-stage objective domain")
        << ": fit/tail/offset weights = "
        << kFitRangeWeight << "/" << kTailValidationWeight << "/" << kOffsetPlausibilityPenaltyWeight
        << ", clusters = " << domain.cluster_by_key.size()
        << ", active atoms = " << domain.active_atom_count
        << ", unique fit/tail samples = " << domain.fit_sample_count << "/" << domain.tail_sample_count
        << ", fixed fit scale median/p99/max = ";
    AppendObjectiveScaleSummary(message, fit_scale_list);
    message << ", fixed tail scale median/p99/max = ";
    AppendObjectiveScaleSummary(message, tail_scale_list);
    message << ".";
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info, message.str());
}

static std::optional<double> BuildFixedObjectiveScale(
    const std::vector<double> & residual_list,
    const std::vector<double> & adjusted_response_list)
{
    if (residual_list.empty() || residual_list.size() != adjusted_response_list.size())
    {
        return std::nullopt;
    }
    const auto scale{
        std::max({
            array_helper::ComputeMedianAbsoluteDeviationScale(residual_list),
            kObjectiveResidualScaleFloorRatio *
                array_helper::ComputeMedianAbsoluteDeviationScale(adjusted_response_list),
            kObjectiveResidualScaleMin
        })
    };
    return numeric_validation::IsFinitePositive(scale) ? std::optional<double>{ scale } : std::nullopt;
}

ObjectiveDomain BuildObjectiveDomain(
    const SecondStageContext & context,
    const SecondStageModelSnapshot & model_snapshot,
    const std::vector<ClusterKey> & cluster_key_list,
    double distance_min,
    double distance_max)
{
    ObjectiveDomain domain;
    domain.owner_key_by_atom_index.resize(context.size());
    domain.fit_sample_mask_by_atom.resize(context.size());
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        domain.fit_sample_mask_by_atom.at(atom_index).resize(
            context.at(atom_index).raw_sampling_entries.size(), 0);
    }
    for (const auto & key : cluster_key_list)
    {
        auto & cluster_domain{ domain.cluster_by_key[key] };
        domain.active_atom_count += key.size();
        std::vector<double> fit_residual_list;
        std::vector<double> fit_response_list;
        std::vector<double> tail_residual_list;
        std::vector<double> tail_response_list;
        for (const auto atom_index : key)
        {
            domain.owner_key_by_atom_index.at(atom_index) = key;
            const auto & raw_sampling_entries{ context.at(atom_index).raw_sampling_entries };
            for (std::size_t sample_index = 0; sample_index < raw_sampling_entries.size(); sample_index++)
            {
                const SampleRef sample_ref{ atom_index, sample_index };
                const auto residual_sample{
                    EvaluateResidualSample(context, sample_ref, model_snapshot)
                };
                const auto distance{
                    raw_sampling_entries.at(sample_index).point.distance
                };
                const auto is_fit_range{ distance >= distance_min && distance <= distance_max };
                domain.fit_sample_mask_by_atom.at(atom_index).at(sample_index) =
                    is_fit_range ? 1 : 0;
                auto & sample_ref_list{ is_fit_range ?
                    cluster_domain.fit_sample_ref_list : cluster_domain.tail_sample_ref_list
                };
                sample_ref_list.emplace_back(sample_ref);
                if (!residual_sample.has_value()) continue;
                auto & residual_list{
                    is_fit_range ? fit_residual_list : tail_residual_list
                };
                auto & response_list{
                    is_fit_range ? fit_response_list : tail_response_list
                };
                residual_list.emplace_back(residual_sample->residual);
                response_list.emplace_back(residual_sample->adjusted_response);
            }
        }
        domain.fit_sample_count += cluster_domain.fit_sample_ref_list.size();
        domain.tail_sample_count += cluster_domain.tail_sample_ref_list.size();
        const auto fit_scale{
            BuildFixedObjectiveScale(fit_residual_list, fit_response_list)
        };
        if (!fit_scale.has_value() ||
            fit_residual_list.size() != cluster_domain.fit_sample_ref_list.size())
        {
            continue;
        }
        ObjectiveScale scale;
        scale.fit = *fit_scale;
        if (!cluster_domain.tail_sample_ref_list.empty())
        {
            const auto tail_scale{
                BuildFixedObjectiveScale(tail_residual_list, tail_response_list)
            };
            if (!tail_scale.has_value() ||
                tail_residual_list.size() != cluster_domain.tail_sample_ref_list.size())
            {
                continue;
            }
            scale.tail = *tail_scale;
        }
        cluster_domain.scale = scale;
    }
    return domain;
}

std::optional<ObjectiveBreakdown> EvaluateAuditObjective(
    const ObjectiveDomain & domain,
    const ResidualBaseline & evaluator)
{
    return EvaluateAuditObjectiveImpl(domain, evaluator);
}

std::optional<ObjectiveBreakdown> EvaluateAuditObjective(
    const ObjectiveDomain & domain,
    const SnapshotResidualEvaluator & evaluator)
{
    return EvaluateAuditObjectiveImpl(domain, evaluator);
}

ObjectiveByKey BuildObjectiveByKey(
    const CouplingGraphPartition & partition,
    const ObjectiveDomain & domain,
    const ResidualBaseline & evaluator)
{
    return BuildObjectiveByKeyImpl(partition, domain, evaluator);
}

ObjectiveByKey BuildObjectiveByKey(
    const CouplingGraphPartition & partition,
    const ObjectiveDomain & domain,
    const SnapshotResidualEvaluator & evaluator)
{
    return BuildObjectiveByKeyImpl(partition, domain, evaluator);
}

BacktrackingWorkspace::BacktrackingWorkspace(
    const SecondStageContext & context,
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
    std::vector<std::size_t> group_id_by_atom_position;
    group_id_by_atom_position.reserve(m_candidate_patch.atom_index_list.size());
    m_previous_model_list.reserve(m_candidate_patch.atom_index_list.size());
    m_endpoint_model_list.reserve(m_candidate_patch.atom_index_list.size());
    for (std::size_t atom_position = 0;
        atom_position < m_candidate_patch.atom_index_list.size();
        atom_position++)
    {
        const auto atom_index{ m_candidate_patch.atom_index_list.at(atom_position) };
        group_id_by_atom_position.emplace_back(context.at(atom_index).group_id);
        m_previous_model_list.emplace_back(previous_state.at(atom_index).mdpde.GetModel());
        const auto & endpoint_mdpde{
            m_candidate_patch.mdpde_list.at(atom_position)
        };
        m_endpoint_model_list.emplace_back(endpoint_mdpde.GetModel());
    }
    m_previous_shared_offset_list = BuildGroupMedianOffsetList(
        group_id_by_atom_position,
        m_previous_model_list);
    m_endpoint_shared_offset_list = BuildGroupMedianOffsetList(
        group_id_by_atom_position,
        m_endpoint_model_list);
}

static std::optional<ObjectiveBreakdown> EvaluateObjectiveDelta(
    const CandidateEvaluationOverlay & candidate_overlay,
    const std::vector<SampleRef> & affected_sample_ref_list,
    const ObjectiveDomain & domain,
    const ObjectiveBreakdown & baseline,
    PerformanceCounters & performance_counters)
{
    const auto & changed_key{
        candidate_overlay.GetState().GetOverrideAtomIndexList()
    };
    const auto unique_sample_count{ domain.fit_sample_count + domain.tail_sample_count };
    performance_counters.RecordObjectiveSampleEvaluation(
        affected_sample_ref_list.size(),
        unique_sample_count);
    const auto candidate_changed{
        EvaluateObjectiveContribution(
            candidate_overlay,
            changed_key,
            affected_sample_ref_list,
            domain)
    };
    const auto previous_changed{
        EvaluateObjectiveContribution(
            candidate_overlay.GetBaseline(),
            changed_key,
            affected_sample_ref_list,
            domain)
    };
    if (!candidate_changed.has_value() || !previous_changed.has_value()) return std::nullopt;
    return BuildObjectiveBreakdown(
        baseline.fit_range_residual_objective +
            candidate_changed->fit_range_residual_objective -
            previous_changed->fit_range_residual_objective,
        baseline.tail_validation_loss +
            candidate_changed->tail_validation_loss -
            previous_changed->tail_validation_loss,
        baseline.offset_plausibility_penalty +
            candidate_changed->offset_plausibility_penalty -
            previous_changed->offset_plausibility_penalty);
}

static std::optional<ObjectiveBreakdown> EvaluateCombinedObjective(
    const CandidateEvaluationOverlay & candidate_overlay,
    const std::vector<SampleRef> & affected_sample_ref_list,
    const ObjectiveDomain & domain,
    const ObjectiveBreakdown * best_objective,
    const ObjectiveBreakdown * previous_objective,
    PerformanceCounters & performance_counters)
{
    if (previous_objective == nullptr) return std::nullopt;
    const auto candidate_objective{
        EvaluateObjectiveDelta(
            candidate_overlay,
            affected_sample_ref_list,
            domain,
            *previous_objective,
            performance_counters)
    };
    if (!candidate_objective.has_value() ||
        !IsAuditObjectiveAcceptableForProgress(
            candidate_objective->GetTotalObjective(),
            previous_objective->GetTotalObjective(),
            best_objective,
            kObjectiveProgressTolerance))
    {
        return std::nullopt;
    }
    return candidate_objective;
}

bool TryUpdateBestAuditState(
    const FitState & candidate_state,
    bool candidate_uses_polish,
    std::size_t source_iteration,
    const ObjectiveBreakdown & candidate_objective,
    BestAuditState & audit_state)
{
    if (audit_state.has_value() &&
        !IsBetterAuditObjective(
            candidate_objective.GetTotalObjective(),
            audit_state->objective.GetTotalObjective(),
            kObjectiveStrictTolerance))
    {
        return false;
    }
    audit_state = AuditedState{
        candidate_objective,
        candidate_state,
        candidate_uses_polish,
        source_iteration
    };
    return true;
}

void ReconcileClusterObjectiveState(
    const ObjectiveByKey & previous_objective_by_key,
    ClusterObjectiveStateMap & state_by_key)
{
    ClusterObjectiveStateMap next_state_by_key;
    for (const auto & [key, previous_objective] : previous_objective_by_key)
    {
        auto state_iter{ state_by_key.find(key) };
        if (state_iter != state_by_key.end())
        {
            next_state_by_key.emplace(key, std::move(state_iter->second));
            continue;
        }
        next_state_by_key.emplace(
            key,
            ClusterObjectiveState{ .best_objective = previous_objective });
    }
    state_by_key = std::move(next_state_by_key);
}

static bool TryCommitClusterCandidate(
    const CandidateEvaluationOverlay & candidate_overlay,
    const ClusterKey & key,
    const std::vector<SampleRef> & objective_sample_ref_list,
    const ObjectiveBreakdown * previous_objective,
    bool requires_strict_improvement,
    const ObjectiveDomain & domain,
    ClusterObjectiveState & objective_state,
    ObjectiveAttemptDiagnostic & diagnostic,
    PerformanceCounters & performance_counters)
{
    const auto unique_sample_count{
        domain.fit_sample_count + domain.tail_sample_count
    };
    performance_counters.RecordObjectiveSampleEvaluation(
        objective_sample_ref_list.size(),
        unique_sample_count);
    const auto transformed_change_summary{
        SummarizeTransformedChanges(
            candidate_overlay.GetState(),
            candidate_overlay.GetBaseline().model_snapshot.selected,
            key)
    };
    const auto maximum_transformed_change{ std::ranges::max(transformed_change_summary.maximum_list) };
    const auto domain_iter{ domain.cluster_by_key.find(key) };
    diagnostic.scale.reset();
    if (domain_iter != domain.cluster_by_key.end())
    {
        diagnostic.fit_sample_count = domain_iter->second.fit_sample_ref_list.size();
        diagnostic.tail_sample_count = domain_iter->second.tail_sample_ref_list.size();
        diagnostic.scale = domain_iter->second.scale;
    }
    diagnostic.candidate_objective =
        EvaluateObjectiveContribution(
            candidate_overlay,
            key,
            objective_sample_ref_list,
            domain);
    diagnostic.previous_objective.reset();
    if (previous_objective != nullptr)
    {
        diagnostic.previous_objective = *previous_objective;
    }
    diagnostic.best_objective = objective_state.best_objective;

    if (!diagnostic.candidate_objective.has_value() || previous_objective == nullptr)
    {
        return false;
    }
    const auto candidate_objective_value{ diagnostic.candidate_objective->GetTotalObjective() };
    const auto previous_objective_value{ previous_objective->GetTotalObjective() };
    diagnostic.rejected_by_previous = IsObjectiveDeteriorated(
        candidate_objective_value,
        previous_objective_value,
        kObjectiveProgressTolerance);
    diagnostic.rejected_by_best = objective_state.best_objective.has_value() &&
        IsObjectiveDeteriorated(
            candidate_objective_value,
            objective_state.best_objective->GetTotalObjective(),
            kObjectiveProgressTolerance);
    if (diagnostic.rejected_by_previous || diagnostic.rejected_by_best) return false;
    if (requires_strict_improvement &&
        !IsBetterAuditObjective(
            candidate_objective_value,
            previous_objective_value,
            kObjectiveStrictTolerance))
    {
        return false;
    }

    auto is_better_than_best{ !objective_state.best_objective.has_value() };
    if (objective_state.best_objective.has_value())
    {
        const auto best_objective_value{ objective_state.best_objective->GetTotalObjective() };
        if (IsBetterAuditObjective(
                candidate_objective_value,
                best_objective_value,
                kObjectiveStrictTolerance))
        {
            is_better_than_best = true;
        }
        else if (IsBetterAuditObjective(
                     best_objective_value,
                     candidate_objective_value,
                     kObjectiveStrictTolerance))
        {
            is_better_than_best = false;
        }
        else
        {
            is_better_than_best = maximum_transformed_change < objective_state.best_maximum_transformed_change;
        }
    }
    if (is_better_than_best)
    {
        objective_state.best_objective = diagnostic.candidate_objective;
        objective_state.best_maximum_transformed_change = maximum_transformed_change;
    }
    return true;
}

static std::optional<FitStateProposal> BuildSharedOffsetProposal(
    const SecondStageContext & context,
    const FitState & outer_previous_state,
    const FitState & operator_proposal_state,
    const SuspiciousBlockActivity & block_activity,
    const ClusterKey & key,
    double factor)
{
    if (key.empty())
    {
        return std::nullopt;
    }

    std::vector<std::size_t> group_id_by_atom_position;
    std::vector<GaussianModel3D> previous_model_list;
    std::vector<GaussianModel3D> raw_model_list;
    group_id_by_atom_position.reserve(key.size());
    previous_model_list.reserve(key.size());
    raw_model_list.reserve(key.size());
    std::set<std::size_t> inactive_offset_group_id_set;
    for (const auto atom_index : key)
    {
        if (!block_activity.HasActiveOffset(atom_index))
        {
            inactive_offset_group_id_set.emplace(context.at(atom_index).group_id);
        }
    }
    for (const auto atom_index : key)
    {
        group_id_by_atom_position.emplace_back(context.at(atom_index).group_id);
        const auto & previous_model{
            outer_previous_state.at(atom_index).mdpde.GetModel()
        };
        auto operator_model{
            operator_proposal_state.at(atom_index).mdpde.GetModel()
        };
        if (!block_activity.HasActiveShape(atom_index))
        {
            operator_model = previous_model.WithOffset(operator_model.GetOffset());
        }
        if (inactive_offset_group_id_set.contains(context.at(atom_index).group_id))
        {
            operator_model = operator_model.WithOffset(previous_model.GetOffset());
        }
        previous_model_list.emplace_back(previous_model);
        raw_model_list.emplace_back(operator_model);
    }
    const auto previous_shared_offset_list{
        BuildGroupMedianOffsetList(group_id_by_atom_position, previous_model_list)
    };
    const auto raw_shared_offset_list{
        BuildGroupMedianOffsetList(group_id_by_atom_position, raw_model_list)
    };

    const auto seed_model_list{
        BuildSharedOffsetDampedModelList(
            previous_model_list,
            raw_model_list,
            previous_shared_offset_list,
            raw_shared_offset_list,
            0.0)
    };
    if (!seed_model_list.has_value())
    {
        return std::nullopt;
    }
    const auto seed_step_norm{
        CalculateModelTrustRegionStepNorm(previous_model_list, *seed_model_list)
    };
    if (!seed_step_norm.has_value())
    {
        return std::nullopt;
    }
    const auto candidate_model_list{
        BuildSharedOffsetDampedModelList(
            previous_model_list,
            raw_model_list,
            previous_shared_offset_list,
            raw_shared_offset_list,
            factor)
    };
    if (!candidate_model_list.has_value())
    {
        return std::nullopt;
    }
    const auto step_norm{
        CalculateModelTrustRegionStepNorm(previous_model_list, *candidate_model_list)
    };
    if (!step_norm.has_value())
    {
        return std::nullopt;
    }
    FitStateProposal proposal{
        .patch{ .atom_index_list = key },
        .effective_damping = factor,
        .step_norm = *step_norm
    };
    proposal.patch.mdpde_list.reserve(key.size());
    for (std::size_t atom_position = 0; atom_position < key.size(); atom_position++)
    {
        const auto atom_index{ key.at(atom_position) };
        proposal.patch.mdpde_list.emplace_back(
            GaussianModel3DWithUncertainty{
                candidate_model_list->at(atom_position),
                operator_proposal_state.at(atom_index).mdpde
                    .GetStandardDeviationModel()
            });
    }
    return proposal;
}

static bool ContainsClusterKey(const std::vector<ClusterKey> & key_list, const ClusterKey & key)
{
    return std::ranges::find(key_list, key) != key_list.end();
}

static void EraseClusterKey(std::vector<ClusterKey> & key_list, const ClusterKey & key)
{
    const auto removed{ std::ranges::remove(key_list, key) };
    key_list.erase(removed.begin(), removed.end());
}

static ClusterKey FlattenClusterKeyList(const std::vector<ClusterKey> & key_list)
{
    ClusterKey atom_index_list;
    for (const auto & key : key_list)
    {
        atom_index_list.insert(
            atom_index_list.end(),
            key.begin(),
            key.end());
    }
    std::ranges::sort(atom_index_list);
    atom_index_list.erase(
        std::ranges::unique(atom_index_list).begin(),
        atom_index_list.end());
    return atom_index_list;
}

static void RejectSelectionKeys(
    const CandidateSelectionInputs & inputs,
    const std::vector<ClusterKey> & key_list,
    bool exhausted,
    ClusterObjectiveStateMap & working_objective_state,
    CandidateSelection & selection)
{
    for (const auto & key : key_list)
    {
        if (!ContainsClusterKey(selection.accepted_key_list, key)) continue;
        for (const auto atom_index : key)
        {
            selection.assembled_state.at(atom_index) = inputs.previous_state.at(atom_index);
            selection.assembled_polish_provenance.at(atom_index) = inputs.previous_polish_provenance.at(atom_index);
        }
        EraseClusterKey(selection.accepted_key_list, key);
        EraseClusterKey(selection.grow_trust_region_key_list, key);
        selection.rejected_key_list.emplace_back(key);
        if (exhausted) selection.exhausted_key_list.emplace_back(key);
        working_objective_state.at(key) = inputs.cluster_objective_state.at(key);

        const auto diagnostic_iter{
            std::ranges::find(
                selection.accepted_cluster_diagnostic_list,
                key,
                &ClusterCandidateDiagnostic::key)
        };
        if (diagnostic_iter != selection.accepted_cluster_diagnostic_list.end())
        {
            selection.rejected_cluster_diagnostic_list.emplace_back(std::move(*diagnostic_iter));
            selection.accepted_cluster_diagnostic_list.erase(diagnostic_iter);
        }
    }
    std::ranges::sort(selection.accepted_key_list);
    std::ranges::sort(selection.rejected_key_list);
    selection.rejected_key_list.erase(
        std::ranges::unique(selection.rejected_key_list).begin(),
        selection.rejected_key_list.end());
    std::ranges::sort(selection.exhausted_key_list);
    selection.exhausted_key_list.erase(
        std::ranges::unique(selection.exhausted_key_list).begin(),
        selection.exhausted_key_list.end());
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
        candidate_snapshot = BuildSecondStageModelSnapshot(
            context,
            BuildFittedGaussianSnapshot(candidate_state));
    }
    catch (const std::exception &)
    {
        result.status = TrustModelPredictionStatus::ModelUnavailable;
        result.shadow_action = DetermineTrustModelShadowAction(result);
        return result;
    }

    double predicted_residual_reduction{ 0.0 };
    std::set<std::size_t> changed_unselected_dependency_set;
    for (const auto & sample_ref : objective_sample_ref_list)
    {
        const auto & owner_key{
            objective_domain.owner_key_by_atom_index.at(sample_ref.atom_index)
        };
        if (owner_key.empty()) continue;
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
        const auto is_fit_range{
            objective_domain.fit_sample_mask_by_atom.at(sample_ref.atom_index)
                .at(sample_ref.sample_index) != 0
        };
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

        const auto & atom_context{ context.at(sample_ref.atom_index) };
        const auto target_direction{
            EvaluateTrustModelResponseDirection(
                residual_baseline.model_snapshot.selected.at(sample_ref.atom_index),
                candidate_snapshot.selected.at(sample_ref.atom_index),
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
                ResolveNeighborAtomModel(neighbor, residual_baseline.model_snapshot)
            };
            const auto & candidate_neighbor{
                ResolveNeighborAtomModel(neighbor, candidate_snapshot)
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
            if (!neighbor.is_selected && *neighbor_direction != 0.0)
            {
                changed_unselected_dependency_set.emplace(neighbor.atom_index);
            }
        }
        const auto linearized_residual{
            previous_residual->residual + residual_direction
        };
        const auto weight{
            algorithm::CalculateCauchyWeight(
                previous_residual->residual,
                scale,
                kObjectiveRobustLossCutoffMultiplier)
        };
        const auto coefficient{
            CalculateClusterAtomWeight(owner_key.size(), objective_domain.active_atom_count) /
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
    result.unselected_dependency_count = changed_unselected_dependency_set.size();
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

static std::optional<StabilizationTerminalDiagnostic>
EvaluateClusterCandidateGuard(
    const CandidateSelectionInputs & inputs,
    const ClusterKey & key,
    const FitStateView & candidate_state,
    const SuspiciousBlockActivity & block_activity)
{
    const auto & previous_snapshot{ inputs.residual_baseline.model_snapshot };
    const auto candidate_snapshot{
        BuildSecondStageModelSnapshot(
            inputs.context,
            BuildFittedGaussianSnapshot(candidate_state))
    };
    for (const auto atom_index : key)
    {
        const auto & atom_context{ inputs.context.at(atom_index) };
        const auto & previous_model{
            inputs.previous_state.at(atom_index).mdpde.GetModel()
        };
        const auto & candidate_model{ candidate_state.GetModel(atom_index) };
        if (block_activity.HasActiveOffset(atom_index))
        {
            const auto offset_assessment{
                AssessSuspiciousGaussianUpdate(
                    atom_context.raw_sampling_entries,
                    candidate_model,
                    inputs.options,
                    BuildPreviousSuspiciousProfileBaseline(
                        atom_context.raw_sampling_entries,
                        previous_model,
                        inputs.options),
                    SuspiciousUpdateMode::OffsetOnly)
            };
            if (offset_assessment.IsSuspicious())
            {
                return StabilizationTerminalDiagnostic{
                    StabilizationTerminalReason::GuardInfeasible,
                    atom_index,
                    SuspiciousUpdateMode::OffsetOnly,
                    offset_assessment.reason
                };
            }
        }
        if (!block_activity.HasActiveShape(atom_index)) continue;
        const auto previous_samples{
            BuildSecondStageAdjustedSamples(atom_context, previous_snapshot)
        };
        const auto candidate_samples{
            BuildSecondStageAdjustedSamples(atom_context, candidate_snapshot)
        };
        const auto shape_assessment{
            AssessSuspiciousGaussianUpdate(
                candidate_samples,
                candidate_model,
                inputs.options,
                BuildPreviousSuspiciousProfileBaseline(
                    previous_samples,
                    previous_model,
                    inputs.options),
                SuspiciousUpdateMode::PostRefit)
        };
        if (shape_assessment.IsSuspicious())
        {
            return StabilizationTerminalDiagnostic{
                StabilizationTerminalReason::GuardInfeasible,
                atom_index,
                SuspiciousUpdateMode::PostRefit,
                shape_assessment.reason
            };
        }
    }
    return std::nullopt;
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
    const auto & previous_objective_state{ inputs.cluster_objective_state.at(key) };
    const auto trust_region_radius{ inputs.trust_region_state.GetRadius(key) };
    auto & performance_counters{ inputs.performance_counters };
    ClusterCandidateResult result;
    result.polish_provenance.reserve(key.size());
    for (const auto atom_index : key)
    {
        result.polish_provenance.emplace_back(previous_polish_provenance.at(atom_index));
    }
    auto search_endpoint_state{ inputs.operator_proposal_state };
    auto search_block_activity{ inputs.block_activity };
    std::vector<StabilizationTerminalDiagnostic> terminal_diagnostic_list;
    std::optional<FitStatePatch> accepted_patch;
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
                        search_block_activity.HasActiveOffset(atom_index);
                })
        };
        if (!has_active_parameter)
        {
            accepted_patch = FitStatePatch::FromState(previous_state, key);
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
                BuildSharedOffsetProposal(
                    context,
                    previous_state,
                    search_endpoint_state,
                    search_block_activity,
                    key,
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
            const FitStateView candidate_state_view{ previous_state, proposal.patch };
            const auto maximum_change{
                std::ranges::max(
                    SummarizeTransformedChanges(
                        candidate_state_view,
                        residual_baseline.model_snapshot.selected,
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
                    accepted_patch = std::move(proposal.patch);
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
                    inputs,
                    key,
                    candidate_state_view,
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
            const CandidateEvaluationOverlay candidate_overlay{
                context,
                residual_baseline,
                candidate_state_view
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
                accepted_patch = std::move(proposal.patch);
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
        if (accepted_patch.has_value()) break;

        if (result.diagnostic.objective_rejected_trial_count == 0 &&
            result.diagnostic.guard_rejected_trial_count != 0 &&
            last_guard_failure.has_value())
        {
            const auto atom_index{ *last_guard_failure->guard_atom_index };
            const auto mode{ *last_guard_failure->guard_mode };
            if (mode == SuspiciousUpdateMode::OffsetOnly)
            {
                const auto group_id{ context.at(atom_index).group_id };
                for (const auto member_index : key)
                {
                    if (context.at(member_index).group_id != group_id) continue;
                    const auto previous_offset{
                        previous_state.at(member_index).mdpde.GetModel().GetOffset()
                    };
                    search_endpoint_state.at(member_index).mdpde =
                        WithPreservedUncertaintyOffset(
                            search_endpoint_state.at(member_index).mdpde,
                            previous_offset);
                    search_block_activity.offset_fixed_atom_mask.at(member_index) = 1;
                }
            }
            else
            {
                const auto retained_offset{
                    search_endpoint_state.at(atom_index).mdpde.GetModel().GetOffset()
                };
                search_endpoint_state.at(atom_index).mdpde =
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

    auto base_patch{ std::move(*accepted_patch) };
    const FitStateView base_state_view{ previous_state, base_patch };
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
            const FitStateView polished_state_view{
                previous_state,
                polished_candidate->patch
            };
            const CandidateEvaluationOverlay polished_overlay{
                context,
                residual_baseline,
                polished_state_view
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
    if (!result.accepted_patch.has_value())
    {
        result.accepted_patch = std::move(base_patch);
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

struct BoundaryCandidateEvaluation
{
    ClusterObjectiveStateMap objective_state_by_key{};
    ObjectiveBreakdown audit_objective{};
    std::size_t locally_deteriorated_member_count{ 0 };
    double maximum_local_deterioration{ 0.0 };
};

static FitStatePatch BuildSelectionPatch(
    const CandidateSelection & selection,
    const std::vector<ClusterKey> & key_list)
{
    return FitStatePatch::FromState(
        selection.assembled_state,
        FlattenClusterKeyList(key_list));
}

static std::optional<BoundaryCandidateEvaluation>
EvaluateBoundaryComponentCandidate(
    const CandidateSelectionInputs & inputs,
    const BoundaryReconciliationComponent & component,
    const CandidateEvaluationOverlay & candidate_overlay,
    const ObjectiveBreakdown * previous_audit_objective,
    bool cooperative = false)
{
    BoundaryCandidateEvaluation evaluation;
    for (const auto & key : component.key_list)
    {
        auto objective_state{ inputs.cluster_objective_state.at(key) };
        ObjectiveAttemptDiagnostic diagnostic;
        const auto & previous_objective{ inputs.previous_objective_by_key.at(key) };
        if (!cooperative)
        {
            if (!TryCommitClusterCandidate(
                    candidate_overlay,
                    key,
                    inputs.partition.sample_id_list_by_key.at(key),
                    previous_objective.has_value() ? &*previous_objective : nullptr,
                    false,
                    inputs.objective_domain,
                    objective_state,
                    diagnostic,
                    inputs.performance_counters))
            {
                return std::nullopt;
            }
        }
        else
        {
            diagnostic.candidate_objective = EvaluateObjectiveContribution(
                candidate_overlay,
                key,
                inputs.partition.sample_id_list_by_key.at(key),
                inputs.objective_domain);
            if (!diagnostic.candidate_objective.has_value() || !previous_objective.has_value())
            {
                return std::nullopt;
            }
            const auto candidate_value{
                diagnostic.candidate_objective->GetTotalObjective()
            };
            const auto previous_value{ previous_objective->GetTotalObjective() };
            if (!std::isfinite(candidate_value) ||
                IsObjectiveDeteriorated(
                    candidate_value,
                    previous_value,
                    kObjectiveProgressTolerance))
            {
                return std::nullopt;
            }
            if (candidate_value > previous_value)
            {
                evaluation.locally_deteriorated_member_count++;
                evaluation.maximum_local_deterioration = std::max(
                    evaluation.maximum_local_deterioration,
                    candidate_value - previous_value);
            }
            const auto improves_member{
                IsBetterAuditObjective(
                    candidate_value,
                    previous_value,
                    kObjectiveStrictTolerance)
            };
            if (improves_member &&
                (!objective_state.best_objective.has_value() ||
                    IsBetterAuditObjective(
                        candidate_value,
                        objective_state.best_objective->GetTotalObjective(),
                        kObjectiveStrictTolerance)))
            {
                objective_state.best_objective = diagnostic.candidate_objective;
                objective_state.best_maximum_transformed_change =
                    std::ranges::max(
                        SummarizeTransformedChanges(
                            candidate_overlay.GetState(),
                            candidate_overlay.GetBaseline().model_snapshot.selected,
                            key).maximum_list);
            }
        }
        evaluation.objective_state_by_key.emplace(key, std::move(objective_state));
    }
    const auto * best_audit_objective{
        cooperative && inputs.best_audit_state.has_value() ?
            &inputs.best_audit_state->objective : nullptr
    };
    const auto audit_objective{ EvaluateCombinedObjective(
        candidate_overlay,
        component.affected_sample_ref_list,
        inputs.objective_domain,
        best_audit_objective,
        previous_audit_objective,
        inputs.performance_counters) };
    if (!audit_objective.has_value()) return std::nullopt;
    if (cooperative &&
        !IsBetterAuditObjective(
            audit_objective->GetTotalObjective(),
            previous_audit_objective->GetTotalObjective(),
            kObjectiveStrictTolerance))
    {
        return std::nullopt;
    }
    evaluation.audit_objective = *audit_objective;
    return evaluation;
}

static void CommitBoundaryObjectiveState(
    const BoundaryCandidateEvaluation & evaluation,
    ClusterObjectiveStateMap & working_objective_state)
{
    for (const auto & [key, objective_state] : evaluation.objective_state_by_key)
    {
        working_objective_state.at(key) = objective_state;
    }
}

static void RemoveTrustGrowthForKeys(
    const std::vector<ClusterKey> & key_list,
    CandidateSelection & selection)
{
    for (const auto & key : key_list)
    {
        EraseClusterKey(selection.grow_trust_region_key_list, key);
    }
}

static bool OverlayFitStatePatch(FitStatePatch & base_patch, const FitStatePatch & overlay_patch)
{
    for (std::size_t position = 0; position < overlay_patch.atom_index_list.size(); position++)
    {
        const auto atom_index{ overlay_patch.atom_index_list.at(position) };
        const auto iter{
            std::ranges::lower_bound(base_patch.atom_index_list, atom_index)
        };
        if (iter == base_patch.atom_index_list.end() || *iter != atom_index)
        {
            return false;
        }
        base_patch.mdpde_list.at(static_cast<std::size_t>(std::distance(
            base_patch.atom_index_list.begin(),
            iter))) = overlay_patch.mdpde_list.at(position);
    }
    return true;
}

static std::size_t CountSuspiciousPolishAtoms(
    const SecondStageContext & context,
    const FitOptions & options,
    const std::vector<std::size_t> & atom_index_list,
    const FitStateView & endpoint_state,
    const FitStateView & candidate_state)
{
    const auto endpoint_snapshot{
        BuildSecondStageModelSnapshot(
            context,
            BuildFittedGaussianSnapshot(endpoint_state))
    };
    const auto candidate_snapshot{
        BuildSecondStageModelSnapshot(
            context,
            BuildFittedGaussianSnapshot(candidate_state))
    };
    std::size_t suspicious_atom_count{ 0 };
    for (const auto atom_index : atom_index_list)
    {
        const auto change{
            CalculateTransformedChange(
                candidate_state.GetModel(atom_index),
                endpoint_state.GetModel(atom_index))
        };
        if (!IsTransformedChangeMaterial(change, kTransformedChangeTolerance))
        {
            continue;
        }
        const auto endpoint_samples{
            BuildSecondStageAdjustedSamples(
                context.at(atom_index),
                endpoint_snapshot)
        };
        const auto candidate_samples{
            BuildSecondStageAdjustedSamples(
                context.at(atom_index),
                candidate_snapshot)
        };
        const auto baseline{
            BuildPreviousSuspiciousProfileBaseline(
                endpoint_samples,
                endpoint_state.GetModel(atom_index),
                options)
        };
        if (AssessSuspiciousGaussianUpdate(
                candidate_samples,
                candidate_state.GetModel(atom_index),
                options,
                baseline,
                SuspiciousUpdateMode::PostRefit).reason !=
            SuspiciousGaussianReason::None)
        {
            suspicious_atom_count++;
        }
    }
    return suspicious_atom_count;
}

static bool TryBoundaryJointCorrection(
    const CandidateSelectionInputs & inputs,
    const BoundaryReconciliationComponent & component,
    const ObjectiveBreakdown & previous_audit_objective,
    const ObjectiveBreakdown & improvement_reference_objective,
    const FitStatePatch & endpoint_patch,
    ClusterObjectiveStateMap & working_objective_state,
    CandidateSelection & selection,
    BoundaryComponentReconciliationDiagnostic & diagnostic)
{
    if (component.offset_closure_atom_index_list.empty()) return false;

    std::vector<std::size_t> shape_active_atom_index_list;
    for (const auto atom_index : component.shape_active_atom_index_list)
    {
        if (inputs.block_activity.HasActiveShape(atom_index))
        {
            shape_active_atom_index_list.emplace_back(atom_index);
        }
    }
    std::set<std::size_t> fixed_offset_group_id_set;
    for (const auto atom_index : component.offset_closure_atom_index_list)
    {
        if (!inputs.block_activity.HasActiveOffset(atom_index))
        {
            fixed_offset_group_id_set.emplace(
                inputs.context.at(atom_index).group_id);
        }
    }
    std::vector<std::size_t> offset_active_atom_index_list;
    for (const auto atom_index : component.offset_closure_atom_index_list)
    {
        if (!fixed_offset_group_id_set.contains(
                inputs.context.at(atom_index).group_id))
        {
            offset_active_atom_index_list.emplace_back(atom_index);
        }
    }
    if (shape_active_atom_index_list.empty() &&
        offset_active_atom_index_list.empty())
    {
        return false;
    }
    diagnostic.shape_active_atom_count = shape_active_atom_index_list.size();
    diagnostic.offset_active_atom_count = offset_active_atom_index_list.size();
    const FitStateView endpoint_state_view{
        inputs.previous_state,
        endpoint_patch
    };
    std::vector<BoundaryJointTrustRegion> trust_region_list;
    trust_region_list.reserve(component.key_list.size());
    for (const auto & key : component.key_list)
    {
        trust_region_list.emplace_back(BoundaryJointTrustRegion{
            key,
            inputs.trust_region_state.GetRadius(key)
        });
    }
    const BoundaryJointCorrectionWorkspaceKey workspace_key{
        shape_active_atom_index_list,
        offset_active_atom_index_list,
        component.offset_closure_atom_index_list,
        component.affected_sample_ref_list
    };
    auto & solver{
        inputs.boundary_joint_correction_workspace_by_key.try_emplace(workspace_key).first->second
    };
    const auto start_time{ std::chrono::steady_clock::now() };
    const auto record_performance = [&](bool accepted)
    {
        inputs.performance_counters.RecordBoundaryJointCorrection(
            accepted,
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - start_time).count());
    };
    diagnostic.joint_reference_component_objective = improvement_reference_objective.GetTotalObjective();
    auto correction_result{
        BuildBoundaryJointCorrection(
            inputs.context,
            endpoint_state_view,
            shape_active_atom_index_list,
            offset_active_atom_index_list,
            component.offset_closure_atom_index_list,
            component.affected_sample_ref_list,
            inputs.ridge_multiplier_list,
            trust_region_list,
            solver)
    };
    diagnostic.joint_correction_status = correction_result.status;
    diagnostic.joint_parameter_count = correction_result.parameter_count;
    if (correction_result.status == BoundaryJointCorrectionStatus::CandidateReady)
    {
        diagnostic.joint_damping = correction_result.damping;
        diagnostic.maximum_normalized_trust_step = correction_result.maximum_normalized_trust_step;
    }
    if (correction_result.status != BoundaryJointCorrectionStatus::CandidateReady ||
        !correction_result.patch.has_value())
    {
        record_performance(false);
        return false;
    }

    auto corrected_component_patch{ endpoint_patch };
    if (!OverlayFitStatePatch(
            corrected_component_patch,
            *correction_result.patch))
    {
        record_performance(false);
        return false;
    }
    const FitStateView corrected_state_view{
        inputs.previous_state,
        corrected_component_patch
    };
    diagnostic.suspicious_candidate_atom_count =
        CountSuspiciousPolishAtoms(
            inputs.context,
            inputs.options,
            component.offset_closure_atom_index_list,
            endpoint_state_view,
            corrected_state_view);
    if (diagnostic.suspicious_candidate_atom_count != 0)
    {
        record_performance(false);
        return false;
    }
    const CandidateEvaluationOverlay corrected_overlay{
        inputs.context,
        inputs.residual_baseline,
        corrected_state_view
    };
    const auto raw_candidate_objective{
        EvaluateObjectiveDelta(
            corrected_overlay,
            component.affected_sample_ref_list,
            inputs.objective_domain,
            previous_audit_objective,
            inputs.performance_counters)
    };
    if (raw_candidate_objective.has_value())
    {
        diagnostic.joint_candidate_component_objective = raw_candidate_objective->GetTotalObjective();
    }
    const auto candidate_evaluation{
        EvaluateBoundaryComponentCandidate(
            inputs,
            component,
            corrected_overlay,
            &previous_audit_objective,
            diagnostic.is_rescue_attempt)
    };
    const auto is_strict_improvement{
        candidate_evaluation.has_value() &&
        IsBetterAuditObjective(
            candidate_evaluation->audit_objective.GetTotalObjective(),
            improvement_reference_objective.GetTotalObjective(),
            kObjectiveStrictTolerance)
    };
    if (!is_strict_improvement)
    {
        record_performance(false);
        return false;
    }

    corrected_component_patch.ApplyTo(selection.assembled_state);
    for (const auto atom_index : component.offset_closure_atom_index_list)
    {
        const auto change{
            CalculateTransformedChange(
                selection.assembled_state.at(atom_index).mdpde.GetModel(),
                endpoint_state_view.GetModel(atom_index))
        };
        if (IsTransformedChangeMaterial(change, kTransformedChangeTolerance))
        {
            selection.assembled_polish_provenance.at(atom_index) = 1;
        }
    }
    CommitBoundaryObjectiveState(*candidate_evaluation, working_objective_state);
    RemoveTrustGrowthForKeys(component.key_list, selection);
    diagnostic.accepted_source = BoundaryComponentAcceptedSource::JointCorrection;
    diagnostic.candidate_component_objective = candidate_evaluation->audit_objective.GetTotalObjective();
    diagnostic.locally_deteriorated_member_count = candidate_evaluation->locally_deteriorated_member_count;
    diagnostic.maximum_local_deterioration = candidate_evaluation->maximum_local_deterioration;
    record_performance(true);
    return true;
}

static void ReconcileBoundaryComponent(
    const CandidateSelectionInputs & inputs,
    const BoundaryReconciliationComponent & component,
    const ObjectiveBreakdown * previous_audit_objective,
    ClusterObjectiveStateMap & working_objective_state,
    CandidateSelection & selection)
{
    BoundaryComponentReconciliationDiagnostic diagnostic;
    diagnostic.key_list = component.key_list;
    diagnostic.atom_count = FlattenClusterKeyList(component.key_list).size();
    diagnostic.boundary_sample_count = component.boundary_sample_count;
    diagnostic.interface_atom_count = component.interface_atom_index_list.size();
    diagnostic.shape_active_atom_count = component.shape_active_atom_index_list.size();
    diagnostic.offset_closure_atom_count = component.offset_closure_atom_index_list.size();
    if (previous_audit_objective != nullptr)
    {
        diagnostic.previous_component_objective = previous_audit_objective->GetTotalObjective();
    }
    const auto endpoint_patch{ BuildSelectionPatch(selection, component.key_list) };
    const FitStateView endpoint_state_view{
        inputs.previous_state,
        endpoint_patch
    };
    const CandidateEvaluationOverlay endpoint_overlay{
        inputs.context,
        inputs.residual_baseline,
        endpoint_state_view
    };
    const auto endpoint_evaluation{
        EvaluateBoundaryComponentCandidate(
            inputs,
            component,
            endpoint_overlay,
            previous_audit_objective)
    };
    if (endpoint_evaluation.has_value())
    {
        diagnostic.endpoint_component_objective = endpoint_evaluation->audit_objective.GetTotalObjective();
        if (previous_audit_objective != nullptr &&
            TryBoundaryJointCorrection(
                inputs,
                component,
                *previous_audit_objective,
                endpoint_evaluation->audit_objective,
                endpoint_patch,
                working_objective_state,
                selection,
                diagnostic))
        {
            selection.boundary_reconciliation_diagnostic_list.emplace_back(std::move(diagnostic));
            return;
        }
        CommitBoundaryObjectiveState(*endpoint_evaluation, working_objective_state);
        diagnostic.accepted_factor = 1.0;
        diagnostic.accepted_source = BoundaryComponentAcceptedSource::Endpoint;
        diagnostic.candidate_component_objective = endpoint_evaluation->audit_objective.GetTotalObjective();
        selection.boundary_reconciliation_diagnostic_list.emplace_back(std::move(diagnostic));
        return;
    }

    if (previous_audit_objective != nullptr &&
        TryBoundaryJointCorrection(
            inputs,
            component,
            *previous_audit_objective,
            *previous_audit_objective,
            endpoint_patch,
            working_objective_state,
            selection,
            diagnostic))
    {
        selection.boundary_reconciliation_diagnostic_list.emplace_back(std::move(diagnostic));
        return;
    }

    BacktrackingWorkspace backtracking_workspace{
        inputs.context,
        inputs.previous_state,
        endpoint_patch,
        kTransformedChangeTolerance
    };
    BacktrackingStep step;
    std::optional<BoundaryCandidateEvaluation> accepted_evaluation;
    for (step = backtracking_workspace.BuildNextCandidate();
        step.status == BacktrackingStepStatus::CandidateReady;
        step = backtracking_workspace.BuildNextCandidate())
    {
        diagnostic.trial_count = step.trial_number;
        const FitStateView candidate_state_view{
            inputs.previous_state,
            backtracking_workspace.GetCandidatePatch()
        };
        const CandidateEvaluationOverlay candidate_overlay{
            inputs.context,
            inputs.residual_baseline,
            candidate_state_view
        };
        accepted_evaluation = EvaluateBoundaryComponentCandidate(
            inputs,
            component,
            candidate_overlay,
            previous_audit_objective);
        if (accepted_evaluation.has_value()) break;
    }
    if (!accepted_evaluation.has_value())
    {
        diagnostic.exhausted = step.status == BacktrackingStepStatus::Exhausted;
        RejectSelectionKeys(
            inputs,
            component.key_list,
            diagnostic.exhausted,
            working_objective_state,
            selection);
        selection.boundary_reconciliation_diagnostic_list.emplace_back(std::move(diagnostic));
        return;
    }

    backtracking_workspace.GetCandidatePatch().ApplyTo(selection.assembled_state);
    const auto reconciled_provenance{
        backtracking_workspace.BuildCandidatePolishProvenance(
            inputs.previous_polish_provenance,
            selection.assembled_polish_provenance)
    };
    for (const auto atom_index : FlattenClusterKeyList(component.key_list))
    {
        selection.assembled_polish_provenance.at(atom_index) = reconciled_provenance.at(atom_index);
    }
    CommitBoundaryObjectiveState(*accepted_evaluation, working_objective_state);
    RemoveTrustGrowthForKeys(component.key_list, selection);
    diagnostic.accepted_factor = step.factor;
    diagnostic.accepted_source = BoundaryComponentAcceptedSource::Backtracking;
    diagnostic.candidate_component_objective = accepted_evaluation->audit_objective.GetTotalObjective();
    selection.boundary_reconciliation_diagnostic_list.emplace_back(std::move(diagnostic));
}

static void PromoteBoundaryRescueKeys(
    const CandidateSelectionInputs & inputs,
    const std::vector<ClusterKey> & rescue_key_list,
    const FitStatePatch & endpoint_patch,
    BoundaryComponentAcceptedSource accepted_source,
    CandidateSelection & selection)
{
    const FitStateView endpoint_state{ inputs.previous_state, endpoint_patch };
    for (const auto & key : rescue_key_list)
    {
        EraseClusterKey(selection.rejected_key_list, key);
        EraseClusterKey(selection.exhausted_key_list, key);
        EraseClusterKey(selection.grow_trust_region_key_list, key);
        if (!ContainsClusterKey(selection.accepted_key_list, key))
        {
            selection.accepted_key_list.emplace_back(key);
        }
        for (const auto atom_index : key)
        {
            const auto changed_from_previous{
                IsTransformedChangeMaterial(
                    CalculateTransformedChange(
                        selection.assembled_state.at(atom_index).mdpde.GetModel(),
                        inputs.previous_state.at(atom_index).mdpde.GetModel()),
                    kTransformedChangeTolerance)
            };
            if (!changed_from_previous)
            {
                selection.assembled_polish_provenance.at(atom_index) = inputs.previous_polish_provenance.at(atom_index);
                continue;
            }
            const auto correction_changed_endpoint{
                accepted_source == BoundaryComponentAcceptedSource::JointCorrection &&
                IsTransformedChangeMaterial(
                    CalculateTransformedChange(
                        selection.assembled_state.at(atom_index).mdpde.GetModel(),
                        endpoint_state.GetModel(atom_index)),
                    kTransformedChangeTolerance)
            };
            selection.assembled_polish_provenance.at(atom_index) = correction_changed_endpoint ? 1 : 0;
        }

        const auto diagnostic_iter{
            std::ranges::find(
                selection.rejected_cluster_diagnostic_list,
                key,
                &ClusterCandidateDiagnostic::key)
        };
        if (diagnostic_iter != selection.rejected_cluster_diagnostic_list.end())
        {
            diagnostic_iter->boundary_rescued = true;
            selection.accepted_cluster_diagnostic_list.emplace_back(std::move(*diagnostic_iter));
            selection.rejected_cluster_diagnostic_list.erase(diagnostic_iter);
        }
    }
    std::ranges::sort(selection.accepted_key_list);
    std::ranges::sort(selection.rejected_key_list);
    std::ranges::sort(
        selection.accepted_cluster_diagnostic_list,
        {},
        &ClusterCandidateDiagnostic::key);
}

static bool TryRescueBoundaryComponent(
    const CandidateSelectionInputs & inputs,
    const BoundaryReconciliationComponent & component,
    const ObjectiveBreakdown & previous_audit_objective,
    const std::map<ClusterKey, FitStatePatch> & rescue_patch_by_key,
    ClusterObjectiveStateMap & working_objective_state,
    CandidateSelection & selection)
{
    std::vector<ClusterKey> rescue_key_list;
    for (const auto & key : component.key_list)
    {
        if (ContainsClusterKey(selection.rejected_key_list, key) &&
            rescue_patch_by_key.contains(key))
        {
            rescue_key_list.emplace_back(key);
        }
    }
    if (rescue_key_list.empty()) return false;

    BoundaryComponentReconciliationDiagnostic diagnostic;
    diagnostic.key_list = component.key_list;
    diagnostic.atom_count = FlattenClusterKeyList(component.key_list).size();
    diagnostic.boundary_sample_count = component.boundary_sample_count;
    diagnostic.interface_atom_count = component.interface_atom_index_list.size();
    diagnostic.shape_active_atom_count = component.shape_active_atom_index_list.size();
    diagnostic.offset_closure_atom_count = component.offset_closure_atom_index_list.size();
    diagnostic.accepted_cluster_count = component.key_list.size() - rescue_key_list.size();
    diagnostic.rescue_candidate_cluster_count = rescue_key_list.size();
    diagnostic.is_rescue_attempt = true;
    diagnostic.previous_component_objective = previous_audit_objective.GetTotalObjective();

    auto endpoint_patch{ BuildSelectionPatch(selection, component.key_list) };
    for (const auto & key : rescue_key_list)
    {
        if (!OverlayFitStatePatch(endpoint_patch, rescue_patch_by_key.at(key)))
        {
            throw std::logic_error(
                "Boundary rescue candidate patch does not match its component.");
        }
    }
    const FitStateView endpoint_state_view{ inputs.previous_state, endpoint_patch };
    const CandidateEvaluationOverlay endpoint_overlay{
        inputs.context,
        inputs.residual_baseline,
        endpoint_state_view
    };
    const auto endpoint_evaluation{
        EvaluateBoundaryComponentCandidate(
            inputs,
            component,
            endpoint_overlay,
            &previous_audit_objective,
            true)
    };
    if (endpoint_evaluation.has_value())
    {
        diagnostic.endpoint_component_objective =
            endpoint_evaluation->audit_objective.GetTotalObjective();
        if (!TryBoundaryJointCorrection(
                inputs,
                component,
                previous_audit_objective,
                endpoint_evaluation->audit_objective,
                endpoint_patch,
                working_objective_state,
                selection,
                diagnostic))
        {
            endpoint_patch.ApplyTo(selection.assembled_state);
            CommitBoundaryObjectiveState(*endpoint_evaluation, working_objective_state);
            diagnostic.accepted_factor = 1.0;
            diagnostic.accepted_source = BoundaryComponentAcceptedSource::Endpoint;
            diagnostic.candidate_component_objective =
                endpoint_evaluation->audit_objective.GetTotalObjective();
            diagnostic.locally_deteriorated_member_count =
                endpoint_evaluation->locally_deteriorated_member_count;
            diagnostic.maximum_local_deterioration =
                endpoint_evaluation->maximum_local_deterioration;
        }
    }
    else if (!TryBoundaryJointCorrection(
            inputs,
            component,
            previous_audit_objective,
            previous_audit_objective,
            endpoint_patch,
            working_objective_state,
            selection,
            diagnostic))
    {
        BacktrackingWorkspace backtracking_workspace{
            inputs.context,
            inputs.previous_state,
            endpoint_patch,
            kTransformedChangeTolerance
        };
        BacktrackingStep step;
        std::optional<BoundaryCandidateEvaluation> accepted_evaluation;
        for (step = backtracking_workspace.BuildNextCandidate();
            step.status == BacktrackingStepStatus::CandidateReady;
            step = backtracking_workspace.BuildNextCandidate())
        {
            diagnostic.trial_count = step.trial_number;
            const FitStateView candidate_state_view{
                inputs.previous_state,
                backtracking_workspace.GetCandidatePatch()
            };
            const CandidateEvaluationOverlay candidate_overlay{
                inputs.context,
                inputs.residual_baseline,
                candidate_state_view
            };
            accepted_evaluation = EvaluateBoundaryComponentCandidate(
                inputs,
                component,
                candidate_overlay,
                &previous_audit_objective,
                true);
            if (accepted_evaluation.has_value()) break;
        }
        if (accepted_evaluation.has_value())
        {
            backtracking_workspace.GetCandidatePatch().ApplyTo(
                selection.assembled_state);
            const auto reconciled_provenance{
                backtracking_workspace.BuildCandidatePolishProvenance(
                    inputs.previous_polish_provenance,
                    selection.assembled_polish_provenance)
            };
            for (const auto atom_index : FlattenClusterKeyList(component.key_list))
            {
                selection.assembled_polish_provenance.at(atom_index) = reconciled_provenance.at(atom_index);
            }
            CommitBoundaryObjectiveState(*accepted_evaluation, working_objective_state);
            diagnostic.accepted_factor = step.factor;
            diagnostic.accepted_source = BoundaryComponentAcceptedSource::Backtracking;
            diagnostic.candidate_component_objective =
                accepted_evaluation->audit_objective.GetTotalObjective();
            diagnostic.locally_deteriorated_member_count =
                accepted_evaluation->locally_deteriorated_member_count;
            diagnostic.maximum_local_deterioration =
                accepted_evaluation->maximum_local_deterioration;
        }
        else
        {
            diagnostic.exhausted = step.status == BacktrackingStepStatus::Exhausted;
        }
    }

    if (diagnostic.accepted_source != BoundaryComponentAcceptedSource::None)
    {
        if (diagnostic.previous_component_objective.has_value() &&
            diagnostic.candidate_component_objective.has_value())
        {
            diagnostic.component_improvement =
                *diagnostic.previous_component_objective -
                *diagnostic.candidate_component_objective;
        }
        diagnostic.rescued_cluster_count = rescue_key_list.size();
        PromoteBoundaryRescueKeys(
            inputs,
            rescue_key_list,
            endpoint_patch,
            diagnostic.accepted_source,
            selection);
    }
    inputs.performance_counters.RecordBoundaryRescue(
        diagnostic.accepted_source != BoundaryComponentAcceptedSource::None,
        diagnostic.accepted_source != BoundaryComponentAcceptedSource::Endpoint);
    selection.boundary_reconciliation_diagnostic_list.emplace_back(std::move(diagnostic));
    return selection.boundary_reconciliation_diagnostic_list.back()
            .accepted_source != BoundaryComponentAcceptedSource::None;
}

static std::vector<BoundaryReconciliationComponent>
BuildExpandedBoundaryReconciliationComponents(
    const CandidateSelectionInputs & inputs,
    const std::vector<ClusterKey> & key_list)
{
    auto component_list{
        BuildBoundaryReconciliationComponents(
            inputs.context,
            inputs.partition,
            key_list)
    };
    for (auto & component : component_list)
    {
        component = ExpandBoundaryReconciliationHalo(
            inputs.context,
            std::move(component),
            inputs.options.second_stage_boundary_halo_depth);
    }
    return component_list;
}

static bool RescueRejectedBoundaryClusters(
    const CandidateSelectionInputs & inputs,
    const ObjectiveBreakdown & previous_audit_objective,
    const std::map<ClusterKey, FitStatePatch> & rescue_patch_by_key,
    ClusterObjectiveStateMap & working_objective_state,
    CandidateSelection & selection)
{
    std::vector<ClusterKey> eligible_key_list;
    for (const auto & key : selection.accepted_key_list)
    {
        eligible_key_list.emplace_back(key);
    }
    for (const auto & key : selection.rejected_key_list)
    {
        if (rescue_patch_by_key.contains(key))
        {
            eligible_key_list.emplace_back(key);
        }
    }
    std::ranges::sort(eligible_key_list);
    eligible_key_list.erase(
        std::ranges::unique(eligible_key_list).begin(),
        eligible_key_list.end());

    bool rescued_any{ false };
    for (const auto & component : BuildExpandedBoundaryReconciliationComponents(
        inputs,
        eligible_key_list))
    {
        rescued_any = TryRescueBoundaryComponent(
            inputs,
            component,
            previous_audit_objective,
            rescue_patch_by_key,
            working_objective_state,
            selection) || rescued_any;
    }
    return rescued_any;
}

static std::optional<ObjectiveBreakdown> EvaluateFinalSelectionAudit(
    const CandidateSelectionInputs & inputs,
    const ObjectiveBreakdown & previous_audit_objective,
    const CandidateSelection & selection)
{
    if (selection.accepted_key_list.empty()) return std::nullopt;
    const auto candidate_patch{
        BuildSelectionPatch(selection, selection.accepted_key_list)
    };
    const FitStateView candidate_state_view{
        inputs.previous_state,
        candidate_patch
    };
    const CandidateEvaluationOverlay candidate_overlay{
        inputs.context,
        inputs.residual_baseline,
        candidate_state_view
    };
    const auto affected_sample_ref_list{
        BuildGraphAffectedSampleUnion(inputs.partition, selection.accepted_key_list)
    };
    const auto * best_audit_objective{
        inputs.best_audit_state.has_value() ? &inputs.best_audit_state->objective : nullptr
    };
    return EvaluateCombinedObjective(
        candidate_overlay,
        affected_sample_ref_list,
        inputs.objective_domain,
        best_audit_objective,
        &previous_audit_objective,
        inputs.performance_counters);
}

static void MarkBoundaryDiagnosticRejected(
    const std::vector<ClusterKey> & key_list,
    bool exhausted,
    CandidateSelection & selection)
{
    auto iter{ std::ranges::find(
        selection.boundary_reconciliation_diagnostic_list | std::views::reverse,
        key_list,
        &BoundaryComponentReconciliationDiagnostic::key_list) };
    if (iter == selection.boundary_reconciliation_diagnostic_list.rend()) return;
    iter->accepted_factor.reset();
    iter->accepted_source = BoundaryComponentAcceptedSource::None;
    iter->exhausted = exhausted;
}

static void AuditAndSalvageFinalSelection(
    const CandidateSelectionInputs & inputs,
    const ObjectiveBreakdown & previous_audit_objective,
    ClusterObjectiveStateMap & working_objective_state,
    CandidateSelection & selection)
{
    selection.final_audit_objective = EvaluateFinalSelectionAudit(
        inputs,
        previous_audit_objective,
        selection);
    if (selection.final_audit_objective.has_value() ||
        selection.accepted_key_list.empty())
    {
        return;
    }

    auto component_list{
        BuildExpandedBoundaryReconciliationComponents(
            inputs,
            selection.accepted_key_list)
    };
    for (const auto & key : selection.accepted_key_list)
    {
        if (std::ranges::any_of(
                component_list,
                [&](const auto & component)
                {
                    return ContainsClusterKey(component.key_list, key);
                }))
        {
            continue;
        }
        component_list.emplace_back(BoundaryReconciliationComponent{
            .key_list = { key },
            .affected_sample_ref_list =
                inputs.partition.sample_id_list_by_key.at(key)
        });
    }

    std::vector<std::pair<double, std::vector<ClusterKey>>>
        rejection_candidate_list;
    for (const auto & component : component_list)
    {
        const auto candidate_patch{
            BuildSelectionPatch(selection, component.key_list)
        };
        const FitStateView candidate_state_view{
            inputs.previous_state,
            candidate_patch
        };
        const CandidateEvaluationOverlay candidate_overlay{
            inputs.context,
            inputs.residual_baseline,
            candidate_state_view
        };
        const auto candidate_audit_objective{
            EvaluateObjectiveDelta(
                candidate_overlay,
                component.affected_sample_ref_list,
                inputs.objective_domain,
                previous_audit_objective,
                inputs.performance_counters)
        };
        if (candidate_audit_objective.has_value() &&
            IsBetterAuditObjective(
                candidate_audit_objective->GetTotalObjective(),
                previous_audit_objective.GetTotalObjective(),
                kObjectiveStrictTolerance))
        {
            continue;
        }
        rejection_candidate_list.emplace_back(
            candidate_audit_objective.has_value() ?
                candidate_audit_objective->GetTotalObjective() :
                std::numeric_limits<double>::infinity(),
            component.key_list);
    }
    std::ranges::sort(
        rejection_candidate_list,
        [](const auto & lhs, const auto & rhs)
        {
            if (lhs.first != rhs.first)
            {
                return lhs.first > rhs.first;
            }
            return lhs.second < rhs.second;
        });

    for (const auto & rejection_candidate : rejection_candidate_list)
    {
        const auto & key_list{ rejection_candidate.second };
        MarkBoundaryDiagnosticRejected(key_list, false, selection);
        RejectSelectionKeys(
            inputs,
            key_list,
            false,
            working_objective_state,
            selection);
        selection.final_audit_objective = EvaluateFinalSelectionAudit(
            inputs,
            previous_audit_objective,
            selection);
        if (selection.final_audit_objective.has_value()) return;
    }

    const auto remaining_key_list{ selection.accepted_key_list };
    for (const auto & component : BuildBoundaryReconciliationComponents(
        inputs.context,
        inputs.partition,
        remaining_key_list))
    {
        MarkBoundaryDiagnosticRejected(component.key_list, true, selection);
    }
    RejectSelectionKeys(
        inputs,
        remaining_key_list,
        true,
        working_objective_state,
        selection);
    selection.final_audit_objective.reset();
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
    auto working_objective_state{ inputs.cluster_objective_state };
    const auto candidate_phase_start{ std::chrono::steady_clock::now() };
    const auto & partition{ inputs.partition };
    auto & solver_workspace_by_key{ inputs.solver_workspace_by_key };
    std::vector<CandidateWork> candidate_work_list;
    candidate_work_list.reserve(partition.sample_id_list_by_key.size());
    for (const auto & entry : partition.sample_id_list_by_key)
    {
        candidate_work_list.emplace_back(
            &entry,
            &solver_workspace_by_key.at(entry.first));
    }

    std::vector<ClusterCandidateResult> result_list(candidate_work_list.size());
    std::vector<std::exception_ptr> exception_list(candidate_work_list.size());
    const auto select_candidate = [&](std::size_t position)
    {
        try
        {
            const auto & entry{ *candidate_work_list.at(position).first };
            result_list.at(position) = SelectClusterCandidate(
                inputs,
                entry.first,
                entry.second,
                *candidate_work_list.at(position).second);
        }
        catch (...)
        {
            exception_list.at(position) = std::current_exception();
        }
    };
#ifdef USE_OPENMP
    if (inputs.options.thread_size > 1 && candidate_work_list.size() > 1)
    {
        eigen_helper::ScopedEigenThreadCount eigen_thread_guard{ 1 };
#pragma omp parallel for schedule(dynamic) num_threads(inputs.options.thread_size)
        for (std::size_t position = 0;
            position < candidate_work_list.size(); position++)
        {
            select_candidate(position);
        }
    }
    else
#endif
    {
        for (std::size_t position = 0;
            position < candidate_work_list.size(); position++)
        {
            select_candidate(position);
        }
    }
    for (const auto & exception : exception_list)
    {
        if (exception) std::rethrow_exception(exception);
    }

    CandidateSelection selection{
        .assembled_state = inputs.previous_state,
        .assembled_polish_provenance = inputs.previous_polish_provenance
    };
    std::map<ClusterKey, FitStatePatch> rescue_patch_by_key;
    std::vector<ClusterKey> locally_polished_key_list;
    locally_polished_key_list.reserve(result_list.size());
    for (std::size_t position = 0; position < result_list.size(); position++)
    {
        auto & result{ result_list.at(position) };
        const auto & key{ candidate_work_list.at(position).first->first };
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
                const auto group_id{ inputs.context.at(atom_index).group_id };
                for (const auto member_index : key)
                {
                    if (inputs.context.at(member_index).group_id == group_id)
                    {
                        inputs.block_activity.offset_fixed_atom_mask.at(
                            member_index) = 1;
                    }
                }
            }
            else
            {
                inputs.block_activity.shape_fixed_atom_mask.at(atom_index) = 1;
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
        working_objective_state.at(key) = std::move(result.objective_state);

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

    const auto boundary_component_list{
        BuildExpandedBoundaryReconciliationComponents(
            inputs,
            selection.accepted_key_list)
    };
    const auto previous_audit_objective{
        EvaluateAuditObjective(inputs.objective_domain, inputs.residual_baseline)
    };
    if (!boundary_component_list.empty())
    {
        const auto boundary_reconciliation_start{ std::chrono::steady_clock::now() };
        for (const auto & component : boundary_component_list)
        {
            ReconcileBoundaryComponent(
                inputs,
                component,
                previous_audit_objective.has_value() ?
                    &*previous_audit_objective : nullptr,
                working_objective_state,
                selection);
        }
        if (!previous_audit_objective.has_value())
        {
            const auto remaining_key_list{ selection.accepted_key_list };
            RejectSelectionKeys(
                inputs,
                remaining_key_list,
                true,
                working_objective_state,
                selection);
        }
        else
        {
            AuditAndSalvageFinalSelection(
                inputs,
                *previous_audit_objective,
                working_objective_state,
                selection);
        }
        const auto backtracked_component_count{
            std::ranges::count_if(
                selection.boundary_reconciliation_diagnostic_list,
                [](const auto & diagnostic)
                {
                    return diagnostic.accepted_source !=
                            BoundaryComponentAcceptedSource::None &&
                        diagnostic.accepted_factor.has_value() &&
                        *diagnostic.accepted_factor < 1.0;
                })
        };
        const auto rejected_component_count{
            std::ranges::count_if(
                selection.boundary_reconciliation_diagnostic_list,
                [](const auto & diagnostic)
                {
                    return diagnostic.accepted_source ==
                        BoundaryComponentAcceptedSource::None;
                })
        };
        inputs.performance_counters.RecordBoundaryReconciliation(
            selection.boundary_reconciliation_diagnostic_list.size(),
            static_cast<std::size_t>(backtracked_component_count),
            static_cast<std::size_t>(rejected_component_count),
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - boundary_reconciliation_start).count());
    }
    if (previous_audit_objective.has_value() &&
        RescueRejectedBoundaryClusters(
            inputs,
            *previous_audit_objective,
            rescue_patch_by_key,
            working_objective_state,
            selection))
    {
        AuditAndSalvageFinalSelection(
            inputs,
            *previous_audit_objective,
            working_objective_state,
            selection);
    }
    if (previous_audit_objective.has_value() &&
        selection.final_audit_objective.has_value())
    {
        const auto global_improvement{
            previous_audit_objective->GetTotalObjective() -
            selection.final_audit_objective->GetTotalObjective()
        };
        for (auto & diagnostic : selection.boundary_reconciliation_diagnostic_list)
        {
            if (diagnostic.is_rescue_attempt &&
                diagnostic.accepted_source !=
                    BoundaryComponentAcceptedSource::None)
            {
                diagnostic.global_improvement = global_improvement;
            }
        }
    }
    inputs.cluster_objective_state = std::move(working_objective_state);
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

FinalDependencyPolishResult RunFinalDependencyPolish(
    const SecondStageContext & context,
    const FitOptions & options,
    const GraphTopology & topology,
    const CouplingGraphPartition & partition,
    const ObjectiveDomain & objective_domain,
    const SuspiciousBlockActivity & block_activity,
    const TrustRegionStateSet & trust_region_state,
    const FitState & base_state,
    BoundaryJointCorrectionWorkspaceMap & workspace_by_key,
    PerformanceCounters & performance_counters)
{
    FinalDependencyPolishResult result{ .state = base_state };
    if (!options.enable_second_stage_dependency_polish) return result;

    const auto polish_start{ std::chrono::steady_clock::now() };
    const auto component_list{
        BuildUncutDependencyPolishComponents(
            topology,
            partition,
            objective_domain.owner_key_by_atom_index)
    };
    result.diagnostic.component_count = component_list.size();
    result.diagnostic.component_list.reserve(component_list.size());
    const auto base_baseline{ BuildResidualBaseline(context, base_state) };
    performance_counters.RecordGaussianCacheMisses();
    const auto base_objective{
        EvaluateAuditObjective(objective_domain, base_baseline)
    };
    if (base_objective.has_value())
    {
        result.objective = base_objective;
        result.diagnostic.objective_before = base_objective->GetTotalObjective();
        result.diagnostic.objective_after = base_objective->GetTotalObjective();
    }

    std::vector<std::optional<FitStatePatch>> accepted_patch_by_component(
        component_list.size());
    std::size_t accepted_component_count{ 0 };
    std::vector<double> ridge_multiplier_list(context.size(), 1.0);

    for (std::size_t component_position = 0;
        component_position < component_list.size();
        component_position++)
    {
        const auto & component{ component_list.at(component_position) };
        std::vector<std::size_t> shape_active_index_list;
        std::set<std::size_t> fixed_offset_group_id_set;
        for (const auto atom_index : component.atom_index_list)
        {
            if (block_activity.HasActiveShape(atom_index))
            {
                shape_active_index_list.emplace_back(atom_index);
            }
            if (!block_activity.HasActiveOffset(atom_index))
            {
                fixed_offset_group_id_set.emplace(
                    context.at(atom_index).group_id);
            }
        }
        std::vector<std::size_t> offset_active_index_list;
        std::set<std::size_t> active_offset_group_id_set;
        for (const auto atom_index : component.atom_index_list)
        {
            const auto group_id{ context.at(atom_index).group_id };
            if (!fixed_offset_group_id_set.contains(group_id))
            {
                offset_active_index_list.emplace_back(atom_index);
                active_offset_group_id_set.emplace(group_id);
            }
        }
        const auto parameter_count{
            2 * shape_active_index_list.size() +
            active_offset_group_id_set.size()
        };
        result.diagnostic.atom_count += component.atom_index_list.size();
        result.diagnostic.parameter_count += parameter_count;
        auto & diagnostic{
            result.diagnostic.component_list.emplace_back(
                FinalDependencyPolishDiagnostic::Component{
                    .key_list = component.key_list,
                    .atom_count = component.atom_index_list.size(),
                    .parameter_count = parameter_count
                })
        };
        if (!base_objective.has_value()) continue;

        const auto component_start{ std::chrono::steady_clock::now() };
        if (shape_active_index_list.empty() && offset_active_index_list.empty())
        {
            diagnostic.elapsed_milliseconds =
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - component_start).count();
            continue;
        }
        diagnostic.objective_before = base_objective->GetTotalObjective();
        result.diagnostic.attempted_component_count++;

        try
        {
            FitStatePatch endpoint_patch{
                FitStatePatch::FromState(base_state, component.atom_index_list)
            };
            auto endpoint_objective{ *base_objective };
            const auto maximum_round_count{
                options.second_stage_dependency_polish_max_iterations
            };
            if (maximum_round_count != 0)
            {
                result.diagnostic.round_count++;
                diagnostic.round_count++;
                std::vector<BoundaryJointTrustRegion> trust_region_list;
                trust_region_list.reserve(component.key_list.size());
                for (const auto & key : component.key_list)
                {
                    trust_region_list.emplace_back(BoundaryJointTrustRegion{
                        key,
                        trust_region_state.GetRadius(key)
                    });
                }
                const BoundaryJointCorrectionWorkspaceKey workspace_key{
                    shape_active_index_list,
                    offset_active_index_list,
                    component.atom_index_list,
                    component.affected_sample_ref_list,
                };
                auto & solver{
                    workspace_by_key.try_emplace(workspace_key).first->second
                };
                for (std::size_t round = 0;
                    round < maximum_round_count;
                    round++)
                {
                    if (round != 0)
                    {
                        result.diagnostic.round_count++;
                        diagnostic.round_count++;
                    }
                    const FitStateView endpoint_state_view{
                        base_state,
                        endpoint_patch
                    };
                    const auto symbolic_analysis_count_before{
                        solver.GetSymbolicAnalysisCount()
                    };
                    const auto correction_result{
                        BuildBoundaryJointCorrection(
                            context,
                            endpoint_state_view,
                            shape_active_index_list,
                            offset_active_index_list,
                            component.atom_index_list,
                            component.affected_sample_ref_list,
                            ridge_multiplier_list,
                            trust_region_list,
                            solver)
                    };
                    diagnostic.symbolic_analysis_count +=
                        solver.GetSymbolicAnalysisCount() -
                        symbolic_analysis_count_before;
                    if (correction_result.status !=
                            BoundaryJointCorrectionStatus::CandidateReady ||
                        !correction_result.patch.has_value())
                    {
                        break;
                    }
                    diagnostic.parameter_count = correction_result.parameter_count;
                    const FitStateView candidate_state_view{
                        base_state,
                        *correction_result.patch
                    };
                    const auto has_invalid_model{
                        std::ranges::any_of(
                            component.atom_index_list,
                            [&](const auto atom_index)
                            {
                                return !IsValidSecondStageGaussianModel(
                                    candidate_state_view.GetModel(atom_index));
                            })
                    };
                    if (has_invalid_model) break;

                    const auto suspicious_atom_count{
                        CountSuspiciousPolishAtoms(
                            context,
                            options,
                            component.atom_index_list,
                            endpoint_state_view,
                            candidate_state_view)
                    };
                    diagnostic.suspicious_candidate_atom_count +=
                        suspicious_atom_count;
                    result.diagnostic.suspicious_candidate_atom_count +=
                        suspicious_atom_count;
                    if (suspicious_atom_count != 0) break;

                    const CandidateEvaluationOverlay candidate_overlay{
                        context,
                        base_baseline,
                        candidate_state_view
                    };
                    const auto candidate_objective{
                        EvaluateObjectiveDelta(
                            candidate_overlay,
                            component.affected_sample_ref_list,
                            objective_domain,
                            *base_objective,
                            performance_counters)
                    };
                    if (!candidate_objective.has_value() ||
                        !IsBetterAuditObjective(
                            candidate_objective->GetTotalObjective(),
                            endpoint_objective.GetTotalObjective(),
                            kObjectiveStrictTolerance))
                    {
                        break;
                    }

                    const auto member_guard_passed{
                        std::ranges::all_of(
                            component.key_list,
                            [&](const auto & key)
                            {
                                const auto sample_iter{
                                    partition.sample_id_list_by_key.find(key)
                                };
                                if (sample_iter ==
                                    partition.sample_id_list_by_key.end())
                                {
                                    return false;
                                }
                                auto owned_sample_ref_list{ sample_iter->second };
                                owned_sample_ref_list.erase(
                                    std::remove_if(
                                        owned_sample_ref_list.begin(),
                                        owned_sample_ref_list.end(),
                                        [&](const auto & sample_ref)
                                        {
                                            return sample_ref.atom_index >=
                                                    objective_domain.owner_key_by_atom_index.size() ||
                                                objective_domain.owner_key_by_atom_index.at(
                                                    sample_ref.atom_index).empty();
                                        }),
                                    owned_sample_ref_list.end());
                                const auto base_contribution{
                                    EvaluateObjectiveContribution(
                                        base_baseline,
                                        key,
                                        owned_sample_ref_list,
                                        objective_domain)
                                };
                                const auto candidate_contribution{
                                    EvaluateObjectiveContribution(
                                        candidate_overlay,
                                        key,
                                        owned_sample_ref_list,
                                        objective_domain)
                                };
                                return base_contribution.has_value() &&
                                    candidate_contribution.has_value() &&
                                    !IsObjectiveDeteriorated(
                                        candidate_contribution->GetTotalObjective(),
                                        base_contribution->GetTotalObjective(),
                                        kObjectiveProgressTolerance);
                            })
                    };
                    if (!member_guard_passed) break;

                    endpoint_patch = *correction_result.patch;
                    endpoint_objective = *candidate_objective;
                    diagnostic.objective_after =
                        candidate_objective->GetTotalObjective();
                }
            }

            if (diagnostic.objective_after.has_value() &&
                IsBetterAuditObjective(
                    *diagnostic.objective_after,
                    *diagnostic.objective_before,
                    kObjectiveStrictTolerance))
            {
                diagnostic.accepted = true;
                accepted_patch_by_component.at(component_position) =
                    std::move(endpoint_patch);
                accepted_component_count++;
            }
        }
        catch (const std::exception &)
        {
            diagnostic.accepted = false;
        }
        diagnostic.elapsed_milliseconds =
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - component_start).count();
    }

    auto assembled_state{ base_state };
    for (const auto & accepted_patch : accepted_patch_by_component)
    {
        if (accepted_patch.has_value())
        {
            accepted_patch->ApplyTo(assembled_state);
        }
    }
    const auto evaluate_global_audit = [&](const FitState & state)
    {
        performance_counters.RecordFullStateMaterialization();
        const auto snapshot{ BuildSecondStageModelSnapshot(context, state) };
        const SnapshotResidualEvaluator evaluator{ context, snapshot };
        return EvaluateAuditObjective(objective_domain, evaluator);
    };

    auto assembled_objective{
        accepted_component_count == 0 ?
            base_objective : evaluate_global_audit(assembled_state)
    };
    while (base_objective.has_value() &&
        accepted_component_count != 0 &&
        (!assembled_objective.has_value() ||
            !IsBetterAuditObjective(
                assembled_objective->GetTotalObjective(),
                base_objective->GetTotalObjective(),
                kObjectiveStrictTolerance)))
    {
        std::optional<std::size_t> removal_position;
        std::optional<ObjectiveBreakdown> best_removal_objective;
        FitState best_removal_state;
        for (std::size_t candidate_position = 0;
            candidate_position < accepted_patch_by_component.size();
            candidate_position++)
        {
            if (!accepted_patch_by_component.at(candidate_position).has_value())
            {
                continue;
            }
            auto candidate_state{ base_state };
            for (std::size_t patch_position = 0;
                patch_position < accepted_patch_by_component.size();
                patch_position++)
            {
                const auto & accepted_patch{
                    accepted_patch_by_component.at(patch_position)
                };
                if (patch_position == candidate_position ||
                    !accepted_patch.has_value())
                {
                    continue;
                }
                accepted_patch->ApplyTo(candidate_state);
            }
            const auto candidate_objective{
                evaluate_global_audit(candidate_state)
            };
            if (!candidate_objective.has_value() ||
                (assembled_objective.has_value() &&
                    candidate_objective->GetTotalObjective() >=
                        assembled_objective->GetTotalObjective()))
            {
                continue;
            }
            if (!best_removal_objective.has_value() ||
                candidate_objective->GetTotalObjective() <
                    best_removal_objective->GetTotalObjective())
            {
                removal_position = candidate_position;
                best_removal_objective = candidate_objective;
                best_removal_state = std::move(candidate_state);
            }
        }
        if (!removal_position.has_value()) break;
        auto & removed_diagnostic{
            result.diagnostic.component_list.at(*removal_position)
        };
        removed_diagnostic.accepted = false;
        accepted_patch_by_component.at(*removal_position).reset();
        accepted_component_count--;
        assembled_state = std::move(best_removal_state);
        assembled_objective = best_removal_objective;
    }

    if (base_objective.has_value() &&
        assembled_objective.has_value() &&
        accepted_component_count != 0 &&
        IsBetterAuditObjective(
            assembled_objective->GetTotalObjective(),
            base_objective->GetTotalObjective(),
            kObjectiveStrictTolerance))
    {
        result.state = std::move(assembled_state);
        result.objective = assembled_objective;
        result.accepted = true;
        result.diagnostic.objective_after =
            assembled_objective->GetTotalObjective();
    }
    else
    {
        for (auto & diagnostic : result.diagnostic.component_list)
        {
            diagnostic.accepted = false;
        }
    }
    result.diagnostic.accepted_component_count =
        result.accepted ? accepted_component_count : 0;
    result.diagnostic.elapsed_milliseconds =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - polish_start).count();
    performance_counters.RecordDependencyPolish(
        result.diagnostic.component_count,
        result.diagnostic.attempted_component_count,
        result.diagnostic.accepted_component_count,
        result.diagnostic.component_count -
            result.diagnostic.accepted_component_count,
        result.diagnostic.atom_count,
        result.diagnostic.parameter_count,
        result.diagnostic.round_count,
        result.diagnostic.elapsed_milliseconds);
    return result;
}

} // namespace rhbm_gem::core::detail
