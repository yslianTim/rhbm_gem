#include "core/detail/CandidateSelection.hpp"

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

#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/EigenHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

namespace rhbm_gem::core::detail {

namespace {

constexpr double kSuspiciousJointOffsetRidgeMultiplier{ 10.0 };
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

enum class SuspiciousProfileAnalysisMode
{
    Candidate,
    PreviousBaseline
};

struct BaseProposalBuildResult
{
    std::optional<FitStateProposal> proposal{};
    PreObjectiveFailureReason failure_reason{ PreObjectiveFailureReason::None };
    std::optional<double> attempted_step_norm{};
};

struct ClusterCandidateResult
{
    std::optional<FitStatePatch> accepted_patch{};
    PolishProvenance polish_provenance{};
    ClusterObjectiveState objective_state{};
    ObjectiveAttemptDiagnostic diagnostic{};
    PolishProgress polish_progress{};
    bool grow_trust_region{ false };
};

} // namespace

TrustRegionStateSet::TrustRegionStateSet(TrustRegionOptions options)
    : m_options{ options }
{
    if (!std::isfinite(m_options.initial_radius) ||
        !std::isfinite(m_options.minimum_radius) ||
        !std::isfinite(m_options.maximum_radius) ||
        m_options.minimum_radius <= 0.0 ||
        m_options.initial_radius < m_options.minimum_radius ||
        m_options.maximum_radius < m_options.initial_radius ||
        !std::isfinite(m_options.shrink_factor) ||
        m_options.shrink_factor <= 0.0 ||
        m_options.shrink_factor >= 1.0 ||
        !std::isfinite(m_options.growth_factor) ||
        m_options.growth_factor <= 1.0)
    {
        throw std::invalid_argument(
            "Local fitting trust-region options are invalid.");
    }
}

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
                m_options.initial_radius : iter->second);
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

TrustRegionRadiusUpdate TrustRegionStateSet::Shrink(
    const std::vector<ClusterKey> & key_list)
{
    TrustRegionRadiusUpdate update;
    for (const auto & key : key_list)
    {
        auto iter{ m_radius_by_key.find(key) };
        if (iter == m_radius_by_key.end())
        {
            throw std::invalid_argument(
                "Local fitting trust-region state is missing.");
        }
        if (iter->second <= m_options.minimum_radius)
        {
            update.saturated_key_list.emplace_back(key);
            continue;
        }
        iter->second = std::max(
            m_options.minimum_radius,
            iter->second * m_options.shrink_factor);
        update.changed_key_list.emplace_back(key);
    }
    return update;
}

void TrustRegionStateSet::Grow(const std::vector<ClusterKey> & key_list)
{
    for (const auto & key : key_list)
    {
        auto iter{ m_radius_by_key.find(key) };
        if (iter == m_radius_by_key.end())
        {
            throw std::invalid_argument(
                "Local fitting trust-region state is missing.");
        }
        iter->second = std::min(
            m_options.maximum_radius,
            iter->second * m_options.growth_factor);
    }
}

TrustRegionIterationUpdate TrustRegionStateSet::UpdateAfterIteration(
    const std::vector<ClusterKey> & grow_key_list,
    const std::vector<ClusterKey> & rejected_key_list,
    const std::vector<ClusterKey> & backtracking_exhausted_key_list)
{
    Grow(grow_key_list);
    auto rejected_cluster_partition{
        PartitionRejectedClusters(
            rejected_key_list,
            backtracking_exhausted_key_list)
    };
    auto radius_update{
        Shrink(rejected_cluster_partition.retryable_key_list)
    };
    return TrustRegionIterationUpdate{
        std::move(rejected_cluster_partition),
        std::move(radius_update)
    };
}

PerformanceCounters::PerformanceCounters(
    bool quiet_mode,
    const SecondStageContext & context,
    const ClusterSolverWorkspaceMap & solver_workspace_by_key)
    : m_quiet_mode{ quiet_mode },
      m_solver_workspace_by_key{ solver_workspace_by_key },
      m_start_time{ std::chrono::steady_clock::now() },
      m_cached_sample_count{ CountRawSamplingEntries(context) }
{
}

PerformanceCounters::~PerformanceCounters()
{
    if (m_quiet_mode) return;

    const auto symbolic_analysis_count{
        m_retired_solver_symbolic_analysis_count +
        CountCurrentSolverSymbolicAnalyses()
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

void PerformanceCounters::RecordFullStateMaterialization()
{
    m_full_state_materialization_count.fetch_add(
        1,
        std::memory_order_relaxed);
}

void PerformanceCounters::RecordGaussianCacheMisses()
{
    m_gaussian_cache_miss_count.fetch_add(
        m_cached_sample_count,
        std::memory_order_relaxed);
}

void PerformanceCounters::RecordGaussianCacheHits()
{
    m_gaussian_cache_hit_count.fetch_add(
        m_cached_sample_count,
        std::memory_order_relaxed);
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

std::chrono::steady_clock::time_point
PerformanceCounters::StartIterationPhase() const
{
    return std::chrono::steady_clock::now();
}

void PerformanceCounters::FinishIterationPhase(
    std::chrono::steady_clock::time_point start_time)
{
    m_iteration_phase_milliseconds +=
        CalculateElapsedMilliseconds(start_time);
}

std::chrono::steady_clock::time_point
PerformanceCounters::StartCandidatePhase() const
{
    return std::chrono::steady_clock::now();
}

void PerformanceCounters::FinishCandidatePhase(
    std::chrono::steady_clock::time_point start_time)
{
    m_candidate_phase_milliseconds +=
        CalculateElapsedMilliseconds(start_time);
}

void PerformanceCounters::RecordSolverWorkspaceReset()
{
    m_retired_solver_symbolic_analysis_count +=
        CountCurrentSolverSymbolicAnalyses();
}

double PerformanceCounters::CalculateElapsedMilliseconds(
    std::chrono::steady_clock::time_point start_time)
{
    return std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start_time).count();
}

std::size_t PerformanceCounters::CountRawSamplingEntries(
    const SecondStageContext & context)
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
    return count;
}

CandidateEvaluationOverlay::CandidateEvaluationOverlay(
    const SecondStageContext & context,
    const ResidualBaseline & baseline,
    const FitStateView & candidate_state)
    : m_context{ context },
      m_baseline{ baseline },
      m_candidate_state{ candidate_state },
      m_changed_group_mask(
          context.selected_atom_index_list_by_group.size(),
          0),
      m_changed_group_median(
          context.selected_atom_index_list_by_group.size())
{
    for (const auto atom_index :
        m_candidate_state.GetOverrideAtomIndexList())
    {
        m_changed_group_mask.at(m_context.at(atom_index).group_id) = 1;
    }
    std::vector<GaussianModel3D> model_list;
    for (std::size_t group_id = 0;
        group_id < m_changed_group_mask.size();
        group_id++)
    {
        if (m_changed_group_mask.at(group_id) == 0) continue;
        const auto & atom_index_list{
            m_context.selected_atom_index_list_by_group.at(group_id)
        };
        model_list.clear();
        model_list.reserve(atom_index_list.size());
        for (const auto atom_index : atom_index_list)
        {
            model_list.emplace_back(
                m_candidate_state.GetModel(atom_index));
        }
        m_changed_group_median.at(group_id) =
            BuildGaussianParameterMedian(model_list);
    }
}

std::optional<ResidualSample> CandidateEvaluationOverlay::operator()(
    const SampleRef & sample_ref) const
{
    const auto & baseline{
        m_baseline.sample_list.at(sample_ref.atom_index)
            .at(sample_ref.sample_index)
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
            if (m_candidate_state.FindOverride(
                    neighbor_atom_sample.atom_index) == nullptr)
            {
                continue;
            }
            baseline_model = &GetFitModel(
                m_baseline.model_snapshot.selected,
                neighbor_atom_sample.atom_index);
            candidate_model = &m_candidate_state.GetModel(
                neighbor_atom_sample.atom_index);
        }
        else
        {
            const auto & unselected_atom_contributor{
                m_context.unselected_atom_list.at(
                    neighbor_atom_sample.atom_index)
            };
            if (!unselected_atom_contributor.selected_group_id.has_value() ||
                m_changed_group_mask.at(
                    *unselected_atom_contributor.selected_group_id) == 0)
            {
                continue;
            }
            baseline_model = &GetFitModel(
                m_baseline.model_snapshot.unselected,
                neighbor_atom_sample.atom_index);
            const auto & median{
                m_changed_group_median.at(
                    *unselected_atom_contributor.selected_group_id)
            };
            if (median.has_value())
            {
                candidate_model = &*median;
            }
            else
            {
                candidate_model =
                    &unselected_atom_contributor.initial_seed.GetModel();
            }
        }
        adjusted_response +=
            baseline_model->ResponseAtDistance(
                neighbor_atom_sample.distance) -
            candidate_model->ResponseAtDistance(
                neighbor_atom_sample.distance);
    }
    const auto expected_response{
        m_candidate_state.FindOverride(sample_ref.atom_index) == nullptr ?
            baseline->adjusted_response - baseline->residual :
            m_candidate_state.GetModel(sample_ref.atom_index).ResponseAtDistance(
                static_cast<double>(sample.point.distance))
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
        if (HasMaterialChange(i))
        {
            provenance.at(m_candidate_patch.atom_index_list.at(i)) =
                endpoint_provenance.at(m_candidate_patch.atom_index_list.at(i));
        }
    }
    return provenance;
}

PolishProvenance
BacktrackingWorkspace::BuildActiveCandidatePolishProvenance(
    const PolishProvenance & previous_provenance,
    const PolishProvenance & endpoint_provenance) const
{
    if (previous_provenance.size() != m_candidate_patch.atom_index_list.size() ||
        endpoint_provenance.size() != m_candidate_patch.atom_index_list.size())
    {
        throw std::invalid_argument(
            "Local fitting active backtracking provenance sizes are inconsistent.");
    }
    auto provenance{ previous_provenance };
    for (std::size_t i = 0; i < m_candidate_patch.atom_index_list.size(); i++)
    {
        if (HasMaterialChange(i))
        {
            provenance.at(i) = endpoint_provenance.at(i);
        }
    }
    return provenance;
}

bool BacktrackingWorkspace::BuildCandidate(double factor)
{
    std::vector<GaussianModel3D> candidate_model_list;
    if (!TryBuildSharedOffsetDampedModelList(
            m_previous_model_list,
            m_endpoint_model_list,
            m_previous_shared_offset_list,
            m_endpoint_shared_offset_list,
            factor,
            candidate_model_list))
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
                candidate_model_list.at(i),
                endpoint_uncertainty
            };
    }
    return true;
}

bool BacktrackingWorkspace::HasMaterialChange(std::size_t atom_position) const
{
    const auto change{
        CalculateTransformedChange(
            m_candidate_patch.mdpde_list.at(atom_position).GetModel(),
            m_previous_model_list.at(atom_position))
    };
    return IsTransformedChangeMaterial(change, m_minimum_transformed_change);
}

double BacktrackingWorkspace::GetMaximumTransformedChange() const
{
    double maximum_change{ 0.0 };
    for (std::size_t i = 0; i < m_candidate_patch.atom_index_list.size(); i++)
    {
        maximum_change = std::max(
            maximum_change,
            detail::GetMaximumTransformedChange(
                CalculateTransformedChange(
                    m_candidate_patch.mdpde_list.at(i).GetModel(),
                    m_previous_model_list.at(i)).value_list));
    }
    return maximum_change;
}

RejectedClusterPartition PartitionRejectedClusters(
    const std::vector<ClusterKey> & rejected_key_list,
    const std::vector<ClusterKey> & exhausted_key_list)
{
    RejectedClusterPartition partition;
    for (const auto & key : rejected_key_list)
    {
        if (std::ranges::find(exhausted_key_list, key) != exhausted_key_list.end())
        {
            partition.exhausted_key_list.emplace_back(key);
        }
        else
        {
            partition.retryable_key_list.emplace_back(key);
        }
    }
    return partition;
}

AllRejectedResolution ResolveAllRejected(
    bool maximum_iterations_reached,
    const RejectedClusterPartition & partition,
    const TrustRegionRadiusUpdate & radius_update)
{
    if (maximum_iterations_reached)
    {
        return AllRejectedResolution::MaximumIterations;
    }
    if (partition.exhausted_key_list.empty() && partition.retryable_key_list.empty())
    {
        throw std::invalid_argument(
            "All-rejected local fitting resolution requires rejected clusters.");
    }
    if (partition.retryable_key_list.empty())
    {
        return AllRejectedResolution::BacktrackingExhausted;
    }
    if (!radius_update.changed_key_list.empty())
    {
        return AllRejectedResolution::Retry;
    }

    const auto all_retryable_saturated{
        std::all_of(
            partition.retryable_key_list.begin(),
            partition.retryable_key_list.end(),
            [&](const ClusterKey & key)
            {
                return std::ranges::find(radius_update.saturated_key_list, key) !=
                    radius_update.saturated_key_list.end();
            })
    };
    if (!all_retryable_saturated)
    {
        throw std::logic_error(
            "Trust-region shrink did not classify every retryable rejection.");
    }
    return partition.exhausted_key_list.empty() ?
        AllRejectedResolution::MinimumRadius : AllRejectedResolution::NoRetryProgress;
}

std::size_t CountSuspiciousAtoms(const SuspiciousUpdateMask & suspicious_mask)
{
    return static_cast<std::size_t>(std::ranges::count_if(suspicious_mask, std::identity{}));
}

static bool HasSuspiciousAtom(
    const std::vector<std::size_t> & atom_index_list,
    const SuspiciousUpdateMask & suspicious_mask)
{
    for (const auto atom_index : atom_index_list)
    {
        if (suspicious_mask.at(atom_index) != 0) return true;
    }
    return false;
}

std::vector<std::size_t> CollectSuspiciousAtomIndices(
    const std::vector<std::size_t> & atom_index_list,
    const SuspiciousUpdateMask & suspicious_mask)
{
    std::vector<std::size_t> suspicious_atom_index_list;
    suspicious_atom_index_list.reserve(atom_index_list.size());
    for (const auto atom_index : atom_index_list)
    {
        if (suspicious_mask.at(atom_index) != 0)
        {
            suspicious_atom_index_list.emplace_back(atom_index);
        }
    }
    return suspicious_atom_index_list;
}

std::vector<double> BuildSuspiciousJointOffsetRidgeMultiplierList(const SuspiciousUpdateMask & suspicious_mask)
{
    std::vector<double> ridge_multiplier_list(suspicious_mask.size(), 1.0);
    for (std::size_t atom_index = 0; atom_index < suspicious_mask.size(); atom_index++)
    {
        if (suspicious_mask.at(atom_index) != 0)
        {
            ridge_multiplier_list.at(atom_index) = kSuspiciousJointOffsetRidgeMultiplier;
        }
    }
    return ridge_multiplier_list;
}

void ClearSuspiciousUpdateMaskForClusters(
    const std::vector<std::vector<std::size_t>> & cluster_key_list,
    SuspiciousUpdateMask & suspicious_mask)
{
    for (const auto & key : cluster_key_list)
    {
        for (const auto atom_index : key)
        {
            suspicious_mask.at(atom_index) = 0;
        }
    }
}

static bool HasSuspiciousCenterSignFlip(
    double previous_innermost_response,
    double candidate_innermost_response,
    double previous_residual_scale)
{
    if (!std::isfinite(previous_innermost_response) ||
        !std::isfinite(candidate_innermost_response) ||
        !std::isfinite(previous_residual_scale) ||
        previous_residual_scale < 0.0)
    {
        return false;
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
    return previous_innermost_response > noise_threshold &&
        candidate_innermost_response < -negative_threshold;
}

static double CalculateZeroOffsetResponse(
    const LocalPotentialSample & sample,
    const GaussianModel3D & model)
{
    const auto distance{ static_cast<double>(sample.point.distance) };
    const auto model_offset{ model.ResponseAtDistance(distance) - model.SignalAtDistance(distance) };
    return static_cast<double>(sample.response) - model_offset;
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
        const auto distance{ static_cast<double>(sample.point.distance) };
        const auto response{ CalculateZeroOffsetResponse(sample, model) };
        if (!std::isfinite(response) ||
            std::abs(response) > static_cast<double>(std::numeric_limits<float>::max()))
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
    ZeroOffsetProfileDiagnostics diagnostics;
    std::ranges::sort(
        profile_samples,
        {},
        &std::pair<double, double>::first);

    diagnostics.distance_range = profile_samples.back().first - profile_samples.front().first;
    for (std::size_t i = 0; i < profile_samples.size();)
    {
        const auto radius{ profile_samples.at(i).first };
        std::vector<double> response_list;
        while (i < profile_samples.size() && IsSameSuspiciousProfileRadius(profile_samples.at(i).first, radius))
        {
            const auto response{ profile_samples.at(i).second };
            diagnostics.max_abs_response = std::max(diagnostics.max_abs_response, std::abs(response));
            response_list.emplace_back(response);
            i++;
        }
        diagnostics.radius_response_median_list.emplace_back(array_helper::ComputeMedian(response_list));
    }
    diagnostics.innermost_response = diagnostics.radius_response_median_list.front();
    if (calculate_residual_scale)
    {
        diagnostics.robust_residual_scale = array_helper::ComputeMedianAbsoluteDeviationScale(residual_list);
    }
    analysis.profile = std::move(diagnostics);
    return analysis;
}

static bool HasUsableSuspiciousProfileBaseline(
    const GaussianModel3D & previous_model,
    const ZeroOffsetProfileDiagnostics & previous_profile)
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
        std::max(std::abs(previous_profile.innermost_response), kSuspiciousProfileScaleMin)
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

static bool HasSuspiciousOffsetMagnitude(
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
    if (!std::isfinite(previous_offset_response) || !std::isfinite(candidate_offset_response))
    {
        return true;
    }
    const auto reference_scale{
        std::max({
            std::abs(previous_model.SignalAtDistance(0.0)),
            std::abs(previous_offset_response),
            previous_profile_max_abs_response,
            kSuspiciousProfileScaleMin
        })
    };
    return std::abs(candidate_offset_response) > kSuspiciousCompensationResponseRatio * reference_scale;
}

static bool HasSuspiciousRadialRebound(
    const ZeroOffsetProfileDiagnostics & previous_profile,
    const ZeroOffsetProfileDiagnostics & candidate_profile)
{
    if (candidate_profile.radius_response_median_list.size() < kSuspiciousProfileMinimumRadiusCount)
    {
        return false;
    }
    const auto reference_innermost_scale{
        std::abs(previous_profile.innermost_response)
    };
    const auto noise_threshold{
        std::max(
            kSuspiciousProfileNoiseScaleMultiplier * previous_profile.robust_residual_scale,
            kSuspiciousProfileScaleMin)
    };
    const auto rebound_magnitude_threshold{
        std::max(
            kSuspiciousProfileReboundReferenceRatio * reference_innermost_scale,
            noise_threshold)
    };
    const auto upward_excursion_threshold{
        std::max(
            kSuspiciousProfileUpwardExcursionReferenceRatio * reference_innermost_scale,
            noise_threshold)
    };
    const auto candidate_innermost_scale{
        std::max(
            std::abs(candidate_profile.innermost_response),
            kSuspiciousProfileScaleMin)
    };
    int upward_excursion_count{ 0 };
    auto previous_abs_response{ std::abs(candidate_profile.radius_response_median_list.front()) };
    for (std::size_t i = 1; i < candidate_profile.radius_response_median_list.size(); i++)
    {
        const auto current_abs_response{
            std::abs(candidate_profile.radius_response_median_list.at(i))
        };
        if (current_abs_response > kSuspiciousProfileReboundCenterRatio * candidate_innermost_scale &&
            current_abs_response > rebound_magnitude_threshold)
        {
            return true;
        }
        if (current_abs_response > previous_abs_response + upward_excursion_threshold)
        {
            upward_excursion_count++;
        }
        previous_abs_response = current_abs_response;
    }
    return upward_excursion_count > kSuspiciousProfileMaximumUpwardExcursions;
}

static bool HasSuspiciousWidthGrowth(
    const GaussianModel3D & previous_model,
    const GaussianModel3D & candidate_model,
    const std::optional<ZeroOffsetProfileDiagnostics> & previous_profile)
{
    if (!std::isfinite(candidate_model.GetWidth()) || candidate_model.GetWidth() <= 0.0) return true;
    if (candidate_model.GetWidth() > kSuspiciousWidthGrowthLimit * previous_model.GetWidth()) return true;
    if (!previous_profile.has_value()) return false;
    const auto distance_range{ previous_profile->distance_range };
    return distance_range > 0.0 && candidate_model.GetWidth() > kSuspiciousWidthRangeLimitRatio * distance_range;
}

static bool HasSuspiciousAmplitudeOffsetCompensation(
    const GaussianModel3D & previous_model,
    const GaussianModel3D & candidate_model,
    const std::optional<ZeroOffsetProfileDiagnostics> & previous_profile)
{
    const auto signal_delta{
        candidate_model.SignalAtDistance(0.0) - previous_model.SignalAtDistance(0.0)
    };
    const auto offset_delta_response{
        candidate_model.GetOffset() * candidate_model.OffsetBasisAtDistance(0.0) -
            previous_model.GetOffset() * previous_model.OffsetBasisAtDistance(0.0)
    };
    const auto reference_scale{
        std::max({
            previous_profile.has_value() ? std::abs(previous_profile->innermost_response) : 0.0,
            std::abs(previous_model.SignalAtDistance(0.0)),
            kSuspiciousProfileScaleMin
        })
    };
    return signal_delta * offset_delta_response < 0.0 &&
        std::abs(signal_delta) > kSuspiciousCompensationResponseRatio * reference_scale &&
        std::abs(offset_delta_response) > kSuspiciousCompensationResponseRatio * reference_scale;
}

SuspiciousGaussianReason EvaluateSuspiciousGaussianUpdate(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & candidate_model,
    const FitOptions & options,
    const SuspiciousUpdateBaseline & previous_baseline,
    SuspiciousUpdateMode mode)
{
    const auto & previous_model{ previous_baseline.previous_model };
    const auto & previous_analysis{ previous_baseline.previous_analysis };
    if (mode == SuspiciousUpdateMode::PostRefit && !IsValidSecondStageGaussianModel(candidate_model))
    {
        return SuspiciousGaussianReason::InvalidModel;
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
        return SuspiciousGaussianReason::NonFiniteResponse;
    }
    if (HasSuspiciousOffsetMagnitude(
            previous_model,
            candidate_model,
            previous_analysis.profile.has_value() ?
                previous_analysis.profile->max_abs_response : 0.0))
    {
        return SuspiciousGaussianReason::OffsetMagnitude;
    }
    const auto has_usable_radial_baseline{
        previous_analysis.all_responses_finite &&
        previous_analysis.profile.has_value() &&
        HasUsableSuspiciousProfileBaseline(previous_model, *previous_analysis.profile)
    };
    if (has_usable_radial_baseline)
    {
        if (!candidate_analysis.profile.has_value())
        {
            return SuspiciousGaussianReason::NonFiniteResponse;
        }
        if (HasSuspiciousCenterSignFlip(
                previous_analysis.profile->innermost_response,
                candidate_analysis.profile->innermost_response,
                previous_analysis.profile->robust_residual_scale))
        {
            return SuspiciousGaussianReason::CenterSignFlip;
        }
        if (HasSuspiciousRadialRebound(
                *previous_analysis.profile,
                *candidate_analysis.profile))
        {
            return SuspiciousGaussianReason::RadialRebound;
        }
    }
    if (mode == SuspiciousUpdateMode::OffsetOnly)
    {
        return SuspiciousGaussianReason::None;
    }
    if (HasSuspiciousWidthGrowth(previous_model, candidate_model, previous_analysis.profile))
    {
        return SuspiciousGaussianReason::WidthGrowth;
    }
    if (HasSuspiciousAmplitudeOffsetCompensation(previous_model, candidate_model, previous_analysis.profile))
    {
        return SuspiciousGaussianReason::AmplitudeOffsetCompensation;
    }
    return SuspiciousGaussianReason::None;
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

SuspiciousGaussianReason EvaluateSuspiciousOffsetUpdate(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & previous_model,
    const GaussianModel3D & candidate_model,
    const FitOptions & options)
{
    const auto previous_baseline{
        BuildPreviousSuspiciousProfileBaseline(sample_entries, previous_model, options)
    };
    return EvaluateSuspiciousGaussianUpdate(
        sample_entries,
        candidate_model,
        options,
        previous_baseline,
        SuspiciousUpdateMode::OffsetOnly);
}

SuspiciousUpdateMask ExpandSuspiciousSharedOffsetGroups(
    const std::vector<std::size_t> & group_id_by_position,
    const SuspiciousUpdateMask & suspicious_seed_mask)
{
    if (group_id_by_position.size() != suspicious_seed_mask.size())
    {
        throw std::invalid_argument("Suspicious shared-offset group input sizes are inconsistent.");
    }

    std::set<std::size_t> suspicious_seed_group_id_set;
    for (std::size_t position = 0; position < group_id_by_position.size(); position++)
    {
        if (suspicious_seed_mask.at(position) == 0) continue;
        suspicious_seed_group_id_set.emplace(group_id_by_position.at(position));
    }

    SuspiciousUpdateMask rollback_mask(group_id_by_position.size(), 0);
    for (std::size_t position = 0; position < group_id_by_position.size(); position++)
    {
        if (suspicious_seed_group_id_set.contains(group_id_by_position.at(position)))
        {
            rollback_mask.at(position) = 1;
        }
    }
    return rollback_mask;
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
                    EvaluateResidualSample(context, model_snapshot.selected, sample_ref, model_snapshot)
                };
                const auto distance{
                    static_cast<double>(raw_sampling_entries.at(sample_index).point.distance)
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

namespace {

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

} // namespace

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

CombinedCandidateObjectiveCheck EvaluateCombinedCandidateObjective(
    const SecondStageContext & context,
    const ResidualBaseline & baseline,
    const CouplingGraphPartition & partition,
    const FitState & previous_state,
    const FitState & candidate_state,
    const std::vector<ClusterKey> & accepted_key_list,
    const ObjectiveDomain & objective_domain,
    const ObjectiveBreakdown * best_objective,
    PerformanceCounters & performance_counters)
{
    CombinedCandidateObjectiveCheck result;
    if (partition.boundary_sample_count == 0 || accepted_key_list.empty()) return result;

    result.previous_objective = EvaluateAuditObjective(objective_domain, baseline);

    ClusterKey changed_atom_index_list;
    for (const auto & key : accepted_key_list)
    {
        changed_atom_index_list.insert(changed_atom_index_list.end(), key.begin(), key.end());
    }
    const auto combined_patch{
        FitStatePatch::FromState(candidate_state, std::move(changed_atom_index_list))
    };
    const FitStateView combined_state_view{ previous_state, combined_patch };
    const CandidateEvaluationOverlay combined_overlay{ context, baseline, combined_state_view };
    const auto affected_sample_ref_list{
        BuildGraphAffectedSampleUnion(partition, accepted_key_list)
    };
    const auto combined_check{
        EvaluateCombinedObjective(
            combined_overlay,
            affected_sample_ref_list,
            objective_domain,
            best_objective,
            result.previous_objective.has_value() ? &*result.previous_objective : nullptr,
            performance_counters)
    };
    result.accepted = combined_check.has_value();
    result.candidate_objective = combined_check;
    return result;
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
    const auto maximum_transformed_change{
        GetMaximumTransformedChange(transformed_change_summary)
    };
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

static BaseProposalBuildResult BuildSharedOffsetBaseProposal(
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
            PreObjectiveFailureReason::InvalidModel,
            std::nullopt
        };
    }

    std::vector<std::size_t> group_id_by_atom_position;
    std::vector<GaussianModel3D> previous_model_list;
    std::vector<GaussianModel3D> raw_model_list;
    group_id_by_atom_position.reserve(key.size());
    previous_model_list.reserve(key.size());
    raw_model_list.reserve(key.size());
    for (const auto atom_index : key)
    {
        group_id_by_atom_position.emplace_back(context.at(atom_index).group_id);
        previous_model_list.emplace_back(outer_previous_state.at(atom_index).mdpde.GetModel());
        raw_model_list.emplace_back(raw_state.at(atom_index).mdpde.GetModel());
    }
    const auto previous_shared_offset_list{
        BuildGroupMedianOffsetList(group_id_by_atom_position, previous_model_list)
    };
    const auto raw_shared_offset_list{
        BuildGroupMedianOffsetList(group_id_by_atom_position, raw_model_list)
    };

    std::vector<GaussianModel3D> seed_model_list;
    if (!TryBuildSharedOffsetDampedModelList(
            previous_model_list,
            raw_model_list,
            previous_shared_offset_list,
            raw_shared_offset_list,
            0.0,
            seed_model_list))
    {
        return BaseProposalBuildResult{
            std::nullopt,
            PreObjectiveFailureReason::InvalidModel,
            std::nullopt
        };
    }
    const auto seed_step_norm{
        CalculateModelTrustRegionStepNorm(previous_model_list, seed_model_list)
    };
    if (!seed_step_norm.has_value())
    {
        return BaseProposalBuildResult{
            std::nullopt,
            PreObjectiveFailureReason::InvalidModel,
            std::nullopt
        };
    }
    if (!IsTrustRegionStepWithinRadius(*seed_step_norm, trust_region_radius))
    {
        return BaseProposalBuildResult{
            std::nullopt,
            PreObjectiveFailureReason::PreviousSharedOffsetProjectionOutsideTrustRegion,
            *seed_step_norm
        };
    }

    double damping{ 1.0 };
    std::optional<double> attempted_step_norm;
    std::vector<GaussianModel3D> candidate_model_list;
    while (damping >= std::numeric_limits<double>::epsilon())
    {
        if (TryBuildSharedOffsetDampedModelList(
                previous_model_list,
                raw_model_list,
                previous_shared_offset_list,
                raw_shared_offset_list,
                damping,
                candidate_model_list))
        {
            const auto step_norm{
                CalculateModelTrustRegionStepNorm(previous_model_list, candidate_model_list)
            };
            if (step_norm.has_value() &&
                IsTrustRegionStepWithinRadius(*step_norm, trust_region_radius))
            {
                FitStateProposal proposal{
                    .patch{ .atom_index_list = key },
                    .effective_damping = damping,
                    .step_norm = *step_norm
                };
                proposal.patch.mdpde_list.reserve(key.size());
                for (std::size_t atom_position = 0; atom_position < key.size(); atom_position++)
                {
                    const auto atom_index{ key.at(atom_position) };
                    proposal.patch.mdpde_list.emplace_back(
                        GaussianModel3DWithUncertainty{
                            candidate_model_list.at(atom_position),
                            raw_state.at(atom_index).mdpde
                                .GetStandardDeviationModel()
                        });
                }
                return BaseProposalBuildResult{
                    std::move(proposal),
                    PreObjectiveFailureReason::None,
                    std::nullopt
                };
            }
            if (step_norm.has_value()) attempted_step_norm = *step_norm;
        }
        damping *= 0.5;
    }
    return BaseProposalBuildResult{
        std::nullopt,
        PreObjectiveFailureReason::NoCandidateWithinTrustRegion,
        attempted_step_norm
    };
}

void RejectCombinedCandidate(
    const FitState & previous_state,
    const PolishProvenance & previous_polish_provenance,
    CandidateSelection & selection)
{
    selection.assembled_state = previous_state;
    selection.assembled_polish_provenance = previous_polish_provenance;
    selection.rejected_key_list.insert(
        selection.rejected_key_list.end(),
        selection.accepted_key_list.begin(),
        selection.accepted_key_list.end());
    std::ranges::sort(selection.rejected_key_list);
    selection.rejected_key_list.erase(
        std::ranges::unique(selection.rejected_key_list).begin(),
        selection.rejected_key_list.end());
    selection.accepted_key_list.clear();
    selection.grow_trust_region_key_list.clear();
    selection.combined_backtracking_objective.reset();
    selection.polish_progress.rejected_count += selection.polish_progress.accepted_count;
    selection.polish_progress.accepted_count = 0;
}

static std::optional<FitStatePatch> BuildDampedCandidatePatch(
    const FitState & previous_state,
    const std::vector<Eigen::Vector3d> & previous_transformed_estimation_list,
    const std::vector<Eigen::Vector3d> & endpoint_transformed_estimation_list,
    const FitState & uncertainty_state,
    const std::vector<std::size_t> & active_index_list,
    double damping)
{
    if (previous_transformed_estimation_list.size() != active_index_list.size() ||
        endpoint_transformed_estimation_list.size() != active_index_list.size() ||
        !std::isfinite(damping) || damping < 0.0 || damping > 1.0)
    {
        throw std::invalid_argument(
            "Local fitting candidate transformed interpolation inputs are invalid.");
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
        const auto candidate_transformed_estimation{
            (previous_transformed_estimation +
                damping * (endpoint_transformed_estimation_list.at(local_position) -
                    previous_transformed_estimation)).eval()
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
            kObjectiveStrictTolerance);
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
    const auto & raw_state{ inputs.raw_state };
    const auto & ridge_multiplier_list{ inputs.ridge_multiplier_list };
    const auto & objective_domain{ inputs.objective_domain };
    const auto & previous_objective_entry{ inputs.previous_objective_by_key.at(key) };
    const auto * previous_objective{
        previous_objective_entry.has_value() ? &*previous_objective_entry : nullptr
    };
    const auto & previous_objective_state{ inputs.cluster_objective_state.at(key) };
    const auto trust_region_radius{ inputs.trust_region_state.GetRadius(key) };
    const auto contains_suspicious_atom{
        HasSuspiciousAtom(key, inputs.rollback_atom_mask)
    };
    const auto is_polish_eligible{
        inputs.health_by_key.at(key).IsStationarityEligible() && !contains_suspicious_atom
    };
    const auto is_unchanged_state_exhausted{
        std::ranges::find(inputs.unchanged_state_exhausted_key_list, key) !=
            inputs.unchanged_state_exhausted_key_list.end()
    };
    auto & performance_counters{ inputs.performance_counters };
    ClusterCandidateResult result;
    result.objective_state = previous_objective_state;
    result.polish_provenance.reserve(key.size());
    for (const auto atom_index : key)
    {
        result.polish_provenance.emplace_back(previous_polish_provenance.at(atom_index));
    }
    if (is_polish_eligible) result.polish_progress.eligible_count = 1;
    result.diagnostic.trust_region_radius = trust_region_radius;
    if (is_unchanged_state_exhausted)
    {
        if (is_polish_eligible) result.polish_progress.skipped_count = 1;
        result.diagnostic.backtracking_exhausted = true;
        return result;
    }

    std::optional<FitStateProposal> base_proposal;
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
        result.diagnostic.pre_objective_failure_reason = proposal_result.failure_reason;
        if (proposal_result.failure_reason != PreObjectiveFailureReason::None &&
            proposal_result.attempted_step_norm.has_value())
        {
            result.diagnostic.pre_objective_attempted_step_norm = proposal_result.attempted_step_norm;
            result.diagnostic.trust_region_step_norm = *proposal_result.attempted_step_norm;
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
            const auto previous_estimation{
                EncodeTransformedCoordinates(previous_state.at(atom_index).mdpde.GetModel())
            };
            const auto raw_estimation{
                EncodeTransformedCoordinates(raw_state.at(atom_index).mdpde.GetModel())
            };
            if (!previous_estimation.has_value() || !raw_estimation.has_value())
            {
                throw std::invalid_argument(
                    "Local fitting state has invalid transformed coordinates.");
            }
            previous_cluster_estimation_list.emplace_back(*previous_estimation);
            raw_cluster_estimation_list.emplace_back(*raw_estimation);
        }
        const auto trust_region_damping{
            LimitTrustRegionDamping(
                previous_cluster_estimation_list,
                raw_cluster_estimation_list,
                1.0,
                trust_region_radius)
        };
        auto base_patch{
            BuildDampedCandidatePatch(
                previous_state,
                previous_cluster_estimation_list,
                raw_cluster_estimation_list,
                raw_state,
                key,
                trust_region_damping.effective_damping)
        };
        if (base_patch.has_value())
        {
            base_proposal = FitStateProposal{
                std::move(*base_patch),
                trust_region_damping.effective_damping,
                trust_region_damping.step_norm
            };
        }
        else
        {
            result.diagnostic.pre_objective_failure_reason = PreObjectiveFailureReason::InvalidModel;
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
        base_patch,
        kTransformedChangeTolerance
    };
    CandidateEvaluationOverlay base_overlay{
        context,
        residual_baseline,
        base_state_view
    };
    auto accepted_base_candidate{
        TryCommitClusterCandidate(
            base_overlay,
            key,
            objective_sample_ref_list,
            previous_objective,
            false,
            objective_domain,
            result.objective_state,
            result.diagnostic,
            performance_counters)
    };
    auto accepted_by_backtracking{ false };
    const PolishProvenance non_polished_endpoint_provenance(key.size(), 0);
    if (!accepted_base_candidate)
    {
        result.diagnostic.accepted_backtracking_factor.reset();
        BacktrackingStep step;
        for (step = backtracking_workspace.BuildNextCandidate();
            step.status == BacktrackingStepStatus::CandidateReady;
            step = backtracking_workspace.BuildNextCandidate())
        {
            const auto factor{ step.factor };
            const auto & backtracked_patch{ backtracking_workspace.GetCandidatePatch() };
            ObjectiveAttemptDiagnostic trial_diagnostic;
            trial_diagnostic.effective_damping = base_proposal->effective_damping * factor;
            trial_diagnostic.trust_region_radius = trust_region_radius;
            trial_diagnostic.trust_region_step_norm = base_proposal->step_norm * factor;
            trial_diagnostic.backtracking_trial_count = step.trial_number;
            const FitStateView backtracked_state_view{
                previous_state,
                backtracked_patch
            };
            const CandidateEvaluationOverlay backtracked_overlay{
                context,
                residual_baseline,
                backtracked_state_view
            };
            if (TryCommitClusterCandidate(
                    backtracked_overlay,
                    key,
                    objective_sample_ref_list,
                    previous_objective,
                    false,
                    objective_domain,
                    result.objective_state,
                    trial_diagnostic,
                    performance_counters))
            {
                trial_diagnostic.accepted_backtracking_factor = factor;
                result.diagnostic = std::move(trial_diagnostic);
                accepted_base_candidate = true;
                accepted_by_backtracking = true;
                break;
            }
            result.diagnostic = std::move(trial_diagnostic);
        }
        if (step.status == BacktrackingStepStatus::Exhausted)
        {
            result.diagnostic.backtracking_exhausted = true;
        }
        else if (step.status == BacktrackingStepStatus::InvalidCandidate)
        {
            result.diagnostic.is_invalid_model = true;
        }
        if (accepted_base_candidate)
        {
            result.polish_provenance =
                backtracking_workspace.BuildActiveCandidatePolishProvenance(
                    result.polish_provenance,
                    non_polished_endpoint_provenance);
            base_patch = backtracking_workspace.TakeCandidatePatch();
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
        return result;
    }

    result.grow_trust_region = !accepted_by_backtracking && ShouldGrowTrustRegion(result.diagnostic);
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
            polish_diagnostic.effective_damping = polished_candidate->effective_damping;
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
            if (!TryCommitClusterCandidate(
                    polished_overlay,
                    key,
                    objective_sample_ref_list,
                    result.diagnostic.candidate_objective.has_value() ?
                        &*result.diagnostic.candidate_objective : nullptr,
                    true,
                    objective_domain,
                    result.objective_state,
                    polish_diagnostic,
                    performance_counters))
            {
                result.polish_progress.rejected_count = 1;
            }
            else
            {
                result.polish_progress.accepted_count = 1;
                for (std::size_t position = 0; position < key.size(); position++)
                {
                    const auto base_coordinates{
                        EncodeTransformedCoordinates(base_state_view.GetModel(key.at(position)))
                    };
                    const auto candidate_coordinates{
                        EncodeTransformedCoordinates(polished_candidate->patch.mdpde_list.at(position).GetModel())
                    };
                    if ((base_coordinates->array() != candidate_coordinates->array()).any())
                    {
                        result.polish_provenance.at(position) = 1;
                    }
                }
                result.accepted_patch = std::move(polished_candidate->patch);
                if (!accepted_by_backtracking && ShouldGrowTrustRegion(polish_diagnostic))
                {
                    result.grow_trust_region = true;
                }
            }
        }
    }
    if (!result.accepted_patch.has_value())
    {
        result.accepted_patch = std::move(base_patch);
    }
    return result;
}

CandidateSelection SelectClusterCandidates(const CandidateSelectionInputs & inputs)
{
    const auto & partition{ inputs.partition };
    auto & solver_workspace_by_key{ inputs.solver_workspace_by_key };
    auto & cluster_objective_state{ inputs.cluster_objective_state };
    const auto & previous_state{ inputs.previous_state };
    const auto & previous_polish_provenance{ inputs.previous_polish_provenance };
    using PartitionEntry = std::pair<const ClusterKey, std::vector<SampleRef>>;
    using CandidateWork = std::pair<const PartitionEntry *, ClusterSolverWorkspace *>;
    std::vector<CandidateWork> candidate_work_list;
    candidate_work_list.reserve(partition.sample_id_list_by_key.size());
    for (const auto & entry : partition.sample_id_list_by_key)
    {
        candidate_work_list.emplace_back(&entry, &solver_workspace_by_key.at(entry.first));
    }

    std::vector<ClusterCandidateResult> result_list(candidate_work_list.size());
    std::vector<std::exception_ptr> exception_list(candidate_work_list.size());
    const auto select_candidate = [&](std::size_t position)
    {
        try
        {
            const auto & entry{ *candidate_work_list.at(position).first };
            const auto & key{ entry.first };
            const auto & sample_ref_list{ entry.second };
            result_list.at(position) = SelectClusterCandidate(
                inputs,
                key,
                sample_ref_list,
                *candidate_work_list.at(position).second
            );
        }
        catch (...)
        {
            exception_list.at(position) = std::current_exception();
        }
    };
#ifdef USE_OPENMP
    if (inputs.thread_size > 1 && candidate_work_list.size() > 1)
    {
        eigen_helper::ScopedEigenThreadCount eigen_thread_guard{ 1 };
#pragma omp parallel for schedule(dynamic) num_threads(inputs.thread_size)
        for (std::size_t position = 0; position < candidate_work_list.size(); position++)
        {
            select_candidate(position);
        }
    }
    else
#endif
    {
        for (std::size_t position = 0; position < candidate_work_list.size(); position++)
        {
            select_candidate(position);
        }
    }
    for (const auto & exception : exception_list)
    {
        if (exception) std::rethrow_exception(exception);
    }

    CandidateSelection selection{
        .assembled_state = previous_state,
        .assembled_polish_provenance = previous_polish_provenance
    };
    for (std::size_t position = 0; position < result_list.size(); position++)
    {
        auto & result{ result_list.at(position) };
        const auto & key{ candidate_work_list.at(position).first->first };
        cluster_objective_state.at(key) = std::move(result.objective_state);
        selection.polish_progress.eligible_count += result.polish_progress.eligible_count;
        selection.polish_progress.accepted_count += result.polish_progress.accepted_count;
        selection.polish_progress.rejected_count += result.polish_progress.rejected_count;
        selection.polish_progress.skipped_count += result.polish_progress.skipped_count;
        if (!result.accepted_patch.has_value())
        {
            selection.rejected_key_list.emplace_back(key);
            if (result.diagnostic.backtracking_exhausted)
            {
                selection.backtracking_exhausted_key_list.emplace_back(key);
            }
            selection.rejected_cluster_diagnostic_list.emplace_back(
                ClusterCandidateDiagnostic{
                    key,
                    std::move(result.diagnostic)
                });
            continue;
        }
        selection.accepted_key_list.emplace_back(key);
        selection.accepted_cluster_diagnostic_list.emplace_back(
            ClusterCandidateDiagnostic{
                key,
                std::move(result.diagnostic)
            });
        if (result.grow_trust_region)
        {
            selection.grow_trust_region_key_list.emplace_back(key);
        }
        result.accepted_patch->ApplyTo(selection.assembled_state);
        for (std::size_t key_position = 0; key_position < key.size(); key_position++)
        {
            selection.assembled_polish_provenance.at(key.at(key_position)) =
                result.polish_provenance.at(key_position);
        }
    }
    return selection;
}

bool TryBacktrackCombinedCandidate(
    const CandidateSelectionInputs & inputs,
    const ObjectiveBreakdown * previous_audit_objective,
    const ObjectiveBreakdown * best_audit_objective,
    const ClusterObjectiveStateMap & committed_objective_state,
    CandidateSelection & selection)
{
    const auto & context{ inputs.context };
    const auto & residual_baseline{ inputs.residual_baseline };
    const auto & partition{ inputs.partition };
    const auto & previous_state{ inputs.previous_state };
    const auto & previous_polish_provenance{ inputs.previous_polish_provenance };
    const auto & objective_domain{ inputs.objective_domain };
    const auto & previous_objective_by_key{ inputs.previous_objective_by_key };
    auto & working_objective_state{ inputs.cluster_objective_state };
    auto & performance_counters{ inputs.performance_counters };
    std::vector<std::size_t> changed_atom_index_list;
    for (const auto & key : selection.accepted_key_list)
    {
        changed_atom_index_list.insert(
            changed_atom_index_list.end(),
            key.begin(),
            key.end());
    }
    const auto affected_sample_ref_list{
        BuildGraphAffectedSampleUnion(partition, selection.accepted_key_list)
    };
    const auto endpoint_patch{
        FitStatePatch::FromState(
            selection.assembled_state,
            changed_atom_index_list)
    };
    BacktrackingWorkspace backtracking_workspace{
        context,
        previous_state,
        endpoint_patch,
        kTransformedChangeTolerance
    };
    selection.combined_backtracking_trial_count = 1;
    ClusterObjectiveStateMap accepted_trial_objective_state;
    BacktrackingStep step;
    for (step = backtracking_workspace.BuildNextCandidate();
        step.status == BacktrackingStepStatus::CandidateReady;
        step = backtracking_workspace.BuildNextCandidate())
    {
        const auto factor{ step.factor };
        const auto & candidate_patch{ backtracking_workspace.GetCandidatePatch() };
        const FitStateView candidate_state_view{ previous_state, candidate_patch};
        const CandidateEvaluationOverlay candidate_overlay{
            context,
            residual_baseline,
            candidate_state_view
        };
        selection.combined_backtracking_trial_count = step.trial_number;
        auto trial_objective_state{ committed_objective_state };
        auto local_criteria_accepted{ true };
        for (const auto & key : selection.accepted_key_list)
        {
            ObjectiveAttemptDiagnostic diagnostic;
            diagnostic.backtracking_trial_count = selection.combined_backtracking_trial_count;
            diagnostic.accepted_backtracking_factor = factor;
            const auto & previous_objective{ previous_objective_by_key.at(key) };
            if (!TryCommitClusterCandidate(
                    candidate_overlay,
                    key,
                    partition.sample_id_list_by_key.at(key),
                    previous_objective.has_value() ? &*previous_objective : nullptr,
                    false,
                    objective_domain,
                    trial_objective_state.at(key),
                    diagnostic,
                    performance_counters))
            {
                local_criteria_accepted = false;
                break;
            }
        }
        const auto combined_check{
            local_criteria_accepted ?
                EvaluateCombinedObjective(
                    candidate_overlay,
                    affected_sample_ref_list,
                    objective_domain,
                    best_audit_objective,
                    previous_audit_objective,
                    performance_counters) :
                std::optional<ObjectiveBreakdown>{}
        };
        if (combined_check.has_value())
        {
            selection.combined_backtracking_factor = factor;
            selection.combined_backtracking_objective = combined_check;
            accepted_trial_objective_state = std::move(trial_objective_state);
            break;
        }
    }
    if (step.status == BacktrackingStepStatus::Exhausted)
    {
        selection.combined_backtracking_exhausted = true;
        return false;
    }
    if (step.status == BacktrackingStepStatus::InvalidCandidate)
    {
        return false;
    }

    selection.assembled_state =
        FitStateView{
            previous_state,
            backtracking_workspace.GetCandidatePatch()
        }.Materialize();
    performance_counters.RecordFullStateMaterialization();
    selection.assembled_polish_provenance =
        backtracking_workspace.BuildCandidatePolishProvenance(
            previous_polish_provenance,
            selection.assembled_polish_provenance);
    selection.grow_trust_region_key_list.clear();
    working_objective_state = std::move(accepted_trial_objective_state);
    return true;
}

} // namespace rhbm_gem::core::detail
