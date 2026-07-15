#include <cstddef>

#include "core/detail/GaussianEstimatorStages.hpp"
#include "core/detail/LocalFittingAudit.hpp"
#include "core/detail/LocalFittingCouplingGraph.hpp"
#include "core/detail/LocalFittingHealth.hpp"
#include "core/detail/LocalFittingJointOffsetConditioning.hpp"
#include "core/detail/LocalFittingJointPolish.hpp"
#include "core/detail/LocalFittingSeedRepair.hpp"
#include "core/detail/LocalFittingTrustRegion.hpp"
#include "core/detail/LocalFittingTransformedChange.hpp"
#include "core/detail/PostRefitRollback.hpp"
#include "data/detail/AtomClassifier.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/algorithm/NormalizedChange.hpp>
#include <rhbm_gem/utils/algorithm/ParameterChangeStats.hpp>
#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>
#include <rhbm_gem/utils/algorithm/ScaleReferenceTracker.hpp>
#include <rhbm_gem/utils/algorithm/WeightedRidgeSolver.hpp>
#include <rhbm_gem/utils/algorithm/WeightedRidgeSystem.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem::core {
namespace {
constexpr double kLocalFittingTransformedChangeTolerance{ 1.0e-4 };
constexpr double kLocalFittingTransformedMaximumChangeTolerance{ 1.0e-3 };
constexpr double kNeighborContributionDistanceMax{ 2.5 };
constexpr double kNeighborAtomSearchRange{ 2.0 * kNeighborContributionDistanceMax };
constexpr std::size_t kLocalFittingMaximumIterations{ 50 };
constexpr std::size_t kLocalFittingAuditPatience{ 3 };
constexpr double kLocalFittingChangePercentile{ 0.99 };
constexpr algorithm::RobustLossKind kSecondStageRobustLossKind{ algorithm::RobustLossKind::Cauchy };
constexpr int kRobustLossMaximumIterations{ 50 };
constexpr double kRobustScaleMultiplier{ 1.4826 };
constexpr double kRobustScaleMin{ 1.0e-12 };
constexpr double kRobustLossCutoffMultiplier{ 1.345 };
constexpr double kJointOffsetRidgeRatio{ 1.0e-3 };
constexpr double kSuspiciousJointOffsetRidgeMultiplier{ 10.0 };
constexpr double kJointOffsetCollinearityOverlapThreshold{ 0.98 };
constexpr double kCollinearJointOffsetRidgeMultiplier{ 10.0 };
constexpr double kJointOffsetConditioningPivotRatioThreshold{ 1.0e-8 };
constexpr double kJointOffsetIrlsScaleFloor{ 1.0e-2 };
constexpr double kJointOffsetIrlsNormalizedChangeTolerance{ 1.0e-6 };
constexpr double kJointOffsetIrlsObjectiveRelativeTolerance{ 1.0e-10 };
constexpr double kLocalFittingObjectiveTieRelativeTolerance{ 1.0e-8 };
constexpr double kLocalFittingConvergenceObjectiveRelativeTolerance{ 1.0e-3 };
constexpr double kLocalFittingCouplingMinimumWeight{ 0.05 };
constexpr std::array<double, 6> kLocalFittingCouplingSensitivityMinimumWeightList{
    0.05,
    0.075,
    0.10,
    0.15,
    0.20,
    0.30
};
constexpr std::size_t kLocalFittingObjectiveScaleWarmupCount{ 5 };
constexpr double kLocalFittingObjectiveResidualScaleFloorRatio{ 1.0e-6 };
constexpr double kLocalFittingWidthPriorPenaltyWeight{ 1.0e-2 };
constexpr double kLocalFittingOffsetPlausibilityPenaltyWeight{ 1.0e-2 };
constexpr double kLocalFittingWidthPriorLogScale{ 0.35 };
constexpr std::array<double, detail::kTransformedChangeSize>
kLocalFittingTrustRegionParameterScale{ 0.50, 0.35, 1.0 };
constexpr double kLocalFittingTrustRegionInitialRadius{ 1.0 };
constexpr double kLocalFittingTrustRegionMinimumRadius{ 0.0625 };
constexpr double kLocalFittingTrustRegionMaximumRadius{ 4.0 };
constexpr double kLocalFittingTrustRegionShrinkFactor{ 0.5 };
constexpr double kLocalFittingTrustRegionGrowthFactor{ 2.0 };
constexpr double kLocalFittingTrustRegionBoundaryRatio{ 0.8 };
constexpr double kLocalFittingOffsetPeakRatioMax{ 1.0 };
constexpr std::size_t kSuspiciousOffsetClusterMaxDepth{ 2 };
constexpr double kSuspiciousOffsetClusterMinimumOverlap{ 0.05 };
constexpr std::size_t kSuspiciousProfileMinimumRadiusCount{ 3 };
constexpr double kSuspiciousProfileDistanceTolerance{ 1.0e-6 };
constexpr double kSuspiciousProfileCenterSignFlipRatio{ 0.25 };
constexpr double kSuspiciousProfileReboundCenterRatio{ 1.5 };
constexpr double kSuspiciousProfileReboundReferenceRatio{ 0.25 };
constexpr double kSuspiciousProfileUpwardExcursionReferenceRatio{ 0.20 };
constexpr int kSuspiciousProfileMaximumUpwardExcursions{ 1 };
constexpr double kSuspiciousWidthGrowthLimit{ 1.5 };
constexpr double kSuspiciousWidthRangeLimitRatio{ 1.5 };
constexpr double kSuspiciousCompensationResponseRatio{ 2.0 };
constexpr std::size_t kPersistentTerminalFailureIterationLimit{ 5 };

using LocalFittingState = std::vector<LocalGaussianResult>;

struct ActiveCouplingEdge
{
    std::size_t neighbor_index{ 0 };
    double overlap{ 0.0 };
};

using ActiveCouplingGraph = std::vector<std::vector<ActiveCouplingEdge>>;
using LocalFittingClusterKey = std::vector<std::size_t>;
using detail::JointOffsetSolveStatus;
using detail::IsJointOffsetSolveHardFailure;
using detail::IsJointOffsetSolveStationarityEligible;

const char * GetJointOffsetSolveStatusText(JointOffsetSolveStatus status)
{
    switch (status)
    {
    case JointOffsetSolveStatus::Converged:
        return "converged";
    case JointOffsetSolveStatus::SystemBuildFailed:
        return "system-build-failed";
    case JointOffsetSolveStatus::EmptySystem:
        return "empty-system";
    case JointOffsetSolveStatus::InitialSolveFailed:
        return "initial-solve-failed";
    case JointOffsetSolveStatus::IrlsSolveFailed:
        return "irls-solve-failed";
    case JointOffsetSolveStatus::IrlsObjectiveDeteriorated:
        return "irls-objective-deteriorated";
    case JointOffsetSolveStatus::IrlsMaximumIterationsReached:
        return "irls-maximum-iterations-reached";
    }
    throw std::logic_error("Joint offset solve status is invalid.");
}

struct ZeroOffsetProfileDiagnostics
{
    double distance_min{ 0.0 };
    double distance_max{ 0.0 };
    double center_response{ 0.0 };
    double max_abs_response{ 0.0 };
    std::vector<double> radius_response_median_list{};
};

struct JointOffsetSolveResult
{
    JointOffsetSolveStatus status{ JointOffsetSolveStatus::SystemBuildFailed };
    Eigen::VectorXd offset{};
    ActiveCouplingGraph active_coupling_graph{};
};

struct JointOffsetBuildResult
{
    algorithm::WeightedRidgeSystem system{};
    ActiveCouplingGraph active_coupling_graph{};
};

struct LocalFittingClusterHealth
{
    JointOffsetSolveStatus joint_offset_status{ JointOffsetSolveStatus::SystemBuildFailed };
    std::vector<std::size_t> stationarity_ineligible_refit_atom_index_list{};
};

using LocalFittingClusterHealthMap = std::map<LocalFittingClusterKey, LocalFittingClusterHealth>;

bool IsLocalFittingClusterHealthy(const LocalFittingClusterHealth & health)
{
    return IsJointOffsetSolveStationarityEligible(health.joint_offset_status) &&
        health.stationarity_ineligible_refit_atom_index_list.empty();
}

void RecordLocalRefitStationarityFailure(
    LocalFittingClusterHealthMap & health_by_key,
    std::size_t atom_index)
{
    for (auto & [key, health] : health_by_key)
    {
        if (!std::binary_search(key.begin(), key.end(), atom_index)) continue;
        health.stationarity_ineligible_refit_atom_index_list.emplace_back(atom_index);
        return;
    }
    throw std::invalid_argument("Unhealthy local refit atom is not in an active cluster.");
}

std::vector<LocalFittingClusterKey> CollectUnhealthyLocalFittingClusterKeys(
    const LocalFittingClusterHealthMap & health_by_key)
{
    std::vector<LocalFittingClusterKey> key_list;
    for (const auto & [key, health] : health_by_key)
    {
        if (!IsLocalFittingClusterHealthy(health))
        {
            key_list.emplace_back(key);
        }
    }
    return key_list;
}

struct ClusteredJointOffsetSolveResult
{
    Eigen::VectorXd offset{};
    ActiveCouplingGraph active_coupling_graph{};
    LocalFittingClusterHealthMap health_by_key{};
};

struct LocalFittingJointPolishStep
{
    std::vector<Eigen::VectorXd> transformed_estimation_list{};
    bool is_stationary{ false };
};

struct LocalFittingIterationResult
{
    LocalFittingState state{};
    std::vector<std::size_t> suspicious_offset_state_index_list{};
    LocalFittingClusterHealthMap health_by_key{};
};

struct LocalAtomRefitResult
{
    LocalGaussianResult result{};
    bool is_stationarity_eligible{ false };
};

struct LocalFittingObjectiveAtomModel
{
    std::size_t atom_index{ 0 };
    GaussianModel3D model{};
};

struct LocalFittingObjectiveSamples
{
    std::vector<double> residual_list{};
    std::vector<LocalFittingObjectiveAtomModel> active_model_list{};
    double scale_sample{ 0.0 };
};

struct LocalFittingTrackedCandidate
{
    algorithm::ParameterChangeStats transformed_change_stats{};
    std::optional<LocalFittingObjectiveSamples> objective_samples{};
};

struct LocalFittingClusterObjectiveState
{
    algorithm::ScaleReferenceTracker scale_tracker;
    LocalFittingTrackedCandidate previous{};
    std::optional<LocalFittingTrackedCandidate> best{};

    LocalFittingClusterObjectiveState(
        std::optional<double> initial_scale_sample,
        LocalFittingTrackedCandidate initial_candidate,
        bool has_initial_objective)
        : scale_tracker{
              kLocalFittingObjectiveScaleWarmupCount,
              initial_scale_sample
          },
          previous{ std::move(initial_candidate) }
    {
        if (has_initial_objective) best = previous;
    }
};

using LocalFittingClusterObjectiveStateMap =
    std::map<LocalFittingClusterKey, LocalFittingClusterObjectiveState>;

using LocalFittingObjectiveSampleRef = detail::LocalFittingCouplingSampleId;

struct LocalFittingParameterPenaltyComponents
{
    double width_prior_penalty_sum{ 0.0 };
    double offset_plausibility_penalty_sum{ 0.0 };
};

struct LocalFittingAuditedState
{
    detail::LocalFittingObjectiveBreakdown objective{};
    LocalFittingState state{};
    std::optional<std::size_t> accepted_iteration{};
};

struct LocalFittingBestAuditState
{
    std::vector<LocalFittingObjectiveSampleRef> sample_ref_list{};
    std::vector<std::size_t> atom_index_list{};
    std::optional<double> fixed_objective_scale{};
    std::optional<LocalFittingAuditedState> best{};
};

const LocalFittingState & SelectLocalFittingFallbackState(
    const LocalFittingState & fallback_state,
    const std::optional<LocalFittingAuditedState> & audited_state)
{
    return audited_state.has_value() ? audited_state->state : fallback_state;
}

struct LocalFittingOffsetStats
{
    std::size_t atom_count{ 0 };
    std::size_t finite_count{ 0 };
    double median_absolute_offset{ 0.0 };
    double percentile_absolute_offset{ 0.0 };
    double maximum_absolute_offset{ 0.0 };
};

LocalFittingOffsetStats SummarizeLocalFittingOffsetValues(
    const std::vector<double> & offset_list);

void AppendLocalFittingOffsetSummary(
    std::ostringstream & stream,
    const LocalFittingOffsetStats & stats);

using PersistentSuspiciousRollbackReason = std::vector<std::size_t>;
using PersistentTerminalFailureReason =
    std::variant<PersistentSuspiciousRollbackReason, JointOffsetSolveStatus>;

struct PersistentTerminalFailureState
{
    PersistentTerminalFailureReason reason{};
    std::size_t stable_iteration_count{ 0 };
};

using PersistentTerminalFailureStateMap =
    std::map<LocalFittingClusterKey, PersistentTerminalFailureState>;
using TerminalPersistentFailureMap =
    std::map<LocalFittingClusterKey, PersistentTerminalFailureReason>;

struct LocalFittingTerminalSummary
{
    std::size_t suspicious_cluster_count{ 0 };
    std::size_t suspicious_atom_count{ 0 };
    std::size_t joint_offset_failure_cluster_count{ 0 };
    std::size_t joint_offset_failure_atom_count{ 0 };
    std::map<JointOffsetSolveStatus, std::size_t> joint_offset_failure_status_count{};

    std::size_t AtomCount() const
    {
        return suspicious_atom_count + joint_offset_failure_atom_count;
    }
};

enum class LocalFittingObjectiveAttemptDiagnosticStatus
{
    Scored,
    InvalidModel,
    ObjectiveUnavailable
};

struct LocalFittingObjectiveAttemptDiagnostic
{
    double effective_damping{ 1.0 };
    LocalFittingObjectiveAttemptDiagnosticStatus status{
        LocalFittingObjectiveAttemptDiagnosticStatus::ObjectiveUnavailable
    };
    std::optional<double> objective_scale{};
    std::optional<detail::LocalFittingObjectiveBreakdown> candidate_objective{};
    std::optional<detail::LocalFittingObjectiveBreakdown> previous_objective{};
    std::optional<detail::LocalFittingObjectiveBreakdown> best_objective{};
    double trust_region_radius{ 0.0 };
    double trust_region_step_norm{ 0.0 };
    bool rejected_by_previous{ false };
    bool rejected_by_best{ false };
};

struct LocalFittingRejectedClusterDiagnostic
{
    LocalFittingClusterKey key{};
    LocalFittingObjectiveAttemptDiagnostic attempt{};
};

struct LocalFittingCandidateSelection
{
    LocalFittingState assembled_state{};
    std::vector<LocalFittingClusterKey> accepted_key_list{};
    std::vector<LocalFittingClusterKey> rejected_key_list{};
    std::vector<LocalFittingRejectedClusterDiagnostic> rejected_cluster_diagnostic_list{};
    std::vector<LocalFittingClusterKey> grow_trust_region_key_list{};
    std::size_t accepted_polish_cluster_count{ 0 };
};

struct LocalFittingIterationProgressSnapshot
{
    std::size_t attempt_number{ 0 };
    std::size_t accepted_iteration_count{ 0 };
    std::size_t active_atom_count{ 0 };
    std::size_t terminal_atom_count{ 0 };
    std::size_t component_count{ 0 };
    std::size_t maximum_component_atom_count{ 0 };
    double maximum_component_active_ratio{ 0.0 };
    std::size_t accepted_cluster_count{ 0 };
    std::size_t rejected_cluster_count{ 0 };
    std::size_t accepted_atom_count{ 0 };
    std::size_t rejected_atom_count{ 0 };
    std::size_t accepted_polish_cluster_count{ 0 };
    std::size_t trust_region_grow_count{ 0 };
    std::size_t trust_region_shrink_count{ 0 };
    std::size_t trust_region_minimum_count{ 0 };
    std::size_t suspicious_atom_count{ 0 };
    std::size_t stationarity_ineligible_cluster_count{ 0 };
    bool combined_objective_rejected{ false };
    std::optional<double> accepted_maximum_transformed_change{};
    std::optional<double> raw_maximum_transformed_change{};
    std::string best_audit_iteration{};
    std::size_t audit_patience_count{ 0 };
};

struct LocalFittingProgressTableLayout
{
    std::array<std::size_t, 9> column_width_list{};
};

struct SecondStageNeighborSample
{
    std::size_t atom_index{ 0 };
    double distance{ 0.0 };
};

struct SecondStageAtomContext
{
    AtomObject * atom{ nullptr };
    LocalPotentialSampleList sample_entries{};
    std::vector<std::vector<SecondStageNeighborSample>> sample_neighbor_list{};
    double alpha_r{ 0.0 };
    double prior_width{ 1.0 };
};

struct SecondStageLocalFittingContext
{
    std::vector<SecondStageAtomContext> atom_context_list{};

    std::size_t AtomSize() const { return atom_context_list.size(); }
};

struct GaussianModelParameterSamples
{
    std::vector<double> amplitude_list{};
    std::vector<double> width_list{};
    std::vector<double> offset_list{};
};

struct SecondStageSeedRepairRecord
{
    std::size_t atom_index{ 0 };
    detail::SecondStageSeedRepairSource source{
        detail::SecondStageSeedRepairSource::GlobalMedian
    };
    GaussianModel3D original_model{};
    GaussianModel3D repaired_model{};
};

struct SecondStageInitialStateBuildResult
{
    std::optional<LocalFittingState> state{};
    std::vector<SecondStageSeedRepairRecord> repair_record_list{};
};

const char * GetSecondStageSeedRepairSourceText(
    detail::SecondStageSeedRepairSource source)
{
    switch (source)
    {
    case detail::SecondStageSeedRepairSource::GroupPosterior:
        return "group-posterior";
    case detail::SecondStageSeedRepairSource::GroupPrior:
        return "group-prior";
    case detail::SecondStageSeedRepairSource::LocalOls:
        return "local-ols";
    case detail::SecondStageSeedRepairSource::GroupMedian:
        return "group-median";
    case detail::SecondStageSeedRepairSource::GlobalMedian:
        return "global-median";
    }
    throw std::logic_error("Unknown second-stage seed repair source.");
}

double CalculateZeroOffsetResponse(
    const LocalPotentialSample & sample,
    const GaussianModel3D & model)
{
    const auto distance{ static_cast<double>(sample.point.distance) };
    const auto model_offset{ model.ResponseAtDistance(distance) - model.SignalAtDistance(distance) };
    return static_cast<double>(sample.response) - model_offset;
}

bool CanBuildFiniteZeroOffsetSamples(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & model)
{
    for (const auto & sample : sample_entries)
    {
        const auto response{ CalculateZeroOffsetResponse(sample, model) };
        if (!std::isfinite(response)) return false;
        if (std::abs(response) > static_cast<double>(std::numeric_limits<float>::max())) return false;
    }
    return true;
}

bool IsSameSuspiciousProfileRadius(double lhs, double rhs)
{
    const auto scale{ std::max({ std::abs(lhs), std::abs(rhs), 1.0 }) };
    return std::abs(lhs - rhs) <= kSuspiciousProfileDistanceTolerance * scale;
}

std::optional<ZeroOffsetProfileDiagnostics> BuildZeroOffsetProfileDiagnostics(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & model,
    const FitOptions & options)
{
    std::vector<std::pair<double, double>> profile_samples;
    profile_samples.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        const auto distance{ static_cast<double>(sample.point.distance) };
        if (distance < options.distance_min || distance > options.distance_max) continue;
        const auto response{ CalculateZeroOffsetResponse(sample, model) };
        if (!std::isfinite(response)) continue;
        if (std::abs(response) > static_cast<double>(std::numeric_limits<float>::max())) continue;
        profile_samples.emplace_back(distance, response);
    }

    if (profile_samples.empty()) return std::nullopt;
    ZeroOffsetProfileDiagnostics diagnostics;
    std::sort(
        profile_samples.begin(),
        profile_samples.end(),
        [](const auto & lhs, const auto & rhs)
        {
            return lhs.first < rhs.first;
        });

    diagnostics.distance_min = profile_samples.front().first;
    diagnostics.distance_max = profile_samples.back().first;
    for (std::size_t i = 0; i < profile_samples.size();)
    {
        const auto radius{ profile_samples.at(i).first };
        std::vector<double> response_list;
        while (i < profile_samples.size() &&
            IsSameSuspiciousProfileRadius(profile_samples.at(i).first, radius))
        {
            const auto response{ profile_samples.at(i).second };
            diagnostics.max_abs_response = std::max(diagnostics.max_abs_response, std::abs(response));
            response_list.emplace_back(response);
            i++;
        }
        diagnostics.radius_response_median_list.emplace_back(array_helper::ComputeMedian(response_list));
    }
    diagnostics.center_response = diagnostics.radius_response_median_list.front();
    return diagnostics;
}

bool HasUsableSuspiciousProfileBaseline(
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
        previous_profile.max_abs_response <= kRobustScaleMin)
    {
        return false;
    }
    const auto center_scale{
        std::max(std::abs(previous_profile.center_response), kRobustScaleMin)
    };
    for (std::size_t i = 1; i < previous_profile.radius_response_median_list.size(); i++)
    {
        const auto current_scale{
            std::abs(previous_profile.radius_response_median_list.at(i))
        };
        if (current_scale > kSuspiciousProfileReboundCenterRatio * center_scale) return false;
    }
    return true;
}

bool HasSuspiciousOffsetMagnitude(
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
            kRobustScaleMin
        })
    };
    return std::abs(candidate_offset_response) > kSuspiciousCompensationResponseRatio * reference_scale;
}

bool HasSuspiciousCenterSignFlip(
    const ZeroOffsetProfileDiagnostics & previous_profile,
    const ZeroOffsetProfileDiagnostics & candidate_profile)
{
    const auto previous_center{ previous_profile.center_response };
    const auto candidate_center{ candidate_profile.center_response };
    const auto reference_scale{ std::max(std::abs(previous_center), kRobustScaleMin) };
    return previous_center * candidate_center < 0.0 &&
        std::abs(candidate_center) > kSuspiciousProfileCenterSignFlipRatio * reference_scale;
}

bool HasSuspiciousRadialRebound(
    const ZeroOffsetProfileDiagnostics & previous_profile,
    const ZeroOffsetProfileDiagnostics & candidate_profile)
{
    if (candidate_profile.radius_response_median_list.size() < kSuspiciousProfileMinimumRadiusCount)
    {
        return false;
    }
    const auto reference_center_scale{
        std::max(std::abs(previous_profile.center_response), kRobustScaleMin)
    };
    const auto candidate_center_scale{
        std::max(std::abs(candidate_profile.center_response), kRobustScaleMin)
    };
    int upward_excursion_count{ 0 };
    auto previous_abs_response{ std::abs(candidate_profile.radius_response_median_list.front()) };
    for (std::size_t i = 1; i < candidate_profile.radius_response_median_list.size(); i++)
    {
        const auto current_abs_response{
            std::abs(candidate_profile.radius_response_median_list.at(i))
        };
        if (current_abs_response > kSuspiciousProfileReboundCenterRatio * candidate_center_scale &&
            current_abs_response > kSuspiciousProfileReboundReferenceRatio * reference_center_scale)
        {
            return true;
        }
        if (current_abs_response >
            previous_abs_response + kSuspiciousProfileUpwardExcursionReferenceRatio * reference_center_scale)
        {
            upward_excursion_count++;
        }
        previous_abs_response = current_abs_response;
    }
    return upward_excursion_count > kSuspiciousProfileMaximumUpwardExcursions;
}

bool HasSuspiciousWidthGrowth(
    const GaussianModel3D & previous_model,
    const GaussianModel3D & candidate_model,
    const ZeroOffsetProfileDiagnostics & previous_profile)
{
    if (!std::isfinite(candidate_model.GetWidth()) || candidate_model.GetWidth() <= 0.0) return true;
    if (candidate_model.GetWidth() > kSuspiciousWidthGrowthLimit * previous_model.GetWidth()) return true;
    const auto distance_range{ previous_profile.distance_max - previous_profile.distance_min };
    return distance_range > 0.0 && candidate_model.GetWidth() > kSuspiciousWidthRangeLimitRatio * distance_range;
}

bool HasSuspiciousAmplitudeOffsetCompensation(
    const GaussianModel3D & previous_model,
    const GaussianModel3D & candidate_model,
    const ZeroOffsetProfileDiagnostics & previous_profile)
{
    const auto signal_delta{
        candidate_model.SignalAtDistance(0.0) - previous_model.SignalAtDistance(0.0)
    };
    const auto offset_delta_response{
        (candidate_model.GetOffset() - previous_model.GetOffset()) * previous_model.OffsetBasisAtDistance(0.0)
    };
    const auto reference_scale{
        std::max({
            std::abs(previous_profile.center_response),
            std::abs(previous_model.SignalAtDistance(0.0)),
            kRobustScaleMin
        })
    };
    return signal_delta * offset_delta_response < 0.0 &&
        std::abs(signal_delta) > kSuspiciousCompensationResponseRatio * reference_scale &&
        std::abs(offset_delta_response) > kSuspiciousCompensationResponseRatio * reference_scale;
}

bool IsSuspiciousJointOffset(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & previous_model,
    const GaussianModel3D & offset_model,
    const FitOptions & options)
{
    if (!CanBuildFiniteZeroOffsetSamples(sample_entries, previous_model)) return false;
    if (!CanBuildFiniteZeroOffsetSamples(sample_entries, offset_model)) return true;
    const auto previous_profile{ BuildZeroOffsetProfileDiagnostics(sample_entries, previous_model, options) };
    if (HasSuspiciousOffsetMagnitude(
            previous_model,
            offset_model,
            previous_profile.has_value() ? previous_profile->max_abs_response : 0.0))
    {
        return true;
    }
    if (!previous_profile.has_value() ||
        !HasUsableSuspiciousProfileBaseline(previous_model, *previous_profile))
    {
        return false;
    }
    const auto candidate_profile{ BuildZeroOffsetProfileDiagnostics(sample_entries, offset_model, options) };
    if (!candidate_profile.has_value()) return true;
    return HasSuspiciousCenterSignFlip(*previous_profile, *candidate_profile) ||
        HasSuspiciousRadialRebound(*previous_profile, *candidate_profile) ||
        HasSuspiciousWidthGrowth(previous_model, offset_model, *previous_profile) ||
        HasSuspiciousAmplitudeOffsetCompensation(previous_model, offset_model, *previous_profile);
}

using FittedGaussianSnapshot = std::vector<GaussianModel3D>;
SecondStageLocalFittingContext BuildSecondStageLocalFittingContext(ModelObject & model_object)
{
    SecondStageLocalFittingContext context;
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    context.atom_context_list.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        context.atom_context_list.emplace_back(SecondStageAtomContext{ atom });
    }
    std::unordered_map<const AtomObject *, std::size_t> atom_index_map;
    atom_index_map.reserve(context.AtomSize());
    for (std::size_t i = 0; i < context.AtomSize(); i++)
    {
        atom_index_map.emplace(context.atom_context_list.at(i).atom, i);
    }
    const auto analysis_view{ model_object.GetAnalysisView() };

    std::unordered_map<GroupKey, std::vector<double>> width_samples_by_group;
    width_samples_by_group.reserve(context.AtomSize());
    for (const auto & atom_context : context.atom_context_list)
    {
        const auto * atom{ atom_context.atom };
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
        const auto width{ local_view.GetEstimateMDPDE().GetWidth() };
        if (numeric_validation::IsFinitePositive(width))
        {
            width_samples_by_group[data_internal::GetGroupKey(atom)].emplace_back(width);
        }
    }

    std::unordered_map<GroupKey, double> median_width_by_group;
    median_width_by_group.reserve(width_samples_by_group.size());
    for (const auto & [group_key, width_samples] : width_samples_by_group)
    {
        if (width_samples.empty()) continue;
        const auto median_width{ array_helper::ComputeMedian(width_samples) };
        if (numeric_validation::IsFinitePositive(median_width))
        {
            median_width_by_group.emplace(group_key, median_width);
        }
    }

    for (std::size_t atom_index = 0; atom_index < context.AtomSize(); atom_index++)
    {
        auto & atom_context{ context.atom_context_list.at(atom_index) };
        const auto * atom{ atom_context.atom };
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
        atom_context.sample_entries = local_view.GetSamplingEntries(false);
        atom_context.alpha_r = local_view.GetAlphaR();
        const auto group_key{ data_internal::GetGroupKey(atom) };
        const auto current_width{ local_view.GetEstimateMDPDE().GetWidth() };
        atom_context.prior_width = numeric_validation::IsFinitePositive(current_width) ?
            current_width : 1.0;
        const auto median_width_iter{ median_width_by_group.find(group_key) };
        if (median_width_iter != median_width_by_group.end())
        {
            atom_context.prior_width = median_width_iter->second;
        }
        if (analysis_view.HasAtomGroup(group_key))
        {
            const auto group_prior_width{ analysis_view.GetAtomGroupPrior(group_key).GetWidth() };
            if (numeric_validation::IsFinitePositive(group_prior_width))
            {
                atom_context.prior_width = group_prior_width;
            }
        }
    }

    for (std::size_t atom_index = 0; atom_index < context.AtomSize(); atom_index++)
    {
        auto & atom_context{ context.atom_context_list.at(atom_index) };
        const auto * atom{ atom_context.atom };
        const auto neighbor_atom_list{ atom->FindNeighborAtoms(kNeighborAtomSearchRange) };
        std::vector<std::size_t> selected_neighbor_index_list;
        selected_neighbor_index_list.reserve(neighbor_atom_list.size());
        for (const auto * neighbor_atom : neighbor_atom_list)
        {
            const auto neighbor_iter{ atom_index_map.find(neighbor_atom) };
            if (neighbor_iter == atom_index_map.end()) continue;

            selected_neighbor_index_list.emplace_back(neighbor_iter->second);
        }

        atom_context.sample_neighbor_list.resize(atom_context.sample_entries.size());
        for (std::size_t sample_index = 0; sample_index < atom_context.sample_entries.size(); sample_index++)
        {
            const auto & sample{ atom_context.sample_entries.at(sample_index) };
            auto & sample_neighbor_list{ atom_context.sample_neighbor_list.at(sample_index) };
            sample_neighbor_list.reserve(selected_neighbor_index_list.size());
            for (const auto neighbor_index : selected_neighbor_index_list)
            {
                const auto distance{
                    static_cast<double>(
                        array_helper::ComputeNorm<float>(
                            sample.point.position,
                            context.atom_context_list.at(neighbor_index).atom->GetPositionRef()))
                };
                if (distance > kNeighborContributionDistanceMax) continue;

                sample_neighbor_list.emplace_back(SecondStageNeighborSample{ neighbor_index, distance });
            }
        }
    }

    return context;
}

detail::LocalFittingCouplingTopology BuildLocalFittingCouplingTopology(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & initial_state)
{
    if (initial_state.size() != context.AtomSize())
    {
        throw std::invalid_argument(
            "Local fitting coupling topology state size is inconsistent.");
    }

    detail::LocalFittingCouplingGraphBuilder builder{ context.AtomSize() };
    const auto invalid_jacobian{
        Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN())
    };
    for (std::size_t atom_index = 0; atom_index < context.AtomSize(); atom_index++)
    {
        const auto & atom_context{ context.atom_context_list.at(atom_index) };
        for (std::size_t sample_index = 0;
            sample_index < atom_context.sample_entries.size();
            sample_index++)
        {
            const auto & sample{ atom_context.sample_entries.at(sample_index) };
            std::vector<detail::LocalFittingCouplingParticipant> participant_list;
            participant_list.reserve(
                atom_context.sample_neighbor_list.at(sample_index).size() + 1);
            const auto target_evaluation{
                detail::EvaluateLocalFittingTransformedResponse(
                    initial_state.at(atom_index).mdpde.GetModel(),
                    static_cast<double>(sample.point.distance))
            };
            participant_list.emplace_back(
                detail::LocalFittingCouplingParticipant{
                    atom_index,
                    target_evaluation.has_value() ?
                        target_evaluation->jacobian : invalid_jacobian
                });
            for (const auto & neighbor_sample :
                atom_context.sample_neighbor_list.at(sample_index))
            {
                const auto neighbor_evaluation{
                    detail::EvaluateLocalFittingTransformedResponse(
                        initial_state.at(neighbor_sample.atom_index).mdpde.GetModel(),
                        neighbor_sample.distance)
                };
                participant_list.emplace_back(
                    detail::LocalFittingCouplingParticipant{
                        neighbor_sample.atom_index,
                        neighbor_evaluation.has_value() ?
                            neighbor_evaluation->jacobian : invalid_jacobian
                    });
            }
            builder.AddSample(
                LocalFittingObjectiveSampleRef{ atom_index, sample_index },
                std::move(participant_list));
        }
    }

    auto weighted_topology{
        builder.BuildWeighted(
            kLocalFittingCouplingMinimumWeight,
            std::vector<double>{
                kLocalFittingCouplingSensitivityMinimumWeightList.begin(),
                kLocalFittingCouplingSensitivityMinimumWeightList.end()
            })
    };
    return weighted_topology.has_value() ? std::move(*weighted_topology) : builder.BuildBinary();
}

void AddGaussianModelParameterSample(
    GaussianModelParameterSamples & samples,
    const GaussianModel3D & model)
{
    if (!detail::IsValidSecondStageGaussianModel(model)) return;
    samples.amplitude_list.emplace_back(model.GetAmplitude());
    samples.width_list.emplace_back(model.GetWidth());
    samples.offset_list.emplace_back(model.GetOffset());
}

std::optional<GaussianModel3DWithUncertainty> BuildValidGaussianParameterMedian(
    const GaussianModelParameterSamples & samples)
{
    if (samples.amplitude_list.empty() ||
        samples.width_list.empty() ||
        samples.offset_list.empty())
    {
        return std::nullopt;
    }
    const GaussianModel3D median_model{
        array_helper::ComputeMedian(samples.amplitude_list),
        array_helper::ComputeMedian(samples.width_list),
        array_helper::ComputeMedian(samples.offset_list)
    };
    if (!detail::IsValidSecondStageGaussianModel(median_model))
    {
        return std::nullopt;
    }
    return GaussianModel3DWithUncertainty{
        median_model,
        GaussianModel3DUncertainty{}
    };
}

SecondStageInitialStateBuildResult BuildInitialLocalFittingState(
    const SecondStageLocalFittingContext & context,
    const ModelAnalysisView & analysis_view)
{
    SecondStageInitialStateBuildResult build_result;
    LocalFittingState state(context.AtomSize());
    std::vector<std::optional<GaussianModel3DWithUncertainty>> group_prior_list(
        context.AtomSize());
    std::unordered_map<GroupKey, GaussianModelParameterSamples> samples_by_group;
    GaussianModelParameterSamples global_samples;

    for (std::size_t i = 0; i < context.AtomSize(); i++)
    {
        const auto * atom{ context.atom_context_list.at(i).atom };
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
        state.at(i) = local_view.GetGaussianResult();
        const auto group_key{ data_internal::GetGroupKey(atom) };
        if (analysis_view.HasAtomGroup(group_key))
        {
            group_prior_list.at(i) = analysis_view.GetAtomGroupPriorWithUncertainty(group_key);
        }

        const auto & result{ state.at(i) };
        const auto & mdpde_model{ result.mdpde.GetModel() };
        std::optional<GaussianModel3DWithUncertainty> preferred_model;
        if (detail::IsValidSecondStageGaussianModel(mdpde_model))
        {
            preferred_model = result.mdpde;
        }
        else
        {
            const auto direct_selection{
                detail::SelectSecondStageSeedRepair(
                    detail::SecondStageSeedRepairCandidates{
                        result.posterior,
                        group_prior_list.at(i),
                        result.ols,
                        std::nullopt,
                        std::nullopt
                    })
            };
            if (direct_selection.has_value())
            {
                preferred_model = direct_selection->model;
            }
        }
        if (!preferred_model.has_value()) continue;

        AddGaussianModelParameterSample(
            samples_by_group[group_key],
            preferred_model->GetModel());
        AddGaussianModelParameterSample(
            global_samples,
            preferred_model->GetModel());
    }

    std::unordered_map<GroupKey, GaussianModel3DWithUncertainty> median_by_group;
    median_by_group.reserve(samples_by_group.size());
    for (const auto & [group_key, samples] : samples_by_group)
    {
        const auto median_model{ BuildValidGaussianParameterMedian(samples) };
        if (median_model.has_value())
        {
            median_by_group.emplace(group_key, *median_model);
        }
    }
    const auto global_median{ BuildValidGaussianParameterMedian(global_samples) };

    for (std::size_t i = 0; i < context.AtomSize(); i++)
    {
        auto & result{ state.at(i) };
        const auto original_model{ result.mdpde.GetModel() };
        if (!detail::IsValidSecondStageGaussianModel(original_model))
        {
            const auto group_key{
                data_internal::GetGroupKey(context.atom_context_list.at(i).atom)
            };
            std::optional<GaussianModel3DWithUncertainty> group_median;
            const auto group_median_iter{ median_by_group.find(group_key) };
            if (group_median_iter != median_by_group.end())
            {
                group_median = group_median_iter->second;
            }
            const auto selection{
                detail::SelectSecondStageSeedRepair(
                    detail::SecondStageSeedRepairCandidates{
                        result.posterior,
                        group_prior_list.at(i),
                        result.ols,
                        group_median,
                        global_median
                    })
            };
            if (!selection.has_value())
            {
                return build_result;
            }

            const auto repaired_seed{
                detail::BuildRepairedSecondStageSeed(original_model, *selection)
            };
            const auto repaired_model{ repaired_seed.GetModel() };
            if (!detail::IsValidSecondStageGaussianModel(repaired_model))
            {
                return build_result;
            }
            result.mdpde = repaired_seed;
            build_result.repair_record_list.emplace_back(
                SecondStageSeedRepairRecord{
                    i,
                    selection->source,
                    original_model,
                    repaired_model
                });
        }
    }
    build_result.state = std::move(state);
    return build_result;
}

void LogSecondStageSeedRepairs(
    const std::vector<SecondStageSeedRepairRecord> & repair_record_list,
    const FitOptions & options)
{
    if (options.quiet_mode || repair_record_list.empty()) return;

    constexpr std::array<detail::SecondStageSeedRepairSource, 5> source_list{
        detail::SecondStageSeedRepairSource::GroupPosterior,
        detail::SecondStageSeedRepairSource::GroupPrior,
        detail::SecondStageSeedRepairSource::LocalOls,
        detail::SecondStageSeedRepairSource::GroupMedian,
        detail::SecondStageSeedRepairSource::GlobalMedian
    };
    std::array<std::size_t, source_list.size()> source_count{};
    for (const auto & record : repair_record_list)
    {
        source_count.at(static_cast<std::size_t>(record.source))++;
    }

    std::ostringstream summary;
    summary << "Repaired invalid second-stage seed atoms = "
        << repair_record_list.size() << ", sources = ";
    for (std::size_t i = 0; i < source_list.size(); i++)
    {
        if (i != 0) summary << ", ";
        summary << GetSecondStageSeedRepairSourceText(source_list.at(i))
            << ":" << source_count.at(i);
    }
    summary << ".";
    Logger::Log(LogLevel::Info, summary.str());

    for (const auto & record : repair_record_list)
    {
        std::ostringstream detail_message;
        detail_message << "Second-stage seed repair: atom index = "
            << record.atom_index
            << ", source = " << GetSecondStageSeedRepairSourceText(record.source)
            << std::scientific << std::setprecision(2)
            << ", original A/B/C = "
            << record.original_model.GetAmplitude() << "/"
            << record.original_model.GetWidth() << "/"
            << record.original_model.GetOffset()
            << ", repaired A/B/C = "
            << record.repaired_model.GetAmplitude() << "/"
            << record.repaired_model.GetWidth() << "/"
            << record.repaired_model.GetOffset() << ".";
        Logger::Log(LogLevel::Debug, detail_message.str());
    }
}

void LogLocalFittingCouplingTopology(
    const detail::LocalFittingCouplingTopology & topology,
    const FitOptions & options)
{
    if (options.quiet_mode) return;

    std::vector<std::size_t> atom_index_list(topology.adjacency_list.size());
    for (std::size_t i = 0; i < atom_index_list.size(); i++) atom_index_list.at(i) = i;
    const auto partition{
        detail::BuildLocalFittingCouplingPartition(topology, atom_index_list)
    };
    std::size_t maximum_component_size{ 0 };
    for (const auto & [key, sample_list] : partition.sample_id_list_by_key)
    {
        static_cast<void>(sample_list);
        maximum_component_size = std::max(maximum_component_size, key.size());
    }
    const auto maximum_component_ratio{
        atom_index_list.empty() ? 0.0 :
            static_cast<double>(maximum_component_size) /
                static_cast<double>(atom_index_list.size())
    };

    const auto & summary{ topology.summary };
    if (!summary.uses_weighted_graph)
    {
        Logger::Log(
            LogLevel::Warning,
            "Weighted local-fitting coupling graph is unavailable; using binary connectivity.");
        if (Logger::GetLogLevel() >= LogLevel::Debug)
        {
            Logger::Log(
                LogLevel::Debug,
                "Local-fitting weighted threshold sensitivity is unavailable in binary fallback mode.");
        }
    }
    std::ostringstream message;
    message << "Local-fitting coupling graph mode = "
        << (summary.uses_weighted_graph ? "weighted" : "binary-fallback")
        << std::scientific << std::setprecision(2)
        << ", minimum weight = " << kLocalFittingCouplingMinimumWeight
        << ", candidate/retained/cut edges = "
        << summary.candidate_edge_count << "/"
        << summary.retained_edge_count << "/"
        << summary.cut_edge_count
        << ", weight p50/p95/max = "
        << summary.weight_median << "/"
        << summary.weight_percentile_95 << "/"
        << summary.weight_maximum
        << ", initial components/max atoms/ratio = "
        << partition.sample_id_list_by_key.size() << "/"
        << maximum_component_size << "/"
        << std::fixed << std::setprecision(2)
        << maximum_component_ratio << ".";
    Logger::Log(LogLevel::Info, message.str());

    for (const auto & sensitivity : summary.threshold_sensitivity_list)
    {
        std::ostringstream sensitivity_message;
        sensitivity_message
            << std::scientific << std::setprecision(2)
            << "Coupling sensitivity: threshold=" << sensitivity.minimum_weight
            << ", retained/cut="
            << sensitivity.retained_edge_count << "/"
            << sensitivity.cut_edge_count
            << ", components/max-atoms/ratio="
            << sensitivity.component_count << "/"
            << sensitivity.maximum_component_size << "/"
            << std::fixed << std::setprecision(2)
            << sensitivity.maximum_component_ratio << ".";
        Logger::Log(LogLevel::Info, sensitivity_message.str());
    }
}

FittedGaussianSnapshot BuildFittedGaussianSnapshot(
    const LocalFittingState & state)
{
    FittedGaussianSnapshot snapshot;
    snapshot.reserve(state.size());
    for (const auto & result : state)
    {
        snapshot.emplace_back(result.mdpde.GetModel());
    }
    return snapshot;
}

JointOffsetBuildResult BuildJointOffsetSystem(
    const SecondStageLocalFittingContext & context,
    const std::vector<std::size_t> & active_index_list,
    const FittedGaussianSnapshot & snapshot,
    double ridge_ratio,
    const std::vector<double> & ridge_multiplier_list,
    bool log_debug_diagnostics)
{
    const auto atom_size{ context.AtomSize() };
    if (snapshot.size() != atom_size || ridge_multiplier_list.size() != atom_size)
    {
        throw std::invalid_argument("Joint offset input sizes are inconsistent.");
    }

    std::vector<int> active_column_by_atom_index(atom_size, -1);
    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto atom_index{ active_index_list.at(i) };
        if (atom_index >= atom_size)
        {
            throw std::invalid_argument("Joint offset active atom index is out of range.");
        }
        active_column_by_atom_index.at(atom_index) = static_cast<int>(i);
    }

    const auto column_count{ static_cast<Eigen::Index>(active_index_list.size()) };
    std::vector<Eigen::Triplet<double>> triplet_list;
    std::vector<double> response_list;
    Eigen::VectorXd column_square_sum{ Eigen::VectorXd::Zero(column_count) };
    std::map<std::pair<Eigen::Index, Eigen::Index>, double> column_cross_sum_map;
    ActiveCouplingGraph active_coupling_graph(active_index_list.size());
    std::vector<std::pair<Eigen::Index, double>> row_basis_entries;
    for (const auto active_index : active_index_list)
    {
        if (active_index >= atom_size)
        {
            throw std::invalid_argument("Joint offset active atom index is out of range.");
        }
        const auto target_column{ active_column_by_atom_index.at(active_index) };
        if (target_column < 0)
        {
            throw std::invalid_argument("Joint offset active column is missing.");
        }

        const auto & atom_context{ context.atom_context_list.at(active_index) };
        const auto & target_model{ snapshot.at(active_index) };
        row_basis_entries.reserve(atom_size);
        for (std::size_t sample_index = 0; sample_index < atom_context.sample_entries.size(); sample_index++)
        {
            const auto & sample{ atom_context.sample_entries.at(sample_index) };
            if (!std::isfinite(static_cast<double>(sample.response)))
            {
                throw std::runtime_error("Joint offset sample response is not finite.");
            }
            const auto target_distance{ static_cast<double>(sample.point.distance) };
            const auto target_signal{ target_model.SignalAtDistance(target_distance) };
            const auto target_basis{ target_model.OffsetBasisAtDistance(target_distance) };
            if (!std::isfinite(target_signal) || !std::isfinite(target_basis))
            {
                throw std::runtime_error("Joint offset target model evaluation is not finite.");
            }
            auto residual{ static_cast<double>(sample.response) - target_signal };
            row_basis_entries.clear();
            if (std::abs(target_basis) > std::numeric_limits<double>::epsilon())
            {
                row_basis_entries.emplace_back(static_cast<Eigen::Index>(target_column), target_basis);
            }

            for (const auto & neighbor_sample : atom_context.sample_neighbor_list.at(sample_index))
            {
                const auto & neighbor_model{ snapshot.at(neighbor_sample.atom_index) };
                const auto neighbor_column{ active_column_by_atom_index.at(neighbor_sample.atom_index) };
                if (neighbor_column < 0)
                {
                    const auto response{ neighbor_model.ResponseAtDistance(neighbor_sample.distance) };
                    if (!std::isfinite(response))
                    {
                        throw std::runtime_error("Joint offset fixed neighbor model evaluation is not finite.");
                    }
                    residual -= response;
                    continue;
                }

                const auto signal{ neighbor_model.SignalAtDistance(neighbor_sample.distance) };
                const auto basis{ neighbor_model.OffsetBasisAtDistance(neighbor_sample.distance) };
                if (!std::isfinite(signal) || !std::isfinite(basis))
                {
                    throw std::runtime_error("Joint offset active neighbor model evaluation is not finite.");
                }
                residual -= signal;
                if (std::abs(basis) > std::numeric_limits<double>::epsilon())
                {
                    row_basis_entries.emplace_back(static_cast<Eigen::Index>(neighbor_column), basis);
                }
            }
            if (!std::isfinite(residual))
            {
                throw std::runtime_error("Joint offset residual is not finite.");
            }
            if (row_basis_entries.empty()) continue;

            const auto row_index{ static_cast<Eigen::Index>(response_list.size()) };
            response_list.emplace_back(residual);
            for (const auto & [column_index, basis] : row_basis_entries)
            {
                triplet_list.emplace_back(row_index, column_index, basis);
                column_square_sum(column_index) += basis * basis;
            }
            for (std::size_t i = 0; i < row_basis_entries.size(); i++)
            {
                const auto [left_column, left_basis]{ row_basis_entries.at(i) };
                for (std::size_t j = i + 1; j < row_basis_entries.size(); j++)
                {
                    const auto [right_column, right_basis]{ row_basis_entries.at(j) };
                    if (left_column == right_column) continue;
                    const auto column_pair{ std::minmax(left_column, right_column) };
                    column_cross_sum_map[column_pair] += left_basis * right_basis;
                }
            }
        }
    }

    Eigen::VectorXd proactive_ridge_multiplier{ Eigen::VectorXd::Ones(column_count) };
    for (const auto & [column_pair, cross_sum] : column_cross_sum_map)
    {
        const auto left_column{ column_pair.first };
        const auto right_column{ column_pair.second };
        const auto left_square_sum{ column_square_sum(left_column) };
        const auto right_square_sum{ column_square_sum(right_column) };
        if (left_square_sum <= std::numeric_limits<double>::epsilon() ||
            right_square_sum <= std::numeric_limits<double>::epsilon())
        {
            continue;
        }
        const auto overlap{ std::abs(cross_sum) / std::sqrt(left_square_sum * right_square_sum) };
        if (!std::isfinite(overlap))
        {
            continue;
        }
        active_coupling_graph.at(static_cast<std::size_t>(left_column)).emplace_back(
            ActiveCouplingEdge{ static_cast<std::size_t>(right_column), overlap });
        active_coupling_graph.at(static_cast<std::size_t>(right_column)).emplace_back(
            ActiveCouplingEdge{ static_cast<std::size_t>(left_column), overlap });
        if (overlap < kJointOffsetCollinearityOverlapThreshold) continue;

        proactive_ridge_multiplier(left_column) = std::max(
            proactive_ridge_multiplier(left_column),
            kCollinearJointOffsetRidgeMultiplier);
        proactive_ridge_multiplier(right_column) = std::max(
            proactive_ridge_multiplier(right_column),
            kCollinearJointOffsetRidgeMultiplier);
    }

    const auto row_count{ static_cast<Eigen::Index>(response_list.size()) };
    Eigen::VectorXd response{ Eigen::VectorXd::Zero(row_count) };
    for (Eigen::Index row_index = 0; row_index < row_count; row_index++)
    {
        response(row_index) = response_list.at(static_cast<std::size_t>(row_index));
    }

    algorithm::WeightedRidgeSystem system;
    system.design_matrix.resize(row_count, column_count);
    system.design_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
    const auto conditioning{
        detail::EvaluateLocalFittingJointOffsetConditioning(
            system.design_matrix,
            kJointOffsetConditioningPivotRatioThreshold)
    };
    if (conditioning.guard_required)
    {
        proactive_ridge_multiplier.array() = proactive_ridge_multiplier.array().max(
            kCollinearJointOffsetRidgeMultiplier);
        if (log_debug_diagnostics)
        {
            std::ostringstream message;
            message
                << std::scientific << std::setprecision(2)
                << "Joint offset conditioning guard: columns = " << column_count
                << ", normalized LDLT pivot ratio = " << conditioning.pivot_ratio
                << ", proactive ridge multiplier = "
                << kCollinearJointOffsetRidgeMultiplier << ".";
            Logger::Log(LogLevel::Debug, message.str());
        }
    }
    system.response = std::move(response);
    system.previous_parameter = Eigen::VectorXd::Zero(column_count);
    system.ridge_diagonal = Eigen::VectorXd::Zero(column_count);
    for (Eigen::Index column_index = 0; column_index < column_count; column_index++)
    {
        const auto atom_index{ active_index_list.at(static_cast<std::size_t>(column_index)) };
        const auto & model{ snapshot.at(atom_index) };
        system.previous_parameter(column_index) = model.GetOffset();
        const auto square_sum{ column_square_sum(column_index) };
        const auto multiplier{ ridge_multiplier_list.at(atom_index) };
        if (!std::isfinite(multiplier) || multiplier <= 0.0)
        {
            throw std::invalid_argument("Joint offset ridge multiplier must be positive and finite.");
        }
        const auto combined_multiplier{
            std::max(multiplier, proactive_ridge_multiplier(column_index))
        };
        const auto base_ridge{
            square_sum > std::numeric_limits<double>::epsilon()
                ? ridge_ratio * square_sum
                : ridge_ratio / kJointOffsetRidgeRatio
        };
        system.ridge_diagonal(column_index) = combined_multiplier * base_ridge;
    }
    return JointOffsetBuildResult{
        std::move(system),
        std::move(active_coupling_graph)
    };
}

double CalculateWeightedRidgeSurrogateObjective(
    const algorithm::WeightedRidgeSystem & system,
    const Eigen::VectorXd & weight,
    const Eigen::VectorXd & offset)
{
    if (system.response.size() != weight.size())
    {
        throw std::invalid_argument("Weighted ridge objective input sizes are inconsistent.");
    }
    if (system.previous_parameter.size() != offset.size() || system.ridge_diagonal.size() != offset.size())
    {
        throw std::invalid_argument("Weighted ridge objective parameter sizes are inconsistent.");
    }
    if (system.response.size() == 0)
    {
        return std::numeric_limits<double>::infinity();
    }

    const Eigen::VectorXd residual{ system.response - system.design_matrix * offset };
    const auto weighted_residual_loss{ weight.cwiseProduct(residual.cwiseAbs2()).sum() };
    const Eigen::VectorXd offset_delta{ offset - system.previous_parameter };
    const auto ridge_loss{ system.ridge_diagonal.cwiseProduct(offset_delta.cwiseAbs2()).sum() };
    const auto objective{ (weighted_residual_loss + ridge_loss) / static_cast<double>(system.response.size()) };
    return std::isfinite(objective) ? objective : std::numeric_limits<double>::infinity();
}

bool IsJointOffsetObjectiveDeteriorated(double updated_objective, double current_objective)
{
    if (!std::isfinite(updated_objective)) return true;
    if (!std::isfinite(current_objective)) return false;
    const auto scale{
        std::max({ std::abs(updated_objective), std::abs(current_objective), 1.0 })
    };
    return updated_objective > current_objective + kJointOffsetIrlsObjectiveRelativeTolerance * scale;
}

double CalculateMedianAbsoluteDeviationScale(const std::vector<double> & value_list)
{
    if (value_list.empty())
    {
        return std::numeric_limits<double>::infinity();
    }

    const auto median_value{ array_helper::ComputeMedian(value_list) };
    std::vector<double> deviation_list;
    deviation_list.reserve(value_list.size());
    for (const auto value : value_list)
    {
        deviation_list.emplace_back(std::abs(value - median_value));
    }
    return kRobustScaleMultiplier * array_helper::ComputeMedian(deviation_list);
}

JointOffsetSolveResult EstimateJointOffsets(
    const SecondStageLocalFittingContext & context,
    const std::vector<std::size_t> & active_index_list,
    const FittedGaussianSnapshot & snapshot,
    double ridge_ratio,
    const std::vector<double> & ridge_multiplier_list,
    bool log_debug_diagnostics)
{
    Eigen::VectorXd previous_offset{
        Eigen::VectorXd::Zero(static_cast<Eigen::Index>(active_index_list.size()))
    };
    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto atom_index{ active_index_list.at(i) };
        if (atom_index >= context.AtomSize())
        {
            throw std::invalid_argument("Joint offset active atom index is out of range.");
        }
        previous_offset(static_cast<Eigen::Index>(i)) = snapshot.at(atom_index).GetOffset();
    }
    JointOffsetBuildResult build_result;
    try
    {
        build_result = BuildJointOffsetSystem(
            context,
            active_index_list,
            snapshot,
            ridge_ratio,
            ridge_multiplier_list,
            log_debug_diagnostics);
    }
    catch (const std::runtime_error &)
    {
        return JointOffsetSolveResult{
            JointOffsetSolveStatus::SystemBuildFailed,
            previous_offset,
            {}
        };
    }
    auto system{ std::move(build_result.system) };
    auto active_coupling_graph{ std::move(build_result.active_coupling_graph) };
    if (system.response.size() == 0 || system.previous_parameter.size() == 0)
    {
        return JointOffsetSolveResult{
            JointOffsetSolveStatus::EmptySystem,
            previous_offset,
            std::move(active_coupling_graph)
        };
    }

    Eigen::VectorXd weight{ Eigen::VectorXd::Ones(system.response.size()) };
    algorithm::WeightedRidgeSolver solver{ system };
    Eigen::VectorXd offset;
    if (!solver.Solve(system, weight, offset))
    {
        return JointOffsetSolveResult{
            JointOffsetSolveStatus::InitialSolveFailed,
            previous_offset,
            std::move(active_coupling_graph)
        };
    }

    for (int iteration = 0; iteration < kRobustLossMaximumIterations; iteration++)
    {
        const Eigen::VectorXd residual{ system.response - system.design_matrix * offset };
        std::vector<double> residual_list(residual.data(), residual.data() + residual.size());
        const auto residual_scale{
            std::max(CalculateMedianAbsoluteDeviationScale(residual_list), kRobustScaleMin)
        };
        for (Eigen::Index i = 0; i < residual.size(); i++)
        {
            weight(i) = algorithm::CalculateRobustWeight(
                kSecondStageRobustLossKind,
                residual(i),
                residual_scale,
                kRobustLossCutoffMultiplier);
        }

        Eigen::VectorXd updated_offset;
        if (!solver.Solve(system, weight, updated_offset))
        {
            return JointOffsetSolveResult{
                JointOffsetSolveStatus::IrlsSolveFailed,
                previous_offset,
                std::move(active_coupling_graph)
            };
        }
        const auto current_objective{
            CalculateWeightedRidgeSurrogateObjective(system, weight, offset)
        };
        const auto updated_objective{
            CalculateWeightedRidgeSurrogateObjective(system, weight, updated_offset)
        };
        if (IsJointOffsetObjectiveDeteriorated(updated_objective, current_objective))
        {
            return JointOffsetSolveResult{
                JointOffsetSolveStatus::IrlsObjectiveDeteriorated,
                offset,
                std::move(active_coupling_graph)
            };
        }
        const auto maximum_change{
            algorithm::CalculateMaximumNormalizedVectorChange(updated_offset, offset, kJointOffsetIrlsScaleFloor)
        };
        offset = std::move(updated_offset);
        if (maximum_change < kJointOffsetIrlsNormalizedChangeTolerance)
        {
            return JointOffsetSolveResult{
                JointOffsetSolveStatus::Converged,
                offset,
                std::move(active_coupling_graph)
            };
        }
    }

    return JointOffsetSolveResult{
        JointOffsetSolveStatus::IrlsMaximumIterationsReached,
        offset,
        std::move(active_coupling_graph)
    };
}

ClusteredJointOffsetSolveResult EstimateClusteredJointOffsets(
    const SecondStageLocalFittingContext & context,
    const std::vector<std::size_t> & active_index_list,
    const std::vector<LocalFittingClusterKey> & cluster_key_list,
    const FittedGaussianSnapshot & snapshot,
    double ridge_ratio,
    const std::vector<double> & ridge_multiplier_list,
    bool log_debug_diagnostics)
{
    std::unordered_map<std::size_t, std::size_t> active_position_by_atom_index;
    active_position_by_atom_index.reserve(active_index_list.size());
    for (std::size_t active_position = 0;
        active_position < active_index_list.size();
        active_position++)
    {
        if (!active_position_by_atom_index.emplace(
                active_index_list.at(active_position), active_position).second)
        {
            throw std::invalid_argument(
                "Clustered joint offset active atom indexes must be unique.");
        }
    }

    ClusteredJointOffsetSolveResult clustered_result;
    clustered_result.offset = Eigen::VectorXd::Zero(
        static_cast<Eigen::Index>(active_index_list.size()));
    clustered_result.active_coupling_graph.resize(active_index_list.size());
    std::vector<char> covered_active_position_list(active_index_list.size(), 0);
    for (const auto & key : cluster_key_list)
    {
        if (key.empty() || !std::is_sorted(key.begin(), key.end()) ||
            std::adjacent_find(key.begin(), key.end()) != key.end())
        {
            throw std::invalid_argument(
                "Clustered joint offset key must be non-empty and canonical.");
        }

        const auto cluster_result{
            EstimateJointOffsets(
                context,
                key,
                snapshot,
                ridge_ratio,
                ridge_multiplier_list,
                log_debug_diagnostics)
        };
        if (cluster_result.offset.size() != static_cast<Eigen::Index>(key.size()))
        {
            throw std::logic_error("Clustered joint offset result size is inconsistent.");
        }
        if (!cluster_result.active_coupling_graph.empty() &&
            cluster_result.active_coupling_graph.size() != key.size())
        {
            throw std::logic_error("Clustered joint offset graph size is inconsistent.");
        }
        if (!clustered_result.health_by_key.emplace(
                key,
                LocalFittingClusterHealth{ cluster_result.status }).second)
        {
            throw std::invalid_argument("Clustered joint offset keys must be unique.");
        }

        for (std::size_t local_position = 0; local_position < key.size(); local_position++)
        {
            const auto atom_index{ key.at(local_position) };
            const auto active_iter{ active_position_by_atom_index.find(atom_index) };
            if (active_iter == active_position_by_atom_index.end())
            {
                throw std::invalid_argument("Clustered joint offset atom is not active.");
            }
            const auto active_position{ active_iter->second };
            if (covered_active_position_list.at(active_position) != 0)
            {
                throw std::invalid_argument(
                    "Clustered joint offset keys must not share active atoms.");
            }
            covered_active_position_list.at(active_position) = 1;
            clustered_result.offset(static_cast<Eigen::Index>(active_position)) =
                cluster_result.offset(static_cast<Eigen::Index>(local_position));

            if (cluster_result.active_coupling_graph.empty()) continue;
            for (const auto & edge : cluster_result.active_coupling_graph.at(local_position))
            {
                if (edge.neighbor_index >= key.size())
                {
                    throw std::logic_error("Clustered joint offset graph neighbor is out of range.");
                }
                const auto neighbor_iter{
                    active_position_by_atom_index.find(key.at(edge.neighbor_index))
                };
                if (neighbor_iter == active_position_by_atom_index.end())
                {
                    throw std::logic_error("Clustered joint offset graph neighbor is not active.");
                }
                clustered_result.active_coupling_graph.at(active_position).emplace_back(
                    ActiveCouplingEdge{ neighbor_iter->second, edge.overlap });
            }
        }
    }

    if (std::find(
            covered_active_position_list.begin(),
            covered_active_position_list.end(),
            0) != covered_active_position_list.end())
    {
        throw std::invalid_argument("Clustered joint offset keys must cover active atoms exactly once.");
    }
    return clustered_result;
}

void ApplyJointOffsetsToSnapshot(
    const std::vector<std::size_t> & active_index_list,
    const Eigen::VectorXd & offset,
    FittedGaussianSnapshot & snapshot)
{
    if (active_index_list.size() != static_cast<std::size_t>(offset.size()))
    {
        throw std::invalid_argument("Joint offset result size is inconsistent.");
    }
    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto atom_index{ active_index_list.at(i) };
        if (atom_index >= snapshot.size())
        {
            throw std::invalid_argument("Joint offset active atom index is out of range.");
        }
        snapshot.at(atom_index) = snapshot.at(atom_index).WithOffset(offset(static_cast<Eigen::Index>(i)));
    }
}

double CalculateSecondStageAdjustedResponse(
    const SecondStageLocalFittingContext & context,
    std::size_t atom_index,
    std::size_t sample_index,
    const FittedGaussianSnapshot & snapshot)
{
    const auto & atom_context{ context.atom_context_list.at(atom_index) };
    auto response_value{ static_cast<double>(atom_context.sample_entries.at(sample_index).response) };
    for (const auto & neighbor_sample : atom_context.sample_neighbor_list.at(sample_index))
    {
        response_value -= snapshot.at(neighbor_sample.atom_index).ResponseAtDistance(neighbor_sample.distance);
    }
    return response_value;
}

LocalPotentialSampleList BuildSecondStageAdjustedSamples(
    const SecondStageLocalFittingContext & context,
    std::size_t atom_index,
    const FittedGaussianSnapshot & snapshot)
{
    const auto & atom_context{ context.atom_context_list.at(atom_index) };
    LocalPotentialSampleList updated_list;
    updated_list.reserve(atom_context.sample_entries.size());
    for (std::size_t sample_index = 0; sample_index < atom_context.sample_entries.size(); sample_index++)
    {
        auto sample{ atom_context.sample_entries.at(sample_index) };
        sample.response = static_cast<float>(CalculateSecondStageAdjustedResponse(context, atom_index, sample_index, snapshot));
        updated_list.emplace_back(sample);
    }
    return updated_list;
}
std::optional<LocalFittingObjectiveSamples> CollectLocalFittingObjectiveSamples(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & state,
    const std::vector<LocalFittingObjectiveSampleRef> & sample_ref_list,
    const std::vector<std::size_t> & active_index_list)
{
    if (state.size() != context.AtomSize())
    {
        throw std::invalid_argument("Local fitting objective state size is inconsistent.");
    }
    const auto snapshot{ BuildFittedGaussianSnapshot(state) };

    LocalFittingObjectiveSamples objective_samples;
    objective_samples.residual_list.reserve(sample_ref_list.size());
    objective_samples.active_model_list.reserve(active_index_list.size());
    for (const auto active_index : active_index_list)
    {
        if (active_index >= snapshot.size())
        {
            throw std::invalid_argument("Local fitting objective active index is out of range.");
        }
        objective_samples.active_model_list.emplace_back(LocalFittingObjectiveAtomModel{
            active_index,
            snapshot.at(active_index)
        });
    }
    std::vector<double> response_list;
    response_list.reserve(sample_ref_list.size());
    for (const auto & sample_ref : sample_ref_list)
    {
        if (sample_ref.atom_index >= context.AtomSize())
        {
            throw std::invalid_argument("Local fitting objective sample atom index is out of range.");
        }
        const auto & atom_context{ context.atom_context_list.at(sample_ref.atom_index) };
        if (sample_ref.sample_index >= atom_context.sample_entries.size())
        {
            throw std::invalid_argument("Local fitting objective sample index is out of range.");
        }

        const auto & sample{ atom_context.sample_entries.at(sample_ref.sample_index) };
        const auto & target_model{ snapshot.at(sample_ref.atom_index) };
        const auto distance{ static_cast<double>(sample.point.distance) };
        const auto expected_response{ target_model.ResponseAtDistance(distance) };
        const auto response{
            CalculateSecondStageAdjustedResponse(
                context, sample_ref.atom_index, sample_ref.sample_index, snapshot)
        };
        const auto residual{ response - expected_response };
        if (!std::isfinite(response) || !std::isfinite(residual)) return std::nullopt;
        objective_samples.residual_list.emplace_back(residual);
        response_list.emplace_back(response);
    }
    if (objective_samples.residual_list.empty())
    {
        return std::nullopt;
    }

    const auto residual_scale{
        CalculateMedianAbsoluteDeviationScale(objective_samples.residual_list)
    };
    const auto response_scale_floor{
        kLocalFittingObjectiveResidualScaleFloorRatio *
        CalculateMedianAbsoluteDeviationScale(response_list)
    };
    objective_samples.scale_sample =
        std::max({ residual_scale, response_scale_floor, kRobustScaleMin });
    if (!std::isfinite(objective_samples.scale_sample)) return std::nullopt;
    return std::move(objective_samples);
}

std::optional<LocalFittingParameterPenaltyComponents>
CalculateLocalFittingParameterPenaltyComponents(
    const SecondStageLocalFittingContext & context,
    const LocalFittingObjectiveSamples & objective_samples,
    double objective_scale)
{
    if (!numeric_validation::IsFinitePositive(objective_scale))
    {
        return std::nullopt;
    }
    LocalFittingParameterPenaltyComponents components;
    const auto residual_scale_floor{ std::max(objective_scale, kRobustScaleMin) };
    for (std::size_t model_position = 0;
        model_position < objective_samples.active_model_list.size();
        model_position++)
    {
        const auto & atom_model{
            objective_samples.active_model_list.at(model_position)
        };
        const auto active_index{ atom_model.atom_index };
        if (active_index >= context.AtomSize()) return std::nullopt;

        const auto & model{ atom_model.model };
        if (!detail::IsValidSecondStageGaussianModel(model)) return std::nullopt;

        const auto prior_width{ context.atom_context_list.at(active_index).prior_width };
        if (!numeric_validation::IsFinitePositive(prior_width)) return std::nullopt;
        const auto normalized_width_difference{
            (std::log(model.GetWidth()) - std::log(prior_width)) /
            kLocalFittingWidthPriorLogScale
        };
        components.width_prior_penalty_sum +=
            normalized_width_difference * normalized_width_difference;

        const auto peak_signal{ model.SignalAtDistance(0.0) };
        const auto offset_peak{ model.GetOffset() * model.OffsetBasisAtDistance(0.0) };
        if (!std::isfinite(peak_signal) || !std::isfinite(offset_peak))
        {
            return std::nullopt;
        }
        const auto offset_ratio{
            std::abs(offset_peak) /
            std::max({ std::abs(peak_signal), residual_scale_floor, kRobustScaleMin })
        };
        if (!std::isfinite(offset_ratio)) return std::nullopt;
        const auto offset_excess{
            std::max(0.0, offset_ratio - kLocalFittingOffsetPeakRatioMax)
        };
        components.offset_plausibility_penalty_sum += offset_excess * offset_excess;
    }

    if (!std::isfinite(components.width_prior_penalty_sum) ||
        !std::isfinite(components.offset_plausibility_penalty_sum))
    {
        return std::nullopt;
    }
    return components;
}

std::optional<double> CalculateLocalFittingResidualObjective(
    const LocalFittingObjectiveSamples & objective_samples,
    double objective_scale)
{
    if (!numeric_validation::IsFinitePositive(objective_scale) ||
        objective_samples.residual_list.empty())
    {
        return std::nullopt;
    }

    double loss_sum{ 0.0 };
    for (const auto residual : objective_samples.residual_list)
    {
        const auto normalized_residual{ residual / objective_scale };
        loss_sum += algorithm::CalculateRobustLoss(
            kSecondStageRobustLossKind,
            normalized_residual,
            kRobustLossCutoffMultiplier);
    }
    const auto residual_objective{
        loss_sum / static_cast<double>(objective_samples.residual_list.size())
    };
    return std::isfinite(residual_objective) ?
        std::optional<double>{ residual_objective } :
        std::nullopt;
}

std::optional<detail::LocalFittingObjectiveBreakdown>
CalculateLocalFittingObjectiveBreakdown(
    const SecondStageLocalFittingContext & context,
    const LocalFittingObjectiveSamples & objective_samples,
    double objective_scale)
{
    const auto residual_objective{
        CalculateLocalFittingResidualObjective(
            objective_samples,
            objective_scale)
    };
    if (!residual_objective.has_value()) return std::nullopt;

    const auto parameter_penalty_components{
        CalculateLocalFittingParameterPenaltyComponents(
            context,
            objective_samples,
            objective_scale)
    };
    if (!parameter_penalty_components.has_value()) return std::nullopt;

    return detail::BuildLocalFittingMeanObjectiveBreakdown(
            *residual_objective,
            parameter_penalty_components->width_prior_penalty_sum,
            parameter_penalty_components->offset_plausibility_penalty_sum,
            objective_samples.active_model_list.size(),
            kLocalFittingWidthPriorPenaltyWeight,
            kLocalFittingOffsetPlausibilityPenaltyWeight);
}

std::optional<detail::LocalFittingObjectiveBreakdown>
EvaluateLocalFittingAuditObjective(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & state,
    LocalFittingBestAuditState & audit_state)
{
    auto objective_samples{
        CollectLocalFittingObjectiveSamples(
            context,
            state,
            audit_state.sample_ref_list,
            audit_state.atom_index_list)
    };
    if (!objective_samples.has_value()) return std::nullopt;

    const auto objective_scale{
        audit_state.fixed_objective_scale.value_or(
            objective_samples->scale_sample)
    };
    const auto objective{
        CalculateLocalFittingObjectiveBreakdown(
            context,
            *objective_samples,
            objective_scale)
    };
    if (!objective.has_value()) return std::nullopt;
    if (!audit_state.fixed_objective_scale.has_value())
    {
        audit_state.fixed_objective_scale = objective_scale;
    }
    return objective;
}

bool IsLocalFittingCombinedObjectiveAcceptable(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & candidate_state,
    const LocalFittingState & previous_state,
    LocalFittingBestAuditState & audit_state)
{
    if (!audit_state.fixed_objective_scale.has_value()) return false;

    const auto candidate_objective{
        EvaluateLocalFittingAuditObjective(context, candidate_state, audit_state)
    };
    const auto previous_objective{
        EvaluateLocalFittingAuditObjective(context, previous_state, audit_state)
    };
    return detail::IsLocalFittingAuditObjectiveAcceptableForProgress(
        candidate_objective.has_value() ?
            std::optional<double>{ candidate_objective->total_objective } :
            std::nullopt,
        previous_objective.has_value() ?
            std::optional<double>{ previous_objective->total_objective } :
            std::nullopt,
        audit_state.best.has_value() ?
            std::optional<double>{ audit_state.best->objective.total_objective } :
            std::nullopt,
        kLocalFittingConvergenceObjectiveRelativeTolerance);
}

void RejectLocalFittingCombinedCandidate(
    const LocalFittingState & previous_state,
    const std::vector<LocalFittingClusterKey> & cluster_key_list,
    LocalFittingCandidateSelection & selection)
{
    selection.assembled_state = previous_state;
    selection.accepted_key_list.clear();
    selection.rejected_key_list = cluster_key_list;
    selection.grow_trust_region_key_list.clear();
    selection.accepted_polish_cluster_count = 0;
}

bool TryUpdateLocalFittingBestAuditState(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & candidate_state,
    const std::optional<std::size_t> & accepted_iteration,
    LocalFittingBestAuditState & audit_state)
{
    const auto candidate_objective{
        EvaluateLocalFittingAuditObjective(context, candidate_state, audit_state)
    };
    if (!candidate_objective.has_value()) return false;
    LocalFittingAuditedState candidate{
        *candidate_objective,
        candidate_state,
        accepted_iteration
    };
    if (audit_state.best.has_value() &&
        !detail::IsBetterLocalFittingAuditObjective(
            candidate_objective->total_objective,
            audit_state.best->objective.total_objective,
            kLocalFittingObjectiveTieRelativeTolerance))
    {
        return false;
    }
    audit_state.best = std::move(candidate);
    return true;
}

LocalFittingBestAuditState BuildInitialLocalFittingBestAuditState(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & initial_state)
{
    LocalFittingBestAuditState audit_state;
    audit_state.atom_index_list.reserve(context.AtomSize());
    for (std::size_t atom_index = 0; atom_index < context.AtomSize(); atom_index++)
    {
        audit_state.atom_index_list.emplace_back(atom_index);
        const auto sample_size{
            context.atom_context_list.at(atom_index).sample_entries.size()
        };
        for (std::size_t sample_index = 0; sample_index < sample_size; sample_index++)
        {
            audit_state.sample_ref_list.emplace_back(
                LocalFittingObjectiveSampleRef{ atom_index, sample_index });
        }
    }
    static_cast<void>(TryUpdateLocalFittingBestAuditState(
        context,
        initial_state,
        std::nullopt,
        audit_state));
    return audit_state;
}

void ReconcileLocalFittingBestAuditTerminalFallback(
    const SecondStageLocalFittingContext & context,
    const std::vector<LocalFittingClusterKey> & terminal_key_list,
    const LocalFittingState & terminal_fallback_state,
    std::size_t accepted_iteration,
    LocalFittingBestAuditState & audit_state)
{
    if (terminal_key_list.empty()) return;
    if (terminal_fallback_state.size() != context.AtomSize())
    {
        throw std::invalid_argument(
            "Local fitting audit terminal fallback state sizes are inconsistent.");
    }

    for (const auto & key : terminal_key_list)
    {
        for (const auto atom_index : key)
        {
            if (atom_index >= context.AtomSize())
            {
                throw std::invalid_argument(
                    "Local fitting audit terminal atom index is out of range.");
            }
        }
    }

    const auto is_terminal_atom = [&](std::size_t atom_index)
    {
        return std::any_of(
            terminal_key_list.begin(),
            terminal_key_list.end(),
            [&](const LocalFittingClusterKey & key)
            {
                return std::binary_search(key.begin(), key.end(), atom_index);
            });
    };
    auto remaining_sample_ref_list{ audit_state.sample_ref_list };
    remaining_sample_ref_list.erase(
        std::remove_if(
            remaining_sample_ref_list.begin(),
            remaining_sample_ref_list.end(),
            [&](const LocalFittingObjectiveSampleRef & sample_ref)
            {
                return is_terminal_atom(sample_ref.atom_index);
            }),
        remaining_sample_ref_list.end());
    auto remaining_atom_index_list{ audit_state.atom_index_list };
    remaining_atom_index_list.erase(
        std::remove_if(
            remaining_atom_index_list.begin(),
            remaining_atom_index_list.end(),
            is_terminal_atom),
        remaining_atom_index_list.end());
    if (!remaining_sample_ref_list.empty() && !remaining_atom_index_list.empty())
    {
        audit_state.sample_ref_list = std::move(remaining_sample_ref_list);
        audit_state.atom_index_list = std::move(remaining_atom_index_list);
    }

    if (!audit_state.best.has_value()) return;
    auto reconciled_state{ terminal_fallback_state };
    const auto reconciled_objective{
        EvaluateLocalFittingAuditObjective(context, reconciled_state, audit_state)
    };
    if (!reconciled_objective.has_value())
    {
        audit_state.best.reset();
        return;
    }
    audit_state.best->objective = *reconciled_objective;
    audit_state.best->state = std::move(reconciled_state);
    audit_state.best->accepted_iteration = accepted_iteration;
}

std::vector<algorithm::ParameterChange> CalculateLocalFittingTransformedChanges(
    const LocalFittingState & current_state,
    const LocalFittingState & previous_state)
{
    if (current_state.size() != previous_state.size())
    {
        throw std::invalid_argument(
            "Local fitting transformed change input sizes are inconsistent.");
    }

    std::vector<algorithm::ParameterChange> change_list(current_state.size());
    for (size_t i = 0; i < current_state.size(); i++)
    {
        change_list.at(i) = detail::CalculateLocalFittingTransformedChange(
            current_state.at(i).mdpde.GetModel(),
            previous_state.at(i).mdpde.GetModel());
    }
    return change_list;
}

bool ContainsLocalFittingClusterKey(
    const std::vector<LocalFittingClusterKey> & key_list,
    const LocalFittingClusterKey & key)
{
    return std::find(key_list.begin(), key_list.end(), key) != key_list.end();
}

algorithm::ParameterChangeStats SummarizeLocalFittingTransformedChanges(
    const LocalFittingState & current_state,
    const LocalFittingState & previous_state,
    const std::vector<std::size_t> & index_list)
{
    if (current_state.size() != previous_state.size())
    {
        throw std::invalid_argument(
            "Local fitting transformed change input sizes are inconsistent.");
    }

    std::vector<algorithm::ParameterChange> change_list;
    change_list.reserve(index_list.size());
    for (const auto i : index_list)
    {
        if (i >= current_state.size())
        {
            throw std::invalid_argument(
                "Local fitting transformed change index is out of range.");
        }
        change_list.emplace_back(detail::CalculateLocalFittingTransformedChange(
            current_state.at(i).mdpde.GetModel(),
            previous_state.at(i).mdpde.GetModel()));
    }

    std::vector<std::size_t> local_index_list(change_list.size());
    for (std::size_t i = 0; i < local_index_list.size(); i++)
    {
        local_index_list.at(i) = i;
    }
    return algorithm::SummarizeParameterChangeStats(
        change_list,
        local_index_list,
        kLocalFittingChangePercentile);
}

bool IsLocalFittingTransformedChangeConverged(
    const algorithm::ParameterChangeStats & percentile_stats,
    const std::vector<double> & maximum_list)
{
    return detail::IsLocalFittingTransformedChangeConverged(
        percentile_stats,
        maximum_list,
        kLocalFittingTransformedChangeTolerance,
        kLocalFittingTransformedMaximumChangeTolerance);
}

bool IsLocalFittingTransformedPercentileConverged(
    const algorithm::ParameterChangeStats & stats)
{
    return std::all_of(
        stats.percentile_list.begin(),
        stats.percentile_list.end(),
        [](double value)
        {
            return std::isfinite(value) && value < kLocalFittingTransformedChangeTolerance;
        });
}

std::vector<std::size_t> CollectClusterSuspiciousAtomIndexes(
    const LocalFittingClusterKey & key,
    const std::vector<std::size_t> & suspicious_atom_index_list)
{
    std::vector<std::size_t> cluster_suspicious_atom_index_list;
    for (const auto atom_index : key)
    {
        if (std::binary_search(
                suspicious_atom_index_list.begin(),
                suspicious_atom_index_list.end(),
                atom_index))
        {
            cluster_suspicious_atom_index_list.emplace_back(atom_index);
        }
    }
    return cluster_suspicious_atom_index_list;
}

TerminalPersistentFailureMap UpdatePersistentTerminalFailureState(
    const std::vector<LocalFittingClusterKey> & cluster_key_list,
    const std::vector<LocalFittingClusterKey> & accepted_key_list,
    const std::vector<std::size_t> & suspicious_atom_index_list,
    const LocalFittingClusterHealthMap & health_by_key,
    const LocalFittingState & assembled_state,
    const LocalFittingState & previous_state,
    PersistentTerminalFailureStateMap & state_by_key)
{
    PersistentTerminalFailureStateMap next_state_by_key;
    TerminalPersistentFailureMap terminal_failure_by_key;
    for (const auto & key : cluster_key_list)
    {
        if (std::find(accepted_key_list.begin(), accepted_key_list.end(), key) ==
            accepted_key_list.end())
        {
            continue;
        }

        auto cluster_suspicious_atom_index_list{
            CollectClusterSuspiciousAtomIndexes(key, suspicious_atom_index_list)
        };
        std::optional<PersistentTerminalFailureReason> reason;
        if (!cluster_suspicious_atom_index_list.empty())
        {
            reason.emplace(std::move(cluster_suspicious_atom_index_list));
        }
        else
        {
            const auto health_iter{ health_by_key.find(key) };
            if (health_iter == health_by_key.end())
            {
                throw std::invalid_argument("Persistent joint-offset failure cluster health is missing.");
            }
            const auto status{ health_iter->second.joint_offset_status };
            if (!IsJointOffsetSolveHardFailure(status)) continue;
            reason.emplace(status);
        }

        const auto transformed_change_stats{
            SummarizeLocalFittingTransformedChanges(
                assembled_state,
                previous_state,
                key)
        };
        if (!IsLocalFittingTransformedPercentileConverged(transformed_change_stats))
        {
            continue;
        }

        PersistentTerminalFailureState next_state{ std::move(*reason), 1 };
        const auto previous_iter{ state_by_key.find(key) };
        if (previous_iter != state_by_key.end() && previous_iter->second.reason == next_state.reason)
        {
            next_state.stable_iteration_count = previous_iter->second.stable_iteration_count + 1;
        }

        if (next_state.stable_iteration_count >= kPersistentTerminalFailureIterationLimit)
        {
            terminal_failure_by_key.emplace(key, std::move(next_state.reason));
            continue;
        }
        next_state_by_key.emplace(key, std::move(next_state));
    }
    state_by_key = std::move(next_state_by_key);
    return terminal_failure_by_key;
}

void ApplyTerminalFallbackClusters(
    const std::vector<LocalFittingClusterKey> & terminal_key_list,
    const LocalFittingState & previous_state,
    std::vector<char> & terminal_atom_mask,
    LocalFittingState & assembled_state)
{
    if (previous_state.size() != terminal_atom_mask.size() ||
        assembled_state.size() != terminal_atom_mask.size())
    {
        throw std::invalid_argument("Local fitting terminal fallback state sizes are inconsistent.");
    }
    for (const auto & key : terminal_key_list)
    {
        for (const auto atom_index : key)
        {
            if (atom_index >= terminal_atom_mask.size())
            {
                throw std::invalid_argument(
                    "Local fitting terminal fallback atom index is out of range.");
            }
            terminal_atom_mask.at(atom_index) = 1;
            assembled_state.at(atom_index) = previous_state.at(atom_index);
        }
    }
}

std::vector<std::size_t> BuildEligibleLocalFittingActiveIndexList(
    std::size_t atom_size,
    const std::vector<char> & terminal_atom_mask)
{
    if (atom_size != terminal_atom_mask.size())
    {
        throw std::invalid_argument(
            "Local fitting terminal fallback active index size is inconsistent.");
    }
    std::vector<std::size_t> active_index_list;
    active_index_list.reserve(atom_size);
    for (std::size_t atom_index = 0; atom_index < atom_size; atom_index++)
    {
        if (terminal_atom_mask.at(atom_index) == 0)
        {
            active_index_list.emplace_back(atom_index);
        }
    }
    return active_index_list;
}

std::vector<Eigen::VectorXd> BuildLocalFittingTransformedEstimationList(
    const LocalFittingState & state)
{
    std::vector<Eigen::VectorXd> transformed_estimation_list;
    transformed_estimation_list.reserve(state.size());
    for (const auto & result : state)
    {
        const auto transformed{
            detail::EncodeLocalFittingTransformedCoordinates(result.mdpde.GetModel())
        };
        if (!transformed.has_value())
        {
            throw std::invalid_argument("Local fitting state has invalid transformed coordinates.");
        }
        transformed_estimation_list.emplace_back(*transformed);
    }
    return transformed_estimation_list;
}

std::optional<LocalFittingJointPolishStep> BuildLocalFittingJointPolishStep(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & candidate_state,
    const LocalFittingClusterKey & key,
    const std::vector<LocalFittingObjectiveSampleRef> & sample_ref_list,
    double ridge_ratio,
    const std::vector<double> & ridge_multiplier_list)
{
    if (candidate_state.size() != context.AtomSize() ||
        ridge_multiplier_list.size() != context.AtomSize())
    {
        throw std::invalid_argument("Local fitting joint polish input sizes are inconsistent.");
    }
    if (key.empty() || sample_ref_list.empty() || !std::isfinite(ridge_ratio) || ridge_ratio <= 0.0)
    {
        return std::nullopt;
    }

    constexpr auto parameter_size{
        static_cast<std::size_t>(detail::kTransformedChangeSize)
    };
    const auto column_count{
        static_cast<Eigen::Index>(key.size() * parameter_size)
    };
    std::vector<int> column_base_by_atom_index(context.AtomSize(), -1);
    for (std::size_t local_position = 0; local_position < key.size(); local_position++)
    {
        const auto atom_index{ key.at(local_position) };
        if (atom_index >= context.AtomSize())
        {
            throw std::invalid_argument("Local fitting joint polish atom index is out of range.");
        }
        column_base_by_atom_index.at(atom_index) = static_cast<int>(local_position * parameter_size);
    }

    const auto snapshot{ BuildFittedGaussianSnapshot(candidate_state) };
    std::vector<Eigen::Triplet<double>> triplet_list;
    std::vector<double> residual_list;
    residual_list.reserve(sample_ref_list.size());
    Eigen::VectorXd column_square_sum{ Eigen::VectorXd::Zero(column_count) };
    for (const auto & sample_ref : sample_ref_list)
    {
        if (sample_ref.atom_index >= context.AtomSize())
        {
            throw std::invalid_argument("Local fitting joint polish sample atom index is out of range.");
        }
        const auto & atom_context{
            context.atom_context_list.at(sample_ref.atom_index)
        };
        if (sample_ref.sample_index >= atom_context.sample_entries.size())
        {
            throw std::invalid_argument("Local fitting joint polish sample index is out of range.");
        }
        const auto & sample{ atom_context.sample_entries.at(sample_ref.sample_index) };
        if (!std::isfinite(static_cast<double>(sample.response))) return std::nullopt;

        const auto row_index{ static_cast<Eigen::Index>(residual_list.size()) };
        double predicted_response{ 0.0 };
        const auto append_model = [&](std::size_t atom_index, double distance) -> bool
        {
            const auto evaluation{
                detail::EvaluateLocalFittingTransformedResponse(snapshot.at(atom_index), distance)
            };
            if (!evaluation.has_value()) return false;
            predicted_response += evaluation->response;

            const auto column_base{ column_base_by_atom_index.at(atom_index) };
            if (column_base < 0) return true;
            for (std::size_t parameter_index = 0; parameter_index < parameter_size; parameter_index++)
            {
                const auto column_index{
                    static_cast<Eigen::Index>(column_base + static_cast<int>(parameter_index))
                };
                const auto derivative{
                    evaluation->jacobian(static_cast<Eigen::Index>(parameter_index))
                };
                if (std::abs(derivative) <= std::numeric_limits<double>::epsilon()) continue;
                triplet_list.emplace_back(row_index, column_index, derivative);
                column_square_sum(column_index) += derivative * derivative;
            }
            return true;
        };

        if (!append_model(
                sample_ref.atom_index,
                static_cast<double>(sample.point.distance)))
        {
            return std::nullopt;
        }
        for (const auto & neighbor_sample :
            atom_context.sample_neighbor_list.at(sample_ref.sample_index))
        {
            if (!append_model(neighbor_sample.atom_index, neighbor_sample.distance))
            {
                return std::nullopt;
            }
        }

        const auto residual{ static_cast<double>(sample.response) - predicted_response };
        if (!std::isfinite(residual)) return std::nullopt;
        residual_list.emplace_back(residual);
    }

    const auto row_count{ static_cast<Eigen::Index>(residual_list.size()) };
    algorithm::WeightedRidgeSystem system;
    system.design_matrix.resize(row_count, column_count);
    system.design_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
    system.response = Eigen::VectorXd::Zero(row_count);
    for (Eigen::Index row_index = 0; row_index < row_count; row_index++)
    {
        system.response(row_index) = residual_list.at(static_cast<std::size_t>(row_index));
    }
    system.previous_parameter = Eigen::VectorXd::Zero(column_count);
    system.ridge_diagonal = Eigen::VectorXd::Zero(column_count);

    const auto conditioning{
        detail::EvaluateLocalFittingJointOffsetConditioning(
            system.design_matrix,
            kJointOffsetConditioningPivotRatioThreshold)
    };
    const auto conditioning_multiplier{
        conditioning.guard_required ? kCollinearJointOffsetRidgeMultiplier : 1.0
    };
    for (Eigen::Index column_index = 0; column_index < column_count; column_index++)
    {
        const auto local_position{
            static_cast<std::size_t>(column_index) / parameter_size
        };
        const auto atom_index{ key.at(local_position) };
        const auto atom_multiplier{ ridge_multiplier_list.at(atom_index) };
        if (!std::isfinite(atom_multiplier) || atom_multiplier <= 0.0)
        {
            throw std::invalid_argument(
                "Local fitting joint polish ridge multiplier must be positive and finite.");
        }
        const auto square_sum{ column_square_sum(column_index) };
        const auto base_ridge{
            square_sum > std::numeric_limits<double>::epsilon()
                ? ridge_ratio * square_sum
                : ridge_ratio / kJointOffsetRidgeRatio
        };
        system.ridge_diagonal(column_index) =
            std::max(atom_multiplier, conditioning_multiplier) * base_ridge;
    }

    const auto residual_scale{
        std::max(CalculateMedianAbsoluteDeviationScale(residual_list), kRobustScaleMin)
    };
    if (!std::isfinite(residual_scale)) return std::nullopt;
    Eigen::VectorXd weight{ Eigen::VectorXd::Ones(row_count) };
    for (Eigen::Index row_index = 0; row_index < row_count; row_index++)
    {
        weight(row_index) = algorithm::CalculateRobustWeight(
            kSecondStageRobustLossKind,
            system.response(row_index),
            residual_scale,
            kRobustLossCutoffMultiplier);
    }

    algorithm::WeightedRidgeSolver solver{ system };
    Eigen::VectorXd direction;
    if (!solver.Solve(system, weight, direction) || !direction.allFinite())
    {
        return std::nullopt;
    }

    auto transformed_estimation_list{
        BuildLocalFittingTransformedEstimationList(candidate_state)
    };
    for (std::size_t local_position = 0; local_position < key.size(); local_position++)
    {
        const auto direction_offset{
            static_cast<Eigen::Index>(local_position * parameter_size)
        };
        transformed_estimation_list.at(key.at(local_position)) +=
            direction.segment(direction_offset, static_cast<Eigen::Index>(parameter_size));
    }
    return LocalFittingJointPolishStep{
        std::move(transformed_estimation_list),
        direction.cwiseAbs().maxCoeff() < kLocalFittingTransformedChangeTolerance
    };
}

std::optional<LocalFittingState> BuildLocalFittingCandidateState(
    const LocalFittingState & previous_state,
    const std::vector<Eigen::VectorXd> & candidate_transformed_estimation_list,
    const LocalFittingState & uncertainty_state,
    const std::vector<std::size_t> & active_index_list,
    double damping)
{
    if (!std::isfinite(damping) || damping <= 0.0 || damping > 1.0)
    {
        throw std::invalid_argument("Local fitting candidate damping is out of range.");
    }
    if (candidate_transformed_estimation_list.size() != previous_state.size() ||
        uncertainty_state.size() != previous_state.size())
    {
        throw std::invalid_argument("Local fitting candidate input sizes are inconsistent.");
    }

    auto candidate_state{ previous_state };
    for (const auto active_index : active_index_list)
    {
        if (active_index >= candidate_state.size())
        {
            throw std::invalid_argument("Local fitting candidate active index is out of range.");
        }
        const auto previous_transformed_estimation{
            detail::EncodeLocalFittingTransformedCoordinates(previous_state.at(active_index).mdpde.GetModel())
        };
        if (!previous_transformed_estimation.has_value())
        {
            throw std::invalid_argument(
                "Local fitting previous state has invalid transformed coordinates.");
        }
        const auto & candidate_transformed_estimation{
            candidate_transformed_estimation_list.at(active_index)
        };
        if (candidate_transformed_estimation.size() != static_cast<Eigen::Index>(detail::kTransformedChangeSize))
        {
            return std::nullopt;
        }
        if (!candidate_transformed_estimation.allFinite())
        {
            return std::nullopt;
        }
        if ((candidate_transformed_estimation.array() == previous_transformed_estimation->array()).all())
        {
            auto & result{ candidate_state.at(active_index) };
            result.mdpde = GaussianModel3DWithUncertainty{
                previous_state.at(active_index).mdpde.GetModel(),
                uncertainty_state.at(active_index).mdpde.GetStandardDeviationModel()
            };
            continue;
        }
        const auto damped_transformed_estimation{
            (*previous_transformed_estimation +
                damping * (candidate_transformed_estimation - *previous_transformed_estimation)).eval()
        };
        const auto damped_model{
            detail::DecodeLocalFittingTransformedCoordinates(damped_transformed_estimation)
        };
        if (!damped_model.has_value())
        {
            return std::nullopt;
        }
        if (!detail::IsValidSecondStageGaussianModel(*damped_model))
        {
            return std::nullopt;
        }
        auto & result{ candidate_state.at(active_index) };
        result.mdpde = GaussianModel3DWithUncertainty{
            *damped_model,
            uncertainty_state.at(active_index).mdpde.GetStandardDeviationModel()
        };
    }
    return std::move(candidate_state);
}

LocalFittingClusterObjectiveState
BuildInitialLocalFittingClusterQualityState(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & previous_state,
    const LocalFittingClusterKey & key,
    const std::vector<LocalFittingObjectiveSampleRef> & objective_sample_ref_list)
{
    auto initial_objective_samples{
        CollectLocalFittingObjectiveSamples(
            context,
            previous_state,
            objective_sample_ref_list,
            key)
    };
    std::optional<double> initial_objective;
    if (initial_objective_samples.has_value())
    {
        const auto objective{
            CalculateLocalFittingObjectiveBreakdown(
                context,
                *initial_objective_samples,
                initial_objective_samples->scale_sample)
        };
        if (objective.has_value())
        {
            initial_objective = objective->total_objective;
        }
    }
    algorithm::ParameterChangeStats initial_change_stats{
        std::vector<double>(
            static_cast<std::size_t>(GaussianModel3D::ParameterSize()),
            0.0)
    };
    return LocalFittingClusterObjectiveState{
        initial_objective_samples.has_value() ?
            std::optional<double>{ initial_objective_samples->scale_sample } :
            std::nullopt,
        LocalFittingTrackedCandidate{
            std::move(initial_change_stats),
            std::move(initial_objective_samples)
        },
        initial_objective.has_value()
    };
}

void ReconcileLocalFittingClusterObjectiveState(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & previous_state,
    const detail::LocalFittingCouplingPartition & partition,
    const std::vector<LocalFittingClusterKey> & key_list,
    LocalFittingClusterObjectiveStateMap & state_by_key)
{
    LocalFittingClusterObjectiveStateMap next_state_by_key;
    for (const auto & key : key_list)
    {
        auto state_iter{ state_by_key.find(key) };
        if (state_iter != state_by_key.end())
        {
            next_state_by_key.emplace(key, std::move(state_iter->second));
            continue;
        }
        const auto sample_iter{ partition.sample_id_list_by_key.find(key) };
        if (sample_iter == partition.sample_id_list_by_key.end())
        {
            throw std::invalid_argument("Local fitting cluster objective samples are missing.");
        }
        next_state_by_key.emplace(
            key,
            BuildInitialLocalFittingClusterQualityState(
                context,
                previous_state,
                key,
                sample_iter->second));
    }
    state_by_key = std::move(next_state_by_key);
}

bool AreLocalFittingObjectiveReferencesLocked(
    const LocalFittingClusterObjectiveStateMap & state_by_key,
    const std::vector<LocalFittingClusterKey> & key_list)
{
    for (const auto & key : key_list)
    {
        const auto iter{ state_by_key.find(key) };
        if (iter == state_by_key.end())
        {
            throw std::invalid_argument("Local fitting cluster objective state is missing.");
        }
        if (iter->second.scale_tracker.GetCommittedReference().has_value() &&
            !iter->second.scale_tracker.IsLocked())
        {
            return false;
        }
    }
    return true;
}

bool TryCommitLocalFittingClusterCandidate(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & candidate_state,
    const LocalFittingState & previous_state,
    const LocalFittingClusterKey & key,
    const std::vector<LocalFittingObjectiveSampleRef> & objective_sample_ref_list,
    bool requires_strict_improvement,
    LocalFittingClusterObjectiveStateMap & cluster_quality_state,
    LocalFittingObjectiveAttemptDiagnostic & diagnostic)
{
    const auto transformed_change_stats{
        SummarizeLocalFittingTransformedChanges(candidate_state, previous_state, key)
    };
    auto state_iter{ cluster_quality_state.find(key) };
    if (state_iter == cluster_quality_state.end())
    {
        throw std::invalid_argument("Local fitting cluster objective state is missing.");
    }
    auto & state{ state_iter->second };
    auto objective_scale{ state.scale_tracker.GetCommittedReference() };
    auto has_objective_reference{ objective_scale.has_value() };
    std::optional<LocalFittingObjectiveSamples> candidate_objective_samples;
    std::optional<double> objective_scale_sample;
    if (objective_scale.has_value())
    {
        candidate_objective_samples = CollectLocalFittingObjectiveSamples(
            context,
            candidate_state,
            objective_sample_ref_list,
            key);
        if (candidate_objective_samples.has_value())
        {
            objective_scale_sample = candidate_objective_samples->scale_sample;
            objective_scale = state.scale_tracker.GetProvisionalReference(
                *objective_scale_sample);
            has_objective_reference = objective_scale.has_value();
        }
    }

    std::optional<double> candidate_quality_objective;
    std::optional<double> previous_quality_objective;
    std::optional<double> best_quality_objective;
    if (objective_scale.has_value())
    {
        diagnostic.objective_scale = *objective_scale;
        if (candidate_objective_samples.has_value() &&
            state.previous.objective_samples.has_value())
        {
            diagnostic.candidate_objective =
                CalculateLocalFittingObjectiveBreakdown(
                    context,
                    *candidate_objective_samples,
                    *objective_scale);
            if (diagnostic.candidate_objective.has_value())
            {
                diagnostic.status = LocalFittingObjectiveAttemptDiagnosticStatus::Scored;
                candidate_quality_objective =
                    diagnostic.candidate_objective->total_objective;
                diagnostic.previous_objective =
                    CalculateLocalFittingObjectiveBreakdown(
                        context,
                        *state.previous.objective_samples,
                        *objective_scale);
                if (diagnostic.previous_objective.has_value())
                {
                    previous_quality_objective =
                        diagnostic.previous_objective->total_objective;
                }
                if (state.best.has_value() &&
                    state.best->objective_samples.has_value())
                {
                    diagnostic.best_objective =
                        CalculateLocalFittingObjectiveBreakdown(
                            context,
                            *state.best->objective_samples,
                            *objective_scale);
                    if (diagnostic.best_objective.has_value())
                    {
                        best_quality_objective =
                            diagnostic.best_objective->total_objective;
                    }
                }
            }
        }

        const auto evaluate_tracked_objective = [&](const LocalFittingTrackedCandidate & tracked)
            -> std::optional<double>
        {
            if (!tracked.objective_samples.has_value()) return std::nullopt;
            const auto objective{
                CalculateLocalFittingObjectiveBreakdown(
                    context,
                    *tracked.objective_samples,
                    *objective_scale)
            };
            return objective.has_value() ?
                std::optional<double>{ objective->total_objective } :
                std::nullopt;
        };
        if (!previous_quality_objective.has_value())
        {
            previous_quality_objective = evaluate_tracked_objective(state.previous);
        }
        if (state.best.has_value() && !best_quality_objective.has_value())
        {
            best_quality_objective = evaluate_tracked_objective(*state.best);
        }
    }

    const auto is_objective_deteriorated = [](
        const std::optional<double> & candidate,
        const std::optional<double> & reference)
    {
        if (!reference.has_value()) return false;
        if (!candidate.has_value()) return true;
        const auto scale{ std::max(std::abs(*reference), 1.0) };
        return *candidate > *reference +
            kLocalFittingConvergenceObjectiveRelativeTolerance * scale;
    };
    if (has_objective_reference)
    {
        diagnostic.rejected_by_previous = is_objective_deteriorated(
            candidate_quality_objective,
            previous_quality_objective);
        diagnostic.rejected_by_best = is_objective_deteriorated(
            candidate_quality_objective,
            best_quality_objective);
        if (diagnostic.rejected_by_previous || diagnostic.rejected_by_best)
        {
            return false;
        }
    }
    if (requires_strict_improvement &&
        (!diagnostic.candidate_objective.has_value() ||
            !diagnostic.previous_objective.has_value() ||
            !detail::IsBetterLocalFittingAuditObjective(
                diagnostic.candidate_objective->total_objective,
                diagnostic.previous_objective->total_objective,
                kLocalFittingObjectiveTieRelativeTolerance)))
    {
        return false;
    }

    if (candidate_quality_objective.has_value() && objective_scale_sample.has_value())
    {
        state.scale_tracker.CommitScaleSample(*objective_scale_sample);
    }
    LocalFittingTrackedCandidate candidate{
        transformed_change_stats,
        std::move(candidate_objective_samples)
    };
    auto is_better_than_best{ !state.best.has_value() };
    if (state.best.has_value())
    {
        if (candidate_quality_objective.has_value() != best_quality_objective.has_value())
        {
            is_better_than_best = candidate_quality_objective.has_value();
        }
        else if (candidate_quality_objective.has_value())
        {
            if (detail::IsBetterLocalFittingAuditObjective(
                    *candidate_quality_objective,
                    *best_quality_objective,
                    kLocalFittingObjectiveTieRelativeTolerance))
            {
                is_better_than_best = true;
            }
            else if (detail::IsBetterLocalFittingAuditObjective(
                         *best_quality_objective,
                         *candidate_quality_objective,
                         kLocalFittingObjectiveTieRelativeTolerance))
            {
                is_better_than_best = false;
            }
            else
            {
                is_better_than_best =
                    algorithm::GetMaximumParameterChange(
                        candidate.transformed_change_stats) <
                    algorithm::GetMaximumParameterChange(
                        state.best->transformed_change_stats);
            }
        }
        else
        {
            is_better_than_best =
                algorithm::GetMaximumParameterChange(
                    candidate.transformed_change_stats) <
                algorithm::GetMaximumParameterChange(
                    state.best->transformed_change_stats);
        }
    }
    state.previous = candidate;
    if (candidate_quality_objective.has_value() && is_better_than_best)
    {
        state.best = std::move(candidate);
    }
    return true;
}

bool ShouldGrowLocalFittingTrustRegion(
    const LocalFittingObjectiveAttemptDiagnostic & diagnostic)
{
    return diagnostic.candidate_objective.has_value() &&
        diagnostic.previous_objective.has_value() &&
        diagnostic.trust_region_step_norm >=
            kLocalFittingTrustRegionBoundaryRatio * diagnostic.trust_region_radius &&
        detail::IsBetterLocalFittingAuditObjective(
            diagnostic.candidate_objective->total_objective,
            diagnostic.previous_objective->total_objective,
            kLocalFittingObjectiveTieRelativeTolerance);
}

LocalFittingCandidateSelection SelectLocalFittingClusterCandidates(
    const SecondStageLocalFittingContext & context,
    const detail::LocalFittingCouplingPartition & partition,
    const std::vector<LocalFittingClusterKey> & key_list,
    const std::vector<LocalFittingClusterKey> & polish_eligible_key_list,
    const LocalFittingState & previous_state,
    const LocalFittingState & raw_state,
    const std::vector<Eigen::VectorXd> & previous_transformed_estimation_list,
    const std::vector<Eigen::VectorXd> & raw_transformed_estimation_list,
    const std::vector<double> & ridge_multiplier_list,
    LocalFittingClusterObjectiveStateMap & cluster_quality_state,
    const detail::LocalFittingTrustRegionStateSet & trust_region_state)
{
    LocalFittingCandidateSelection selection;
    selection.assembled_state = previous_state;
    for (const auto & key : key_list)
    {
        const auto sample_iter{ partition.sample_id_list_by_key.find(key) };
        if (sample_iter == partition.sample_id_list_by_key.end())
        {
            throw std::invalid_argument("Local fitting cluster objective samples are missing.");
        }
        const auto trust_region_radius{ trust_region_state.GetRadius(key) };
        const auto trust_region_damping{
            detail::LimitLocalFittingTrustRegionDamping(
                previous_transformed_estimation_list,
                raw_transformed_estimation_list,
                key,
                kLocalFittingTrustRegionParameterScale,
                1.0,
                trust_region_radius)
        };
        LocalFittingObjectiveAttemptDiagnostic base_diagnostic;
        base_diagnostic.effective_damping = trust_region_damping.effective_damping;
        base_diagnostic.trust_region_radius = trust_region_radius;
        base_diagnostic.trust_region_step_norm = trust_region_damping.step_norm;
        auto base_state{
            BuildLocalFittingCandidateState(
                previous_state,
                raw_transformed_estimation_list,
                raw_state,
                key,
                trust_region_damping.effective_damping)
        };
        if (!base_state.has_value())
        {
            base_diagnostic.status = LocalFittingObjectiveAttemptDiagnosticStatus::InvalidModel;
            selection.rejected_key_list.emplace_back(key);
            selection.rejected_cluster_diagnostic_list.emplace_back(
                LocalFittingRejectedClusterDiagnostic{
                    key,
                    std::move(base_diagnostic)
                });
            continue;
        }
        if (!TryCommitLocalFittingClusterCandidate(
                context,
                *base_state,
                previous_state,
                key,
                sample_iter->second,
                false,
                cluster_quality_state,
                base_diagnostic))
        {
            selection.rejected_key_list.emplace_back(key);
            selection.rejected_cluster_diagnostic_list.emplace_back(
                LocalFittingRejectedClusterDiagnostic{
                    key,
                    std::move(base_diagnostic)
                });
            continue;
        }
        selection.accepted_key_list.emplace_back(key);
        if (ShouldGrowLocalFittingTrustRegion(base_diagnostic))
        {
            selection.grow_trust_region_key_list.emplace_back(key);
        }
        for (const auto active_index : key)
        {
            selection.assembled_state.at(active_index) = base_state->at(active_index);
        }

        if (!ContainsLocalFittingClusterKey(polish_eligible_key_list, key)) continue;
        const auto polish_step{
            BuildLocalFittingJointPolishStep(
                context,
                *base_state,
                key,
                sample_iter->second,
                kJointOffsetRidgeRatio,
                ridge_multiplier_list)
        };
        if (!polish_step.has_value()) continue;
        if (polish_step->is_stationary) continue;
        const auto base_transformed_estimation_list{
            BuildLocalFittingTransformedEstimationList(*base_state)
        };
        const auto polish_damping{
            detail::LimitLocalFittingTrustRegionSubstepDamping(
                previous_transformed_estimation_list,
                base_transformed_estimation_list,
                polish_step->transformed_estimation_list,
                key,
                kLocalFittingTrustRegionParameterScale,
                1.0,
                trust_region_radius)
        };
        if (polish_damping.effective_damping <= 0.0) continue;
        auto polished_state{
            BuildLocalFittingCandidateState(
                *base_state,
                polish_step->transformed_estimation_list,
                *base_state,
                key,
                polish_damping.effective_damping)
        };
        if (!polished_state.has_value()) continue;
        LocalFittingObjectiveAttemptDiagnostic polish_diagnostic;
        polish_diagnostic.trust_region_radius = trust_region_radius;
        polish_diagnostic.trust_region_step_norm = polish_damping.step_norm;
        if (!TryCommitLocalFittingClusterCandidate(
                context,
                *polished_state,
                *base_state,
                key,
                sample_iter->second,
                true,
                cluster_quality_state,
                polish_diagnostic))
        {
            continue;
        }
        for (const auto active_index : key)
        {
            selection.assembled_state.at(active_index) =
                polished_state->at(active_index);
        }
        selection.accepted_polish_cluster_count++;
        if (ShouldGrowLocalFittingTrustRegion(polish_diagnostic) &&
            !ContainsLocalFittingClusterKey(
                selection.grow_trust_region_key_list,
                key))
        {
            selection.grow_trust_region_key_list.emplace_back(key);
        }
    }
    return selection;
}

std::optional<LocalAtomRefitResult> FitAtomWithJointOffsetFallback(
    const SecondStageLocalFittingContext & context,
    std::size_t atom_index,
    const LocalGaussianResult & previous_result,
    const FittedGaussianSnapshot & offset_snapshot,
    const FitOptions & options)
{
    auto sample_entries{ BuildSecondStageAdjustedSamples(context, atom_index, offset_snapshot) };
    const auto & offset_model{ offset_snapshot.at(atom_index) };
    const auto & previous_model{ previous_result.mdpde.GetModel() };
    const auto is_acceptable = [&](const GaussianModel3D & model)
    {
        return detail::IsValidSecondStageGaussianModel(model) &&
            CanBuildFiniteZeroOffsetSamples(sample_entries, model) &&
            !IsSuspiciousJointOffset(
                sample_entries,
                previous_model,
                model,
                options);
    };
    try
    {
        auto candidate_result{
            EstimateLocalGaussian(
                sample_entries,
                context.atom_context_list.at(atom_index).alpha_r,
                options,
                offset_model)
        };
        const auto candidate_model{ candidate_result.mdpde.GetModel() };
        if (is_acceptable(candidate_model))
        {
            const auto is_stationarity_eligible{
                candidate_result.fit_result.has_value() &&
                detail::IsLocalGaussianRefitStatusStationarityEligible(candidate_result.fit_result->status)
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
    if (!is_acceptable(result.mdpde.GetModel()))
    {
        return std::nullopt;
    }
    return LocalAtomRefitResult{ std::move(result), false };
}

std::vector<std::size_t> ExpandSuspiciousOffsetClusters(
    const ActiveCouplingGraph & active_coupling_graph,
    const std::vector<std::size_t> & suspicious_offset_seed_position_list,
    std::vector<char> & suspicious_offset_mask)
{
    if (!active_coupling_graph.empty() &&
        active_coupling_graph.size() != suspicious_offset_mask.size())
    {
        throw std::invalid_argument("Suspicious offset cluster graph size is inconsistent.");
    }

    std::vector<char> visited(suspicious_offset_mask.size(), 0);
    std::vector<std::pair<std::size_t, std::size_t>> queue;
    queue.reserve(suspicious_offset_mask.size());
    for (const auto seed_position : suspicious_offset_seed_position_list)
    {
        if (seed_position >= suspicious_offset_mask.size())
        {
            throw std::invalid_argument("Suspicious offset seed position is out of range.");
        }
        if (visited.at(seed_position) != 0) continue;
        visited.at(seed_position) = 1;
        queue.emplace_back(seed_position, 0);
    }

    std::vector<std::size_t> added_position_list;
    for (std::size_t queue_index = 0; queue_index < queue.size(); queue_index++)
    {
        const auto [active_position, depth]{ queue.at(queue_index) };
        if (suspicious_offset_mask.at(active_position) == 0)
        {
            suspicious_offset_mask.at(active_position) = 1;
            added_position_list.emplace_back(active_position);
        }
        if (active_coupling_graph.empty() || depth >= kSuspiciousOffsetClusterMaxDepth) continue;

        for (const auto & edge : active_coupling_graph.at(active_position))
        {
            if (edge.neighbor_index >= active_coupling_graph.size())
            {
                throw std::invalid_argument("Suspicious offset cluster adjacency index is out of range.");
            }
            if (!std::isfinite(edge.overlap) || edge.overlap < kSuspiciousOffsetClusterMinimumOverlap) continue;
            if (visited.at(edge.neighbor_index) != 0) continue;
            visited.at(edge.neighbor_index) = 1;
            queue.emplace_back(edge.neighbor_index, depth + 1);
        }
    }
    return added_position_list;
}

void RollBackSuspiciousOffsetClusters(
    const SecondStageLocalFittingContext & context,
    const std::vector<std::size_t> & active_index_list,
    const LocalFittingState & previous_state,
    const std::vector<std::size_t> & suspicious_active_position_list,
    FittedGaussianSnapshot & current_snapshot,
    LocalFittingState & iteration_state)
{
    for (const auto active_position : suspicious_active_position_list)
    {
        if (active_position >= active_index_list.size())
        {
            throw std::invalid_argument("Suspicious offset rollback position is out of range.");
        }

        const auto state_index{ active_index_list.at(active_position) };
        if (state_index >= context.AtomSize() ||
            state_index >= previous_state.size() ||
            state_index >= iteration_state.size())
        {
            throw std::invalid_argument("Suspicious offset rollback state index is out of range.");
        }
        const auto & previous_model{
            previous_state.at(state_index).mdpde.GetModel()
        };
        current_snapshot.at(state_index) = previous_model;
        iteration_state.at(state_index) = previous_state.at(state_index);
    }
}

LocalFittingIterationResult RunLocalFittingIteration(
    const SecondStageLocalFittingContext & context,
    const std::vector<std::size_t> & active_index_list,
    const std::vector<LocalFittingClusterKey> & cluster_key_list,
    const LocalFittingState & previous_state,
    const FitOptions & options,
    double ridge_ratio,
    const std::vector<double> & ridge_multiplier_list)
{
    const auto selected_atom_size{ context.AtomSize() };
    if (previous_state.size() != selected_atom_size ||
        ridge_multiplier_list.size() != selected_atom_size)
    {
        throw std::invalid_argument("Local fitting iteration input sizes are inconsistent.");
    }
    auto current_snapshot{ BuildFittedGaussianSnapshot(previous_state) };
    auto clustered_joint_offset_result{
        EstimateClusteredJointOffsets(
            context,
            active_index_list,
            cluster_key_list,
            current_snapshot,
            ridge_ratio,
            ridge_multiplier_list,
            !options.quiet_mode && Logger::GetLogLevel() >= LogLevel::Debug)
    };
    ApplyJointOffsetsToSnapshot(
        active_index_list,
        clustered_joint_offset_result.offset,
        current_snapshot);
    auto iteration_state{ previous_state };
    std::vector<char> suspicious_offset_mask(active_index_list.size(), 0);
    std::vector<std::size_t> pre_refit_suspicious_seed_position_list;

    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto state_index{ active_index_list.at(i) };
        const auto & previous_model{
            previous_state.at(state_index).mdpde.GetModel()
        };
        const auto & offset_model{ current_snapshot.at(state_index) };
        if (!IsSuspiciousJointOffset(
                context.atom_context_list.at(state_index).sample_entries,
                previous_model,
                offset_model,
                options))
        {
            continue;
        }
        pre_refit_suspicious_seed_position_list.emplace_back(i);
    }
    const auto pre_refit_suspicious_position_list{
        ExpandSuspiciousOffsetClusters(
            clustered_joint_offset_result.active_coupling_graph,
            pre_refit_suspicious_seed_position_list,
            suspicious_offset_mask)
    };
    RollBackSuspiciousOffsetClusters(
        context,
        active_index_list,
        previous_state,
        pre_refit_suspicious_position_list,
        current_snapshot,
        iteration_state);

    const auto refit_snapshot{ current_snapshot };
    std::vector<std::size_t> post_refit_suspicious_seed_position_list;
    for (size_t i = 0; i < active_index_list.size(); i++)
    {
        if (suspicious_offset_mask.at(i) != 0) continue;

        const auto state_index{ active_index_list.at(i) };
        auto refit_result{
            FitAtomWithJointOffsetFallback(
                context,
                state_index,
                previous_state.at(state_index),
                refit_snapshot,
                options)
        };
        if (!refit_result.has_value())
        {
            RecordLocalRefitStationarityFailure(
                clustered_joint_offset_result.health_by_key,
                state_index);
            post_refit_suspicious_seed_position_list.emplace_back(i);
            continue;
        }
        if (!refit_result->is_stationarity_eligible)
        {
            RecordLocalRefitStationarityFailure(
                clustered_joint_offset_result.health_by_key,
                state_index);
        }
        auto result{ std::move(refit_result->result) };
        iteration_state.at(state_index) = std::move(result);
    }
    const auto post_refit_suspicious_position_list{
        detail::ExpandPostRefitRollbackClusters(
            active_index_list,
            cluster_key_list,
            post_refit_suspicious_seed_position_list,
            suspicious_offset_mask)
    };
    RollBackSuspiciousOffsetClusters(
        context,
        active_index_list,
        previous_state,
        post_refit_suspicious_position_list,
        current_snapshot,
        iteration_state);

    LocalFittingIterationResult iteration_result;
    iteration_result.state = std::move(iteration_state);
    iteration_result.health_by_key = std::move(clustered_joint_offset_result.health_by_key);
    for (std::size_t i = 0; i < suspicious_offset_mask.size(); i++)
    {
        if (suspicious_offset_mask.at(i) == 0) continue;
        iteration_result.suspicious_offset_state_index_list.emplace_back(active_index_list.at(i));
    }
    return iteration_result;
}

void ApplyLocalFittingState(
    ModelObject & model_object,
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & iteration_state)
{
    if (context.AtomSize() != iteration_state.size())
    {
        throw std::invalid_argument("Local fitting context and state sizes are inconsistent.");
    }

    auto analysis{ model_object.EditAnalysis() };
    for (std::size_t i = 0; i < context.AtomSize(); i++)
    {
        auto local_editor{
            analysis.EnsureAtomLocalPotential(*context.atom_context_list.at(i).atom)
        };
        local_editor.SetGaussianResult(iteration_state.at(i));
    }
}

LocalFittingOffsetStats SummarizeLocalFittingOffsetValues(
    const std::vector<double> & offset_list)
{
    LocalFittingOffsetStats stats;
    stats.atom_count = offset_list.size();
    std::vector<double> absolute_offset_list;
    absolute_offset_list.reserve(offset_list.size());
    for (const auto offset : offset_list)
    {
        if (!std::isfinite(offset)) continue;
        stats.finite_count++;
        absolute_offset_list.emplace_back(std::abs(offset));
    }
    if (absolute_offset_list.empty()) return stats;

    stats.median_absolute_offset =
        array_helper::ComputeMedian(absolute_offset_list);
    stats.percentile_absolute_offset =
        array_helper::ComputePercentile(
            absolute_offset_list,
            kLocalFittingChangePercentile);
    stats.maximum_absolute_offset =
        *std::max_element(
            absolute_offset_list.begin(),
            absolute_offset_list.end());
    return stats;
}

LocalFittingOffsetStats SummarizeLocalFittingOffsets(
    const LocalFittingState & state)
{
    std::vector<double> offset_list;
    offset_list.reserve(state.size());
    for (const auto & result : state)
    {
        offset_list.emplace_back(result.mdpde.GetModel().GetOffset());
    }
    return SummarizeLocalFittingOffsetValues(offset_list);
}

void AppendLocalFittingOffsetSummary(
    std::ostringstream & stream,
    const LocalFittingOffsetStats & stats)
{
    stream
        << std::scientific << std::setprecision(2)
        << "; offsets finite = " << stats.finite_count
        << " of " << stats.atom_count
        << ", |C| median/p99/max = "
        << stats.median_absolute_offset << "/"
        << stats.percentile_absolute_offset << "/"
        << stats.maximum_absolute_offset;
}

void AppendLocalFittingAuditSummary(
    std::ostringstream & stream,
    const LocalFittingAuditedState & audited_state)
{
    const auto & objective{ audited_state.objective };
    stream
        << "; audit best source = ";
    if (audited_state.accepted_iteration.has_value())
    {
        stream << "accepted iteration " << *audited_state.accepted_iteration;
    }
    else
    {
        stream << "initial";
    }
    stream
        << std::scientific << std::setprecision(2)
        << ", fixed audit objective residual/width/offset/total = "
        << objective.residual_objective << "/"
        << objective.width_prior_penalty << "/"
        << objective.offset_plausibility_penalty << "/"
        << objective.total_objective;
}

void AppendLocalFittingTerminalSummary(
    std::ostringstream & stream,
    const LocalFittingTerminalSummary & summary)
{
    if (summary.suspicious_atom_count > 0)
    {
        stream << "; terminal suspicious rollback fallback clusters/atoms = "
            << summary.suspicious_cluster_count
            << "/" << summary.suspicious_atom_count;
    }
    if (summary.joint_offset_failure_atom_count > 0)
    {
        stream << "; terminal joint-offset failure fallback clusters/atoms = "
            << summary.joint_offset_failure_cluster_count
            << "/" << summary.joint_offset_failure_atom_count;
        if (!summary.joint_offset_failure_status_count.empty())
        {
            stream << ", statuses = ";
            bool is_first_status{ true };
            for (const auto & [status, count] :
                summary.joint_offset_failure_status_count)
            {
                if (!is_first_status) stream << ",";
                stream << GetJointOffsetSolveStatusText(status) << ":" << count;
                is_first_status = false;
            }
        }
    }
}

void LogLocalFittingTerminalFallback(
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    const LocalFittingTerminalSummary & terminal_summary,
    const LocalFittingOffsetStats & offset_stats)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message
        << "Completed local fitting after " << accepted_iteration_count
        << " accepted iterations with last validated states retained";
    AppendLocalFittingTerminalSummary(warning_message, terminal_summary);
    AppendLocalFittingOffsetSummary(warning_message, offset_stats);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void FinishLocalFittingWithNoActiveAtoms(
    ModelObject & model_object,
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & state,
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    const LocalFittingTerminalSummary & terminal_summary)
{
    ApplyLocalFittingState(model_object, context, state);
    const auto offset_stats{ SummarizeLocalFittingOffsets(state) };
    if (terminal_summary.AtomCount() > 0)
    {
        LogLocalFittingTerminalFallback(
            options,
            accepted_iteration_count,
            terminal_summary,
            offset_stats);
        return;
    }
    if (!options.quiet_mode)
    {
        Logger::FinishProgressLine();
        Logger::Log(
            LogLevel::Info,
            "Skip 2nd-stage local atom fitting because no atoms are selected.");
    }
}

void AppendLocalFittingObjectiveBreakdown(
    std::ostringstream & stream,
    const std::optional<detail::LocalFittingObjectiveBreakdown> & breakdown)
{
    if (!breakdown.has_value())
    {
        stream << "unavailable";
        return;
    }
    stream
        << breakdown->residual_objective << "/"
        << breakdown->width_prior_penalty << "/"
        << breakdown->offset_plausibility_penalty << "/"
        << breakdown->total_objective;
}

void LogRejectedLocalFittingClusterDiagnostics(
    const FitOptions & options,
    const std::vector<LocalFittingRejectedClusterDiagnostic> & diagnostic_list)
{
    if (options.quiet_mode || Logger::GetLogLevel() < LogLevel::Debug || diagnostic_list.empty())
    {
        return;
    }

    Logger::FinishProgressLine();
    for (const auto & cluster_diagnostic : diagnostic_list)
    {
        if (cluster_diagnostic.key.empty())
        {
            throw std::logic_error("Rejected local fitting cluster diagnostic key is empty.");
        }
        std::ostringstream header;
        header
            << "Rejected local fitting cluster objective diagnostics: atoms = "
            << cluster_diagnostic.key.size()
            << ", key first/last = "
            << cluster_diagnostic.key.front() << "/" << cluster_diagnostic.key.back()
            << ", breakdown order = residual/width/offset/total";
        Logger::Log(LogLevel::Debug, header.str());

        const auto & diagnostic{ cluster_diagnostic.attempt };
        std::ostringstream message;
        message
            << std::scientific << std::setprecision(2)
            << "  fixed-point effective damping = "
            << diagnostic.effective_damping
            << ", trust radius/step norm = "
            << diagnostic.trust_region_radius << "/"
            << diagnostic.trust_region_step_norm;

        if (diagnostic.status == LocalFittingObjectiveAttemptDiagnosticStatus::InvalidModel)
        {
            message << ", status = invalid-model";
            Logger::Log(LogLevel::Debug, message.str());
            continue;
        }

        message << ", objective scale = ";
        if (diagnostic.objective_scale.has_value())
        {
            message << *diagnostic.objective_scale;
        }
        else
        {
            message << "unavailable";
        }
        message << ", candidate = ";
        AppendLocalFittingObjectiveBreakdown(message, diagnostic.candidate_objective);
        message << ", previous = ";
        AppendLocalFittingObjectiveBreakdown(message, diagnostic.previous_objective);
        message << ", best = ";
        AppendLocalFittingObjectiveBreakdown(message, diagnostic.best_objective);
        message << ", rejected-by = ";
        if (diagnostic.status == LocalFittingObjectiveAttemptDiagnosticStatus::ObjectiveUnavailable)
        {
            message << "objective-unavailable";
        }
        else if (diagnostic.rejected_by_previous && diagnostic.rejected_by_best)
        {
            message << "previous+best";
        }
        else if (diagnostic.rejected_by_previous)
        {
            message << "previous";
        }
        else if (diagnostic.rejected_by_best)
        {
            message << "best";
        }
        else
        {
            message << "none";
        }
        Logger::Log(LogLevel::Debug, message.str());
    }
}

constexpr std::array<std::string_view, 9> kLocalFittingProgressHeaderList{
    "Try/Acc",
    "Atom A/T",
    "Cmp/Max/R",
    "Cand C;A A/R",
    "Pol",
    "TR G/S/M",
    "Guard S/U/C",
    "dMax A/R",
    "Audit B/P"
};

std::string FormatLocalFittingProgressTable(
    const LocalFittingProgressTableLayout & layout,
    const std::array<std::string, 9> & cell_list)
{
    std::ostringstream stream;
    for (std::size_t i = 0; i < cell_list.size(); i++)
    {
        if (i > 0) stream << " | ";
        stream
            << std::left
            << std::setw(static_cast<int>(layout.column_width_list.at(i)))
            << cell_list.at(i);
    }
    return stream.str();
}

LocalFittingProgressTableLayout BuildLocalFittingProgressTableLayout(
    std::size_t atom_size)
{
    const auto maximum_iteration_text{
        std::to_string(kLocalFittingMaximumIterations)
    };
    const auto maximum_atom_text{ std::to_string(atom_size) };
    const std::string scientific_placeholder{ "1.00e+308" };
    const std::array<std::string, 9> maximum_cell_list{
        maximum_iteration_text + "/" + maximum_iteration_text,
        maximum_atom_text + "/" + maximum_atom_text,
        maximum_atom_text + "/" + maximum_atom_text + "/1.00",
        maximum_atom_text + "/" + maximum_atom_text + ";" +
            maximum_atom_text + "/" + maximum_atom_text,
        maximum_atom_text,
        maximum_atom_text + "/" + maximum_atom_text + "/" + maximum_atom_text,
        maximum_atom_text + "/" + maximum_atom_text + "/1",
        scientific_placeholder + "/" + scientific_placeholder,
        maximum_iteration_text + "/" + maximum_iteration_text
    };

    LocalFittingProgressTableLayout layout;
    for (std::size_t i = 0; i < layout.column_width_list.size(); i++)
    {
        layout.column_width_list.at(i) = std::max(
            kLocalFittingProgressHeaderList.at(i).size(),
            maximum_cell_list.at(i).size());
    }
    return layout;
}

std::string FormatLocalFittingProgressMaximum(
    const std::optional<double> & value)
{
    if (!value.has_value()) return "-";
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(2) << *value;
    return stream.str();
}

std::optional<double> SummarizeLocalFittingProgressMaximum(
    const std::vector<double> & maximum_list)
{
    if (maximum_list.empty()) return std::nullopt;
    return *std::max_element(maximum_list.begin(), maximum_list.end());
}

std::size_t CountLocalFittingClusterAtoms(
    const std::vector<LocalFittingClusterKey> & key_list)
{
    std::size_t atom_count{ 0 };
    for (const auto & key : key_list) atom_count += key.size();
    return atom_count;
}

std::string FormatLocalFittingBestAuditIteration(
    const LocalFittingBestAuditState & best_audit_state)
{
    if (!best_audit_state.best.has_value()) return "-";
    if (!best_audit_state.best->accepted_iteration.has_value()) return "I";
    return std::to_string(*best_audit_state.best->accepted_iteration);
}

void LogLocalFittingProgressHeader(
    const FitOptions & options,
    const LocalFittingProgressTableLayout & layout)
{
    if (options.quiet_mode) return;
    std::array<std::string, 9> header_list;
    for (std::size_t i = 0; i < header_list.size(); i++)
    {
        header_list.at(i) = kLocalFittingProgressHeaderList.at(i);
    }
    Logger::Log(
        LogLevel::Info,
        FormatLocalFittingProgressTable(layout, header_list));
}

void LogLocalFittingIterationProgress(
    const FitOptions & options,
    const LocalFittingProgressTableLayout & layout,
    const LocalFittingIterationProgressSnapshot & snapshot)
{
    if (options.quiet_mode) return;

    std::ostringstream component_cell;
    component_cell
        << snapshot.component_count << "/"
        << snapshot.maximum_component_atom_count << "/"
        << std::fixed << std::setprecision(2)
        << snapshot.maximum_component_active_ratio;
    const std::array<std::string, 9> cell_list{
        std::to_string(snapshot.attempt_number) + "/" +
            std::to_string(snapshot.accepted_iteration_count),
        std::to_string(snapshot.active_atom_count) + "/" +
            std::to_string(snapshot.terminal_atom_count),
        component_cell.str(),
        std::to_string(snapshot.accepted_cluster_count) + "/" +
            std::to_string(snapshot.rejected_cluster_count) + ";" +
            std::to_string(snapshot.accepted_atom_count) + "/" +
            std::to_string(snapshot.rejected_atom_count),
        std::to_string(snapshot.accepted_polish_cluster_count),
        std::to_string(snapshot.trust_region_grow_count) + "/" +
            std::to_string(snapshot.trust_region_shrink_count) + "/" +
            std::to_string(snapshot.trust_region_minimum_count),
        std::to_string(snapshot.suspicious_atom_count) + "/" +
            std::to_string(snapshot.stationarity_ineligible_cluster_count) + "/" +
            (snapshot.combined_objective_rejected ? "1" : "0"),
        FormatLocalFittingProgressMaximum(
            snapshot.accepted_maximum_transformed_change) + "/" +
            FormatLocalFittingProgressMaximum(
                snapshot.raw_maximum_transformed_change),
        snapshot.best_audit_iteration + "/" +
            std::to_string(snapshot.audit_patience_count)
    };
    Logger::ProgressLine(
        FormatLocalFittingProgressTable(layout, cell_list));
}

void LogLocalFittingConverged(
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    const algorithm::ParameterChangeStats & transformed_change_stats,
    const LocalFittingOffsetStats & offset_stats)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream message;
    message
        << "Converged after " << accepted_iteration_count
        << " iterations with percentile log-peak-height change = "
        << transformed_change_stats.percentile_list.at(detail::kLogPeakHeightChangeIndex)
        << ", percentile log-width change = "
        << transformed_change_stats.percentile_list.at(detail::kLogWidthChangeIndex)
        << ", and percentile offset-to-peak-ratio change = "
        << transformed_change_stats.percentile_list.at(detail::kOffsetToPeakRatioChangeIndex);
    AppendLocalFittingOffsetSummary(message, offset_stats);
    message << ".";
    Logger::Log(LogLevel::Info, message.str());
}

void LogLocalFittingMaximumIterations(
    const FitOptions & options,
    const LocalFittingAuditedState * applied_audit_state,
    const LocalFittingTerminalSummary & terminal_summary,
    const LocalFittingOffsetStats & applied_offset_stats)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message << "Reached maximum iteration size";
    AppendLocalFittingTerminalSummary(warning_message, terminal_summary);
    if (applied_audit_state != nullptr)
    {
        warning_message << "; applying best validated audit state";
        AppendLocalFittingAuditSummary(warning_message, *applied_audit_state);
    }
    else
    {
        warning_message
            << "; applying current accepted candidate because no finite fixed audit state is available";
    }
    AppendLocalFittingOffsetSummary(warning_message, applied_offset_stats);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void LogSecondStageLocalFittingSummary(
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    std::string_view stop_reason,
    const LocalFittingBestAuditState & best_audit_state)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream message;
    message
        << "Second-stage local fitting summary: accepted_iterations="
        << accepted_iteration_count << ", best_iteration=";
    if (!best_audit_state.best.has_value())
    {
        message << "unavailable";
    }
    else if (best_audit_state.best->accepted_iteration.has_value())
    {
        message << *best_audit_state.best->accepted_iteration;
    }
    else
    {
        message << "initial";
    }
    message << ", stop_reason=" << stop_reason << ", best_audit_objective=";
    if (best_audit_state.best.has_value())
    {
        message << std::scientific << std::setprecision(8)
            << best_audit_state.best->objective.total_objective;
    }
    else
    {
        message << "unavailable";
    }
    message << ".";
    Logger::Log(LogLevel::Info, message.str());
}

} // namespace

void RunSecondStageLocalFitting(
    ModelObject & model_object,
    const FitOptions & options)
{
    const auto context{ BuildSecondStageLocalFittingContext(model_object) };
    const auto atom_size{ context.AtomSize() };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run 2nd-stage local atom fitting with iterations...");
    }

    auto initial_state_build_result{
        BuildInitialLocalFittingState(context, model_object.GetAnalysisView())
    };
    if (!initial_state_build_result.state.has_value())
    {
        if (!options.quiet_mode)
        {
            Logger::Log(
                LogLevel::Warning,
                "Skip 2nd-stage local atom fitting because no valid Gaussian seed "
                "is available for every selected atom.");
            Logger::Log(
                LogLevel::Info,
                "Second-stage local fitting summary: accepted_iterations=0, "
                "best_iteration=unavailable, stop_reason=no-valid-seed, "
                "best_audit_objective=unavailable.");
        }
        return;
    }
    LogSecondStageSeedRepairs(initial_state_build_result.repair_record_list, options);
    auto previous_state{ std::move(*initial_state_build_result.state) };
    const auto coupling_topology{
        BuildLocalFittingCouplingTopology(context, previous_state)
    };
    LogLocalFittingCouplingTopology(coupling_topology, options);
    auto best_audit_state{
        BuildInitialLocalFittingBestAuditState(context, previous_state)
    };

    std::vector<std::size_t> suspicious_offset_state_index_list;
    PersistentTerminalFailureStateMap persistent_terminal_failure_state_by_key;
    std::vector<char> terminal_fallback_atom_mask(atom_size, 0);
    LocalFittingTerminalSummary terminal_summary;
    LocalFittingClusterObjectiveStateMap cluster_quality_state;
    detail::LocalFittingTrustRegionStateSet trust_region_state{
        detail::LocalFittingTrustRegionOptions{
            kLocalFittingTrustRegionInitialRadius,
            kLocalFittingTrustRegionMinimumRadius,
            kLocalFittingTrustRegionMaximumRadius,
            kLocalFittingTrustRegionShrinkFactor,
            kLocalFittingTrustRegionGrowthFactor
        }
    };
    const auto progress_table_layout{
        BuildLocalFittingProgressTableLayout(atom_size)
    };
    LogLocalFittingProgressHeader(options, progress_table_layout);

    std::size_t accepted_iteration_count{ 0 };
    std::size_t audit_patience_count{ 0 };
    for (std::size_t iter = 0; iter < kLocalFittingMaximumIterations; iter++)
    {
        const auto active_index_list{
            BuildEligibleLocalFittingActiveIndexList(
                atom_size,
                terminal_fallback_atom_mask)
        };
        if (active_index_list.empty())
        {
            FinishLocalFittingWithNoActiveAtoms(
                model_object,
                context,
                previous_state,
                options,
                accepted_iteration_count,
                terminal_summary);
            LogSecondStageLocalFittingSummary(
                options,
                accepted_iteration_count,
                "terminal-isolation",
                best_audit_state);
            return;
        }

        const auto cluster_partition{
            detail::BuildLocalFittingCouplingPartition(
                coupling_topology,
                active_index_list)
        };
        std::vector<LocalFittingClusterKey> cluster_key_list;
        cluster_key_list.reserve(cluster_partition.sample_id_list_by_key.size());
        std::size_t maximum_cluster_atom_count{ 0 };
        for (const auto & [key, sample_id_list] :
            cluster_partition.sample_id_list_by_key)
        {
            static_cast<void>(sample_id_list);
            cluster_key_list.emplace_back(key);
            maximum_cluster_atom_count = std::max(maximum_cluster_atom_count, key.size());
        }

        ReconcileLocalFittingClusterObjectiveState(
            context,
            previous_state,
            cluster_partition,
            cluster_key_list,
            cluster_quality_state);
        trust_region_state.Reconcile(cluster_key_list);

        std::vector<double> joint_offset_ridge_multiplier_list(atom_size, 1.0);
        for (const auto state_index : suspicious_offset_state_index_list)
        {
            if (state_index >= atom_size)
            {
                throw std::invalid_argument(
                    "Local fitting suspicious offset atom index is out of range.");
            }
            joint_offset_ridge_multiplier_list.at(state_index) = kSuspiciousJointOffsetRidgeMultiplier;
        }

        auto iteration_result{
            RunLocalFittingIteration(
                context,
                active_index_list,
                cluster_key_list,
                previous_state,
                options,
                kJointOffsetRidgeRatio,
                joint_offset_ridge_multiplier_list)
        };
        const auto current_health_by_key{ std::move(iteration_result.health_by_key) };
        const auto stationarity_ineligible_key_list{
            CollectUnhealthyLocalFittingClusterKeys(current_health_by_key)
        };
        suspicious_offset_state_index_list = std::move(iteration_result.suspicious_offset_state_index_list);
        const auto iteration_suspicious_atom_count{ suspicious_offset_state_index_list.size() };
        const auto has_suspicious_offset_fallback{ !suspicious_offset_state_index_list.empty() };

        std::vector<LocalFittingClusterKey> polish_eligible_key_list;
        for (const auto & key : cluster_key_list)
        {
            if (ContainsLocalFittingClusterKey(stationarity_ineligible_key_list, key))
            {
                continue;
            }
            const auto contains_suspicious_atom{
                std::any_of(
                    key.begin(),
                    key.end(),
                    [&](std::size_t atom_index)
                    {
                        return std::binary_search(
                            suspicious_offset_state_index_list.begin(),
                            suspicious_offset_state_index_list.end(),
                            atom_index);
                    })
            };
            if (!contains_suspicious_atom)
            {
                polish_eligible_key_list.emplace_back(key);
            }
        }

        const auto raw_state{ std::move(iteration_result.state) };
        const auto raw_fixed_point_change_list{
            CalculateLocalFittingTransformedChanges(raw_state, previous_state)
        };
        const auto raw_fixed_point_residual_stats{
            algorithm::SummarizeParameterChangeStats(
                raw_fixed_point_change_list,
                active_index_list,
                kLocalFittingChangePercentile)
        };
        const auto raw_fixed_point_residual_maximum_list{
            detail::SummarizeLocalFittingMaximumTransformedChanges(
                raw_fixed_point_change_list,
                active_index_list)
        };
        const auto previous_transformed_estimation_list{
            BuildLocalFittingTransformedEstimationList(previous_state)
        };
        const auto raw_transformed_estimation_list{
            BuildLocalFittingTransformedEstimationList(raw_state)
        };

        auto working_cluster_quality_state{ cluster_quality_state };
        auto selection{
            SelectLocalFittingClusterCandidates(
                context,
                cluster_partition,
                cluster_key_list,
                polish_eligible_key_list,
                previous_state,
                raw_state,
                previous_transformed_estimation_list,
                raw_transformed_estimation_list,
                joint_offset_ridge_multiplier_list,
                working_cluster_quality_state,
                trust_region_state)
        };

        const auto needs_combined_objective_guard{
            cluster_partition.boundary_sample_count > 0 && !selection.accepted_key_list.empty()
        };
        const auto combined_objective_accepted{
            !needs_combined_objective_guard ||
            IsLocalFittingCombinedObjectiveAcceptable(
                context,
                selection.assembled_state,
                previous_state,
                best_audit_state)
        };
        if (!combined_objective_accepted)
        {
            RejectLocalFittingCombinedCandidate(previous_state, cluster_key_list, selection);
        }
        else
        {
            cluster_quality_state = std::move(working_cluster_quality_state);
        }

        auto assembled_state{ std::move(selection.assembled_state) };
        const auto terminal_failure_by_key{
            UpdatePersistentTerminalFailureState(
                cluster_key_list,
                selection.accepted_key_list,
                suspicious_offset_state_index_list,
                current_health_by_key,
                assembled_state,
                previous_state,
                persistent_terminal_failure_state_by_key)
        };

        std::vector<LocalFittingClusterKey> terminal_key_list;
        for (const auto & [key, reason] : terminal_failure_by_key)
        {
            terminal_key_list.emplace_back(key);
            if (std::holds_alternative<PersistentSuspiciousRollbackReason>(reason))
            {
                terminal_summary.suspicious_cluster_count++;
                terminal_summary.suspicious_atom_count += key.size();
                continue;
            }
            const auto status{ std::get<JointOffsetSolveStatus>(reason) };
            terminal_summary.joint_offset_failure_cluster_count++;
            terminal_summary.joint_offset_failure_atom_count += key.size();
            terminal_summary.joint_offset_failure_status_count[status]++;
        }
        ApplyTerminalFallbackClusters(
            terminal_key_list,
            previous_state,
            terminal_fallback_atom_mask,
            assembled_state);
        ReconcileLocalFittingBestAuditTerminalFallback(
            context,
            terminal_key_list,
            assembled_state,
            accepted_iteration_count + 1,
            best_audit_state);
        if (!terminal_key_list.empty())
        {
            suspicious_offset_state_index_list.erase(
                std::remove_if(
                    suspicious_offset_state_index_list.begin(),
                    suspicious_offset_state_index_list.end(),
                    [&](std::size_t atom_index)
                    {
                        return terminal_fallback_atom_mask.at(atom_index) != 0;
                    }),
                suspicious_offset_state_index_list.end());
        }

        trust_region_state.Grow(selection.grow_trust_region_key_list);
        const auto trust_region_radius_update{
            trust_region_state.Shrink(selection.rejected_key_list)
        };
        const auto terminal_atom_count{ terminal_summary.AtomCount() };
        if (terminal_atom_count > atom_size)
        {
            throw std::logic_error(
                "Local fitting terminal atom count exceeds selected atom count.");
        }
        LocalFittingIterationProgressSnapshot progress_snapshot{
            iter + 1,
            accepted_iteration_count,
            atom_size - terminal_atom_count,
            terminal_atom_count,
            cluster_key_list.size(),
            maximum_cluster_atom_count,
            static_cast<double>(maximum_cluster_atom_count) /
                static_cast<double>(active_index_list.size()),
            selection.accepted_key_list.size(),
            selection.rejected_key_list.size(),
            CountLocalFittingClusterAtoms(selection.accepted_key_list),
            CountLocalFittingClusterAtoms(selection.rejected_key_list),
            selection.accepted_polish_cluster_count,
            selection.grow_trust_region_key_list.size(),
            trust_region_radius_update.changed_key_list.size(),
            trust_region_radius_update.saturated_key_list.size(),
            iteration_suspicious_atom_count,
            stationarity_ineligible_key_list.size(),
            !combined_objective_accepted,
            std::nullopt,
            SummarizeLocalFittingProgressMaximum(
                raw_fixed_point_residual_maximum_list),
            FormatLocalFittingBestAuditIteration(best_audit_state),
            audit_patience_count
        };
        if (selection.accepted_key_list.empty())
        {
            LogRejectedLocalFittingClusterDiagnostics(
                options,
                selection.rejected_cluster_diagnostic_list);
            LogLocalFittingIterationProgress(
                options,
                progress_table_layout,
                progress_snapshot);
            if (!trust_region_radius_update.changed_key_list.empty() && iter + 1 < kLocalFittingMaximumIterations)
            {
                continue;
            }

            const auto fallback{
                SelectLocalFittingFallbackState(previous_state, best_audit_state.best)
            };
            ApplyLocalFittingState(model_object, context, fallback);
            LogSecondStageLocalFittingSummary(
                options,
                accepted_iteration_count,
                iter + 1 >= kLocalFittingMaximumIterations ?
                    "maximum-iterations" :
                    "all-rejected-minimum-radius",
                best_audit_state);
            return;
        }

        const auto change_list{
            CalculateLocalFittingTransformedChanges(assembled_state, previous_state)
        };
        const auto transformed_change_stats{
            algorithm::SummarizeParameterChangeStats(
                change_list,
                active_index_list,
                kLocalFittingChangePercentile)
        };
        const auto transformed_change_maximum_list{
            detail::SummarizeLocalFittingMaximumTransformedChanges(change_list, active_index_list)
        };
        const auto accepted_offset_stats{
            SummarizeLocalFittingOffsets(assembled_state)
        };

        accepted_iteration_count++;
        const auto improved_best_audit{
            TryUpdateLocalFittingBestAuditState(
                context,
                assembled_state,
                accepted_iteration_count,
                best_audit_state)
        };
        audit_patience_count = detail::AdvanceLocalFittingAuditPatience(
            audit_patience_count,
            improved_best_audit,
            !selection.rejected_key_list.empty() &&
                !trust_region_radius_update.changed_key_list.empty());
        progress_snapshot.accepted_iteration_count = accepted_iteration_count;
        progress_snapshot.accepted_maximum_transformed_change =
            SummarizeLocalFittingProgressMaximum(
                transformed_change_maximum_list);
        progress_snapshot.best_audit_iteration =
            FormatLocalFittingBestAuditIteration(best_audit_state);
        progress_snapshot.audit_patience_count = audit_patience_count;
        LogRejectedLocalFittingClusterDiagnostics(
            options,
            selection.rejected_cluster_diagnostic_list);
        LogLocalFittingIterationProgress(
            options,
            progress_table_layout,
            progress_snapshot);

        if (audit_patience_count >= kLocalFittingAuditPatience)
        {
            const auto fallback{
                SelectLocalFittingFallbackState(
                    assembled_state,
                    best_audit_state.best)
            };
            ApplyLocalFittingState(model_object, context, fallback);
            LogSecondStageLocalFittingSummary(
                options,
                accepted_iteration_count,
                "audit-patience",
                best_audit_state);
            return;
        }

        const auto converged{
            stationarity_ineligible_key_list.empty() &&
            !has_suspicious_offset_fallback &&
            selection.rejected_key_list.empty() &&
            AreLocalFittingObjectiveReferencesLocked(
                cluster_quality_state,
                cluster_key_list) &&
            IsLocalFittingTransformedChangeConverged(
                transformed_change_stats,
                transformed_change_maximum_list) &&
            IsLocalFittingTransformedChangeConverged(
                raw_fixed_point_residual_stats,
                raw_fixed_point_residual_maximum_list)
        };
        if (converged)
        {
            ApplyLocalFittingState(model_object, context, assembled_state);
            if (terminal_summary.AtomCount() > 0)
            {
                LogLocalFittingTerminalFallback(
                    options,
                    accepted_iteration_count,
                    terminal_summary,
                    accepted_offset_stats);
            }
            else
            {
                LogLocalFittingConverged(
                    options,
                    accepted_iteration_count,
                    transformed_change_stats,
                    accepted_offset_stats);
            }
            LogSecondStageLocalFittingSummary(
                options,
                accepted_iteration_count,
                "converged",
                best_audit_state);
            return;
        }

        if (iter + 1 == kLocalFittingMaximumIterations)
        {
            const auto fallback{
                SelectLocalFittingFallbackState(
                    assembled_state,
                    best_audit_state.best)
            };
            ApplyLocalFittingState(model_object, context, fallback);
            LogLocalFittingMaximumIterations(
                options,
                best_audit_state.best.has_value() ? &*best_audit_state.best : nullptr,
                terminal_summary,
                SummarizeLocalFittingOffsets(fallback));
            LogSecondStageLocalFittingSummary(
                options,
                accepted_iteration_count,
                "maximum-iterations",
                best_audit_state);
            return;
        }
        previous_state = std::move(assembled_state);
    }
}

} // namespace rhbm_gem::core
