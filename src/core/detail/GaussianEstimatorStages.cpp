#include <cstddef>

#include "core/detail/GaussianEstimatorStages.hpp"
#include "core/detail/LocalFittingAudit.hpp"
#include "core/detail/LocalFittingCouplingGraph.hpp"
#include "core/detail/LocalFittingGroupMedian.hpp"
#include "core/detail/LocalFittingHealth.hpp"
#include "core/detail/LocalFittingJointOffset.hpp"
#include "core/detail/LocalFittingJointOffsetConditioning.hpp"
#include "core/detail/LocalFittingJointPolish.hpp"
#include "core/detail/LocalFittingSeedRepair.hpp"
#include "core/detail/LocalFittingTrustRegion.hpp"
#include "core/detail/LocalFittingTransformedChange.hpp"
#include "data/detail/AtomClassifier.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/algorithm/Convergence.hpp>
#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>
#include <rhbm_gem/utils/algorithm/ScaleReferenceTracker.hpp>
#include <rhbm_gem/utils/algorithm/WeightedRidgeSolver.hpp>
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
using LocalFittingPolishProvenance = std::vector<char>;

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
    detail::LocalFittingJointOffsetParameterization parameterization{};
};

struct LocalFittingClusterHealth
{
    JointOffsetSolveStatus joint_offset_status{ JointOffsetSolveStatus::SystemBuildFailed };
    bool is_refit_stationarity_eligible{ true };
};

using LocalFittingClusterHealthMap = std::map<LocalFittingClusterKey, LocalFittingClusterHealth>;

struct LocalFittingIterationResult
{
    LocalFittingState state{};
    std::vector<char> suspicious_offset_atom_mask{};
    LocalFittingClusterHealthMap health_by_key{};
};

struct LocalAtomRefitResult
{
    LocalGaussianResult result{};
    bool is_stationarity_eligible{ false };
};

struct LocalFittingTransformedChangeSummary
{
    algorithm::ParameterChangeStats percentile_stats{};
    std::vector<double> maximum_list{};
};

struct LocalFittingObjectiveAtomSample
{
    GaussianModel3D model{};
    double prior_width{ 1.0 };
};

struct LocalFittingObjectiveSamples
{
    std::vector<double> residual_list{};
    std::vector<LocalFittingObjectiveAtomSample> atom_sample_list{};
    double scale_sample{ 0.0 };
};

struct LocalFittingClusterObjectiveState
{
    algorithm::ScaleReferenceTracker scale_tracker;
    std::optional<LocalFittingObjectiveSamples> previous_objective_samples{};
    std::optional<LocalFittingObjectiveSamples> best_objective_samples{};
    double best_maximum_transformed_change{ 0.0 };

    LocalFittingClusterObjectiveState(
        std::optional<LocalFittingObjectiveSamples> initial_objective_samples,
        bool has_initial_objective)
        : scale_tracker{
              kLocalFittingObjectiveScaleWarmupCount,
              initial_objective_samples.has_value() ?
                  std::optional<double>{ initial_objective_samples->scale_sample } :
                  std::nullopt
          },
          previous_objective_samples{ std::move(initial_objective_samples) }
    {
        if (has_initial_objective)
        {
            best_objective_samples = previous_objective_samples;
        }
    }
};

using LocalFittingClusterObjectiveStateMap =
    std::map<LocalFittingClusterKey, LocalFittingClusterObjectiveState>;

using LocalFittingObjectiveSampleRef = detail::LocalFittingCouplingSampleId;

struct LocalFittingAuditedState
{
    detail::LocalFittingObjectiveBreakdown objective{};
    LocalFittingState state{};
    LocalFittingPolishProvenance polish_provenance{};
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

bool UsesLocalFittingPolish(const LocalFittingPolishProvenance & provenance)
{
    return std::any_of(
        provenance.begin(),
        provenance.end(),
        [](char is_polished)
        {
            return is_polished != 0;
        });
}

bool SelectLocalFittingFallbackUsesPolish(
    const LocalFittingPolishProvenance & fallback_provenance,
    const std::optional<LocalFittingAuditedState> & audited_state)
{
    return UsesLocalFittingPolish(
        audited_state.has_value() ? audited_state->polish_provenance : fallback_provenance);
}

struct LocalFittingOffsetStats
{
    std::size_t atom_count{ 0 };
    std::size_t finite_count{ 0 };
    double median_absolute_offset{ 0.0 };
    double percentile_absolute_offset{ 0.0 };
    double maximum_absolute_offset{ 0.0 };
};

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

struct LocalFittingObjectiveAttemptDiagnostic
{
    double effective_damping{ 1.0 };
    bool is_invalid_model{ false };
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

struct LocalFittingPolishProgress
{
    std::size_t eligible_count{ 0 };
    std::size_t accepted_count{ 0 };
    std::size_t rejected_count{ 0 };
    std::size_t skipped_count{ 0 };
};

struct LocalFittingCandidateSelection
{
    LocalFittingState assembled_state{};
    LocalFittingPolishProvenance assembled_polish_provenance{};
    std::vector<LocalFittingClusterKey> accepted_key_list{};
    std::vector<LocalFittingClusterKey> rejected_key_list{};
    std::vector<LocalFittingRejectedClusterDiagnostic> rejected_cluster_diagnostic_list{};
    std::vector<LocalFittingClusterKey> grow_trust_region_key_list{};
    LocalFittingPolishProgress polish_progress{};
};

struct LocalFittingIterationProgress
{
    std::size_t attempt_number{ 0 };
    std::size_t accepted_iteration_count{ 0 };
    std::size_t active_atom_count{ 0 };
    std::size_t terminal_atom_count{ 0 };
    std::size_t accepted_cluster_count{ 0 };
    std::size_t rejected_cluster_count{ 0 };
    LocalFittingPolishProgress polish_progress{};
    std::size_t suspicious_atom_count{ 0 };
    std::optional<double> accepted_maximum_transformed_change{};
    std::optional<double> raw_maximum_transformed_change{};
};

using LocalFittingProgressColumnWidths = std::array<std::size_t, 6>;

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

using SecondStageLocalFittingContext = std::vector<SecondStageAtomContext>;

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
    LocalFittingState state{};
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
    if (!CanBuildFiniteZeroOffsetSamples(sample_entries, offset_model)) return true;
    if (!CanBuildFiniteZeroOffsetSamples(sample_entries, previous_model)) return false;
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
    context.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        context.emplace_back(SecondStageAtomContext{ atom });
    }
    std::unordered_map<const AtomObject *, std::size_t> atom_index_map;
    atom_index_map.reserve(context.size());
    for (std::size_t i = 0; i < context.size(); i++)
    {
        atom_index_map.emplace(context.at(i).atom, i);
    }
    const auto analysis_view{ model_object.GetAnalysisView() };

    std::unordered_map<GroupKey, std::vector<double>> width_samples_by_group;
    width_samples_by_group.reserve(context.size());
    for (const auto & atom_context : context)
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

    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        auto & atom_context{ context.at(atom_index) };
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

    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        auto & atom_context{ context.at(atom_index) };
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
                            context.at(neighbor_index).atom->GetPositionRef()))
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
    detail::LocalFittingCouplingGraphBuilder builder{ context.size() };
    const auto invalid_jacobian{
        Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN())
    };
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        const auto & atom_context{ context.at(atom_index) };
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

std::optional<GaussianModel3DWithUncertainty> BuildValidGaussianParameterMedian(
    const std::vector<GaussianModel3D> & model_list)
{
    const auto median_model{
        detail::BuildLocalFittingGaussianParameterMedian(model_list)
    };
    if (!median_model.has_value()) return std::nullopt;
    return GaussianModel3DWithUncertainty{
        *median_model,
        GaussianModel3DUncertainty{}
    };
}

std::optional<SecondStageInitialStateBuildResult> BuildInitialLocalFittingState(
    const SecondStageLocalFittingContext & context,
    const ModelAnalysisView & analysis_view)
{
    SecondStageInitialStateBuildResult build_result;
    auto & state{ build_result.state };
    state.resize(context.size());
    std::vector<std::optional<GaussianModel3DWithUncertainty>> group_prior_list(
        context.size());
    std::unordered_map<GroupKey, std::vector<GaussianModel3D>> models_by_group;
    std::vector<GaussianModel3D> global_models;

    for (std::size_t i = 0; i < context.size(); i++)
    {
        const auto * atom{ context.at(i).atom };
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

        models_by_group[group_key].emplace_back(preferred_model->GetModel());
        global_models.emplace_back(preferred_model->GetModel());
    }

    std::unordered_map<GroupKey, GaussianModel3DWithUncertainty> median_by_group;
    median_by_group.reserve(models_by_group.size());
    for (const auto & [group_key, models] : models_by_group)
    {
        const auto median_model{ BuildValidGaussianParameterMedian(models) };
        if (median_model.has_value())
        {
            median_by_group.emplace(group_key, *median_model);
        }
    }
    const auto global_median{ BuildValidGaussianParameterMedian(global_models) };

    for (std::size_t i = 0; i < context.size(); i++)
    {
        auto & result{ state.at(i) };
        const auto original_model{ result.mdpde.GetModel() };
        if (!detail::IsValidSecondStageGaussianModel(original_model))
        {
            const auto group_key{
                data_internal::GetGroupKey(context.at(i).atom)
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
                return std::nullopt;
            }

            const auto repaired_seed{
                detail::BuildRepairedSecondStageSeed(original_model, *selection)
            };
            const auto repaired_model{ repaired_seed.GetModel() };
            if (!detail::IsValidSecondStageGaussianModel(repaired_model))
            {
                return std::nullopt;
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

double CalculateLocalFittingRidgeDiagonal(
    double column_square_sum,
    double multiplier)
{
    if (!std::isfinite(multiplier) || multiplier <= 0.0)
    {
        throw std::invalid_argument(
            "Local fitting ridge multiplier must be positive and finite.");
    }
    const auto base_ridge{
        column_square_sum > std::numeric_limits<double>::epsilon() ?
            kJointOffsetRidgeRatio * column_square_sum : 1.0
    };
    return multiplier * base_ridge;
}

JointOffsetBuildResult BuildJointOffsetSystem(
    const SecondStageLocalFittingContext & context,
    const std::vector<std::size_t> & active_index_list,
    const FittedGaussianSnapshot & snapshot,
    const std::vector<double> & ridge_multiplier_list,
    bool log_debug_diagnostics)
{
    const auto atom_size{ context.size() };
    std::vector<GroupKey> group_key_by_atom_position;
    std::vector<GaussianModel3D> active_model_list;
    group_key_by_atom_position.reserve(active_index_list.size());
    active_model_list.reserve(active_index_list.size());
    for (const auto atom_index : active_index_list)
    {
        group_key_by_atom_position.emplace_back(
            data_internal::GetGroupKey(context.at(atom_index).atom));
        active_model_list.emplace_back(snapshot.at(atom_index));
    }
    auto parameterization{
        detail::BuildLocalFittingJointOffsetParameterization(
            group_key_by_atom_position,
            active_model_list)
    };
    if (!parameterization.has_value())
    {
        throw std::runtime_error(
            "Joint offset group parameterization is invalid.");
    }

    std::vector<int> active_position_by_atom_index(atom_size, -1);
    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto atom_index{ active_index_list.at(i) };
        active_position_by_atom_index.at(atom_index) = static_cast<int>(i);
    }

    const auto column_count{ parameterization->ParameterCount() };
    std::vector<Eigen::Triplet<double>> triplet_list;
    std::vector<double> response_list;
    Eigen::VectorXd group_column_square_sum{
        Eigen::VectorXd::Zero(column_count)
    };
    std::map<std::pair<Eigen::Index, Eigen::Index>, double>
        group_column_cross_sum_map;
    Eigen::VectorXd atom_column_square_sum{
        Eigen::VectorXd::Zero(
            static_cast<Eigen::Index>(active_index_list.size()))
    };
    std::map<std::pair<std::size_t, std::size_t>, double>
        atom_column_cross_sum_map;
    ActiveCouplingGraph active_coupling_graph(active_index_list.size());
    std::vector<std::pair<std::size_t, double>> atom_row_basis_entries;
    std::vector<std::pair<Eigen::Index, double>> group_row_basis_entries;
    for (const auto active_index : active_index_list)
    {
        const auto target_position{
            active_position_by_atom_index.at(active_index)
        };
        const auto & atom_context{ context.at(active_index) };
        const auto & target_model{ snapshot.at(active_index) };
        atom_row_basis_entries.reserve(atom_size);
        group_row_basis_entries.reserve(parameterization->GroupCount());
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
            atom_row_basis_entries.clear();
            if (std::abs(target_basis) > std::numeric_limits<double>::epsilon())
            {
                atom_row_basis_entries.emplace_back(
                    static_cast<std::size_t>(target_position),
                    target_basis);
            }

            for (const auto & neighbor_sample : atom_context.sample_neighbor_list.at(sample_index))
            {
                const auto & neighbor_model{ snapshot.at(neighbor_sample.atom_index) };
                const auto neighbor_position{
                    active_position_by_atom_index.at(neighbor_sample.atom_index)
                };
                if (neighbor_position < 0)
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
                    atom_row_basis_entries.emplace_back(
                        static_cast<std::size_t>(neighbor_position),
                        basis);
                }
            }
            if (!std::isfinite(residual))
            {
                throw std::runtime_error("Joint offset residual is not finite.");
            }
            if (atom_row_basis_entries.empty()) continue;

            for (const auto & [atom_position, basis] : atom_row_basis_entries)
            {
                atom_column_square_sum(
                    static_cast<Eigen::Index>(atom_position)) += basis * basis;
            }
            for (std::size_t i = 0; i < atom_row_basis_entries.size(); i++)
            {
                const auto [left_position, left_basis]{
                    atom_row_basis_entries.at(i)
                };
                for (std::size_t j = i + 1;
                    j < atom_row_basis_entries.size();
                    j++)
                {
                    const auto [right_position, right_basis]{
                        atom_row_basis_entries.at(j)
                    };
                    if (left_position == right_position) continue;
                    const auto position_pair{
                        std::minmax(left_position, right_position)
                    };
                    atom_column_cross_sum_map[position_pair] +=
                        left_basis * right_basis;
                }
            }

            const auto group_basis{
                parameterization->AggregateBasis(atom_row_basis_entries)
            };
            if (!group_basis.has_value())
            {
                throw std::runtime_error(
                    "Joint offset group basis is invalid.");
            }
            group_row_basis_entries.clear();
            for (Eigen::Index column_index = 0;
                column_index < group_basis->size();
                column_index++)
            {
                const auto basis{ (*group_basis)(column_index) };
                if (std::abs(basis) <= std::numeric_limits<double>::epsilon())
                {
                    continue;
                }
                group_row_basis_entries.emplace_back(column_index, basis);
            }
            if (group_row_basis_entries.empty()) continue;

            const auto row_index{
                static_cast<Eigen::Index>(response_list.size())
            };
            response_list.emplace_back(residual);
            for (const auto & [column_index, basis] : group_row_basis_entries)
            {
                triplet_list.emplace_back(row_index, column_index, basis);
                group_column_square_sum(column_index) += basis * basis;
            }
            for (std::size_t i = 0; i < group_row_basis_entries.size(); i++)
            {
                const auto [left_column, left_basis]{
                    group_row_basis_entries.at(i)
                };
                for (std::size_t j = i + 1;
                    j < group_row_basis_entries.size();
                    j++)
                {
                    const auto [right_column, right_basis]{
                        group_row_basis_entries.at(j)
                    };
                    const auto column_pair{
                        std::minmax(left_column, right_column)
                    };
                    group_column_cross_sum_map[column_pair] +=
                        left_basis * right_basis;
                }
            }
        }
    }

    for (const auto & [position_pair, cross_sum] :
        atom_column_cross_sum_map)
    {
        const auto left_position{ position_pair.first };
        const auto right_position{ position_pair.second };
        const auto left_square_sum{ atom_column_square_sum(
            static_cast<Eigen::Index>(left_position)) };
        const auto right_square_sum{ atom_column_square_sum(
            static_cast<Eigen::Index>(right_position)) };
        if (left_square_sum <= std::numeric_limits<double>::epsilon() ||
            right_square_sum <= std::numeric_limits<double>::epsilon())
        {
            continue;
        }
        const auto overlap{
            std::abs(cross_sum) /
                std::sqrt(left_square_sum * right_square_sum)
        };
        if (!std::isfinite(overlap)) continue;
        active_coupling_graph.at(left_position).emplace_back(
            ActiveCouplingEdge{ right_position, overlap });
        active_coupling_graph.at(right_position).emplace_back(
            ActiveCouplingEdge{ left_position, overlap });
    }

    Eigen::VectorXd proactive_ridge_multiplier{ Eigen::VectorXd::Ones(column_count) };
    for (const auto & [column_pair, cross_sum] :
        group_column_cross_sum_map)
    {
        const auto left_column{ column_pair.first };
        const auto right_column{ column_pair.second };
        const auto left_square_sum{ group_column_square_sum(left_column) };
        const auto right_square_sum{ group_column_square_sum(right_column) };
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
    system.previous_parameter = parameterization->seed_offset;
    system.ridge_diagonal = Eigen::VectorXd::Zero(column_count);
    for (Eigen::Index column_index = 0; column_index < column_count; column_index++)
    {
        double multiplier{ 1.0 };
        for (const auto atom_position :
            parameterization->atom_position_list_by_group.at(
                static_cast<std::size_t>(column_index)))
        {
            const auto atom_index{ active_index_list.at(atom_position) };
            multiplier = std::max(
                multiplier,
                ridge_multiplier_list.at(atom_index));
        }
        const auto square_sum{ group_column_square_sum(column_index) };
        const auto combined_multiplier{
            std::max(multiplier, proactive_ridge_multiplier(column_index))
        };
        system.ridge_diagonal(column_index) =
            CalculateLocalFittingRidgeDiagonal(square_sum, combined_multiplier);
    }
    return JointOffsetBuildResult{
        std::move(system),
        std::move(active_coupling_graph),
        std::move(*parameterization)
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
    const std::vector<double> & ridge_multiplier_list,
    bool log_debug_diagnostics)
{
    Eigen::VectorXd previous_offset{
        Eigen::VectorXd::Zero(static_cast<Eigen::Index>(active_index_list.size()))
    };
    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto atom_index{ active_index_list.at(i) };
        previous_offset(static_cast<Eigen::Index>(i)) = snapshot.at(atom_index).GetOffset();
    }
    JointOffsetBuildResult build_result;
    try
    {
        build_result = BuildJointOffsetSystem(
            context,
            active_index_list,
            snapshot,
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
    auto parameterization{ std::move(build_result.parameterization) };
    const auto make_progress_result = [&](
        JointOffsetSolveStatus status,
        const Eigen::VectorXd & group_offset)
    {
        auto atom_offset{ parameterization.ExpandOffsets(group_offset) };
        if (!atom_offset.has_value())
        {
            return JointOffsetSolveResult{
                JointOffsetSolveStatus::IrlsSolveFailed,
                previous_offset,
                std::move(active_coupling_graph)
            };
        }
        return JointOffsetSolveResult{
            status,
            std::move(*atom_offset),
            std::move(active_coupling_graph)
        };
    };
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
            weight(i) = algorithm::CalculateCauchyWeight(
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
            return make_progress_result(
                JointOffsetSolveStatus::IrlsObjectiveDeteriorated,
                offset);
        }
        const auto maximum_change{
            algorithm::CalculateMaximumNormalizedVectorChange(updated_offset, offset, kJointOffsetIrlsScaleFloor)
        };
        offset = std::move(updated_offset);
        if (maximum_change < kJointOffsetIrlsNormalizedChangeTolerance)
        {
            return make_progress_result(
                JointOffsetSolveStatus::Converged,
                offset);
        }
    }

    return make_progress_result(
        JointOffsetSolveStatus::IrlsMaximumIterationsReached,
        offset);
}

double CalculateSecondStageAdjustedResponse(
    const SecondStageLocalFittingContext & context,
    std::size_t atom_index,
    std::size_t sample_index,
    const FittedGaussianSnapshot & snapshot)
{
    const auto & atom_context{ context.at(atom_index) };
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
    const auto & atom_context{ context.at(atom_index) };
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
    LocalFittingObjectiveSamples objective_samples;
    objective_samples.residual_list.reserve(sample_ref_list.size());
    objective_samples.atom_sample_list.reserve(active_index_list.size());
    for (const auto active_index : active_index_list)
    {
        objective_samples.atom_sample_list.emplace_back(LocalFittingObjectiveAtomSample{
            state.at(active_index).mdpde.GetModel(),
            context.at(active_index).prior_width
        });
    }
    std::vector<double> response_list;
    response_list.reserve(sample_ref_list.size());
    for (const auto & sample_ref : sample_ref_list)
    {
        const auto & atom_context{ context.at(sample_ref.atom_index) };
        const auto & sample{ atom_context.sample_entries.at(sample_ref.sample_index) };
        const auto & target_model{
            state.at(sample_ref.atom_index).mdpde.GetModel()
        };
        const auto distance{ static_cast<double>(sample.point.distance) };
        const auto expected_response{ target_model.ResponseAtDistance(distance) };
        auto response{ static_cast<double>(sample.response) };
        for (const auto & neighbor_sample :
            atom_context.sample_neighbor_list.at(sample_ref.sample_index))
        {
            response -= state.at(neighbor_sample.atom_index).mdpde.GetModel()
                .ResponseAtDistance(neighbor_sample.distance);
        }
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
    return objective_samples;
}

std::optional<detail::LocalFittingObjectiveBreakdown>
CalculateLocalFittingObjectiveBreakdown(
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
        loss_sum += algorithm::CalculateCauchyLoss(
            normalized_residual,
            kRobustLossCutoffMultiplier);
    }
    const auto residual_objective{
        loss_sum / static_cast<double>(objective_samples.residual_list.size())
    };
    if (!std::isfinite(residual_objective)) return std::nullopt;

    double width_prior_penalty_sum{ 0.0 };
    double offset_plausibility_penalty_sum{ 0.0 };
    const auto residual_scale_floor{ std::max(objective_scale, kRobustScaleMin) };
    for (const auto & atom_sample : objective_samples.atom_sample_list)
    {
        const auto & model{ atom_sample.model };
        if (!detail::IsValidSecondStageGaussianModel(model)) return std::nullopt;

        const auto prior_width{ atom_sample.prior_width };
        if (!numeric_validation::IsFinitePositive(prior_width)) return std::nullopt;
        const auto normalized_width_difference{
            (std::log(model.GetWidth()) - std::log(prior_width)) /
            kLocalFittingWidthPriorLogScale
        };
        width_prior_penalty_sum +=
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
        offset_plausibility_penalty_sum += offset_excess * offset_excess;
    }

    if (!std::isfinite(width_prior_penalty_sum) ||
        !std::isfinite(offset_plausibility_penalty_sum))
    {
        return std::nullopt;
    }

    return detail::BuildLocalFittingMeanObjectiveBreakdown(
            residual_objective,
            width_prior_penalty_sum,
            offset_plausibility_penalty_sum,
            objective_samples.atom_sample_list.size(),
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
    const LocalFittingPolishProvenance & previous_polish_provenance,
    const std::vector<LocalFittingClusterKey> & cluster_key_list,
    LocalFittingCandidateSelection & selection)
{
    selection.assembled_state = previous_state;
    selection.assembled_polish_provenance = previous_polish_provenance;
    selection.accepted_key_list.clear();
    selection.rejected_key_list = cluster_key_list;
    selection.grow_trust_region_key_list.clear();
    selection.polish_progress.rejected_count += selection.polish_progress.accepted_count;
    selection.polish_progress.accepted_count = 0;
}

bool TryUpdateLocalFittingBestAuditState(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & candidate_state,
    const LocalFittingPolishProvenance & candidate_polish_provenance,
    const std::optional<std::size_t> & accepted_iteration,
    LocalFittingBestAuditState & audit_state)
{
    const auto candidate_objective{
        EvaluateLocalFittingAuditObjective(context, candidate_state, audit_state)
    };
    if (!candidate_objective.has_value()) return false;
    if (audit_state.best.has_value() &&
        !detail::IsBetterLocalFittingAuditObjective(
            candidate_objective->total_objective,
            audit_state.best->objective.total_objective,
            kLocalFittingObjectiveTieRelativeTolerance))
    {
        return false;
    }
    audit_state.best = LocalFittingAuditedState{
        *candidate_objective,
        candidate_state,
        candidate_polish_provenance,
        accepted_iteration
    };
    return true;
}

LocalFittingBestAuditState BuildInitialLocalFittingBestAuditState(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & initial_state)
{
    LocalFittingBestAuditState audit_state;
    audit_state.atom_index_list.reserve(context.size());
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        audit_state.atom_index_list.emplace_back(atom_index);
        const auto sample_size{
            context.at(atom_index).sample_entries.size()
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
        LocalFittingPolishProvenance(initial_state.size(), 0),
        std::nullopt,
        audit_state));
    return audit_state;
}

void ReconcileLocalFittingBestAuditTerminalFallback(
    const SecondStageLocalFittingContext & context,
    const std::vector<LocalFittingClusterKey> & terminal_key_list,
    const LocalFittingState & terminal_fallback_state,
    const LocalFittingPolishProvenance & terminal_fallback_polish_provenance,
    std::size_t accepted_iteration,
    LocalFittingBestAuditState & audit_state)
{
    if (terminal_key_list.empty()) return;
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
    audit_state.best->polish_provenance = terminal_fallback_polish_provenance;
    audit_state.best->accepted_iteration = accepted_iteration;
}

LocalFittingTransformedChangeSummary SummarizeLocalFittingTransformedChanges(
    const LocalFittingState & current_state,
    const LocalFittingState & previous_state,
    const std::vector<std::size_t> & index_list)
{
    std::vector<algorithm::ParameterChange> change_list;
    change_list.reserve(index_list.size());
    for (const auto i : index_list)
    {
        change_list.emplace_back(detail::CalculateLocalFittingTransformedChange(
            current_state.at(i).mdpde.GetModel(),
            previous_state.at(i).mdpde.GetModel()));
    }

    std::vector<std::size_t> local_index_list(change_list.size());
    for (std::size_t i = 0; i < local_index_list.size(); i++)
    {
        local_index_list.at(i) = i;
    }
    return LocalFittingTransformedChangeSummary{
        algorithm::SummarizeParameterChangeStats(
            change_list,
            local_index_list,
            kLocalFittingChangePercentile),
        detail::SummarizeLocalFittingMaximumTransformedChanges(
            change_list,
            local_index_list)
    };
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

TerminalPersistentFailureMap UpdatePersistentTerminalFailureState(
    const std::vector<LocalFittingClusterKey> & accepted_key_list,
    const std::vector<char> & suspicious_atom_mask,
    const LocalFittingClusterHealthMap & health_by_key,
    const LocalFittingState & assembled_state,
    const LocalFittingState & previous_state,
    PersistentTerminalFailureStateMap & state_by_key)
{
    PersistentTerminalFailureStateMap next_state_by_key;
    TerminalPersistentFailureMap terminal_failure_by_key;
    for (const auto & [key, health] : health_by_key)
    {
        if (std::find(accepted_key_list.begin(), accepted_key_list.end(), key) ==
            accepted_key_list.end())
        {
            continue;
        }

        PersistentSuspiciousRollbackReason cluster_suspicious_atom_index_list;
        for (const auto atom_index : key)
        {
            if (suspicious_atom_mask.at(atom_index) != 0)
            {
                cluster_suspicious_atom_index_list.emplace_back(atom_index);
            }
        }
        PersistentTerminalFailureReason reason;
        if (!cluster_suspicious_atom_index_list.empty())
        {
            reason = std::move(cluster_suspicious_atom_index_list);
        }
        else
        {
            const auto status{ health.joint_offset_status };
            if (!IsJointOffsetSolveHardFailure(status)) continue;
            reason = status;
        }

        const auto transformed_change_summary{
            SummarizeLocalFittingTransformedChanges(
                assembled_state,
                previous_state,
                key)
        };
        if (!IsLocalFittingTransformedPercentileConverged(
                transformed_change_summary.percentile_stats))
        {
            continue;
        }

        PersistentTerminalFailureState next_state{ std::move(reason), 1 };
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
    const LocalFittingPolishProvenance & previous_polish_provenance,
    std::vector<char> & terminal_atom_mask,
    LocalFittingState & assembled_state,
    LocalFittingPolishProvenance & assembled_polish_provenance)
{
    for (const auto & key : terminal_key_list)
    {
        for (const auto atom_index : key)
        {
            terminal_atom_mask.at(atom_index) = 1;
            assembled_state.at(atom_index) = previous_state.at(atom_index);
            assembled_polish_provenance.at(atom_index) = previous_polish_provenance.at(atom_index);
        }
    }
}

std::vector<std::size_t> BuildEligibleLocalFittingActiveIndexList(
    const std::vector<char> & terminal_atom_mask)
{
    const auto atom_size{ terminal_atom_mask.size() };
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

std::vector<Eigen::Vector3d> BuildLocalFittingTransformedEstimationList(
    const LocalFittingState & state)
{
    std::vector<Eigen::Vector3d> transformed_estimation_list;
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

std::optional<Eigen::VectorXd> BuildLocalFittingJointPolishDirection(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & base_state,
    const std::vector<GaussianModel3D> & seed_model_list,
    const LocalFittingClusterKey & key,
    const std::vector<LocalFittingObjectiveSampleRef> & sample_ref_list,
    const std::vector<double> & ridge_multiplier_list,
    const detail::LocalFittingJointPolishParameterization & parameterization)
{
    if (key.empty() || sample_ref_list.empty() ||
        seed_model_list.size() != key.size() ||
        parameterization.AtomCount() != key.size())
    {
        return std::nullopt;
    }

    const auto column_count{ parameterization.ParameterCount() };
    std::vector<int> local_position_by_atom_index(context.size(), -1);
    for (std::size_t local_position = 0; local_position < key.size(); local_position++)
    {
        const auto atom_index{ key.at(local_position) };
        local_position_by_atom_index.at(atom_index) =
            static_cast<int>(local_position);
    }

    std::vector<Eigen::Triplet<double>> triplet_list;
    std::vector<double> residual_list;
    residual_list.reserve(sample_ref_list.size());
    for (const auto & sample_ref : sample_ref_list)
    {
        const auto & atom_context{
            context.at(sample_ref.atom_index)
        };
        const auto & sample{ atom_context.sample_entries.at(sample_ref.sample_index) };
        if (!std::isfinite(static_cast<double>(sample.response))) return std::nullopt;

        const auto row_index{ static_cast<Eigen::Index>(residual_list.size()) };
        double predicted_response{ 0.0 };
        const auto append_model = [&](std::size_t atom_index, double distance) -> bool
        {
            const auto local_position_value{
                local_position_by_atom_index.at(atom_index)
            };
            const auto & model{
                local_position_value >= 0 ?
                    seed_model_list.at(static_cast<std::size_t>(
                        local_position_value)) :
                    base_state.at(atom_index).mdpde.GetModel()
            };
            const auto evaluation{
                detail::EvaluateLocalFittingSharedOffsetResponse(
                    model,
                    distance)
            };
            if (!evaluation.has_value()) return false;
            predicted_response += evaluation->response;

            if (local_position_value < 0) return true;
            const auto local_position{
                static_cast<std::size_t>(local_position_value)
            };
            for (std::size_t parameter_index = 0;
                parameter_index < detail::kLocalFittingJointPolishShapeParameterSize;
                parameter_index++)
            {
                const auto column_index{
                    parameterization.ShapeColumn(
                        local_position,
                        parameter_index)
                };
                const auto derivative{
                    evaluation->shape_jacobian(
                        static_cast<Eigen::Index>(parameter_index))
                };
                if (std::abs(derivative) <= std::numeric_limits<double>::epsilon()) continue;
                triplet_list.emplace_back(row_index, column_index, derivative);
            }
            if (std::abs(evaluation->offset_jacobian) >
                std::numeric_limits<double>::epsilon())
            {
                triplet_list.emplace_back(
                    row_index,
                    parameterization.OffsetColumn(local_position),
                    evaluation->offset_jacobian);
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
    Eigen::VectorXd column_square_sum{ Eigen::VectorXd::Zero(column_count) };
    for (Eigen::Index column_index = 0;
        column_index < system.design_matrix.outerSize();
        column_index++)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator iter(
                system.design_matrix,
                column_index);
            iter;
            ++iter)
        {
            column_square_sum(column_index) += iter.value() * iter.value();
        }
    }
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
        double parameter_multiplier{ 1.0 };
        const auto offset_column_base{
            static_cast<Eigen::Index>(
                key.size() * detail::kLocalFittingJointPolishShapeParameterSize)
        };
        if (column_index < offset_column_base)
        {
            const auto local_position{
                static_cast<std::size_t>(column_index) /
                detail::kLocalFittingJointPolishShapeParameterSize
            };
            parameter_multiplier = ridge_multiplier_list.at(
                key.at(local_position));
        }
        else
        {
            const auto group_position{
                static_cast<std::size_t>(column_index - offset_column_base)
            };
            for (const auto local_position :
                parameterization.atom_position_list_by_group.at(group_position))
            {
                parameter_multiplier = std::max(
                    parameter_multiplier,
                    ridge_multiplier_list.at(key.at(local_position)));
            }
        }
        const auto square_sum{ column_square_sum(column_index) };
        system.ridge_diagonal(column_index) =
            CalculateLocalFittingRidgeDiagonal(
                square_sum,
                std::max(parameter_multiplier, conditioning_multiplier));
    }

    const auto residual_scale{
        std::max(CalculateMedianAbsoluteDeviationScale(residual_list), kRobustScaleMin)
    };
    if (!std::isfinite(residual_scale)) return std::nullopt;
    Eigen::VectorXd weight{ Eigen::VectorXd::Ones(row_count) };
    for (Eigen::Index row_index = 0; row_index < row_count; row_index++)
    {
        weight(row_index) = algorithm::CalculateCauchyWeight(
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

    return direction;
}

std::optional<double> CalculateLocalFittingClusterModelTrustRegionStepNorm(
    const LocalFittingState & outer_previous_state,
    const LocalFittingClusterKey & key,
    const std::vector<GaussianModel3D> & candidate_model_list)
{
    if (candidate_model_list.size() != key.size()) return std::nullopt;
    double step_norm{ 0.0 };
    for (std::size_t atom_position = 0;
        atom_position < key.size();
        atom_position++)
    {
        const auto atom_index{ key.at(atom_position) };
        const auto previous{
            detail::EncodeLocalFittingTransformedCoordinates(
                outer_previous_state.at(atom_index).mdpde.GetModel())
        };
        const auto candidate{
            detail::EncodeLocalFittingTransformedCoordinates(
                candidate_model_list.at(atom_position))
        };
        if (!previous.has_value() || !candidate.has_value())
        {
            return std::nullopt;
        }
        for (std::size_t parameter_index = 0;
            parameter_index < detail::kTransformedChangeSize;
            parameter_index++)
        {
            const auto eigen_index{ static_cast<Eigen::Index>(parameter_index) };
            step_norm = std::max(
                step_norm,
                std::abs((*candidate)(eigen_index) - (*previous)(eigen_index)) /
                    kLocalFittingTrustRegionParameterScale.at(parameter_index));
        }
    }
    return std::isfinite(step_norm) ?
        std::optional<double>{ step_norm } : std::nullopt;
}

struct LocalFittingBaseProposal
{
    LocalFittingState state{};
    double effective_damping{ 0.0 };
    double step_norm{ 0.0 };
};

std::optional<LocalFittingBaseProposal>
BuildLocalFittingSharedOffsetBaseProposal(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & outer_previous_state,
    const LocalFittingState & raw_state,
    const LocalFittingClusterKey & key,
    double trust_region_radius)
{
    constexpr double trust_region_tolerance{ 1.0e-12 };
    if (key.empty()) return std::nullopt;

    std::vector<GroupKey> group_key_by_atom_position;
    std::vector<GaussianModel3D> previous_model_list;
    std::vector<GaussianModel3D> raw_model_list;
    group_key_by_atom_position.reserve(key.size());
    previous_model_list.reserve(key.size());
    raw_model_list.reserve(key.size());
    for (const auto atom_index : key)
    {
        group_key_by_atom_position.emplace_back(
            data_internal::GetGroupKey(context.at(atom_index).atom));
        previous_model_list.emplace_back(
            outer_previous_state.at(atom_index).mdpde.GetModel());
        raw_model_list.emplace_back(
            raw_state.at(atom_index).mdpde.GetModel());
    }
    const auto shared_offset_model_list{
        detail::BuildLocalFittingGroupMedianModelList(
            group_key_by_atom_position,
            raw_model_list)
    };

    const auto seed_model_list{
        detail::BuildLocalFittingSharedOffsetDampedModelList(
            previous_model_list,
            raw_model_list,
            shared_offset_model_list,
            0.0)
    };
    if (!seed_model_list.has_value()) return std::nullopt;
    const auto seed_step_norm{
        CalculateLocalFittingClusterModelTrustRegionStepNorm(
            outer_previous_state,
            key,
            *seed_model_list)
    };
    if (!seed_step_norm.has_value() ||
        *seed_step_norm > trust_region_radius + trust_region_tolerance)
    {
        return std::nullopt;
    }

    double damping{ 1.0 };
    while (damping >= std::numeric_limits<double>::epsilon())
    {
        auto candidate_model_list{
            detail::BuildLocalFittingSharedOffsetDampedModelList(
                previous_model_list,
                raw_model_list,
                shared_offset_model_list,
                damping)
        };
        if (candidate_model_list.has_value())
        {
            const auto step_norm{
                CalculateLocalFittingClusterModelTrustRegionStepNorm(
                    outer_previous_state,
                    key,
                    *candidate_model_list)
            };
            if (step_norm.has_value() &&
                *step_norm <= trust_region_radius + trust_region_tolerance)
            {
                LocalFittingBaseProposal proposal;
                proposal.state = outer_previous_state;
                proposal.effective_damping = damping;
                proposal.step_norm = *step_norm;
                for (std::size_t atom_position = 0;
                    atom_position < key.size();
                    atom_position++)
                {
                    const auto atom_index{ key.at(atom_position) };
                    proposal.state.at(atom_index).mdpde =
                        GaussianModel3DWithUncertainty{
                            candidate_model_list->at(atom_position),
                            raw_state.at(atom_index).mdpde
                                .GetStandardDeviationModel()
                        };
                }
                return proposal;
            }
        }
        damping *= 0.5;
    }
    return std::nullopt;
}

bool HasMaterialLocalFittingJointPolishChange(
    const std::vector<GaussianModel3D> & candidate_model_list,
    const std::vector<GaussianModel3D> & seed_model_list)
{
    if (candidate_model_list.size() != seed_model_list.size()) return false;
    for (std::size_t atom_position = 0;
        atom_position < candidate_model_list.size();
        atom_position++)
    {
        const auto change{
            detail::CalculateLocalFittingTransformedChange(
                candidate_model_list.at(atom_position),
                seed_model_list.at(atom_position))
        };
        if (std::any_of(
                change.value_list.begin(),
                change.value_list.end(),
                [](double value)
                {
                    return std::isfinite(value) &&
                        value >= kLocalFittingTransformedChangeTolerance;
                }))
        {
            return true;
        }
    }
    return false;
}

struct LocalFittingJointPolishProposal
{
    LocalFittingState state{};
    double effective_damping{ 0.0 };
    double step_norm{ 0.0 };
    std::vector<std::size_t> changed_atom_index_list{};
};

std::optional<LocalFittingJointPolishProposal>
BuildLocalFittingJointPolishProposal(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & outer_previous_state,
    const LocalFittingState & base_state,
    const LocalFittingClusterKey & key,
    const std::vector<LocalFittingObjectiveSampleRef> & sample_ref_list,
    const std::vector<double> & ridge_multiplier_list,
    double trust_region_radius)
{
    constexpr double trust_region_tolerance{ 1.0e-12 };
    std::vector<GroupKey> group_key_by_atom_position;
    std::vector<GaussianModel3D> base_model_list;
    group_key_by_atom_position.reserve(key.size());
    base_model_list.reserve(key.size());
    for (const auto atom_index : key)
    {
        group_key_by_atom_position.emplace_back(
            data_internal::GetGroupKey(context.at(atom_index).atom));
        base_model_list.emplace_back(
            base_state.at(atom_index).mdpde.GetModel());
    }
    const auto parameterization{
        detail::BuildLocalFittingJointPolishParameterization(
            group_key_by_atom_position,
            base_model_list)
    };
    if (!parameterization.has_value()) return std::nullopt;

    const Eigen::VectorXd zero_direction{
        Eigen::VectorXd::Zero(parameterization->ParameterCount())
    };
    const auto seed_model_list{
        parameterization->DecodeModels(zero_direction, 0.0)
    };
    if (!seed_model_list.has_value() ||
        std::any_of(
            seed_model_list->begin(),
            seed_model_list->end(),
            [](const GaussianModel3D & model)
            {
                return !detail::IsValidSecondStageGaussianModel(model);
            }))
    {
        return std::nullopt;
    }
    const auto seed_step_norm{
        CalculateLocalFittingClusterModelTrustRegionStepNorm(
            outer_previous_state,
            key,
            *seed_model_list)
    };
    if (!seed_step_norm.has_value() ||
        *seed_step_norm > trust_region_radius + trust_region_tolerance)
    {
        return std::nullopt;
    }
    const auto direction{
        BuildLocalFittingJointPolishDirection(
            context,
            base_state,
            *seed_model_list,
            key,
            sample_ref_list,
            ridge_multiplier_list,
            *parameterization)
    };
    if (!direction.has_value()) return std::nullopt;

    double damping{ 1.0 };
    while (damping >= std::numeric_limits<double>::epsilon())
    {
        auto candidate_model_list{
            parameterization->DecodeModels(*direction, damping)
        };
        if (candidate_model_list.has_value() &&
            std::none_of(
                candidate_model_list->begin(),
                candidate_model_list->end(),
                [](const GaussianModel3D & model)
                {
                    return !detail::IsValidSecondStageGaussianModel(model);
                }))
        {
            const auto step_norm{
                CalculateLocalFittingClusterModelTrustRegionStepNorm(
                    outer_previous_state,
                    key,
                    *candidate_model_list)
            };
            if (step_norm.has_value() &&
                *step_norm <= trust_region_radius + trust_region_tolerance)
            {
                if (!HasMaterialLocalFittingJointPolishChange(
                        *candidate_model_list,
                        *seed_model_list))
                {
                    return std::nullopt;
                }

                LocalFittingJointPolishProposal proposal;
                proposal.state = base_state;
                proposal.effective_damping = damping;
                proposal.step_norm = *step_norm;
                for (std::size_t atom_position = 0;
                    atom_position < key.size();
                    atom_position++)
                {
                    const auto atom_index{ key.at(atom_position) };
                    const auto & base_result{ base_state.at(atom_index) };
                    const auto & candidate_model{
                        candidate_model_list->at(atom_position)
                    };
                    proposal.state.at(atom_index).mdpde =
                        GaussianModel3DWithUncertainty{
                            candidate_model,
                            base_result.mdpde.GetStandardDeviationModel()
                        };
                    const auto base_coordinates{
                        detail::EncodeLocalFittingTransformedCoordinates(
                            base_result.mdpde.GetModel())
                    };
                    const auto candidate_coordinates{
                        detail::EncodeLocalFittingTransformedCoordinates(
                            candidate_model)
                    };
                    if (!base_coordinates.has_value() ||
                        !candidate_coordinates.has_value())
                    {
                        return std::nullopt;
                    }
                    if ((base_coordinates->array() !=
                            candidate_coordinates->array()).any())
                    {
                        proposal.changed_atom_index_list.emplace_back(
                            atom_index);
                    }
                }
                return proposal;
            }
        }
        damping *= 0.5;
    }
    return std::nullopt;
}

std::vector<Eigen::Vector3d> InterpolateLocalFittingTransformedEstimations(
    const std::vector<Eigen::Vector3d> & previous_estimation_list,
    const std::vector<Eigen::Vector3d> & candidate_estimation_list,
    double damping)
{
    if (previous_estimation_list.size() != candidate_estimation_list.size() ||
        !std::isfinite(damping) || damping < 0.0 || damping > 1.0)
    {
        throw std::invalid_argument(
            "Local fitting transformed interpolation inputs are invalid.");
    }
    std::vector<Eigen::Vector3d> interpolated_list;
    interpolated_list.reserve(previous_estimation_list.size());
    for (std::size_t i = 0; i < previous_estimation_list.size(); i++)
    {
        interpolated_list.emplace_back(
            (previous_estimation_list.at(i) +
                damping * (candidate_estimation_list.at(i) -
                    previous_estimation_list.at(i))).eval());
    }
    return interpolated_list;
}

std::optional<LocalFittingState> BuildLocalFittingCandidateState(
    const LocalFittingState & previous_state,
    const std::vector<Eigen::Vector3d> & previous_transformed_estimation_list,
    const std::vector<Eigen::Vector3d> & candidate_transformed_estimation_list,
    const LocalFittingState & uncertainty_state,
    const std::vector<std::size_t> & active_index_list)
{
    if (previous_transformed_estimation_list.size() != active_index_list.size() ||
        candidate_transformed_estimation_list.size() != active_index_list.size())
    {
        throw std::invalid_argument(
            "Local fitting candidate transformed coordinate count is inconsistent.");
    }
    auto candidate_state{ previous_state };
    for (std::size_t local_position = 0;
        local_position < active_index_list.size();
        local_position++)
    {
        const auto active_index{ active_index_list.at(local_position) };
        const auto & previous_transformed_estimation{
            previous_transformed_estimation_list.at(local_position)
        };
        const auto & candidate_transformed_estimation{
            candidate_transformed_estimation_list.at(local_position)
        };
        if (!previous_transformed_estimation.allFinite() ||
            !candidate_transformed_estimation.allFinite())
        {
            return std::nullopt;
        }
        if ((candidate_transformed_estimation.array() ==
                previous_transformed_estimation.array()).all())
        {
            auto & result{ candidate_state.at(active_index) };
            result.mdpde = GaussianModel3DWithUncertainty{
                previous_state.at(active_index).mdpde.GetModel(),
                uncertainty_state.at(active_index).mdpde.GetStandardDeviationModel()
            };
            continue;
        }
        const auto candidate_model{
            detail::DecodeLocalFittingTransformedCoordinates(
                candidate_transformed_estimation)
        };
        if (!candidate_model.has_value())
        {
            return std::nullopt;
        }
        if (!detail::IsValidSecondStageGaussianModel(*candidate_model))
        {
            return std::nullopt;
        }
        auto & result{ candidate_state.at(active_index) };
        result.mdpde = GaussianModel3DWithUncertainty{
            *candidate_model,
            uncertainty_state.at(active_index).mdpde.GetStandardDeviationModel()
        };
    }
    return candidate_state;
}

LocalFittingClusterObjectiveState
BuildInitialLocalFittingClusterObjectiveState(
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
    bool has_initial_objective{ false };
    if (initial_objective_samples.has_value())
    {
        has_initial_objective = CalculateLocalFittingObjectiveBreakdown(
            *initial_objective_samples,
            initial_objective_samples->scale_sample).has_value();
    }
    return LocalFittingClusterObjectiveState{
        std::move(initial_objective_samples),
        has_initial_objective
    };
}

void ReconcileLocalFittingClusterObjectiveState(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & previous_state,
    const detail::LocalFittingCouplingPartition & partition,
    LocalFittingClusterObjectiveStateMap & state_by_key)
{
    LocalFittingClusterObjectiveStateMap next_state_by_key;
    for (const auto & [key, objective_sample_ref_list] :
        partition.sample_id_list_by_key)
    {
        auto state_iter{ state_by_key.find(key) };
        if (state_iter != state_by_key.end())
        {
            next_state_by_key.emplace(key, std::move(state_iter->second));
            continue;
        }
        next_state_by_key.emplace(
            key,
            BuildInitialLocalFittingClusterObjectiveState(
                context,
                previous_state,
                key,
                objective_sample_ref_list));
    }
    state_by_key = std::move(next_state_by_key);
}

bool AreLocalFittingObjectiveReferencesLocked(
    const LocalFittingClusterObjectiveStateMap & state_by_key)
{
    for (const auto & state_entry : state_by_key)
    {
        const auto & state{ state_entry.second };
        if (state.scale_tracker.GetCommittedReference().has_value() &&
            !state.scale_tracker.IsLocked())
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
    LocalFittingClusterObjectiveStateMap & cluster_objective_state,
    LocalFittingObjectiveAttemptDiagnostic & diagnostic)
{
    const auto maximum_transformed_change{
        algorithm::GetMaximumParameterChange(
            SummarizeLocalFittingTransformedChanges(
                candidate_state,
                previous_state,
                key).percentile_stats)
    };
    auto & state{ cluster_objective_state.at(key) };
    auto objective_scale{ state.scale_tracker.GetCommittedReference() };
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
        }
    }

    std::optional<double> candidate_objective_value;
    std::optional<double> previous_objective_value;
    std::optional<double> best_objective_value;
    if (objective_scale.has_value())
    {
        diagnostic.objective_scale = *objective_scale;
        const auto evaluate_objective = [&](
            const std::optional<LocalFittingObjectiveSamples> & samples)
            -> std::optional<detail::LocalFittingObjectiveBreakdown>
        {
            if (!samples.has_value()) return std::nullopt;
            return CalculateLocalFittingObjectiveBreakdown(
                *samples,
                *objective_scale);
        };

        const auto candidate_objective{
            evaluate_objective(candidate_objective_samples)
        };
        const auto previous_objective{
            evaluate_objective(state.previous_objective_samples)
        };
        const auto best_objective{
            evaluate_objective(state.best_objective_samples)
        };
        if (candidate_objective.has_value() &&
            state.previous_objective_samples.has_value())
        {
            diagnostic.candidate_objective = candidate_objective;
            diagnostic.previous_objective = previous_objective;
            diagnostic.best_objective = best_objective;
            candidate_objective_value = candidate_objective->total_objective;
        }
        if (previous_objective.has_value())
        {
            previous_objective_value = previous_objective->total_objective;
        }
        if (best_objective.has_value())
        {
            best_objective_value = best_objective->total_objective;
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
    if (objective_scale.has_value())
    {
        diagnostic.rejected_by_previous = is_objective_deteriorated(
            candidate_objective_value,
            previous_objective_value);
        diagnostic.rejected_by_best = is_objective_deteriorated(
            candidate_objective_value,
            best_objective_value);
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

    if (candidate_objective_value.has_value() && objective_scale_sample.has_value())
    {
        state.scale_tracker.CommitScaleSample(*objective_scale_sample);
    }
    auto is_better_than_best{ !state.best_objective_samples.has_value() };
    if (state.best_objective_samples.has_value())
    {
        if (candidate_objective_value.has_value() != best_objective_value.has_value())
        {
            is_better_than_best = candidate_objective_value.has_value();
        }
        else if (candidate_objective_value.has_value())
        {
            if (detail::IsBetterLocalFittingAuditObjective(
                    *candidate_objective_value,
                    *best_objective_value,
                    kLocalFittingObjectiveTieRelativeTolerance))
            {
                is_better_than_best = true;
            }
            else if (detail::IsBetterLocalFittingAuditObjective(
                         *best_objective_value,
                         *candidate_objective_value,
                         kLocalFittingObjectiveTieRelativeTolerance))
            {
                is_better_than_best = false;
            }
            else
            {
                is_better_than_best =
                    maximum_transformed_change <
                    state.best_maximum_transformed_change;
            }
        }
        else
        {
            is_better_than_best =
                maximum_transformed_change <
                state.best_maximum_transformed_change;
        }
    }
    if (candidate_objective_value.has_value() && is_better_than_best)
    {
        state.best_objective_samples = candidate_objective_samples;
        state.best_maximum_transformed_change = maximum_transformed_change;
    }
    state.previous_objective_samples = std::move(candidate_objective_samples);
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
    const std::vector<LocalFittingClusterKey> & polish_eligible_key_list,
    const LocalFittingState & previous_state,
    const LocalFittingPolishProvenance & previous_polish_provenance,
    const LocalFittingState & raw_state,
    const std::vector<Eigen::Vector3d> & previous_transformed_estimation_list,
    const std::vector<Eigen::Vector3d> & raw_transformed_estimation_list,
    const std::vector<char> & suspicious_offset_atom_mask,
    const std::vector<double> & ridge_multiplier_list,
    LocalFittingClusterObjectiveStateMap & cluster_objective_state,
    const detail::LocalFittingTrustRegionStateSet & trust_region_state)
{
    LocalFittingCandidateSelection selection;
    selection.assembled_state = previous_state;
    selection.assembled_polish_provenance = previous_polish_provenance;
    for (const auto & [key, objective_sample_ref_list] :
        partition.sample_id_list_by_key)
    {
        const auto is_polish_eligible{
            std::find(
                polish_eligible_key_list.begin(),
                polish_eligible_key_list.end(),
                key) != polish_eligible_key_list.end()
        };
        if (is_polish_eligible) selection.polish_progress.eligible_count++;
        const auto trust_region_radius{ trust_region_state.GetRadius(key) };
        const auto contains_suspicious_atom{
            std::any_of(
                key.begin(),
                key.end(),
                [&](std::size_t atom_index)
                {
                    return suspicious_offset_atom_mask.at(atom_index) != 0;
                })
        };
        LocalFittingObjectiveAttemptDiagnostic base_diagnostic;
        base_diagnostic.trust_region_radius = trust_region_radius;
        std::optional<LocalFittingBaseProposal> base_proposal;
        if (!contains_suspicious_atom)
        {
            base_proposal = BuildLocalFittingSharedOffsetBaseProposal(
                context,
                previous_state,
                raw_state,
                key,
                trust_region_radius);
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
                detail::LimitLocalFittingTrustRegionDamping(
                    previous_cluster_estimation_list,
                    raw_cluster_estimation_list,
                    kLocalFittingTrustRegionParameterScale,
                    1.0,
                    trust_region_radius)
            };
            const auto base_cluster_estimation_list{
                InterpolateLocalFittingTransformedEstimations(
                    previous_cluster_estimation_list,
                    raw_cluster_estimation_list,
                    trust_region_damping.effective_damping)
            };
            auto base_state{
                BuildLocalFittingCandidateState(
                    previous_state,
                    previous_cluster_estimation_list,
                    base_cluster_estimation_list,
                    raw_state,
                    key)
            };
            if (base_state.has_value())
            {
                base_proposal = LocalFittingBaseProposal{
                    std::move(*base_state),
                    trust_region_damping.effective_damping,
                    trust_region_damping.step_norm
                };
            }
            else
            {
                base_diagnostic.is_invalid_model = true;
            }
        }
        if (!base_proposal.has_value())
        {
            if (is_polish_eligible) selection.polish_progress.skipped_count++;
            selection.rejected_key_list.emplace_back(key);
            selection.rejected_cluster_diagnostic_list.emplace_back(
                LocalFittingRejectedClusterDiagnostic{
                    key,
                    std::move(base_diagnostic)
                });
            continue;
        }
        base_diagnostic.effective_damping =
            base_proposal->effective_damping;
        base_diagnostic.trust_region_step_norm =
            base_proposal->step_norm;
        auto & base_state{ base_proposal->state };
        if (!TryCommitLocalFittingClusterCandidate(
                context,
                base_state,
                previous_state,
                key,
                objective_sample_ref_list,
                false,
                cluster_objective_state,
                base_diagnostic))
        {
            if (is_polish_eligible) selection.polish_progress.skipped_count++;
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
            selection.assembled_state.at(active_index) = base_state.at(active_index);
        }

        for (const auto atom_index : key)
        {
            const auto base_transformed{
                detail::EncodeLocalFittingTransformedCoordinates(
                    base_state.at(atom_index).mdpde.GetModel())
            };
            if (base_transformed.has_value() &&
                (base_transformed->array() ==
                    previous_transformed_estimation_list.at(atom_index).array()).all())
            {
                continue;
            }
            selection.assembled_polish_provenance.at(atom_index) = 0;
        }

        if (!is_polish_eligible)
        {
            continue;
        }
        auto polished_candidate{
            BuildLocalFittingJointPolishProposal(
                context,
                previous_state,
                base_state,
                key,
                objective_sample_ref_list,
                ridge_multiplier_list,
                trust_region_radius)
        };
        if (!polished_candidate.has_value())
        {
            selection.polish_progress.skipped_count++;
            continue;
        }
        LocalFittingObjectiveAttemptDiagnostic polish_diagnostic;
        polish_diagnostic.effective_damping =
            polished_candidate->effective_damping;
        polish_diagnostic.trust_region_radius = trust_region_radius;
        polish_diagnostic.trust_region_step_norm =
            polished_candidate->step_norm;
        if (!TryCommitLocalFittingClusterCandidate(
                context,
                polished_candidate->state,
                base_state,
                key,
                objective_sample_ref_list,
                true,
                cluster_objective_state,
                polish_diagnostic))
        {
            selection.polish_progress.rejected_count++;
            continue;
        }
        selection.polish_progress.accepted_count++;
        for (const auto active_index : key)
        {
            selection.assembled_state.at(active_index) =
                polished_candidate->state.at(active_index);
        }
        for (const auto atom_index :
            polished_candidate->changed_atom_index_list)
        {
            selection.assembled_polish_provenance.at(atom_index) = 1;
        }
        if (ShouldGrowLocalFittingTrustRegion(polish_diagnostic) &&
            std::find(
                selection.grow_trust_region_key_list.begin(),
                selection.grow_trust_region_key_list.end(),
                key) == selection.grow_trust_region_key_list.end())
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
    const FittedGaussianSnapshot & refit_model_snapshot,
    const FitOptions & options)
{
    auto sample_entries{
        BuildSecondStageAdjustedSamples(
            context,
            atom_index,
            refit_model_snapshot)
    };
    const auto & offset_model{ refit_model_snapshot.at(atom_index) };
    const auto & previous_model{ previous_result.mdpde.GetModel() };
    const auto is_acceptable = [&](const GaussianModel3D & model)
    {
        return detail::IsValidSecondStageGaussianModel(model) &&
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
                context.at(atom_index).alpha_r,
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
    std::size_t active_position_count)
{
    std::vector<char> visited(active_position_count, 0);
    std::vector<std::pair<std::size_t, std::size_t>> queue;
    queue.reserve(active_position_count);
    for (const auto seed_position : suspicious_offset_seed_position_list)
    {
        if (visited.at(seed_position) != 0) continue;
        visited.at(seed_position) = 1;
        queue.emplace_back(seed_position, 0);
    }

    std::vector<std::size_t> expanded_position_list;
    for (std::size_t queue_index = 0; queue_index < queue.size(); queue_index++)
    {
        const auto [active_position, depth]{ queue.at(queue_index) };
        expanded_position_list.emplace_back(active_position);
        if (active_coupling_graph.empty() || depth >= kSuspiciousOffsetClusterMaxDepth) continue;

        for (const auto & edge : active_coupling_graph.at(active_position))
        {
            if (!std::isfinite(edge.overlap) || edge.overlap < kSuspiciousOffsetClusterMinimumOverlap) continue;
            if (visited.at(edge.neighbor_index) != 0) continue;
            visited.at(edge.neighbor_index) = 1;
            queue.emplace_back(edge.neighbor_index, depth + 1);
        }
    }
    return expanded_position_list;
}

void ExpandPostRefitRollbackClusters(
    const std::vector<LocalFittingClusterKey> & cluster_key_list,
    const std::vector<std::size_t> & seed_atom_index_list,
    std::vector<char> & suspicious_mask)
{
    if (seed_atom_index_list.empty()) return;

    for (const auto & key : cluster_key_list)
    {
        const auto is_affected{
            std::any_of(
                key.begin(),
                key.end(),
                [&](std::size_t atom_index)
                {
                    return std::find(
                        seed_atom_index_list.begin(),
                        seed_atom_index_list.end(),
                        atom_index) != seed_atom_index_list.end();
                })
        };
        if (!is_affected) continue;
        for (const auto atom_index : key)
        {
            suspicious_mask.at(atom_index) = 1;
        }
    }
}

LocalFittingIterationResult RunLocalFittingIteration(
    const SecondStageLocalFittingContext & context,
    const std::vector<LocalFittingClusterKey> & cluster_key_list,
    const LocalFittingState & previous_state,
    const FitOptions & options,
    const std::vector<double> & ridge_multiplier_list)
{
    const auto selected_atom_size{ context.size() };
    auto current_snapshot{ BuildFittedGaussianSnapshot(previous_state) };
    std::map<LocalFittingClusterKey, JointOffsetSolveResult>
        joint_offset_result_by_key;
    for (const auto & key : cluster_key_list)
    {
        joint_offset_result_by_key.emplace(
            key,
            EstimateJointOffsets(
                context,
                key,
                current_snapshot,
                ridge_multiplier_list,
                !options.quiet_mode &&
                    Logger::GetLogLevel() >= LogLevel::Debug));
    }
    LocalFittingClusterHealthMap health_by_key;
    for (const auto & [key, result] : joint_offset_result_by_key)
    {
        for (std::size_t i = 0; i < key.size(); i++)
        {
            const auto atom_index{ key.at(i) };
            current_snapshot.at(atom_index) =
                current_snapshot.at(atom_index).WithOffset(
                    result.offset(static_cast<Eigen::Index>(i)));
        }
        health_by_key.emplace(
            key,
            LocalFittingClusterHealth{ result.status });
    }

    auto iteration_state{ previous_state };
    std::vector<char> suspicious_offset_mask(selected_atom_size, 0);
    for (const auto & [key, result] : joint_offset_result_by_key)
    {
        std::vector<std::size_t> seed_position_list;
        for (std::size_t position = 0; position < key.size(); position++)
        {
            const auto atom_index{ key.at(position) };
            if (IsSuspiciousJointOffset(
                    context.at(atom_index).sample_entries,
                    previous_state.at(atom_index).mdpde.GetModel(),
                    current_snapshot.at(atom_index),
                    options))
            {
                seed_position_list.emplace_back(position);
            }
        }
        const auto suspicious_position_list{
            ExpandSuspiciousOffsetClusters(
                result.active_coupling_graph,
                seed_position_list,
                key.size())
        };
        for (const auto position : suspicious_position_list)
        {
            suspicious_offset_mask.at(key.at(position)) = 1;
        }
    }
    for (std::size_t atom_index = 0;
        atom_index < suspicious_offset_mask.size();
        atom_index++)
    {
        if (suspicious_offset_mask.at(atom_index) == 0) continue;
        current_snapshot.at(atom_index) =
            previous_state.at(atom_index).mdpde.GetModel();
    }

    std::vector<GroupKey> group_key_by_atom_index;
    group_key_by_atom_index.reserve(context.size());
    for (const auto & atom_context : context)
    {
        group_key_by_atom_index.emplace_back(
            data_internal::GetGroupKey(atom_context.atom));
    }
    const auto refit_model_snapshot{
        detail::BuildLocalFittingGroupMedianModelList(
            group_key_by_atom_index,
            current_snapshot)
    };
    std::vector<std::size_t> post_refit_suspicious_seed_atom_index_list;
    for (auto & [key, health] : health_by_key)
    {
        for (const auto atom_index : key)
        {
            if (suspicious_offset_mask.at(atom_index) != 0) continue;

            auto refit_result{
                FitAtomWithJointOffsetFallback(
                    context,
                    atom_index,
                    previous_state.at(atom_index),
                    refit_model_snapshot,
                    options)
            };
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
    ExpandPostRefitRollbackClusters(
        cluster_key_list,
        post_refit_suspicious_seed_atom_index_list,
        suspicious_offset_mask);
    for (std::size_t atom_index = 0;
        atom_index < suspicious_offset_mask.size();
        atom_index++)
    {
        if (suspicious_offset_mask.at(atom_index) == 0) continue;
        iteration_state.at(atom_index) = previous_state.at(atom_index);
    }

    LocalFittingIterationResult iteration_result;
    iteration_result.state = std::move(iteration_state);
    iteration_result.suspicious_offset_atom_mask =
        std::move(suspicious_offset_mask);
    iteration_result.health_by_key = std::move(health_by_key);
    return iteration_result;
}

void ApplyLocalFittingState(
    ModelObject & model_object,
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & iteration_state)
{
    auto analysis{ model_object.EditAnalysis() };
    for (std::size_t i = 0; i < context.size(); i++)
    {
        auto local_editor{
            analysis.EnsureAtomLocalPotential(*context.at(i).atom)
        };
        local_editor.SetGaussianResult(iteration_state.at(i));
    }
}

LocalFittingOffsetStats SummarizeLocalFittingOffsets(
    const LocalFittingState & state)
{
    LocalFittingOffsetStats stats;
    stats.atom_count = state.size();
    std::vector<double> absolute_offset_list;
    absolute_offset_list.reserve(state.size());
    for (const auto & result : state)
    {
        const auto offset{ result.mdpde.GetModel().GetOffset() };
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

        if (diagnostic.is_invalid_model)
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
        if (!diagnostic.candidate_objective.has_value())
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

constexpr std::array<std::string_view, 6> kLocalFittingProgressHeaderList{
    "Try/Acc",
    "Atom A/T",
    "Cluster A/R",
    "Polish E/A/R/S",
    "Suspicious",
    "dMax A/R"
};

std::string FormatLocalFittingProgressRow(
    const LocalFittingProgressColumnWidths & column_widths,
    const std::array<std::string, 6> & cell_list)
{
    std::ostringstream stream;
    for (std::size_t i = 0; i < cell_list.size(); i++)
    {
        if (i > 0) stream << " | ";
        stream
            << std::left
            << std::setw(static_cast<int>(column_widths.at(i)))
            << cell_list.at(i);
    }
    return stream.str();
}

LocalFittingProgressColumnWidths BuildLocalFittingProgressColumnWidths(
    std::size_t atom_size)
{
    const auto maximum_iteration_text{
        std::to_string(kLocalFittingMaximumIterations)
    };
    const auto maximum_atom_text{ std::to_string(atom_size) };
    const auto maximum_change_text{
        FormatLocalFittingProgressMaximum(
            std::numeric_limits<double>::max())
    };
    const std::array<std::string, 6> maximum_cell_list{
        maximum_iteration_text + "/" + maximum_iteration_text,
        maximum_atom_text + "/" + maximum_atom_text,
        maximum_atom_text + "/" + maximum_atom_text,
        maximum_atom_text + "/" + maximum_atom_text + "/" +
            maximum_atom_text + "/" + maximum_atom_text,
        maximum_atom_text,
        maximum_change_text + "/" + maximum_change_text
    };

    LocalFittingProgressColumnWidths column_widths;
    for (std::size_t i = 0; i < column_widths.size(); i++)
    {
        column_widths.at(i) = std::max(
            kLocalFittingProgressHeaderList.at(i).size(),
            maximum_cell_list.at(i).size());
    }
    return column_widths;
}

void LogLocalFittingProgressHeader(
    const FitOptions & options,
    const LocalFittingProgressColumnWidths & column_widths)
{
    if (options.quiet_mode) return;
    std::array<std::string, 6> header_list;
    for (std::size_t i = 0; i < header_list.size(); i++)
    {
        header_list.at(i) = kLocalFittingProgressHeaderList.at(i);
    }
    Logger::Log(
        LogLevel::Info,
        FormatLocalFittingProgressRow(column_widths, header_list));
}

void LogLocalFittingIterationProgress(
    const FitOptions & options,
    const LocalFittingProgressColumnWidths & column_widths,
    const LocalFittingIterationProgress & progress)
{
    if (options.quiet_mode) return;

    const std::array<std::string, 6> cell_list{
        std::to_string(progress.attempt_number) + "/" +
            std::to_string(progress.accepted_iteration_count),
        std::to_string(progress.active_atom_count) + "/" +
            std::to_string(progress.terminal_atom_count),
        std::to_string(progress.accepted_cluster_count) + "/" +
            std::to_string(progress.rejected_cluster_count),
        std::to_string(progress.polish_progress.eligible_count) + "/" +
            std::to_string(progress.polish_progress.accepted_count) + "/" +
            std::to_string(progress.polish_progress.rejected_count) + "/" +
            std::to_string(progress.polish_progress.skipped_count),
        std::to_string(progress.suspicious_atom_count),
        FormatLocalFittingProgressMaximum(
            progress.accepted_maximum_transformed_change) + "/" +
            FormatLocalFittingProgressMaximum(
                progress.raw_maximum_transformed_change)
    };
    Logger::ProgressLine(
        FormatLocalFittingProgressRow(column_widths, cell_list));
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
    const LocalFittingBestAuditState & best_audit_state,
    std::optional<bool> final_uses_polish)
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
    message << ", final_uses_polish=";
    if (!final_uses_polish.has_value())
    {
        message << "unavailable";
    }
    else
    {
        message << (*final_uses_polish ? "yes" : "no");
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
    const auto atom_size{ context.size() };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run 2nd-stage local atom fitting with iterations...");
    }

    auto initial_state_build_result{
        BuildInitialLocalFittingState(context, model_object.GetAnalysisView())
    };
    if (!initial_state_build_result.has_value())
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
                "best_audit_objective=unavailable, final_uses_polish=unavailable.");
        }
        return;
    }
    LogSecondStageSeedRepairs(initial_state_build_result->repair_record_list, options);
    auto previous_state{ std::move(initial_state_build_result->state) };
    LocalFittingPolishProvenance previous_polish_provenance(atom_size, 0);
    const auto coupling_topology{
        BuildLocalFittingCouplingTopology(context, previous_state)
    };
    LogLocalFittingCouplingTopology(coupling_topology, options);
    auto best_audit_state{
        BuildInitialLocalFittingBestAuditState(context, previous_state)
    };

    std::vector<char> suspicious_offset_atom_mask(atom_size, 0);
    PersistentTerminalFailureStateMap persistent_terminal_failure_state_by_key;
    std::vector<char> terminal_fallback_atom_mask(atom_size, 0);
    LocalFittingTerminalSummary terminal_summary;
    LocalFittingClusterObjectiveStateMap cluster_objective_state;
    detail::LocalFittingTrustRegionStateSet trust_region_state{
        detail::LocalFittingTrustRegionOptions{
            kLocalFittingTrustRegionInitialRadius,
            kLocalFittingTrustRegionMinimumRadius,
            kLocalFittingTrustRegionMaximumRadius,
            kLocalFittingTrustRegionShrinkFactor,
            kLocalFittingTrustRegionGrowthFactor
        }
    };
    const auto progress_column_widths{
        BuildLocalFittingProgressColumnWidths(atom_size)
    };
    LogLocalFittingProgressHeader(options, progress_column_widths);

    std::size_t accepted_iteration_count{ 0 };
    std::size_t audit_patience_count{ 0 };
    for (std::size_t iter = 0; iter < kLocalFittingMaximumIterations; iter++)
    {
        const auto active_index_list{
            BuildEligibleLocalFittingActiveIndexList(
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
                best_audit_state,
                UsesLocalFittingPolish(previous_polish_provenance));
            return;
        }

        const auto cluster_partition{
            detail::BuildLocalFittingCouplingPartition(
                coupling_topology,
                active_index_list)
        };
        std::vector<LocalFittingClusterKey> cluster_key_list;
        cluster_key_list.reserve(cluster_partition.sample_id_list_by_key.size());
        for (const auto & [key, sample_id_list] :
            cluster_partition.sample_id_list_by_key)
        {
            static_cast<void>(sample_id_list);
            cluster_key_list.emplace_back(key);
        }

        ReconcileLocalFittingClusterObjectiveState(
            context,
            previous_state,
            cluster_partition,
            cluster_objective_state);
        trust_region_state.Reconcile(cluster_key_list);

        std::vector<double> joint_offset_ridge_multiplier_list(atom_size, 1.0);
        for (std::size_t atom_index = 0; atom_index < atom_size; atom_index++)
        {
            if (suspicious_offset_atom_mask.at(atom_index) == 0) continue;
            joint_offset_ridge_multiplier_list.at(atom_index) =
                kSuspiciousJointOffsetRidgeMultiplier;
        }

        auto iteration_result{
            RunLocalFittingIteration(
                context,
                cluster_key_list,
                previous_state,
                options,
                joint_offset_ridge_multiplier_list)
        };
        const auto current_health_by_key{ std::move(iteration_result.health_by_key) };
        suspicious_offset_atom_mask =
            std::move(iteration_result.suspicious_offset_atom_mask);
        const auto iteration_suspicious_atom_count{
            static_cast<std::size_t>(std::count_if(
                suspicious_offset_atom_mask.begin(),
                suspicious_offset_atom_mask.end(),
                [](char is_suspicious)
                {
                    return is_suspicious != 0;
                }))
        };
        const auto has_suspicious_offset_fallback{
            iteration_suspicious_atom_count > 0
        };

        std::size_t stationarity_ineligible_cluster_count{ 0 };
        std::vector<LocalFittingClusterKey> polish_eligible_key_list;
        for (const auto & [key, health] : current_health_by_key)
        {
            if (!IsJointOffsetSolveStationarityEligible(
                    health.joint_offset_status) ||
                !health.is_refit_stationarity_eligible)
            {
                stationarity_ineligible_cluster_count++;
                continue;
            }
            const auto contains_suspicious_atom{
                std::any_of(
                    key.begin(),
                    key.end(),
                    [&](std::size_t atom_index)
                    {
                        return suspicious_offset_atom_mask.at(atom_index) != 0;
                    })
            };
            if (!contains_suspicious_atom)
            {
                polish_eligible_key_list.emplace_back(key);
            }
        }
        const auto raw_state{ std::move(iteration_result.state) };
        const auto raw_fixed_point_change_summary{
            SummarizeLocalFittingTransformedChanges(
                raw_state,
                previous_state,
                active_index_list)
        };
        const auto previous_transformed_estimation_list{
            BuildLocalFittingTransformedEstimationList(previous_state)
        };
        const auto raw_transformed_estimation_list{
            BuildLocalFittingTransformedEstimationList(raw_state)
        };

        auto working_cluster_objective_state{ cluster_objective_state };
        auto selection{
            SelectLocalFittingClusterCandidates(
                context,
                cluster_partition,
                polish_eligible_key_list,
                previous_state,
                previous_polish_provenance,
                raw_state,
                previous_transformed_estimation_list,
                raw_transformed_estimation_list,
                suspicious_offset_atom_mask,
                joint_offset_ridge_multiplier_list,
                working_cluster_objective_state,
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
            RejectLocalFittingCombinedCandidate(
                previous_state,
                previous_polish_provenance,
                cluster_key_list,
                selection);
        }
        else
        {
            cluster_objective_state = std::move(working_cluster_objective_state);
        }

        auto assembled_state{ std::move(selection.assembled_state) };
        auto assembled_polish_provenance{
            std::move(selection.assembled_polish_provenance)
        };
        const auto terminal_failure_by_key{
            UpdatePersistentTerminalFailureState(
                selection.accepted_key_list,
                suspicious_offset_atom_mask,
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
            previous_polish_provenance,
            terminal_fallback_atom_mask,
            assembled_state,
            assembled_polish_provenance);
        ReconcileLocalFittingBestAuditTerminalFallback(
            context,
            terminal_key_list,
            assembled_state,
            assembled_polish_provenance,
            accepted_iteration_count + 1,
            best_audit_state);
        if (!terminal_key_list.empty())
        {
            for (const auto & key : terminal_key_list)
            {
                for (const auto atom_index : key)
                {
                    suspicious_offset_atom_mask.at(atom_index) = 0;
                }
            }
        }

        trust_region_state.Grow(selection.grow_trust_region_key_list);
        const auto trust_region_radius_update{
            trust_region_state.Shrink(selection.rejected_key_list)
        };
        const auto terminal_atom_count{ terminal_summary.AtomCount() };
        LocalFittingIterationProgress progress{
            iter + 1,
            accepted_iteration_count,
            atom_size - terminal_atom_count,
            terminal_atom_count,
            selection.accepted_key_list.size(),
            selection.rejected_key_list.size(),
            selection.polish_progress,
            iteration_suspicious_atom_count,
            std::nullopt,
            SummarizeLocalFittingProgressMaximum(
                raw_fixed_point_change_summary.maximum_list)
        };
        if (selection.accepted_key_list.empty())
        {
            LogRejectedLocalFittingClusterDiagnostics(
                options,
                selection.rejected_cluster_diagnostic_list);
            LogLocalFittingIterationProgress(
                options,
                progress_column_widths,
                progress);
            if (!trust_region_radius_update.changed_key_list.empty() && iter + 1 < kLocalFittingMaximumIterations)
            {
                continue;
            }

            const auto fallback{
                SelectLocalFittingFallbackState(previous_state, best_audit_state.best)
            };
            const auto final_uses_polish{
                SelectLocalFittingFallbackUsesPolish(
                    previous_polish_provenance,
                    best_audit_state.best)
            };
            ApplyLocalFittingState(model_object, context, fallback);
            LogSecondStageLocalFittingSummary(
                options,
                accepted_iteration_count,
                iter + 1 >= kLocalFittingMaximumIterations ?
                    "maximum-iterations" :
                    "all-rejected-minimum-radius",
                best_audit_state,
                final_uses_polish);
            return;
        }

        const auto transformed_change_summary{
            SummarizeLocalFittingTransformedChanges(
                assembled_state,
                previous_state,
                active_index_list)
        };

        accepted_iteration_count++;
        const auto improved_best_audit{
            TryUpdateLocalFittingBestAuditState(
                context,
                assembled_state,
                assembled_polish_provenance,
                accepted_iteration_count,
                best_audit_state)
        };
        audit_patience_count = detail::AdvanceLocalFittingAuditPatience(
            audit_patience_count,
            improved_best_audit,
            !selection.rejected_key_list.empty() &&
                !trust_region_radius_update.changed_key_list.empty());
        progress.accepted_iteration_count = accepted_iteration_count;
        progress.accepted_maximum_transformed_change =
            SummarizeLocalFittingProgressMaximum(
                transformed_change_summary.maximum_list);
        LogRejectedLocalFittingClusterDiagnostics(
            options,
            selection.rejected_cluster_diagnostic_list);
        LogLocalFittingIterationProgress(
            options,
            progress_column_widths,
            progress);

        if (audit_patience_count >= kLocalFittingAuditPatience)
        {
            const auto fallback{
                SelectLocalFittingFallbackState(
                    assembled_state,
                    best_audit_state.best)
            };
            const auto final_uses_polish{
                SelectLocalFittingFallbackUsesPolish(
                    assembled_polish_provenance,
                    best_audit_state.best)
            };
            ApplyLocalFittingState(model_object, context, fallback);
            LogSecondStageLocalFittingSummary(
                options,
                accepted_iteration_count,
                "audit-patience",
                best_audit_state,
                final_uses_polish);
            return;
        }

        const auto converged{
            stationarity_ineligible_cluster_count == 0 &&
            !has_suspicious_offset_fallback &&
            selection.rejected_key_list.empty() &&
            AreLocalFittingObjectiveReferencesLocked(cluster_objective_state) &&
            IsLocalFittingTransformedChangeConverged(
                transformed_change_summary.percentile_stats,
                transformed_change_summary.maximum_list) &&
            IsLocalFittingTransformedChangeConverged(
                raw_fixed_point_change_summary.percentile_stats,
                raw_fixed_point_change_summary.maximum_list)
        };
        if (converged)
        {
            const auto accepted_offset_stats{
                SummarizeLocalFittingOffsets(assembled_state)
            };
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
                    transformed_change_summary.percentile_stats,
                    accepted_offset_stats);
            }
            LogSecondStageLocalFittingSummary(
                options,
                accepted_iteration_count,
                "converged",
                best_audit_state,
                UsesLocalFittingPolish(assembled_polish_provenance));
            return;
        }

        if (iter + 1 == kLocalFittingMaximumIterations)
        {
            const auto fallback{
                SelectLocalFittingFallbackState(
                    assembled_state,
                    best_audit_state.best)
            };
            const auto final_uses_polish{
                SelectLocalFittingFallbackUsesPolish(
                    assembled_polish_provenance,
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
                best_audit_state,
                final_uses_polish);
            return;
        }
        previous_state = std::move(assembled_state);
        previous_polish_provenance = std::move(assembled_polish_provenance);
    }
}

} // namespace rhbm_gem::core
