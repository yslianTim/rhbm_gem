#include <cstddef>
#include <rhbm_gem/core/GaussianEstimator.hpp>

#include "data/detail/AtomClassifier.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/algorithm/AndersonAcceleration.hpp>
#include <rhbm_gem/utils/algorithm/ClusteredAndersonAcceleration.hpp>
#include <rhbm_gem/utils/algorithm/ClusteredFittingQualityState.hpp>
#include <rhbm_gem/utils/algorithm/ConvergenceFreezeTracker.hpp>
#include <rhbm_gem/utils/algorithm/DependencyThawHysteresisTracker.hpp>
#include <rhbm_gem/utils/algorithm/IterationState.hpp>
#include <rhbm_gem/utils/algorithm/NormalizedChange.hpp>
#include <rhbm_gem/utils/algorithm/ParameterChangeStats.hpp>
#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>
#include <rhbm_gem/utils/algorithm/ScaleReferenceTracker.hpp>
#include <rhbm_gem/utils/algorithm/WeightedRidgeSolver.hpp>
#include <rhbm_gem/utils/algorithm/WeightedRidgeSystem.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/SampleFilter.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem::core {
namespace {
constexpr double kLocalFittingNormalizedChangeTolerance{ 1.0e-4 };
constexpr double kLocalFittingNormalizedChangeScaleFloor{ 1.0 };
constexpr std::size_t kMinimumAlphaRTrainingSampleCount{ 10 };
constexpr std::size_t kMinimumAlphaGTrainingMemberCount{ 10 };
constexpr double kNeighborContributionDistanceMax{ 2.5 };
constexpr double kNeighborAtomSearchRange{ 2.0 * kNeighborContributionDistanceMax };
constexpr std::size_t kLocalFittingMaximumIterations{ 200 };
constexpr double kLocalFittingParameterChangeTolerance{ 1.0e-6 };
constexpr double kLocalFittingChangePercentile{ 0.99 };
constexpr algorithm::RobustLossKind kSecondStageRobustLossKind{ algorithm::RobustLossKind::Cauchy };
constexpr int kRobustLossMaximumIterations{ 50 };
constexpr double kRobustScaleMultiplier{ 1.4826 };
constexpr double kRobustScaleMin{ 1.0e-12 };
constexpr double kRobustLossCutoffMultiplier{ 1.345 };
constexpr double kJointOffsetRidgeRatio{ 1.0e-3 };
constexpr double kJointOffsetRidgeRatioMin{ 1.0e-4 };
constexpr double kJointOffsetRidgeRatioMax{ 1.0 };
constexpr double kJointOffsetRidgeGrowth{ 2.0 };
constexpr double kJointOffsetRidgeShrink{ 0.8 };
constexpr double kSuspiciousJointOffsetRidgeMultiplier{ 10.0 };
constexpr double kJointOffsetCollinearityOverlapThreshold{ 0.98 };
constexpr double kCollinearJointOffsetRidgeMultiplier{ 10.0 };
constexpr double kJointOffsetIrlsScaleFloor{ 1.0e-2 };
constexpr double kJointOffsetIrlsNormalizedChangeTolerance{ 1.0e-6 };
constexpr double kJointOffsetIrlsObjectiveRelativeTolerance{ 1.0e-10 };
constexpr std::size_t kLocalFittingAndersonHistoryDepth{ 5 };
constexpr double kLocalFittingAndersonCoefficientL1Limit{ 10.0 };
constexpr double kLocalFittingAndersonRegularization{ 1.0e-4 };
constexpr double kLocalFittingAndersonCoefficientAbsLimit{ 5.0 };
constexpr std::array<double, 5> kLocalFittingAccelerationDampingList{ 1.0, 0.5, 0.25, 0.125, 0.0625 };
constexpr double kLocalFittingFreezeChangeRatio{ 0.1 };
constexpr int kLocalFittingFreezeStableIterations{ 3 };
constexpr double kLocalFittingDependencyThawHysteresisGrowth{ 2.0 };
constexpr double kLocalFittingDependencyThawHysteresisMax{ 8.0 };
constexpr double kLocalFittingDependencyThawHysteresisFrozenDecay{ 0.9 };
constexpr int kLocalFittingDependencyThawMaximumCount{ 5 };
constexpr double kLocalFittingObjectiveTieRelativeTolerance{ 1.0e-8 };
constexpr double kLocalFittingConvergenceObjectiveRelativeTolerance{ 1.0e-3 };
constexpr std::size_t kLocalFittingObjectiveScaleWarmupCount{ 5 };
constexpr double kLocalFittingObjectiveResidualScaleFloorRatio{ 1.0e-6 };
constexpr double kLocalFittingWidthPriorPenaltyWeight{ 1.0e-2 };
constexpr double kLocalFittingOffsetPlausibilityPenaltyWeight{ 1.0e-2 };
constexpr double kLocalFittingMovementPenaltyWeight{ 1.0e-2 };
constexpr double kLocalFittingWidthPriorLogScale{ 0.35 };
constexpr double kLocalFittingMovementAmplitudeLogScale{ 0.50 };
constexpr double kLocalFittingMovementWidthLogScale{ 0.35 };
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

using GaussianFittingState = algorithm::IterationState<LocalGaussianResult, Eigen::VectorXd>;

enum class LocalFittingCandidateKind
{
    Anderson,
    FixedPoint
};

enum class LocalFittingBacktrackingStopReason
{
    MaximumGlobalRidge,
    MaximumIterationLimit
};

struct LocalFittingAccelerationAttempt
{
    LocalFittingCandidateKind kind{ LocalFittingCandidateKind::FixedPoint };
    double damping{ 1.0 };
};

struct ActiveCouplingEdge
{
    std::size_t neighbor_index{ 0 };
    double overlap{ 0.0 };
};

using ActiveCouplingGraph = std::vector<std::vector<ActiveCouplingEdge>>;

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
    Eigen::VectorXd offset{};
    ActiveCouplingGraph active_coupling_graph{};
};

struct JointOffsetBuildResult
{
    algorithm::WeightedRidgeSystem system{};
    ActiveCouplingGraph active_coupling_graph{};
};

std::vector<double> CombineLocalFittingRidgeMultiplierLists(
    const std::vector<double> & left_multiplier_list,
    const std::vector<double> & right_multiplier_list)
{
    if (left_multiplier_list.size() != right_multiplier_list.size())
    {
        throw std::invalid_argument("Local fitting ridge multiplier sizes are inconsistent.");
    }
    std::vector<double> multiplier_list(left_multiplier_list.size(), 1.0);
    for (std::size_t i = 0; i < multiplier_list.size(); i++)
    {
        multiplier_list.at(i) = std::max(left_multiplier_list.at(i), right_multiplier_list.at(i));
    }
    return multiplier_list;
}

struct LocalFittingIterationResult
{
    GaussianFittingState state{};
    std::vector<std::size_t> suspicious_offset_state_index_list{};
};

struct LocalFittingObjectiveAtomModel
{
    std::size_t atom_index{ 0 };
    GaussianModel3D model{};
};

struct LocalFittingObjectiveSamples
{
    std::vector<double> residual_list{};
    std::vector<double> response_list{};
    std::vector<LocalFittingObjectiveAtomModel> active_model_list{};
};

struct LocalFittingObjectiveSampleRef
{
    std::size_t atom_index{ 0 };
    std::size_t sample_index{ 0 };
};

using LocalFittingClusterKey = algorithm::ClusterKey;

enum class LocalFittingClusterAttemptState
{
    Pending,
    Stopped,
    AcceptedAnderson,
    AcceptedFixedPoint
};

struct LocalFittingClusterWork
{
    std::vector<LocalFittingObjectiveSampleRef> objective_sample_ref_list{};
    LocalFittingClusterAttemptState attempt_state{ LocalFittingClusterAttemptState::Pending };
    bool rejected_anderson{ false };
};

using LocalFittingClusterMap = std::map<LocalFittingClusterKey, LocalFittingClusterWork>;

struct LocalFittingCandidateSelection
{
    GaussianFittingState assembled_state{};
    std::vector<algorithm::ClusteredFittingQualityAcceptedScore<LocalFittingObjectiveSamples>> accepted_score_list{};
    std::vector<LocalFittingClusterKey> accepted_key_list{};
    std::vector<LocalFittingClusterKey> rejected_key_list{};
    LocalFittingAccelerationAttempt accepted_acceleration_attempt{};
    bool has_objective_backtracking_rejection{ false };
};

struct SecondStageNeighborSample
{
    std::size_t atom_index{ 0 };
    double distance{ 0.0 };
};

struct SecondStageAtomContext
{
    LocalPotentialSampleList sample_entries{};
    std::vector<std::size_t> selected_neighbor_index_list{};
    std::vector<std::vector<SecondStageNeighborSample>> sample_neighbor_list{};
    double alpha_r{ 0.0 };
    double prior_width{ 1.0 };
};

struct SecondStageLocalFittingContext
{
    std::vector<AtomObject *> atom_list{};
    std::vector<SecondStageAtomContext> atom_context_list{};

    std::size_t AtomSize() const { return atom_list.size(); }
};

struct ParameterSummaryStats
{
    double mean{ 0.0 };
    double standard_deviation{ 0.0 };
};

struct GaussianModelParameterSamples
{
    std::vector<double> amplitude_list{};
    std::vector<double> width_list{};
    std::vector<double> offset_list{};
};

enum class LocalFittingPass
{
    FirstStage,
    ThirdStage
};

std::vector<AtomLocalPotentialEditor> BuildAtomLocalEditors(
    ModelObject & model_object,
    const std::vector<AtomObject *> & atom_list)
{
    auto analysis{ model_object.EditAnalysis() };
    std::vector<AtomLocalPotentialEditor> local_editor_list;
    local_editor_list.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        local_editor_list.emplace_back(analysis.EnsureAtomLocalPotential(*atom));
    }
    return local_editor_list;
}

bool HasEnoughSamplesInFitRange(
    const LocalPotentialSampleList & sample_entries,
    double fit_range_min,
    double fit_range_max,
    std::size_t minimum_sample_count)
{
    std::size_t count{ 0 };
    for (const auto & sample : sample_entries)
    {
        if (sample.point.distance < fit_range_min || sample.point.distance > fit_range_max) continue;
        count++;
        if (count >= minimum_sample_count) return true;
    }
    return false;
}

RHBMExecutionOptions MakeExecutionOptions(const FitOptions & options)
{
    RHBMExecutionOptions execution_options;
    execution_options.quiet_mode = false;
    execution_options.thread_size = options.thread_size;
    return execution_options;
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
    return std::abs(candidate_offset_response) >
        kSuspiciousCompensationResponseRatio * reference_scale;
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

rhbm_trainer::RHBMTrainingOptions MakeTrainingOptions(const FitOptions & options)
{
    rhbm_trainer::RHBMTrainingOptions training_options;
    training_options.execution_options = MakeExecutionOptions(options);
    return training_options;
}

std::size_t GetMinimumDatasetResponseCount(const std::vector<RHBMMemberDataset> & dataset_list)
{
    std::size_t minimum_response_count{ std::numeric_limits<std::size_t>::max() };
    for (const auto & dataset : dataset_list)
    {
        const auto response_count{ static_cast<std::size_t>(dataset.y.size()) };
        if (response_count < minimum_response_count)
        {
            minimum_response_count = response_count;
        }
    }
    return minimum_response_count;
}

LocalPotentialSampleList BuildSamplesForZeroOffsetGaussianFit(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & model)
{
    LocalPotentialSampleList updated_sample_entries;
    updated_sample_entries.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        const auto response{ static_cast<float>(CalculateZeroOffsetResponse(sample, model)) };
        updated_sample_entries.emplace_back(LocalPotentialSample{ response, sample.point });
    }
    return updated_sample_entries;
}

LocalGaussianResult DecodeLocalGaussianResult(
    double alpha_r,
    const RHBMBetaEstimateResult & fit_result,
    double offset = 0.0)
{
    const auto ols_model{
        linearization_service::DecodeParameterVector(fit_result.beta_ols).WithOffset(offset)
    };
    const auto mdpde_model{
        linearization_service::DecodeParameterVector(fit_result.beta_mdpde).WithOffset(offset)
    };
    return LocalGaussianResult{
        alpha_r,
        GaussianModel3DWithUncertainty{ ols_model, GaussianModel3DUncertainty{} },
        GaussianModel3DWithUncertainty{ mdpde_model, GaussianModel3DUncertainty{} },
        std::nullopt,
        false,
        0.0,
        fit_result
    };
}

GroupGaussianResult DecodeGroupGaussianResult(
    double alpha_g,
    const RHBMGroupEstimationResult & result,
    double offset)
{
    const auto prior{
        linearization_service::DecodeParameterVector(result.mu_prior, result.capital_lambda)
    };
    return GroupGaussianResult{
        alpha_g,
        linearization_service::DecodeParameterVector(result.mu_mean).WithOffset(offset),
        linearization_service::DecodeParameterVector(result.mu_mdpde).WithOffset(offset),
        GaussianModel3DWithUncertainty{
            prior.GetModel().WithOffset(offset),
            prior.GetStandardDeviationModel()
        }
    };
}

std::vector<LocalGaussianResult> DecodeMemberGaussianResults(
    const RHBMGroupEstimationResult & result,
    const std::vector<LocalGaussianResult> & member_result_list)
{
    const auto member_count{ static_cast<std::size_t>(result.beta_posterior_matrix.cols()) };
    if (member_result_list.size() != member_count)
    {
        throw std::invalid_argument("Group Gaussian member result count is inconsistent.");
    }
    if (result.capital_sigma_posterior_list.size() != member_count)
    {
        throw std::invalid_argument("Group Gaussian member result count is inconsistent.");
    }
    eigen_validation::RequireVectorSize(
        result.outlier_flag_array, result.beta_posterior_matrix.cols(),
        "outlier_flag_array", "Group Gaussian member result count is inconsistent.");
    eigen_validation::RequireVectorSize(
        result.statistical_distance_array, result.beta_posterior_matrix.cols(),
        "statistical_distance_array", "Group Gaussian member result count is inconsistent.");

    std::vector<LocalGaussianResult> member_results;
    member_results.reserve(member_count);
    for (Eigen::Index i = 0; i < result.beta_posterior_matrix.cols(); i++)
    {
        const auto member_index{ static_cast<std::size_t>(i) };
        const auto offset{
            member_result_list.at(member_index).mdpde.GetModel().GetOffset()
        };
        const auto gaussian{
            linearization_service::DecodeParameterVector(
                result.beta_posterior_matrix.col(i),
                result.capital_sigma_posterior_list.at(member_index))
        };
        const auto gaussian_with_offset{
            GaussianModel3DWithUncertainty{
                gaussian.GetModel().WithOffset(offset),
                gaussian.GetStandardDeviationModel()
            }
        };
        member_results.emplace_back(LocalGaussianResult{
            0.0,
            gaussian_with_offset,
            gaussian_with_offset,
            gaussian_with_offset,
            static_cast<bool>(result.outlier_flag_array(i)),
            result.statistical_distance_array(i)
        });
    }
    return member_results;
}

using FittedGaussianSnapshot = std::vector<GaussianModel3D>;
using GroupMedianModelMap = std::unordered_map<GroupKey, GaussianModel3D>;

SecondStageLocalFittingContext BuildSecondStageLocalFittingContext(ModelObject & model_object)
{
    SecondStageLocalFittingContext context;
    context.atom_list = model_object.GetSelectedAtoms();
    context.atom_context_list.resize(context.atom_list.size());
    std::unordered_map<const AtomObject *, std::size_t> atom_index_map;
    atom_index_map.reserve(context.atom_list.size());
    for (std::size_t i = 0; i < context.atom_list.size(); i++)
    {
        atom_index_map.emplace(context.atom_list.at(i), i);
    }
    const auto analysis_view{ model_object.GetAnalysisView() };

    std::unordered_map<GroupKey, std::vector<double>> width_samples_by_group;
    width_samples_by_group.reserve(context.atom_list.size());
    for (const auto * atom : context.atom_list)
    {
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

    for (std::size_t atom_index = 0; atom_index < context.atom_list.size(); atom_index++)
    {
        const auto * atom{ context.atom_list.at(atom_index) };
        auto & atom_context{ context.atom_context_list.at(atom_index) };
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

    for (std::size_t atom_index = 0; atom_index < context.atom_list.size(); atom_index++)
    {
        const auto * atom{ context.atom_list.at(atom_index) };
        auto & atom_context{ context.atom_context_list.at(atom_index) };
        const auto neighbor_atom_list{ atom->FindNeighborAtoms(kNeighborAtomSearchRange) };
        atom_context.selected_neighbor_index_list.reserve(neighbor_atom_list.size());
        for (const auto * neighbor_atom : neighbor_atom_list)
        {
            const auto neighbor_iter{ atom_index_map.find(neighbor_atom) };
            if (neighbor_iter == atom_index_map.end()) continue;

            atom_context.selected_neighbor_index_list.emplace_back(neighbor_iter->second);
        }

        atom_context.sample_neighbor_list.resize(atom_context.sample_entries.size());
        for (std::size_t sample_index = 0; sample_index < atom_context.sample_entries.size(); sample_index++)
        {
            const auto & sample{ atom_context.sample_entries.at(sample_index) };
            auto & sample_neighbor_list{ atom_context.sample_neighbor_list.at(sample_index) };
            sample_neighbor_list.reserve(atom_context.selected_neighbor_index_list.size());
            for (const auto neighbor_index : atom_context.selected_neighbor_index_list)
            {
                const auto distance{
                    static_cast<double>(
                        array_helper::ComputeNorm<float>(
                            sample.point.position,
                            context.atom_list.at(neighbor_index)->GetPositionRef()))
                };
                if (distance > kNeighborContributionDistanceMax) continue;

                sample_neighbor_list.emplace_back(SecondStageNeighborSample{ neighbor_index, distance });
            }
        }
    }

    return context;
}

std::size_t FindLocalFittingClusterRoot(std::vector<std::size_t> & parent_list, std::size_t index)
{
    if (index >= parent_list.size())
    {
        throw std::invalid_argument("Local fitting cluster index is out of range.");
    }
    auto root{ index };
    while (parent_list.at(root) != root)
    {
        root = parent_list.at(root);
    }
    while (parent_list.at(index) != index)
    {
        const auto parent{ parent_list.at(index) };
        parent_list.at(index) = root;
        index = parent;
    }
    return root;
}

void MergeLocalFittingClusters(
    std::vector<std::size_t> & parent_list,
    std::size_t left_index,
    std::size_t right_index)
{
    const auto left_root{ FindLocalFittingClusterRoot(parent_list, left_index) };
    const auto right_root{ FindLocalFittingClusterRoot(parent_list, right_index) };
    if (left_root == right_root) return;

    parent_list.at(right_root) = left_root;
}

LocalFittingClusterMap BuildLocalFittingClusters(
    const SecondStageLocalFittingContext & context,
    const std::vector<std::size_t> & active_index_list)
{
    std::vector<int> active_position_by_atom_index(context.AtomSize(), -1);
    for (std::size_t active_position = 0; active_position < active_index_list.size(); active_position++)
    {
        const auto active_index{ active_index_list.at(active_position) };
        if (active_index >= context.AtomSize())
        {
            throw std::invalid_argument("Local fitting active index is out of range.");
        }
        active_position_by_atom_index.at(active_index) = static_cast<int>(active_position);
    }

    std::vector<std::size_t> parent_list(active_index_list.size());
    for (std::size_t i = 0; i < parent_list.size(); i++)
    {
        parent_list.at(i) = i;
    }

    for (std::size_t atom_index = 0; atom_index < context.AtomSize(); atom_index++)
    {
        const auto & atom_context{ context.atom_context_list.at(atom_index) };
        for (std::size_t sample_index = 0; sample_index < atom_context.sample_entries.size(); sample_index++)
        {
            std::vector<std::size_t> contributor_position_list;
            const auto target_active_position{ active_position_by_atom_index.at(atom_index) };
            if (target_active_position >= 0)
            {
                contributor_position_list.emplace_back(static_cast<std::size_t>(target_active_position));
            }
            for (const auto & neighbor_sample : atom_context.sample_neighbor_list.at(sample_index))
            {
                if (neighbor_sample.atom_index >= active_position_by_atom_index.size())
                {
                    throw std::invalid_argument("Local fitting neighbor index is out of range.");
                }
                const auto neighbor_active_position{
                    active_position_by_atom_index.at(neighbor_sample.atom_index)
                };
                if (neighbor_active_position >= 0)
                {
                    contributor_position_list.emplace_back(static_cast<std::size_t>(neighbor_active_position));
                }
            }

            if (contributor_position_list.size() < 2) continue;

            std::sort(contributor_position_list.begin(), contributor_position_list.end());
            contributor_position_list.erase(
                std::unique(contributor_position_list.begin(), contributor_position_list.end()),
                contributor_position_list.end());
            for (std::size_t i = 1; i < contributor_position_list.size(); i++)
            {
                MergeLocalFittingClusters(
                    parent_list,
                    contributor_position_list.front(),
                    contributor_position_list.at(i));
            }
        }
    }

    std::map<std::size_t, LocalFittingClusterKey> active_index_list_by_root;
    for (std::size_t active_position = 0; active_position < active_index_list.size(); active_position++)
    {
        const auto root{ FindLocalFittingClusterRoot(parent_list, active_position) };
        active_index_list_by_root[root].emplace_back(active_index_list.at(active_position));
    }

    std::map<std::size_t, std::vector<LocalFittingObjectiveSampleRef>> sample_ref_list_by_root;
    for (std::size_t atom_index = 0; atom_index < context.AtomSize(); atom_index++)
    {
        const auto & atom_context{ context.atom_context_list.at(atom_index) };
        for (std::size_t sample_index = 0; sample_index < atom_context.sample_entries.size(); sample_index++)
        {
            std::vector<std::size_t> contributor_root_list;
            const auto target_active_position{ active_position_by_atom_index.at(atom_index) };
            if (target_active_position >= 0)
            {
                contributor_root_list.emplace_back(
                    FindLocalFittingClusterRoot(
                        parent_list, static_cast<std::size_t>(target_active_position)));
            }
            for (const auto & neighbor_sample : atom_context.sample_neighbor_list.at(sample_index))
            {
                const auto neighbor_active_position{
                    active_position_by_atom_index.at(neighbor_sample.atom_index)
                };
                if (neighbor_active_position >= 0)
                {
                    contributor_root_list.emplace_back(
                        FindLocalFittingClusterRoot(
                            parent_list,
                            static_cast<std::size_t>(neighbor_active_position)));
                }
            }
            if (contributor_root_list.empty()) continue;

            std::sort(contributor_root_list.begin(), contributor_root_list.end());
            contributor_root_list.erase(
                std::unique(contributor_root_list.begin(), contributor_root_list.end()),
                contributor_root_list.end());
            if (contributor_root_list.size() != 1)
            {
                throw std::logic_error(
                    "Local fitting objective sample spans multiple clusters.");
            }
            sample_ref_list_by_root[contributor_root_list.front()].emplace_back(
                LocalFittingObjectiveSampleRef{ atom_index, sample_index });
        }
    }

    LocalFittingClusterMap cluster_map;
    for (auto & [root, cluster_active_index_list] : active_index_list_by_root)
    {
        std::sort(cluster_active_index_list.begin(), cluster_active_index_list.end());
        const auto inserted{
            cluster_map.emplace(
                std::move(cluster_active_index_list),
                LocalFittingClusterWork{ std::move(sample_ref_list_by_root[root]) })
        };
        if (!inserted.second)
        {
            throw std::logic_error("Local fitting cluster key is duplicated.");
        }
    }
    return cluster_map;
}

GaussianFittingState BuildInitialLocalFittingState(const SecondStageLocalFittingContext & context)
{
    GaussianFittingState state{
        std::vector<LocalGaussianResult>(context.AtomSize()),
        std::vector<Eigen::VectorXd>(context.AtomSize())
    };
    for (std::size_t i = 0; i < context.AtomSize(); i++)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*context.atom_list.at(i)) };
        state.result_list.at(i) = local_view.GetGaussianResult();
        state.estimation_list.at(i) = state.result_list.at(i).mdpde.GetModel().ToVector();
    }
    return state;
}

FittedGaussianSnapshot BuildFittedGaussianSnapshot(
    const SecondStageLocalFittingContext & context,
    const std::vector<Eigen::VectorXd> & estimation_list)
{
    if (context.AtomSize() != estimation_list.size())
    {
        throw std::invalid_argument("atom context and estimation_list sizes are inconsistent.");
    }

    FittedGaussianSnapshot snapshot;
    snapshot.reserve(context.AtomSize());
    for (const auto & estimation : estimation_list)
    {
        snapshot.emplace_back(GaussianModel3D::FromVector(estimation));
    }
    return snapshot;
}

GroupMedianModelMap BuildGroupMedianMDPDEModelMap(const std::vector<AtomObject *> & atom_list)
{
    std::unordered_map<GroupKey, GaussianModelParameterSamples> parameter_samples_by_group;
    parameter_samples_by_group.reserve(atom_list.size());
    for (const auto * atom : atom_list)
    {
        const auto local_view{ AtomLocalPotentialView::For(*atom) };
        if (!local_view.IsAvailable()) continue;

        const auto & model{ local_view.GetEstimateMDPDE() };
        auto & parameter_samples{
            parameter_samples_by_group[data_internal::GetGroupKey(atom)]
        };
        parameter_samples.amplitude_list.emplace_back(model.GetAmplitude());
        parameter_samples.width_list.emplace_back(model.GetWidth());
        parameter_samples.offset_list.emplace_back(model.GetOffset());
    }

    GroupMedianModelMap median_model_by_group;
    median_model_by_group.reserve(parameter_samples_by_group.size());
    for (const auto & [group_key, parameter_samples] : parameter_samples_by_group)
    {
        if (parameter_samples.amplitude_list.empty()) continue;

        median_model_by_group.emplace(
            group_key,
            GaussianModel3D{
                array_helper::ComputeMedian(parameter_samples.amplitude_list),
                array_helper::ComputeMedian(parameter_samples.width_list),
                array_helper::ComputeMedian(parameter_samples.offset_list)
            });
    }
    return median_model_by_group;
}

JointOffsetBuildResult BuildJointOffsetSystem(
    const SecondStageLocalFittingContext & context,
    const std::vector<std::size_t> & active_index_list,
    const FittedGaussianSnapshot & snapshot,
    double ridge_ratio,
    const std::vector<double> & ridge_multiplier_list)
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
        row_basis_entries.reserve(atom_context.selected_neighbor_index_list.size() + 1);
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
    if (!std::isfinite(updated_objective))
    {
        return true;
    }
    if (!std::isfinite(current_objective))
    {
        return false;
    }
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
    const std::vector<double> & ridge_multiplier_list)
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
        build_result = BuildJointOffsetSystem(context, active_index_list, snapshot, ridge_ratio, ridge_multiplier_list);
    }
    catch (const std::exception &)
    {
        return JointOffsetSolveResult{ previous_offset, {} };
    }
    auto system{ std::move(build_result.system) };
    auto active_coupling_graph{ std::move(build_result.active_coupling_graph) };
    if (system.response.size() == 0 || system.previous_parameter.size() == 0)
    {
        return JointOffsetSolveResult{
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
                system.previous_parameter,
                std::move(active_coupling_graph)
            };
        }
        const auto current_objective{
            CalculateWeightedRidgeSurrogateObjective(system, weight, offset)
        };
        const auto updated_objective{
            CalculateWeightedRidgeSurrogateObjective(system, weight, updated_offset)
        };
        if (IsJointOffsetObjectiveDeteriorated(updated_objective, current_objective)) break;
        const auto maximum_change{
            algorithm::CalculateMaximumNormalizedVectorChange(updated_offset, offset, kJointOffsetIrlsScaleFloor)
        };
        offset = std::move(updated_offset);
        if (maximum_change < kJointOffsetIrlsNormalizedChangeTolerance) break;
    }

    return JointOffsetSolveResult{
        offset,
        std::move(active_coupling_graph)
    };
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

template <typename GaussianLookup>
LocalPotentialSampleList UpdateSampleListWithGaussianLookup(
    const AtomObject & atom,
    GaussianLookup lookup_gaussian)
{
    const auto local_view{ AtomLocalPotentialView::RequireFor(atom) };
    const auto sample_entries{ local_view.GetSamplingEntries(false) };
    const auto & neighbor_atom_list{ atom.FindNeighborAtoms(kNeighborAtomSearchRange) };
    LocalPotentialSampleList updated_list;
    updated_list.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        auto sample_position{ sample.point.position };
        auto response_value{ sample.response };
        for (const auto * neighbor_atom : neighbor_atom_list)
        {
            const auto * gaussian{ lookup_gaussian(*neighbor_atom) };
            if (gaussian == nullptr) continue;

            auto neighbor_position{ neighbor_atom->GetPosition() };
            auto distance{
                static_cast<double>(array_helper::ComputeNorm<float>(sample_position, neighbor_position))
            };
            if (distance > kNeighborContributionDistanceMax) continue;
            response_value -= static_cast<float>(gaussian->ResponseAtDistance(distance));
        }
        updated_list.emplace_back(LocalPotentialSample{response_value, sample.point });
    }
    return updated_list;
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

LocalPotentialSampleList UpdateSampleListWithGroupMedianGaussian(
    const AtomObject & atom,
    const GroupMedianModelMap & median_model_by_group)
{
    return UpdateSampleListWithGaussianLookup(
        atom,
        [&median_model_by_group](const AtomObject & neighbor_atom) -> const GaussianModel3D *
        {
            const auto median_model_iter{
                median_model_by_group.find(data_internal::GetGroupKey(&neighbor_atom))
            };
            if (median_model_iter != median_model_by_group.end())
            {
                return &median_model_iter->second;
            }

            const auto local_view{ AtomLocalPotentialView::For(neighbor_atom) };
            return local_view.IsAvailable() ? &local_view.GetEstimateMDPDE() : nullptr;
        });
}

LocalPotentialSampleList UpdateSampleListWithFittedGroupGaussian(
    const AtomObject & atom,
    const ModelAnalysisView & analysis_view)
{
    return UpdateSampleListWithGaussianLookup(
        atom,
        [&analysis_view](const AtomObject & neighbor_atom) -> const GaussianModel3D *
        {
            const auto group_key{ data_internal::GetGroupKey(&neighbor_atom) };
            if (!analysis_view.HasAtomGroup(group_key))
            {
                return nullptr;
            }
            return &analysis_view.GetAtomGroupPrior(group_key);
        });
}

void SetUpdatedSamplingEntriesFromGroupMedianGaussian(ModelObject & model_object)
{
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    auto local_editor_list{ BuildAtomLocalEditors(model_object, atom_list) };
    const auto median_model_by_group{ BuildGroupMedianMDPDEModelMap(atom_list) };
    for (size_t i = 0; i < atom_list.size(); i++)
    {
        local_editor_list[i].SetUpdatedSamplingEntries(
            UpdateSampleListWithGroupMedianGaussian(*atom_list[i], median_model_by_group));
    }
}

void SetUpdatedSamplingEntriesFromFittedGroupGaussian(ModelObject & model_object)
{
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    auto local_editor_list{ BuildAtomLocalEditors(model_object, atom_list) };
    const auto analysis_view{ model_object.GetAnalysisView() };
    for (size_t i = 0; i < atom_list.size(); i++)
    {
        local_editor_list[i].SetUpdatedSamplingEntries(
            UpdateSampleListWithFittedGroupGaussian(*atom_list[i], analysis_view));
    }
}

std::optional<LocalFittingObjectiveSamples> CollectLocalFittingObjectiveSamples(
    const SecondStageLocalFittingContext & context,
    const GaussianFittingState & fitting_state,
    const std::vector<LocalFittingObjectiveSampleRef> & sample_ref_list,
    const std::vector<std::size_t> & active_index_list)
{
    const auto snapshot{
        BuildFittedGaussianSnapshot(context, fitting_state.estimation_list)
    };

    LocalFittingObjectiveSamples objective_samples;
    objective_samples.residual_list.reserve(sample_ref_list.size());
    objective_samples.response_list.reserve(sample_ref_list.size());
    objective_samples.active_model_list.reserve(active_index_list.size());
    for (const auto active_index : active_index_list)
    {
        if (active_index >= fitting_state.estimation_list.size())
        {
            throw std::invalid_argument("Local fitting objective active index is out of range.");
        }
        objective_samples.active_model_list.emplace_back(LocalFittingObjectiveAtomModel{
            active_index,
            GaussianModel3D::FromVector(fitting_state.estimation_list.at(active_index))
        });
    }
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
        objective_samples.response_list.emplace_back(response);
    }
    return objective_samples;
}

std::optional<double> CalculateLocalFittingResidualScaleSample(
    const LocalFittingObjectiveSamples & objective_samples)
{
    if (objective_samples.residual_list.empty() ||
        objective_samples.residual_list.size() != objective_samples.response_list.size())
    {
        return std::nullopt;
    }

    const auto residual_scale{
        CalculateMedianAbsoluteDeviationScale(objective_samples.residual_list)
    };
    const auto response_scale_floor{
        kLocalFittingObjectiveResidualScaleFloorRatio *
        CalculateMedianAbsoluteDeviationScale(objective_samples.response_list)
    };
    const auto scale_sample{
        std::max({ residual_scale, response_scale_floor, kRobustScaleMin })
    };
    if (!std::isfinite(scale_sample)) return std::nullopt;
    return scale_sample;
}

double CalculateSquaredValue(double value)
{
    return value * value;
}

bool IsLocalFittingObjectiveModelValid(const GaussianModel3D & model)
{
    return numeric_validation::IsFinitePositive(model.GetAmplitude()) &&
        numeric_validation::IsFinitePositive(model.GetWidth()) &&
        std::isfinite(model.GetOffset());
}

std::optional<double> CalculateLocalFittingParameterPenalty(
    const SecondStageLocalFittingContext & context,
    const LocalFittingObjectiveSamples & objective_samples,
    const algorithm::ScaleReference & objective_reference,
    const GaussianFittingState * movement_reference_state)
{
    if (!objective_reference.has_reference || !std::isfinite(objective_reference.scale))
    {
        return std::nullopt;
    }
    if (movement_reference_state != nullptr &&
        movement_reference_state->estimation_list.size() != context.AtomSize())
    {
        return std::nullopt;
    }

    double width_prior_penalty{ 0.0 };
    double offset_plausibility_penalty{ 0.0 };
    double movement_penalty{ 0.0 };
    const auto residual_scale_floor{ std::max(objective_reference.scale, kRobustScaleMin) };
    for (const auto & atom_model : objective_samples.active_model_list)
    {
        const auto active_index{ atom_model.atom_index };
        if (active_index >= context.AtomSize()) return std::nullopt;

        const auto & model{ atom_model.model };
        if (!IsLocalFittingObjectiveModelValid(model)) return std::nullopt;

        const auto prior_width{ context.atom_context_list.at(active_index).prior_width };
        if (!numeric_validation::IsFinitePositive(prior_width)) return std::nullopt;
        width_prior_penalty += CalculateSquaredValue(
            (std::log(model.GetWidth()) - std::log(prior_width)) / kLocalFittingWidthPriorLogScale);

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
        offset_plausibility_penalty += CalculateSquaredValue(offset_excess);

        if (movement_reference_state != nullptr)
        {
            const auto previous_model{
                GaussianModel3D::FromVector(
                    movement_reference_state->estimation_list.at(active_index))
            };
            if (!IsLocalFittingObjectiveModelValid(previous_model)) return std::nullopt;

            const auto previous_offset_basis{ previous_model.OffsetBasisAtDistance(0.0) };
            const auto previous_peak_signal{ previous_model.SignalAtDistance(0.0) };
            if (!std::isfinite(previous_offset_basis) ||
                previous_offset_basis == 0.0 ||
                !std::isfinite(previous_peak_signal))
            {
                return std::nullopt;
            }
            const auto offset_scale{
                std::max({
                    std::abs(previous_model.GetOffset()),
                    std::abs(previous_peak_signal / previous_offset_basis),
                    residual_scale_floor,
                    kRobustScaleMin
                })
            };
            movement_penalty += CalculateSquaredValue(
                    (std::log(model.GetAmplitude()) - std::log(previous_model.GetAmplitude())) /
                    kLocalFittingMovementAmplitudeLogScale) +
                CalculateSquaredValue(
                    (std::log(model.GetWidth()) - std::log(previous_model.GetWidth())) /
                    kLocalFittingMovementWidthLogScale) +
                CalculateSquaredValue(
                    (model.GetOffset() - previous_model.GetOffset()) / offset_scale);
        }
    }

    const auto penalty{
        kLocalFittingWidthPriorPenaltyWeight * width_prior_penalty +
        kLocalFittingOffsetPlausibilityPenaltyWeight * offset_plausibility_penalty +
        kLocalFittingMovementPenaltyWeight * movement_penalty
    };
    return std::isfinite(penalty) ? std::optional<double>{ penalty } : std::nullopt;
}

std::optional<double> CalculateLocalFittingObjective(
    const SecondStageLocalFittingContext & context,
    const LocalFittingObjectiveSamples & objective_samples,
    const algorithm::ScaleReference & objective_reference,
    const GaussianFittingState * movement_reference_state)
{
    if (!objective_reference.has_reference ||
        !numeric_validation::IsFinitePositive(objective_reference.scale) ||
        objective_samples.residual_list.empty())
    {
        return std::nullopt;
    }

    double loss_sum{ 0.0 };
    for (const auto residual : objective_samples.residual_list)
    {
        const auto normalized_residual{ residual / objective_reference.scale };
        loss_sum += algorithm::CalculateRobustLoss(
            kSecondStageRobustLossKind,
            normalized_residual,
            kRobustLossCutoffMultiplier);
    }
    const auto residual_objective{
        loss_sum / static_cast<double>(objective_samples.residual_list.size())
    };
    if (!std::isfinite(residual_objective)) return std::nullopt;

    const auto parameter_penalty{
        CalculateLocalFittingParameterPenalty(
            context,
            objective_samples,
            objective_reference,
            movement_reference_state)
    };
    if (!parameter_penalty.has_value()) return std::nullopt;

    const auto objective{ residual_objective + *parameter_penalty };
    return std::isfinite(objective) ?
        std::optional<double>{ objective } :
        std::nullopt;
}

algorithm::ClusteredFittingQualityCandidateScore<LocalFittingObjectiveSamples>
ScoreLocalFittingClusterCandidate(
    const SecondStageLocalFittingContext & context,
    const GaussianFittingState & candidate_state,
    const GaussianFittingState & previous_state,
    const std::vector<std::size_t> & active_index_list,
    const std::vector<LocalFittingObjectiveSampleRef> & sample_ref_list,
    const algorithm::ScaleReferenceTracker & objective_scale_tracker,
    algorithm::FittingQualityCandidateStats & previous_candidate_stats,
    const std::optional<LocalFittingObjectiveSamples> & previous_objective_samples,
    const std::optional<algorithm::ClusteredFittingQualityTrackedCandidate<LocalFittingObjectiveSamples>> & best_candidate,
    const algorithm::ParameterChangeStats & normalized_change_stats)
{
    algorithm::ClusteredFittingQualityCandidateScore<LocalFittingObjectiveSamples> score;
    auto objective_reference{ objective_scale_tracker.GetCommittedReference() };
    score.has_objective_reference = objective_reference.has_reference;
    if (best_candidate.has_value())
    {
        score.best_candidate_stats = best_candidate->candidate_stats;
    }

    if (objective_scale_tracker.HasReference())
    {
        score.objective_samples =
            CollectLocalFittingObjectiveSamples(
                context,
                candidate_state,
                sample_ref_list,
                active_index_list);
        if (score.objective_samples.has_value())
        {
            const auto scale_sample{
                CalculateLocalFittingResidualScaleSample(*score.objective_samples)
            };
            if (scale_sample.has_value())
            {
                objective_reference = objective_scale_tracker.GetProvisionalReference(*scale_sample);
                score.has_objective_reference = objective_reference.has_reference;
                const auto candidate_objective{
                    CalculateLocalFittingObjective(
                        context,
                        *score.objective_samples,
                        objective_reference,
                        &previous_state)
                };
                const auto commit_objective{
                    CalculateLocalFittingObjective(
                        context,
                        *score.objective_samples,
                        objective_reference,
                        nullptr)
                };
                score.objective_scale_sample = *scale_sample;
                if (candidate_objective.has_value())
                {
                    if (previous_objective_samples.has_value())
                    {
                        const auto recalculated_previous_objective{
                            CalculateLocalFittingObjective(
                                context,
                                *previous_objective_samples,
                                objective_reference,
                                nullptr)
                        };
                        previous_candidate_stats.has_quality_objective =
                            recalculated_previous_objective.has_value();
                        previous_candidate_stats.quality_objective =
                            recalculated_previous_objective.value_or(
                                std::numeric_limits<double>::infinity());
                    }
                    if (best_candidate.has_value() &&
                        best_candidate->objective_samples.has_value() &&
                        score.best_candidate_stats.has_value())
                    {
                        const auto recalculated_best_objective{
                            CalculateLocalFittingObjective(
                                context,
                                *best_candidate->objective_samples,
                                objective_reference,
                                nullptr)
                        };
                        score.best_candidate_stats->has_quality_objective =
                            recalculated_best_objective.has_value();
                        score.best_candidate_stats->quality_objective =
                            recalculated_best_objective.value_or(
                                std::numeric_limits<double>::infinity());
                    }
                }
                score.candidate_stats.has_quality_objective = candidate_objective.has_value();
                score.candidate_stats.quality_objective = candidate_objective.value_or(
                    std::numeric_limits<double>::infinity());
                score.commit_candidate_stats = algorithm::FittingQualityCandidateStats{
                    commit_objective.has_value(),
                    commit_objective.value_or(std::numeric_limits<double>::infinity()),
                    normalized_change_stats
                };
            }
        }
    }

    score.candidate_stats.parameter_change_stats = normalized_change_stats;
    return score;
}

ParameterSummaryStats SummarizeParameterValues(const std::vector<double> & value_list)
{
    if (value_list.empty()) return {};

    double sum{ 0.0 };
    for (const auto value : value_list)
    {
        sum += value;
    }
    const auto mean{ sum / static_cast<double>(value_list.size()) };
    if (value_list.size() < 2)
    {
        return ParameterSummaryStats{ mean, 0.0 };
    }

    double squared_error_sum{ 0.0 };
    for (const auto value : value_list)
    {
        const auto error{ value - mean };
        squared_error_sum += error * error;
    }
    return ParameterSummaryStats{
        mean,
        std::sqrt(squared_error_sum / static_cast<double>(value_list.size() - 1))
    };
}

std::vector<std::string> BuildGroupPriorSpotSummaryLines(const ModelObject & model_object)
{
    const auto analysis_view{ model_object.GetAnalysisView() };
    std::map<std::string, GaussianModelParameterSamples> spot_sample_map;
    for (const auto group_key : analysis_view.CollectAtomGroupKeys())
    {
        const auto & atom_list{ analysis_view.GetAtomObjectList(group_key) };
        if (atom_list.empty()) continue;

        const auto & prior{ analysis_view.GetAtomGroupPrior(group_key) };
        auto & sample_list{
            spot_sample_map["Spot::" + ChemicalDataHelper::GetLabel(atom_list.front()->GetSpot())]
        };
        sample_list.amplitude_list.emplace_back(prior.GetAmplitude());
        sample_list.width_list.emplace_back(prior.GetWidth());
        sample_list.offset_list.emplace_back(prior.GetOffset());
    }

    std::vector<std::string> summary_lines;
    summary_lines.reserve(spot_sample_map.size());
    for (const auto & [spot_label, sample_list] : spot_sample_map)
    {
        if (spot_label != "Spot::CA" && spot_label != "Spot::C" && spot_label != "Spot::N" && spot_label != "Spot::O")
        {
            continue;
        }
        const auto amplitude_stats{ SummarizeParameterValues(sample_list.amplitude_list) };
        const auto width_stats{ SummarizeParameterValues(sample_list.width_list) };
        const auto offset_stats{ SummarizeParameterValues(sample_list.offset_list) };

        std::ostringstream stream;
        stream << spot_label << std::fixed << std::setprecision(2)
            << " , amplitude mean = " << amplitude_stats.mean
            << ", amplitude s.d. = " << amplitude_stats.standard_deviation
            << ", width mean = " << width_stats.mean
            << ", width s.d. = " << width_stats.standard_deviation
            << ", offset mean = " << offset_stats.mean
            << ", offset s.d. = " << offset_stats.standard_deviation;
        summary_lines.emplace_back(stream.str());
    }
    return summary_lines;
}

void LogGroupPriorSpotSummary(const ModelObject & model_object)
{
    const auto summary_lines{ BuildGroupPriorSpotSummaryLines(model_object) };
    if (summary_lines.empty())
    {
        Logger::Log(LogLevel::Info, "Group fitting prior summary by Spot: no atom groups available.");
        return;
    }

    Logger::Log(LogLevel::Info, "Group fitting prior summary by Spot:");
    for (const auto & line : summary_lines)
    {
        Logger::Log(LogLevel::Info, line);
    }
}

std::vector<algorithm::ParameterChange> CalculateLocalFittingParameterChanges(
    const std::vector<Eigen::VectorXd> & current_estimation_list,
    const std::vector<Eigen::VectorXd> & previous_estimation_list)
{
    if (current_estimation_list.size() != previous_estimation_list.size())
    {
        throw std::invalid_argument("Local fitting parameter change input sizes are inconsistent.");
    }

    std::vector<algorithm::ParameterChange> change_list(current_estimation_list.size());
    for (size_t i = 0; i < current_estimation_list.size(); i++)
    {
        const auto & current_estimation{ current_estimation_list.at(i) };
        const auto & previous_estimation{ previous_estimation_list.at(i) };
        change_list.at(i).value_list = {
            std::abs(
                current_estimation(GaussianModel3D::AmplitudeIndex()) -
                previous_estimation(GaussianModel3D::AmplitudeIndex())),
            std::abs(
                current_estimation(GaussianModel3D::WidthIndex()) -
                previous_estimation(GaussianModel3D::WidthIndex())),
            std::abs(
                current_estimation(GaussianModel3D::OffsetIndex()) -
                previous_estimation(GaussianModel3D::OffsetIndex()))
        };
    }
    return change_list;
}

std::vector<algorithm::ParameterChange> CalculateLocalFittingNormalizedParameterChanges(
    const std::vector<Eigen::VectorXd> & current_estimation_list,
    const std::vector<Eigen::VectorXd> & previous_estimation_list)
{
    if (current_estimation_list.size() != previous_estimation_list.size())
    {
        throw std::invalid_argument("Local fitting normalized parameter change input sizes are inconsistent.");
    }

    std::vector<algorithm::ParameterChange> change_list(current_estimation_list.size());
    for (size_t i = 0; i < current_estimation_list.size(); i++)
    {
        const auto & current_estimation{ current_estimation_list.at(i) };
        const auto & previous_estimation{ previous_estimation_list.at(i) };
        change_list.at(i).value_list = {
            algorithm::CalculateNormalizedChange(
                current_estimation(GaussianModel3D::AmplitudeIndex()),
                previous_estimation(GaussianModel3D::AmplitudeIndex()),
                kLocalFittingNormalizedChangeScaleFloor),
            algorithm::CalculateNormalizedChange(
                current_estimation(GaussianModel3D::WidthIndex()),
                previous_estimation(GaussianModel3D::WidthIndex()),
                kLocalFittingNormalizedChangeScaleFloor),
            algorithm::CalculateNormalizedChange(
                current_estimation(GaussianModel3D::OffsetIndex()),
                previous_estimation(GaussianModel3D::OffsetIndex()),
                kLocalFittingNormalizedChangeScaleFloor)
        };
    }
    return change_list;
}

bool IsLocalFittingNormalizedParameterChangeConverged(const algorithm::ParameterChangeStats & stats)
{
    for (std::size_t i = 0; i < stats.percentile_list.size(); i++)
    {
        if (stats.percentile_list.at(i) >= kLocalFittingNormalizedChangeTolerance) return false;
    }
    return true;
}

bool TryApplyLocalFittingCandidate(
    GaussianFittingState & current_state,
    const GaussianFittingState & previous_state,
    const std::vector<Eigen::VectorXd> & candidate_estimation_list,
    const std::vector<LocalGaussianResult> & uncertainty_result_list,
    const std::vector<std::size_t> & active_index_list,
    double damping)
{
    if (!std::isfinite(damping) || damping <= 0.0 || damping > 1.0)
    {
        throw std::invalid_argument("Local fitting candidate damping is out of range.");
    }
    if (current_state.estimation_list.size() != previous_state.estimation_list.size() ||
        current_state.result_list.size() != previous_state.result_list.size() ||
        candidate_estimation_list.size() != previous_state.estimation_list.size() ||
        uncertainty_result_list.size() != previous_state.result_list.size())
    {
        throw std::invalid_argument("Local fitting candidate input sizes are inconsistent.");
    }

    std::vector<std::pair<std::size_t, Eigen::VectorXd>> damped_estimation_list;
    damped_estimation_list.reserve(active_index_list.size());
    for (const auto active_index : active_index_list)
    {
        if (active_index >= current_state.estimation_list.size())
        {
            throw std::invalid_argument("Local fitting candidate active index is out of range.");
        }
        const auto damped_estimation{
            (previous_state.estimation_list.at(active_index) +
                damping * (
                    candidate_estimation_list.at(active_index) -
                    previous_state.estimation_list.at(active_index))).eval()
        };
        if (damped_estimation.size() != GaussianModel3D::ParameterSize() ||
            !damped_estimation.allFinite() ||
            damped_estimation(GaussianModel3D::WidthIndex()) <= 0.0)
        {
            return false;
        }
        damped_estimation_list.emplace_back(active_index, damped_estimation);
    }

    for (const auto & [active_index, damped_estimation] : damped_estimation_list)
    {
        const auto damped_model{ GaussianModel3D::FromVector(damped_estimation) };
        auto & result{ current_state.result_list.at(active_index) };
        result.mdpde = GaussianModel3DWithUncertainty{
            damped_model,
            uncertainty_result_list.at(active_index).mdpde.GetStandardDeviationModel()
        };
        current_state.estimation_list.at(active_index) = damped_estimation;
    }
    return true;
}

const char * GetLocalFittingAccelerationText(const LocalFittingAccelerationAttempt & attempt)
{
    if (attempt.kind == LocalFittingCandidateKind::FixedPoint) return "damped-fixed-point";
    return attempt.damping == 1.0 ? "aa" : "damped-aa";
}

std::vector<LocalFittingClusterKey> BuildLocalFittingClusterKeyList(
    const LocalFittingClusterMap & cluster_map)
{
    std::vector<LocalFittingClusterKey> key_list;
    key_list.reserve(cluster_map.size());
    for (const auto & [key, cluster] : cluster_map)
    {
        static_cast<void>(cluster);
        key_list.emplace_back(key);
    }
    return key_list;
}

algorithm::ClusteredFittingQualityInitialState<LocalFittingObjectiveSamples>
BuildInitialLocalFittingClusterQualityState(
    const SecondStageLocalFittingContext & context,
    const GaussianFittingState & previous_state,
    const LocalFittingClusterKey & key,
    const LocalFittingClusterWork & cluster)
{
    auto initial_objective_samples{
        CollectLocalFittingObjectiveSamples(
            context,
            previous_state,
            cluster.objective_sample_ref_list,
            key)
    };
    std::optional<double> initial_objective_scale_sample;
    if (initial_objective_samples.has_value())
    {
        initial_objective_scale_sample =
            CalculateLocalFittingResidualScaleSample(*initial_objective_samples);
    }
    algorithm::ScaleReferenceTracker initial_tracker{
        kLocalFittingObjectiveScaleWarmupCount,
        initial_objective_scale_sample
    };
    std::optional<double> initial_objective;
    if (initial_objective_samples.has_value() &&
        initial_objective_scale_sample.has_value())
    {
        initial_objective = CalculateLocalFittingObjective(
            context,
            *initial_objective_samples,
            initial_tracker.GetCommittedReference(),
            nullptr);
    }
    algorithm::FittingQualityCandidateStats initial_candidate_stats{
        initial_objective.has_value(),
        initial_objective.value_or(std::numeric_limits<double>::infinity()),
        algorithm::ParameterChangeStats{
            std::vector<double>(
                static_cast<std::size_t>(GaussianModel3D::ParameterSize()),
                0.0)
        }
    };
    return algorithm::ClusteredFittingQualityInitialState<LocalFittingObjectiveSamples>{
        initial_objective_scale_sample,
        std::move(initial_candidate_stats),
        std::move(initial_objective_samples)
    };
}

LocalFittingCandidateSelection SelectLocalFittingClusterCandidates(
    const SecondStageLocalFittingContext & context,
    LocalFittingClusterMap cluster_map,
    const std::vector<LocalFittingClusterKey> & key_list,
    const GaussianFittingState & previous_state,
    const GaussianFittingState & raw_state,
    algorithm::ClusteredAndersonAccelerationHistorySet & acceleration_history,
    algorithm::ClusteredFittingQualityStateSet<LocalFittingObjectiveSamples> & cluster_quality_state)
{
    const auto localized_anderson_candidate{
        acceleration_history.BuildCandidate(
            key_list,
            previous_state.estimation_list,
            raw_state.estimation_list)
    };
    LocalFittingCandidateSelection selection;
    selection.assembled_state = previous_state;
    bool has_accepted_acceleration_attempt{ false };

    const auto copy_cluster_state = [](
        GaussianFittingState & target_state,
        const GaussianFittingState & source_state,
        const LocalFittingClusterKey & key)
    {
        for (const auto active_index : key)
        {
            target_state.result_list.at(active_index) = source_state.result_list.at(active_index);
            target_state.estimation_list.at(active_index) = source_state.estimation_list.at(active_index);
        }
    };

    const auto run_attempt_group = [&](LocalFittingCandidateKind candidate_kind)
    {
        const auto attempt_size{ static_cast<int>(kLocalFittingAccelerationDampingList.size()) };
        for (int attempt = 0; attempt < attempt_size; attempt++)
        {
            const auto acceleration_attempt{
                LocalFittingAccelerationAttempt{
                    candidate_kind,
                    kLocalFittingAccelerationDampingList.at(static_cast<std::size_t>(attempt))
                }
            };
            for (auto & entry : cluster_map)
            {
                const auto & key{ entry.first };
                auto & cluster{ entry.second };
                if (cluster.attempt_state != LocalFittingClusterAttemptState::Pending) continue;
                if (candidate_kind == LocalFittingCandidateKind::Anderson &&
                    (!localized_anderson_candidate.has_value() ||
                        std::find(
                            localized_anderson_candidate->used_cluster_key_list.begin(),
                            localized_anderson_candidate->used_cluster_key_list.end(),
                            key) == localized_anderson_candidate->used_cluster_key_list.end()))
                {
                    continue;
                }

                auto attempt_state{ previous_state };
                const auto & candidate_estimation_list{
                    candidate_kind == LocalFittingCandidateKind::Anderson ?
                        localized_anderson_candidate->state_list :
                        raw_state.estimation_list
                };
                const auto & uncertainty_result_list{
                    candidate_kind == LocalFittingCandidateKind::Anderson ?
                        previous_state.result_list :
                        raw_state.result_list
                };
                const auto valid_attempt{
                    TryApplyLocalFittingCandidate(
                        attempt_state,
                        previous_state,
                        candidate_estimation_list,
                        uncertainty_result_list,
                        key,
                        acceleration_attempt.damping)
                };
                if (!valid_attempt)
                {
                    if (candidate_kind == LocalFittingCandidateKind::Anderson)
                    {
                        cluster.rejected_anderson = true;
                    }
                    continue;
                }

                const auto normalized_change_list{
                    CalculateLocalFittingNormalizedParameterChanges(
                        attempt_state.estimation_list,
                        previous_state.estimation_list)
                };
                auto normalized_change_stats{
                    algorithm::SummarizeParameterChangeStats(
                        normalized_change_list,
                        key,
                        kLocalFittingChangePercentile)
                };
                auto evaluated_attempt{
                    cluster_quality_state.EvaluateCandidate(
                        key,
                        attempt,
                        attempt_size,
                        [&](const algorithm::ScaleReferenceTracker & objective_scale_tracker,
                            algorithm::FittingQualityCandidateStats & previous_candidate_stats,
                            const std::optional<LocalFittingObjectiveSamples> & previous_objective_samples,
                            const std::optional<
                                algorithm::ClusteredFittingQualityTrackedCandidate<LocalFittingObjectiveSamples>> & best_candidate)
                        {
                            return ScoreLocalFittingClusterCandidate(
                                context,
                                attempt_state,
                                previous_state,
                                key,
                                cluster.objective_sample_ref_list,
                                objective_scale_tracker,
                                previous_candidate_stats,
                                previous_objective_samples,
                                best_candidate,
                                normalized_change_stats);
                        })
                };
                if (evaluated_attempt.outcome == algorithm::ClusteredFittingQualityAttemptOutcome::Accepted)
                {
                    copy_cluster_state(selection.assembled_state, attempt_state, key);
                    cluster.attempt_state =
                        candidate_kind == LocalFittingCandidateKind::Anderson ?
                            LocalFittingClusterAttemptState::AcceptedAnderson :
                            LocalFittingClusterAttemptState::AcceptedFixedPoint;
                    selection.accepted_score_list.emplace_back(
                        std::move(evaluated_attempt.accepted_score));
                    if (!has_accepted_acceleration_attempt)
                    {
                        selection.accepted_acceleration_attempt = acceleration_attempt;
                        has_accepted_acceleration_attempt = true;
                    }
                    continue;
                }

                selection.has_objective_backtracking_rejection = true;
                if (candidate_kind == LocalFittingCandidateKind::Anderson)
                {
                    cluster.rejected_anderson = true;
                }
                if (evaluated_attempt.outcome == algorithm::ClusteredFittingQualityAttemptOutcome::ObjectiveStop)
                {
                    cluster.attempt_state = LocalFittingClusterAttemptState::Stopped;
                }
            }
        }
    };

    if (localized_anderson_candidate.has_value())
    {
        run_attempt_group(LocalFittingCandidateKind::Anderson);
    }
    std::vector<LocalFittingClusterKey> anderson_failure_keys;
    for (auto & [key, cluster] : cluster_map)
    {
        if (cluster.rejected_anderson)
        {
            anderson_failure_keys.emplace_back(key);
        }
        if (cluster.attempt_state != LocalFittingClusterAttemptState::AcceptedAnderson)
        {
            cluster.attempt_state = LocalFittingClusterAttemptState::Pending;
        }
    }
    if (!anderson_failure_keys.empty())
    {
        acceleration_history.ClearAndSuppress(anderson_failure_keys);
    }
    run_attempt_group(LocalFittingCandidateKind::FixedPoint);

    std::vector<LocalFittingClusterKey> fixed_point_progress_keys;
    for (const auto & [key, cluster] : cluster_map)
    {
        if (cluster.attempt_state == LocalFittingClusterAttemptState::AcceptedAnderson ||
            cluster.attempt_state == LocalFittingClusterAttemptState::AcceptedFixedPoint)
        {
            selection.accepted_key_list.emplace_back(key);
            if (cluster.attempt_state == LocalFittingClusterAttemptState::AcceptedFixedPoint)
            {
                fixed_point_progress_keys.emplace_back(key);
            }
        }
        else
        {
            selection.rejected_key_list.emplace_back(key);
        }
    }
    acceleration_history.ReleaseSuppression(fixed_point_progress_keys);
    return selection;
}

std::optional<LocalGaussianResult> FitAtomWithJointOffsetFallback(
    const SecondStageLocalFittingContext & context,
    std::size_t atom_index,
    const LocalGaussianResult & previous_result,
    const FittedGaussianSnapshot & offset_snapshot,
    const FitOptions & options)
{
    auto sample_entries{ BuildSecondStageAdjustedSamples(context, atom_index, offset_snapshot) };
    const auto & offset_model{ offset_snapshot.at(atom_index) };
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
        if (
            CanBuildFiniteZeroOffsetSamples(sample_entries, candidate_model) &&
            !IsSuspiciousJointOffset(
                sample_entries,
                previous_result.mdpde.GetModel(),
                candidate_model,
                options))
        {
            return candidate_result;
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
    if (
        !CanBuildFiniteZeroOffsetSamples(sample_entries, result.mdpde.GetModel()) ||
        IsSuspiciousJointOffset(
            sample_entries,
            previous_result.mdpde.GetModel(),
            result.mdpde.GetModel(),
            options))
    {
        return std::nullopt;
    }
    return result;
}

void ThawChangedActiveAtomNeighbors(
    const SecondStageLocalFittingContext & context,
    const std::vector<algorithm::ParameterChange> & change_list,
    const std::vector<std::size_t> & active_index_list,
    algorithm::ConvergenceFreezeTracker & freeze_tracker,
    algorithm::DependencyThawHysteresisTracker & thaw_hysteresis_tracker)
{
    if (change_list.size() != context.AtomSize())
    {
        throw std::invalid_argument("Local fitting dependency thaw input size is inconsistent.");
    }

    const double thaw_threshold{ std::sqrt(kLocalFittingParameterChangeTolerance) };
    for (const auto active_index : active_index_list)
    {
        if (active_index >= context.AtomSize())
        {
            throw std::invalid_argument("Local fitting dependency thaw active index is out of range.");
        }
        const auto active_change{ algorithm::GetMaximumParameterChange(change_list.at(active_index)) };
        for (const auto neighbor_index : context.atom_context_list.at(active_index).selected_neighbor_index_list)
        {
            if (!freeze_tracker.IsFrozen(neighbor_index)) continue;
            if (!thaw_hysteresis_tracker.ShouldThaw(neighbor_index, active_change, thaw_threshold)) continue;
            if (!thaw_hysteresis_tracker.CanDependencyThaw(neighbor_index)) continue;
            if (freeze_tracker.Thaw(neighbor_index))
            {
                thaw_hysteresis_tracker.RecordDependencyThaw(neighbor_index);
            }
        }
    }
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
    const GaussianFittingState & previous_state,
    const std::vector<std::size_t> & suspicious_active_position_list,
    FittedGaussianSnapshot & current_snapshot,
    GaussianFittingState & iteration_state)
{
    for (const auto active_position : suspicious_active_position_list)
    {
        if (active_position >= active_index_list.size())
        {
            throw std::invalid_argument("Suspicious offset rollback position is out of range.");
        }

        const auto state_index{ active_index_list.at(active_position) };
        if (state_index >= context.AtomSize() ||
            state_index >= previous_state.estimation_list.size() ||
            state_index >= previous_state.result_list.size() ||
            state_index >= iteration_state.estimation_list.size() ||
            state_index >= iteration_state.result_list.size())
        {
            throw std::invalid_argument("Suspicious offset rollback state index is out of range.");
        }
        const auto previous_model{
            GaussianModel3D::FromVector(previous_state.estimation_list.at(state_index))
        };
        current_snapshot.at(state_index) = previous_model;
        iteration_state.estimation_list.at(state_index) = previous_state.estimation_list.at(state_index);
        iteration_state.result_list.at(state_index) = previous_state.result_list.at(state_index);
    }
}

LocalFittingIterationResult RunLocalFittingIteration(
    const SecondStageLocalFittingContext & context,
    const std::vector<std::size_t> & active_index_list,
    const GaussianFittingState & previous_state,
    const FitOptions & options,
    double ridge_ratio,
    const std::vector<double> & ridge_multiplier_list)
{
    const auto selected_atom_size{ context.AtomSize() };
    if (previous_state.result_list.size() != selected_atom_size ||
        previous_state.estimation_list.size() != selected_atom_size ||
        ridge_multiplier_list.size() != selected_atom_size)
    {
        throw std::invalid_argument("Local fitting iteration input sizes are inconsistent.");
    }
    auto current_snapshot{
        BuildFittedGaussianSnapshot(context, previous_state.estimation_list)
    };
    const auto joint_offset_result{
        EstimateJointOffsets(
            context, active_index_list, current_snapshot, ridge_ratio, ridge_multiplier_list)
    };
    ApplyJointOffsetsToSnapshot(active_index_list, joint_offset_result.offset, current_snapshot);
    auto iteration_state{ previous_state };
    std::vector<char> suspicious_offset_mask(active_index_list.size(), 0);
    std::vector<std::size_t> pre_refit_suspicious_seed_position_list;

    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto state_index{ active_index_list.at(i) };
        const auto previous_model{
            GaussianModel3D::FromVector(previous_state.estimation_list.at(state_index))
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
            joint_offset_result.active_coupling_graph,
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

    std::vector<std::size_t> post_refit_suspicious_seed_position_list;
    for (size_t i = 0; i < active_index_list.size(); i++)
    {
        if (suspicious_offset_mask.at(i) != 0) continue;

        const auto state_index{ active_index_list.at(i) };
        auto refit_result{
            FitAtomWithJointOffsetFallback(
                context,
                state_index,
                previous_state.result_list.at(state_index),
                current_snapshot,
                options)
        };
        if (!refit_result.has_value())
        {
            post_refit_suspicious_seed_position_list.emplace_back(i);
            continue;
        }
        auto result{ std::move(*refit_result) };
        const auto fitted_model{ result.mdpde.GetModel() };
        iteration_state.estimation_list.at(state_index) = fitted_model.ToVector();
        iteration_state.result_list.at(state_index) = std::move(result);
    }
    const auto post_refit_suspicious_position_list{
        ExpandSuspiciousOffsetClusters(
            joint_offset_result.active_coupling_graph,
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
    for (std::size_t i = 0; i < suspicious_offset_mask.size(); i++)
    {
        if (suspicious_offset_mask.at(i) == 0) continue;
        iteration_result.suspicious_offset_state_index_list.emplace_back(active_index_list.at(i));
    }
    return iteration_result;
}

void ApplyLocalFittingState(
    const GaussianFittingState & iteration_state,
    std::vector<AtomLocalPotentialEditor> & local_editor_list)
{
    if (local_editor_list.size() != iteration_state.result_list.size())
    {
        throw std::invalid_argument("local_editor_list and local fitting state sizes are inconsistent.");
    }

    for (std::size_t i = 0; i < local_editor_list.size(); i++)
    {
        local_editor_list.at(i).SetGaussianResult(iteration_state.result_list.at(i));
    }
}

void LogLocalFittingAllAtomsFrozen(const FitOptions & options, std::size_t accepted_iteration_count)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info,
        "Converged after " + std::to_string(accepted_iteration_count) +
        " iterations because all local atoms are frozen.");
}

void LogLocalFittingBacktrackingRetry(
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    double ridge_ratio,
    bool uses_cluster_local_objective_ridge)
{
    if (options.quiet_mode) return;

    std::ostringstream progress_message;
    progress_message
        << "Objective backtracking rejected all attempts; retrying after local fitting iteration "
        << accepted_iteration_count
        << std::fixed << std::setprecision(5)
        << "; acceleration history reset";
    if (uses_cluster_local_objective_ridge)
    {
        progress_message
            << ", next attempt uses increased cluster-local objective ridge"
            << ", global ridge ratio remains = " << ridge_ratio;
    }
    else
    {
        progress_message
            << ", next attempt uses increased global ridge ratio = " << ridge_ratio;
    }
    Logger::ProgressLine(progress_message.str());
}

void LogLocalFittingBacktrackingStop(
    const FitOptions & options,
    LocalFittingBacktrackingStopReason reason)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message
        << "Stopped local fitting because objective backtracking rejected all "
        << "acceleration and fixed-point attempts "
        << (reason == LocalFittingBacktrackingStopReason::MaximumGlobalRidge ?
            "at the maximum joint-offset ridge ratio" : "at the maximum iteration limit")
        << "; applying previous state.";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void LogLocalFittingProgress(
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    const LocalFittingAccelerationAttempt & current_acceleration_attempt,
    const algorithm::ConvergenceFreezeTracker & freeze_tracker)
{
    if (options.quiet_mode) return;

    std::ostringstream progress_message;
    progress_message << "Iter. " << accepted_iteration_count
        << '/' << kLocalFittingMaximumIterations
        << std::fixed << std::setprecision(4)
        << ", acceleration = "<< GetLocalFittingAccelerationText(current_acceleration_attempt)
        << ", damping = "<< current_acceleration_attempt.damping
        << ", active/frozen atoms = "<< freeze_tracker.GetActiveCount()
        << "/" << freeze_tracker.GetFrozenCount();
    Logger::ProgressLine(progress_message.str());
}

void LogLocalFittingConverged(
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    const algorithm::ParameterChangeStats & normalized_change_stats)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info,
        "Converged after " + std::to_string(accepted_iteration_count) +
        " iterations with normalized percentile amplitude change = " +
        std::to_string(normalized_change_stats.percentile_list.at(GaussianModel3D::AmplitudeIndex())) +
        ", normalized percentile width change = " +
        std::to_string(normalized_change_stats.percentile_list.at(GaussianModel3D::WidthIndex())) +
        ", and normalized percentile offset change = " +
        std::to_string(normalized_change_stats.percentile_list.at(GaussianModel3D::OffsetIndex())) +
        ".");
}

void LogLocalFittingMaximumIterations(
    const FitOptions & options,
    const algorithm::ParameterChangeStats & normalized_change_stats)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Warning,
        "Reached maximum iteration size; refitting at current accepted candidate "
        "with normalized percentile amplitude change = " +
        std::to_string(
            normalized_change_stats.percentile_list.at(GaussianModel3D::AmplitudeIndex())) +
        ", normalized percentile width change = " +
        std::to_string(
            normalized_change_stats.percentile_list.at(GaussianModel3D::WidthIndex())) +
        ", and normalized percentile offset change = " +
        std::to_string(
            normalized_change_stats.percentile_list.at(GaussianModel3D::OffsetIndex())));
}

void InitializeLocalFittingSeedModels(ModelObject & model_object)
{
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    auto local_editor_list{ BuildAtomLocalEditors(model_object, atom_list) };
    const auto seed_model{ GaussianModel3D{ 0.0, 1.0, 0.0 } };
    for (size_t i = 0; i < atom_list.size(); i++)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom_list[i]) };
        auto result{ local_view.GetGaussianResult() };
        result.ols = GaussianModel3DWithUncertainty{
            seed_model,
            GaussianModel3DUncertainty{}
        };
        result.mdpde = GaussianModel3DWithUncertainty{
            seed_model,
            GaussianModel3DUncertainty{}
        };
        result.posterior.reset();
        result.is_outlier = false;
        result.statistical_distance = 0.0;
        result.fit_result.reset();
        local_editor_list[i].SetGaussianResult(result);
    }
}

void RunFixedOffsetLocalFittingPass(
    ModelObject & model_object,
    const FitOptions & options,
    LocalFittingPass pass)
{
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    const auto selected_atom_size{ atom_list.size() };
    auto local_editor_list{ BuildAtomLocalEditors(model_object, atom_list) };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto stage_label{
        pass == LocalFittingPass::FirstStage ? "1st-stage" : "3rd-stage"
    };
    std::atomic<size_t> atom_count{ 0 };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info,
            "Run " + std::string{ stage_label } + " local atom fitting for " +
            std::to_string(selected_atom_size) + " atoms.");
    }

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        auto & atom{ *atom_list[i] };
        const auto local_view{ AtomLocalPotentialView::RequireFor(atom) };
        LocalPotentialSampleList sample_entries;
        GaussianModel3D offset_model;
        switch (pass)
        {
        case LocalFittingPass::FirstStage:
            sample_entries = local_view.GetSamplingEntries();
            offset_model = local_view.GetGaussianResult().mdpde.GetModel();
            break;
        case LocalFittingPass::ThirdStage:
            sample_entries = local_view.GetSamplingEntries(false, true);
            offset_model = analysis_view.GetAtomGroupPrior(data_internal::GetGroupKey(&atom));
            break;
        }

        auto result{
            EstimateLocalGaussian(sample_entries, local_view.GetAlphaR(), options, offset_model)
        };
        local_editor_list[i].SetGaussianResult(result);

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            atom_count++;
            if (!options.quiet_mode)
            {
                Logger::ProgressPercent(atom_count, selected_atom_size);
            }
        }
    }
}

} // namespace

double TrainAlphaR(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");

    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(sample_entries_list.size());
    for (const auto & sample_entries : sample_entries_list)
    {
        dataset_list.emplace_back(
            rhbm_helper::BuildMemberDataset(
                sample_entries,
                options.distance_min,
                options.distance_max));
    }
    auto training_options{ MakeTrainingOptions(options) };
    if (!dataset_list.empty())
    {
        const auto minimum_response_count{ GetMinimumDatasetResponseCount(dataset_list) };
        if (minimum_response_count < 2)
        {
            return training_options.alpha_min;
        }
        if (training_options.subset_size > minimum_response_count)
        {
            training_options.subset_size = minimum_response_count;
        }
    }
    return rhbm_trainer::CrossValidationAlphaR(dataset_list, training_options).best_alpha;
}

double TrainAlphaG(
    const std::vector<std::vector<LocalGaussianResult>> & member_result_list,
    const FitOptions & options)
{
    std::vector<std::vector<RHBMParameterVector>> beta_group_list;
    beta_group_list.reserve(member_result_list.size());
    for (const auto & member_results : member_result_list)
    {
        std::vector<RHBMParameterVector> beta_list;
        beta_list.reserve(member_results.size());
        for (const auto & member_result : member_results)
        {
            beta_list.emplace_back(
                linearization_service::EncodeGaussianToParameterVector(member_result.mdpde.GetModel()));
        }
        beta_group_list.emplace_back(std::move(beta_list));
    }

    const auto training_options{ MakeTrainingOptions(options) };
    if (beta_group_list.empty())
    {
        return training_options.alpha_min;
    }

    return rhbm_trainer::CrossValidationAlphaG(beta_group_list, training_options).best_alpha;
}

LocalGaussianResult EstimateLocalGaussian(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options,
    const GaussianModel3D & offset_model)
{
    auto range_min{ options.distance_min };
    auto range_max{ options.distance_max };
    numeric_validation::RequireFiniteNonNegativeRange(range_min, range_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");
    numeric_validation::RequireFinite(offset_model.GetOffset(), "offset");

    auto execution_options{ MakeExecutionOptions(options) };
    const auto updated_sample_entries{
        BuildSamplesForZeroOffsetGaussianFit(sample_entries, offset_model)
    };
    auto dataset{
        rhbm_helper::BuildMemberDataset(updated_sample_entries, range_min, range_max)
    };
    const auto result{ rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset, execution_options) };
    return DecodeLocalGaussianResult(alpha_r, result, offset_model.GetOffset());
}

GroupGaussianResult EstimateGroupGaussian(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const std::vector<LocalGaussianResult> & member_result_list,
    double alpha_g,
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_g, "alpha_g");

    if (sample_entries_list.size() != member_result_list.size())
    {
        throw std::invalid_argument("sample_entries_list and member_result_list sizes are inconsistent.");
    }

    auto execution_options{ MakeExecutionOptions(options) };
    const auto range_min{ options.distance_min };
    const auto range_max{ options.distance_max };
    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(sample_entries_list.size());
    for (std::size_t i = 0; i < sample_entries_list.size(); i++)
    {
        const auto sampling_entries{
            BuildSamplesForZeroOffsetGaussianFit(
                sample_entries_list.at(i),
                member_result_list.at(i).mdpde.GetModel())
        };
        dataset_list.emplace_back(
            rhbm_helper::BuildMemberDataset(sampling_entries, range_min, range_max));
    }
    if (dataset_list.size() != member_result_list.size())
    {
        throw std::invalid_argument("dataset_list and member_result_list sizes are inconsistent.");
    }

    std::vector<RHBMBetaEstimateResult> fit_result_list;
    fit_result_list.reserve(member_result_list.size());
    for (std::size_t i = 0; i < member_result_list.size(); i++)
    {
        fit_result_list.emplace_back(
            rhbm_helper::EstimateBetaMDPDE(
                member_result_list.at(i).alpha_r,
                dataset_list.at(i),
                execution_options));
    }
    const auto group_input{ rhbm_helper::BuildGroupInput(dataset_list, fit_result_list) };
    const auto raw_result{ rhbm_helper::EstimateGroup(alpha_g, group_input, execution_options) };
    std::vector<double> member_offset_list;
    member_offset_list.reserve(member_result_list.size());
    for (const auto & member_result : member_result_list)
    {
        member_offset_list.emplace_back(member_result.mdpde.GetModel().GetOffset());
    }
    const auto group_offset{ array_helper::ComputeMedian(member_offset_list) };
    auto result{ DecodeGroupGaussianResult(alpha_g, raw_result, group_offset) };
    result.member_results = DecodeMemberGaussianResults(raw_result, member_result_list);
    return result;
}

void RunLocalAlphaTraining(
    ModelObject & model_object,
    const FitOptions & options,
    bool apply_selection,
    bool use_updated_sample)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto group_key_list{ analysis_view.CollectAtomGroupKeys() };

    size_t count{ 0 };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run local alpha training for " + std::to_string(group_key_list.size()) + " groups.");
    }
    for (const auto group_key : group_key_list)
    {
        const auto & group_atom_list{ analysis_view.GetAtomObjectList(group_key) };
        std::vector<LocalPotentialSampleList> sample_entries_list;
        sample_entries_list.reserve(group_atom_list.size());
        for (auto * atom : group_atom_list)
        {
            analysis.EnsureAtomLocalPotential(*atom);
            const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
            const auto & sample_entries{ local_view.GetSamplingEntries(apply_selection, use_updated_sample) };
            if (!HasEnoughSamplesInFitRange(
                    sample_entries,
                    options.distance_min,
                    options.distance_max,
                    kMinimumAlphaRTrainingSampleCount)) continue;
            sample_entries_list.emplace_back(sample_entries);
        }
        sample_entries_list.shrink_to_fit();
        if (!sample_entries_list.empty())
        {
            const auto alpha_r{ TrainAlphaR(sample_entries_list, options) };
            for (auto * atom : group_atom_list)
            {
                analysis.EnsureAtomLocalPotential(*atom).SetAlphaR(alpha_r);
            }
        }
        count++;
        if (!options.quiet_mode)
        {
            Logger::ProgressPercent(count, group_key_list.size());
        }
    }
}

void RunGroupAlphaTraining(ModelObject & model_object, const FitOptions & options)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto group_key_list{ analysis_view.CollectAtomGroupKeys() };

    std::vector<std::vector<LocalGaussianResult>> member_result_list;
    member_result_list.reserve(group_key_list.size());
    for (const auto group_key : group_key_list)
    {
        const auto & group_atom_list{ analysis_view.GetAtomObjectList(group_key) };
        if (group_atom_list.size() < kMinimumAlphaGTrainingMemberCount) continue;
        if (group_atom_list.front()->IsMainChainAtom() == false) continue;

        std::vector<LocalGaussianResult> group_member_results;
        group_member_results.reserve(group_atom_list.size());
        for (auto * atom : group_atom_list)
        {
            analysis.EnsureAtomLocalPotential(*atom);
            const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
            group_member_results.emplace_back(local_view.GetGaussianResult());
        }
        member_result_list.emplace_back(std::move(group_member_results));
    }

    const auto alpha_g{ TrainAlphaG(member_result_list, options) };
    for (const auto group_key : analysis_view.CollectAtomGroupKeys())
    {
        analysis.SetAtomGroupAlphaG(group_key, alpha_g);
    }
}

void RunFirstStageLocalFitting(ModelObject & model_object, const FitOptions & options)
{
    RunFixedOffsetLocalFittingPass(model_object, options, LocalFittingPass::FirstStage);
}

void RunSecondStageLocalFitting(ModelObject & model_object, const FitOptions & options)
{
    const auto context{ BuildSecondStageLocalFittingContext(model_object) };
    auto local_editor_list{ BuildAtomLocalEditors(model_object, context.atom_list) };
    const auto atom_size{ context.AtomSize() };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run 2nd-stage local atom fitting with iterations...");
    }

    auto previous_state{ BuildInitialLocalFittingState(context) };
    algorithm::ClusteredAndersonAccelerationHistorySet acceleration_history{
        algorithm::AndersonAccelerationOptions{
            kLocalFittingAndersonHistoryDepth,
            kLocalFittingNormalizedChangeScaleFloor,
            kLocalFittingAndersonCoefficientL1Limit,
            kLocalFittingAndersonRegularization,
            kLocalFittingAndersonCoefficientAbsLimit
        }
    };
    algorithm::ConvergenceFreezeTracker freeze_tracker{
        atom_size,
        kLocalFittingParameterChangeTolerance,
        kLocalFittingFreezeChangeRatio,
        kLocalFittingFreezeStableIterations
    };
    algorithm::DependencyThawHysteresisTracker thaw_hysteresis_tracker{
        atom_size,
        kLocalFittingDependencyThawHysteresisGrowth,
        kLocalFittingDependencyThawHysteresisMax,
        kLocalFittingDependencyThawHysteresisFrozenDecay,
        kLocalFittingDependencyThawMaximumCount
    };
    double ridge_ratio{ kJointOffsetRidgeRatio };
    std::vector<double> suspicious_joint_offset_ridge_multiplier_list(atom_size, 1.0);
    algorithm::ClusteredFittingQualityStateSet<LocalFittingObjectiveSamples> cluster_quality_state{
        algorithm::ClusteredFittingQualityOptions{
            kLocalFittingObjectiveScaleWarmupCount,
            kLocalFittingConvergenceObjectiveRelativeTolerance,
            kLocalFittingObjectiveTieRelativeTolerance,
            1.0,
            kSuspiciousJointOffsetRidgeMultiplier,
            kJointOffsetRidgeGrowth,
            kJointOffsetRidgeShrink
        }
    };
    std::size_t accepted_iteration_count{ 0 };
    for (size_t iter = 0; iter < kLocalFittingMaximumIterations; iter++)
    {
        const auto active_index_list{ freeze_tracker.BuildActiveIndexList() };
        if (active_index_list.empty())
        {
            ApplyLocalFittingState(previous_state, local_editor_list);
            LogLocalFittingAllAtomsFrozen(options, accepted_iteration_count);
            break;
        }

        auto cluster_map{ BuildLocalFittingClusters(context, active_index_list) };
        const auto cluster_key_list{ BuildLocalFittingClusterKeyList(cluster_map) };
        cluster_quality_state.Reconcile(
            cluster_key_list,
            [&](const LocalFittingClusterKey & key)
            {
                return BuildInitialLocalFittingClusterQualityState(
                    context,
                    previous_state,
                    key,
                    cluster_map.at(key));
            });
        acceleration_history.Reconcile(cluster_key_list);
        const auto joint_offset_ridge_multiplier_list{
            CombineLocalFittingRidgeMultiplierLists(
                suspicious_joint_offset_ridge_multiplier_list,
                cluster_quality_state.BuildObjectiveRidgeMultiplierList(atom_size))
        };
        auto iteration_result{
            RunLocalFittingIteration(
                context,
                active_index_list,
                previous_state,
                options,
                ridge_ratio,
                joint_offset_ridge_multiplier_list)
        };
        const auto has_suspicious_offset_fallback{ !iteration_result.suspicious_offset_state_index_list.empty() };
        for (const auto active_index : active_index_list)
        {
            suspicious_joint_offset_ridge_multiplier_list.at(active_index) = 1.0;
        }
        for (const auto state_index : iteration_result.suspicious_offset_state_index_list)
        {
            if (state_index >= atom_size)
            {
                throw std::invalid_argument("Local fitting suspicious offset atom index is out of range.");
            }
            suspicious_joint_offset_ridge_multiplier_list.at(state_index) = kSuspiciousJointOffsetRidgeMultiplier;
        }
        const auto raw_state{ std::move(iteration_result.state) };
        auto selection{
            SelectLocalFittingClusterCandidates(
                context,
                std::move(cluster_map),
                cluster_key_list,
                previous_state,
                raw_state,
                acceleration_history,
                cluster_quality_state)
        };
        auto assembled_state{ std::move(selection.assembled_state) };
        if (selection.accepted_key_list.empty())
        {
            bool increased_cluster_objective_ridge{ false };
            bool increased_global_ridge_ratio{ false };
            if (!selection.rejected_key_list.empty())
            {
                increased_cluster_objective_ridge =
                    cluster_quality_state.IncreaseObjectiveRidge(selection.rejected_key_list);
            }
            if (selection.has_objective_backtracking_rejection && !increased_cluster_objective_ridge)
            {
                const auto previous_ridge_ratio{ ridge_ratio };
                ridge_ratio = std::min(kJointOffsetRidgeRatioMax, ridge_ratio * kJointOffsetRidgeGrowth);
                increased_global_ridge_ratio = ridge_ratio > previous_ridge_ratio;
            }
            if ((increased_cluster_objective_ridge || increased_global_ridge_ratio) &&
                iter + 1 < kLocalFittingMaximumIterations)
            {
                LogLocalFittingBacktrackingRetry(
                    options, accepted_iteration_count, ridge_ratio,
                    increased_cluster_objective_ridge);
                continue;
            }

            const auto stop_reason{
                iter + 1 >= kLocalFittingMaximumIterations ?
                    LocalFittingBacktrackingStopReason::MaximumIterationLimit :
                    LocalFittingBacktrackingStopReason::MaximumGlobalRidge
            };
            ApplyLocalFittingState(previous_state, local_editor_list);
            LogLocalFittingBacktrackingStop(options, stop_reason);
            break;
        }

        auto change_list{
            CalculateLocalFittingParameterChanges(
                assembled_state.estimation_list,
                previous_state.estimation_list)
        };
        const auto normalized_change_list{
            CalculateLocalFittingNormalizedParameterChanges(
                assembled_state.estimation_list,
                previous_state.estimation_list)
        };
        auto normalized_change_stats{
            algorithm::SummarizeParameterChangeStats(
                normalized_change_list,
                active_index_list,
                kLocalFittingChangePercentile)
        };
        accepted_iteration_count++;
        cluster_quality_state.CommitAccepted(selection.accepted_score_list);
        cluster_quality_state.DecreaseObjectiveRidge(selection.accepted_key_list);
        if (!selection.rejected_key_list.empty())
        {
            cluster_quality_state.IncreaseObjectiveRidge(selection.rejected_key_list);
        }

        if (!selection.has_objective_backtracking_rejection)
        {
            ridge_ratio = std::max(kJointOffsetRidgeRatioMin, ridge_ratio * kJointOffsetRidgeShrink);
        }
        if (has_suspicious_offset_fallback)
        {
            acceleration_history.ClearAndSuppressContaining(iteration_result.suspicious_offset_state_index_list);
        }
        acceleration_history.Commit(
            selection.accepted_key_list,
            previous_state.estimation_list,
            raw_state.estimation_list);
        std::vector<std::size_t> accepted_active_index_list;
        for (const auto & key : selection.accepted_key_list)
        {
            accepted_active_index_list.insert(
                accepted_active_index_list.end(), key.begin(), key.end());
        }
        freeze_tracker.Update(change_list, accepted_active_index_list);
        ThawChangedActiveAtomNeighbors(
            context, change_list, accepted_active_index_list,
            freeze_tracker, thaw_hysteresis_tracker);
        for (const auto state_index : iteration_result.suspicious_offset_state_index_list)
        {
            freeze_tracker.Thaw(state_index);
        }
        for (std::size_t state_index = 0; state_index < atom_size; state_index++)
        {
            if (freeze_tracker.IsFrozen(state_index))
            {
                thaw_hysteresis_tracker.DecayFrozen(state_index);
            }
        }

        LogLocalFittingProgress(
            options,
            accepted_iteration_count,
            selection.accepted_acceleration_attempt,
            freeze_tracker);

        if (freeze_tracker.GetActiveCount() == 0)
        {
            ApplyLocalFittingState(assembled_state, local_editor_list);
            LogLocalFittingAllAtomsFrozen(options, accepted_iteration_count);
            break;
        }

        const auto converged{
            !has_suspicious_offset_fallback &&
            !selection.has_objective_backtracking_rejection &&
            cluster_quality_state.AllActiveReferencesLocked(cluster_key_list) &&
            IsLocalFittingNormalizedParameterChangeConverged(normalized_change_stats)
        };
        if (converged)
        {
            ApplyLocalFittingState(assembled_state, local_editor_list);
            LogLocalFittingConverged(options, accepted_iteration_count, normalized_change_stats);
            break;
        }

        if (iter + 1 == kLocalFittingMaximumIterations)
        {
            ApplyLocalFittingState(assembled_state, local_editor_list);
            LogLocalFittingMaximumIterations(options, normalized_change_stats);
        }
        previous_state = std::move(assembled_state);
    }
}

void RunThirdStageLocalFitting(ModelObject & model_object, const FitOptions & options)
{
    RunFixedOffsetLocalFittingPass(model_object, options, LocalFittingPass::ThirdStage);
}

void RunGroupPotentialFitting(ModelObject & model_object, const FitOptions & options)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto & selected_atom_list{ model_object.GetSelectedAtoms() };
    for (auto * atom : selected_atom_list)
    {
        analysis.EnsureAtomLocalPotential(*atom);
    }
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run atom group fitting.");
    }

    auto group_key_list{ analysis_view.CollectAtomGroupKeys() };
    auto group_key_size{ group_key_list.size() };
    std::atomic<size_t> key_count{ 0 };

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t k = 0; k < group_key_size; k++)
    {
        auto group_key{ group_key_list[k] };
        const auto & atom_list{ analysis_view.GetAtomObjectList(group_key) };
        const auto alpha_g{ analysis_view.GetAtomAlphaG(group_key) };
        std::vector<LocalPotentialSampleList> sample_entries_list;
        std::vector<LocalGaussianResult> member_result_list;
        sample_entries_list.reserve(atom_list.size());
        member_result_list.reserve(atom_list.size());
        for (const auto & atom : atom_list)
        {
            const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
            sample_entries_list.emplace_back(local_view.GetSamplingEntries(false, true));
            member_result_list.emplace_back(local_view.GetGaussianResult());
        }
        const auto result{
            EstimateGroupGaussian(sample_entries_list, member_result_list, alpha_g, options)
        };

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            analysis.ApplyAtomGroupGaussianResult(group_key, result);
            key_count++;
            if (!options.quiet_mode)
            {
                Logger::ProgressBar(key_count, group_key_size);
            }
        }
    }
}

void RunPotentialFittingWorkflow(ModelObject & model_object, const FitOptions & options)
{
    RunLocalAlphaTraining(model_object, options, true, false);

    InitializeLocalFittingSeedModels(model_object);
    RunFirstStageLocalFitting(model_object, options);
    RunGroupAlphaTraining(model_object, options);
    RunGroupPotentialFitting(model_object, options);

    RunSecondStageLocalFitting(model_object, options);
    RunGroupAlphaTraining(model_object, options);
    SetUpdatedSamplingEntriesFromGroupMedianGaussian(model_object);
    RunGroupPotentialFitting(model_object, options);
    SetUpdatedSamplingEntriesFromFittedGroupGaussian(model_object);
    RunLocalAlphaTraining(model_object, options, false, true);
    RunThirdStageLocalFitting(model_object, options);

    RunGroupAlphaTraining(model_object, options);
    RunGroupPotentialFitting(model_object, options);
    if (!options.quiet_mode)
    {
        LogGroupPriorSpotSummary(model_object);
    }
}

} // namespace rhbm_gem::core
