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
#include <unordered_set>
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
constexpr double kSuspiciousProfileInnermostSignFlipRatio{ 0.25 };
constexpr double kSuspiciousProfileNoiseScaleMultiplier{ 3.0 };
constexpr double kSuspiciousProfileNoiseScaleMin{ 1.0e-12 };
constexpr double kLocalFittingTransformedChangeTolerance{ 1.0e-4 };
constexpr double kLocalFittingTransformedMaximumChangeTolerance{ 1.0e-3 };
constexpr double kNeighborContributionDistanceMax{ 2.5 };
constexpr double kNeighborAtomSearchRange{ 2.0 * kNeighborContributionDistanceMax };
constexpr std::size_t kLocalFittingMaximumIterations{ 100 };
constexpr std::size_t kLocalFittingAuditPatience{ 3 };
constexpr bool kApplyLocalFittingBestIteration{ true };
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
constexpr detail::LocalFittingObjectiveTolerance
kLocalFittingObjectiveStrictTolerance{ 1.0e-10, 1.0e-8 };
constexpr detail::LocalFittingObjectiveTolerance
kLocalFittingObjectiveProgressTolerance{ 1.0e-8, 1.0e-3 };
constexpr double kLocalFittingCouplingMinimumWeight{ 0.05 };
constexpr std::size_t kLocalFittingCouplingMaximumResidueCount{ 10 };
constexpr std::array<double, 6> kLocalFittingCouplingSensitivityMinimumWeightList{
    0.05,
    0.075,
    0.10,
    0.15,
    0.20,
    0.30
};
constexpr double kLocalFittingObjectiveResidualScaleFloorRatio{ 1.0e-6 };
constexpr double kLocalFittingFitRangeWeight{ 1.0 };
constexpr double kLocalFittingTailValidationWeight{ 0.25 };
constexpr double kLocalFittingOffsetPlausibilityPenaltyWeight{ 1.0e-2 };
constexpr std::array<double, detail::kTransformedChangeSize>
kLocalFittingTrustRegionParameterScale{ 0.50, 0.35, 1.0 };
constexpr double kLocalFittingTrustRegionInitialRadius{ 1.0 };
constexpr double kLocalFittingTrustRegionMinimumRadius{ 0.0625 };
constexpr double kLocalFittingTrustRegionMaximumRadius{ 4.0 };
constexpr double kLocalFittingTrustRegionShrinkFactor{ 0.5 };
constexpr double kLocalFittingTrustRegionGrowthFactor{ 2.0 };
constexpr double kLocalFittingTrustRegionBoundaryRatio{ 0.8 };
constexpr double kLocalFittingOffsetPeakRatioMax{ 1.0 };
constexpr std::size_t kSuspiciousProfileMinimumRadiusCount{ 3 };
constexpr double kSuspiciousProfileDistanceTolerance{ 1.0e-6 };
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
    double innermost_response{ 0.0 };
    double max_abs_response{ 0.0 };
    double robust_residual_scale{ 0.0 };
    std::vector<double> radius_response_median_list{};
};

struct SuspiciousProfileAnalysis
{
    bool all_responses_finite{ true };
    std::optional<ZeroOffsetProfileDiagnostics> profile{};
};

struct JointOffsetSolveResult
{
    JointOffsetSolveStatus status{ JointOffsetSolveStatus::SystemBuildFailed };
    Eigen::VectorXd offset{};
};

struct JointOffsetBuildResult
{
    algorithm::WeightedRidgeSystem system{};
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
    std::vector<char> rollback_atom_mask{};
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

using LocalFittingObjectiveSampleRef = detail::LocalFittingCouplingSampleId;

struct LocalFittingObjectiveScale
{
    double fit{ 0.0 };
    std::optional<double> tail{};
};

struct LocalFittingObjectiveClusterDomain
{
    std::vector<LocalFittingObjectiveSampleRef> fit_sample_ref_list{};
    std::vector<LocalFittingObjectiveSampleRef> tail_sample_ref_list{};
    std::optional<LocalFittingObjectiveScale> scale{};
};

struct LocalFittingObjectiveDomain
{
    std::map<LocalFittingClusterKey, LocalFittingObjectiveClusterDomain>
        cluster_by_key{};
    std::vector<LocalFittingClusterKey> owner_key_by_atom_index{};
    std::vector<std::vector<char>> fit_sample_mask_by_atom{};
    std::size_t active_atom_count{ 0 };
    std::size_t fit_sample_count{ 0 };
    std::size_t tail_sample_count{ 0 };
};

struct LocalFittingClusterObjectiveState
{
    std::optional<detail::LocalFittingObjectiveBreakdown> best_objective{};
    double best_maximum_transformed_change{ 0.0 };
};

using LocalFittingClusterObjectiveStateMap =
    std::map<LocalFittingClusterKey, LocalFittingClusterObjectiveState>;

struct LocalFittingAuditedState
{
    detail::LocalFittingObjectiveBreakdown objective{};
    LocalFittingState state{};
    LocalFittingPolishProvenance polish_provenance{};
    std::optional<std::size_t> accepted_iteration{};
};

struct LocalFittingBestAuditState
{
    std::optional<LocalFittingAuditedState> best{};
};

struct LocalFittingFinalStateSelection
{
    const LocalFittingState * state{ nullptr };
    const LocalFittingPolishProvenance * polish_provenance{ nullptr };
    const LocalFittingAuditedState * audit_state{ nullptr };
    detail::LocalFittingFinalStateSource source{
        detail::LocalFittingFinalStateSource::Unavailable
    };
};

LocalFittingFinalStateSelection SelectLocalFittingFinalState(
    const LocalFittingState & latest_validated_state,
    const LocalFittingPolishProvenance & latest_validated_polish_provenance,
    const std::optional<LocalFittingAuditedState> & audited_state)
{
    const auto source{ detail::SelectLocalFittingFinalStateSource(
        kApplyLocalFittingBestIteration,
        true,
        audited_state.has_value()) };
    if (source == detail::LocalFittingFinalStateSource::BestAudit)
    {
        return LocalFittingFinalStateSelection{
            &audited_state->state,
            &audited_state->polish_provenance,
            &*audited_state,
            source
        };
    }
    return LocalFittingFinalStateSelection{
        &latest_validated_state,
        &latest_validated_polish_provenance,
        nullptr,
        source
    };
}

std::string_view GetLocalFittingFinalStateSourceText(
    detail::LocalFittingFinalStateSource source)
{
    switch (source)
    {
    case detail::LocalFittingFinalStateSource::BestAudit:
        return "best-audit";
    case detail::LocalFittingFinalStateSource::LatestValidated:
        return "latest-validated";
    case detail::LocalFittingFinalStateSource::Unavailable:
        return "unavailable";
    }
    return "unavailable";
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

enum class LocalFittingPreObjectiveFailureReason
{
    None,
    InvalidModel,
    PreviousSharedOffsetProjectionOutsideTrustRegion,
    NoCandidateWithinTrustRegion
};

struct LocalFittingObjectiveAttemptDiagnostic
{
    double effective_damping{ 1.0 };
    bool is_invalid_model{ false };
    LocalFittingPreObjectiveFailureReason pre_objective_failure_reason{
        LocalFittingPreObjectiveFailureReason::None
    };
    std::optional<double> pre_objective_attempted_step_norm{};
    std::optional<double> fit_scale{};
    std::optional<double> tail_scale{};
    std::size_t fit_sample_count{ 0 };
    std::size_t tail_sample_count{ 0 };
    std::optional<detail::LocalFittingObjectiveBreakdown> candidate_objective{};
    std::optional<detail::LocalFittingObjectiveBreakdown> previous_objective{};
    std::optional<detail::LocalFittingObjectiveBreakdown> best_objective{};
    double trust_region_radius{ 0.0 };
    double trust_region_step_norm{ 0.0 };
    bool rejected_by_previous{ false };
    bool rejected_by_best{ false };
    std::size_t backtracking_trial_count{ 0 };
    std::optional<double> accepted_backtracking_factor{};
    bool backtracking_exhausted{ false };
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
    std::vector<LocalFittingRejectedClusterDiagnostic>
        accepted_cluster_diagnostic_list{};
    std::vector<LocalFittingRejectedClusterDiagnostic> rejected_cluster_diagnostic_list{};
    std::vector<LocalFittingClusterKey> grow_trust_region_key_list{};
    std::vector<LocalFittingClusterKey> backtracking_exhausted_key_list{};
    std::size_t combined_backtracking_trial_count{ 0 };
    std::optional<double> combined_backtracking_factor{};
    bool combined_backtracking_exhausted{ false };
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
    bool is_selected{ true };
    std::size_t atom_index{ 0 };
    double distance{ 0.0 };
};

struct SecondStageUnselectedContributor
{
    AtomObject * atom{ nullptr };
    GroupKey group_key{};
    std::optional<GaussianModel3DWithUncertainty> initial_seed{};
    detail::SecondStageSeedSource seed_source{
        detail::SecondStageSeedSource::GlobalMedian
    };
};

struct SecondStageAtomContext
{
    AtomObject * atom{ nullptr };
    LocalPotentialSampleList raw_sampling_entries{};
    std::vector<std::vector<SecondStageNeighborSample>> sample_neighbor_list{};
    double alpha_r{ 0.0 };
};

struct SecondStageLocalFittingContext
{
    std::vector<SecondStageAtomContext> selected_atom_list{};
    std::vector<SecondStageUnselectedContributor> unselected_atom_list{};

    std::size_t size() const { return selected_atom_list.size(); }
    SecondStageAtomContext & at(std::size_t index)
    {
        return selected_atom_list.at(index);
    }
    const SecondStageAtomContext & at(std::size_t index) const
    {
        return selected_atom_list.at(index);
    }
    auto begin() { return selected_atom_list.begin(); }
    auto end() { return selected_atom_list.end(); }
    auto begin() const { return selected_atom_list.begin(); }
    auto end() const { return selected_atom_list.end(); }
};

struct SecondStageSeedSelectionRecord
{
    std::size_t atom_index{ 0 };
    detail::SecondStageSeedSource source{
        detail::SecondStageSeedSource::GlobalMedian
    };
    GaussianModel3D original_model{};
    GaussianModel3D selected_model{};
};

struct SecondStageInitialStateBuildResult
{
    LocalFittingState state{};
    std::vector<SecondStageSeedSelectionRecord> selection_record_list{};
    std::vector<SecondStageSeedSelectionRecord>
        unselected_selection_record_list{};
};

bool HasSuspiciousCenterSignFlip(
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
            kSuspiciousProfileNoiseScaleMin)
    };
    const auto negative_threshold{
        std::max(
            kSuspiciousProfileInnermostSignFlipRatio * previous_innermost_response,
            noise_threshold)
    };
    return previous_innermost_response > noise_threshold &&
        candidate_innermost_response < -negative_threshold;
}

const char * GetSecondStageSeedSourceText(
    detail::SecondStageSeedSource source)
{
    switch (source)
    {
    case detail::SecondStageSeedSource::GroupPosterior:
        return "group-posterior";
    case detail::SecondStageSeedSource::GroupPrior:
        return "group-prior";
    case detail::SecondStageSeedSource::GroupMedian:
        return "group-median";
    case detail::SecondStageSeedSource::GlobalMedian:
        return "global-median";
    }
    throw std::logic_error("Unknown second-stage seed source.");
}

double CalculateZeroOffsetResponse(
    const LocalPotentialSample & sample,
    const GaussianModel3D & model)
{
    const auto distance{ static_cast<double>(sample.point.distance) };
    const auto model_offset{ model.ResponseAtDistance(distance) - model.SignalAtDistance(distance) };
    return static_cast<double>(sample.response) - model_offset;
}

double CalculateMedianAbsoluteDeviationScale(const std::vector<double> & value_list);

bool IsSameSuspiciousProfileRadius(double lhs, double rhs)
{
    const auto scale{ std::max({ std::abs(lhs), std::abs(rhs), 1.0 }) };
    return std::abs(lhs - rhs) <= kSuspiciousProfileDistanceTolerance * scale;
}

SuspiciousProfileAnalysis BuildSuspiciousProfileAnalysis(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & model,
    const FitOptions & options,
    bool calculate_residual_scale)
{
    SuspiciousProfileAnalysis analysis;
    std::vector<std::pair<double, double>> profile_samples;
    std::vector<double> residual_list;
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
    diagnostics.innermost_response = diagnostics.radius_response_median_list.front();
    if (calculate_residual_scale)
    {
        diagnostics.robust_residual_scale =
            CalculateMedianAbsoluteDeviationScale(residual_list);
    }
    analysis.profile = std::move(diagnostics);
    return analysis;
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
        previous_profile.max_abs_response <= kRobustScaleMin ||
        !std::isfinite(previous_profile.robust_residual_scale))
    {
        return false;
    }
    const auto innermost_scale{
        std::max(std::abs(previous_profile.innermost_response), kRobustScaleMin)
    };
    for (std::size_t i = 1; i < previous_profile.radius_response_median_list.size(); i++)
    {
        const auto current_scale{
            std::abs(previous_profile.radius_response_median_list.at(i))
        };
        if (current_scale >
            kSuspiciousProfileReboundCenterRatio * innermost_scale)
        {
            return false;
        }
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

bool HasSuspiciousRadialRebound(
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
            kSuspiciousProfileNoiseScaleMultiplier *
                previous_profile.robust_residual_scale,
            kSuspiciousProfileNoiseScaleMin)
    };
    const auto rebound_magnitude_threshold{
        std::max(
            kSuspiciousProfileReboundReferenceRatio *
                reference_innermost_scale,
            noise_threshold)
    };
    const auto upward_excursion_threshold{
        std::max(
            kSuspiciousProfileUpwardExcursionReferenceRatio *
                reference_innermost_scale,
            noise_threshold)
    };
    const auto candidate_innermost_scale{
        std::max(
            std::abs(candidate_profile.innermost_response),
            kSuspiciousProfileNoiseScaleMin)
    };
    int upward_excursion_count{ 0 };
    auto previous_abs_response{ std::abs(candidate_profile.radius_response_median_list.front()) };
    for (std::size_t i = 1; i < candidate_profile.radius_response_median_list.size(); i++)
    {
        const auto current_abs_response{
            std::abs(candidate_profile.radius_response_median_list.at(i))
        };
        if (current_abs_response >
                kSuspiciousProfileReboundCenterRatio * candidate_innermost_scale &&
            current_abs_response > rebound_magnitude_threshold)
        {
            return true;
        }
        if (current_abs_response >
            previous_abs_response + upward_excursion_threshold)
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
    const std::optional<ZeroOffsetProfileDiagnostics> & previous_profile)
{
    if (!std::isfinite(candidate_model.GetWidth()) || candidate_model.GetWidth() <= 0.0) return true;
    if (candidate_model.GetWidth() > kSuspiciousWidthGrowthLimit * previous_model.GetWidth()) return true;
    if (!previous_profile.has_value()) return false;
    const auto distance_range{ previous_profile->distance_max - previous_profile->distance_min };
    return distance_range > 0.0 && candidate_model.GetWidth() > kSuspiciousWidthRangeLimitRatio * distance_range;
}

bool HasSuspiciousAmplitudeOffsetCompensation(
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
            previous_profile.has_value() ?
                std::abs(previous_profile->innermost_response) : 0.0,
            std::abs(previous_model.SignalAtDistance(0.0)),
            kRobustScaleMin
        })
    };
    return signal_delta * offset_delta_response < 0.0 &&
        std::abs(signal_delta) > kSuspiciousCompensationResponseRatio * reference_scale &&
        std::abs(offset_delta_response) > kSuspiciousCompensationResponseRatio * reference_scale;
}

detail::SuspiciousGaussianReason EvaluateSuspiciousGaussianUpdate(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & previous_model,
    const GaussianModel3D & candidate_model,
    const FitOptions & options,
    const SuspiciousProfileAnalysis & previous_analysis,
    bool is_post_refit)
{
    if (is_post_refit && !detail::IsValidSecondStageGaussianModel(candidate_model))
    {
        return detail::SuspiciousGaussianReason::InvalidModel;
    }
    const auto candidate_analysis{
        BuildSuspiciousProfileAnalysis(
            sample_entries,
            candidate_model,
            options,
            false)
    };
    if (!candidate_analysis.all_responses_finite)
    {
        return detail::SuspiciousGaussianReason::NonFiniteResponse;
    }
    if (HasSuspiciousOffsetMagnitude(
            previous_model,
            candidate_model,
            previous_analysis.profile.has_value() ?
                previous_analysis.profile->max_abs_response : 0.0))
    {
        return detail::SuspiciousGaussianReason::OffsetMagnitude;
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
            return detail::SuspiciousGaussianReason::NonFiniteResponse;
        }
        if (HasSuspiciousCenterSignFlip(
                previous_analysis.profile->innermost_response,
                candidate_analysis.profile->innermost_response,
                previous_analysis.profile->robust_residual_scale))
        {
            return detail::SuspiciousGaussianReason::CenterSignFlip;
        }
        if (HasSuspiciousRadialRebound(
                *previous_analysis.profile,
                *candidate_analysis.profile))
        {
            return detail::SuspiciousGaussianReason::RadialRebound;
        }
    }
    if (!is_post_refit) return detail::SuspiciousGaussianReason::None;
    if (HasSuspiciousWidthGrowth(previous_model, candidate_model, previous_analysis.profile))
    {
        return detail::SuspiciousGaussianReason::WidthGrowth;
    }
    if (HasSuspiciousAmplitudeOffsetCompensation(previous_model, candidate_model, previous_analysis.profile))
    {
        return detail::SuspiciousGaussianReason::AmplitudeOffsetCompensation;
    }
    return detail::SuspiciousGaussianReason::None;
}

SuspiciousProfileAnalysis BuildPreviousSuspiciousProfileBaseline(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & previous_model,
    const FitOptions & options)
{
    return BuildSuspiciousProfileAnalysis(
        sample_entries,
        previous_model,
        options,
        true);
}

using FittedGaussianSnapshot = std::vector<GaussianModel3D>;
FittedGaussianSnapshot BuildFittedGaussianSnapshot(
    const LocalFittingState & state);
FittedGaussianSnapshot BuildUnselectedContributorSnapshot(
    const SecondStageLocalFittingContext & context,
    const FittedGaussianSnapshot & selected_snapshot);
const GaussianModel3D & ResolveSecondStageNeighborModel(
    const SecondStageNeighborSample & neighbor_sample,
    const FittedGaussianSnapshot & selected_snapshot,
    const FittedGaussianSnapshot & unselected_snapshot);

SecondStageLocalFittingContext BuildSecondStageLocalFittingContext(
    ModelObject & model_object,
    const FitOptions & options)
{
    SecondStageLocalFittingContext context;
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    context.selected_atom_list.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        context.selected_atom_list.emplace_back(SecondStageAtomContext{ atom });
    }
    std::unordered_map<const AtomObject *, std::size_t> atom_index_map;
    atom_index_map.reserve(context.size());
    for (std::size_t i = 0; i < context.size(); i++)
    {
        atom_index_map.emplace(context.at(i).atom, i);
    }
    std::unordered_map<const AtomObject *, std::size_t>
        unselected_atom_index_map;
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        auto & atom_context{ context.at(atom_index) };
        const auto * atom{ atom_context.atom };
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
        atom_context.raw_sampling_entries = local_view.GetRawSamplingEntries(false);
        atom_context.alpha_r = local_view.GetAlphaR();
    }

    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        auto & atom_context{ context.at(atom_index) };
        const auto * atom{ atom_context.atom };
        const auto neighbor_atom_list{ atom->FindNeighborAtoms(kNeighborAtomSearchRange) };

        atom_context.sample_neighbor_list.resize(atom_context.raw_sampling_entries.size());
        for (std::size_t sample_index = 0;
            sample_index < atom_context.raw_sampling_entries.size();
            sample_index++)
        {
            const auto & sample{ atom_context.raw_sampling_entries.at(sample_index) };
            auto & sample_neighbor_list{ atom_context.sample_neighbor_list.at(sample_index) };
            sample_neighbor_list.reserve(neighbor_atom_list.size());
            for (auto * neighbor_atom : neighbor_atom_list)
            {
                if (options.exclude_hydrogen &&
                    neighbor_atom->GetElement() == Element::HYDROGEN)
                {
                    continue;
                }
                const auto distance{
                    static_cast<double>(
                        array_helper::ComputeNorm<float>(
                            sample.point.position,
                            neighbor_atom->GetPositionRef()))
                };
                if (distance > kNeighborContributionDistanceMax) continue;

                const auto selected_iter{ atom_index_map.find(neighbor_atom) };
                if (selected_iter != atom_index_map.end())
                {
                    sample_neighbor_list.emplace_back(
                        SecondStageNeighborSample{
                            true,
                            selected_iter->second,
                            distance
                        });
                    continue;
                }

                auto contributor_iter{
                    unselected_atom_index_map.find(neighbor_atom)
                };
                if (contributor_iter == unselected_atom_index_map.end())
                {
                    const auto contributor_index{
                        context.unselected_atom_list.size()
                    };
                    context.unselected_atom_list.emplace_back(
                        SecondStageUnselectedContributor{
                            neighbor_atom,
                            data_internal::GetGroupKey(neighbor_atom)
                        });
                    contributor_iter = unselected_atom_index_map.emplace(
                        neighbor_atom,
                        contributor_index).first;
                }
                sample_neighbor_list.emplace_back(
                    SecondStageNeighborSample{
                        false,
                        contributor_iter->second,
                        distance
                    });
            }
        }
    }

    return context;
}

SecondStageLocalFittingDiagnostics BuildSecondStageLocalFittingDiagnostics(
    const SecondStageLocalFittingContext & context)
{
    SecondStageLocalFittingDiagnostics diagnostics;
    diagnostics.neighbor_count_by_serial_id.reserve(context.size());
    for (const auto & atom_context : context)
    {
        std::unordered_set<const AtomObject *> neighbor_atom_set;
        for (const auto & sample_neighbor_list : atom_context.sample_neighbor_list)
        {
            for (const auto & neighbor_sample : sample_neighbor_list)
            {
                const auto * neighbor_atom{ neighbor_sample.is_selected
                    ? context.selected_atom_list.at(neighbor_sample.atom_index).atom
                    : context.unselected_atom_list.at(neighbor_sample.atom_index).atom
                };
                neighbor_atom_set.emplace(neighbor_atom);
            }
        }
        diagnostics.neighbor_count_by_serial_id.emplace(
            atom_context.atom->GetSerialID(),
            neighbor_atom_set.size());
    }
    return diagnostics;
}

detail::LocalFittingCouplingTopology BuildLocalFittingCouplingTopology(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & initial_state)
{
    detail::LocalFittingCouplingGraphBuilder builder{ context.size() };
    const auto selected_snapshot{ BuildFittedGaussianSnapshot(initial_state) };
    const auto unselected_snapshot{
        BuildUnselectedContributorSnapshot(context, selected_snapshot)
    };
    const auto invalid_jacobian{
        Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN())
    };
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        const auto & atom_context{ context.at(atom_index) };
        for (std::size_t sample_index = 0;
            sample_index < atom_context.raw_sampling_entries.size();
            sample_index++)
        {
            const auto & sample{ atom_context.raw_sampling_entries.at(sample_index) };
            std::map<std::size_t, Eigen::Vector3d> jacobian_by_atom_index;
            const auto add_participant = [&](
                std::size_t participant_atom_index,
                const Eigen::Vector3d & jacobian)
            {
                auto [iter, inserted]{ jacobian_by_atom_index.emplace(
                    participant_atom_index,
                    jacobian) };
                if (!inserted)
                {
                    if (!iter->second.allFinite() || !jacobian.allFinite())
                    {
                        iter->second = invalid_jacobian;
                    }
                    else
                    {
                        iter->second += jacobian;
                    }
                }
            };
            const auto target_evaluation{
                detail::EvaluateLocalFittingTransformedResponse(
                    initial_state.at(atom_index).mdpde.GetModel(),
                    static_cast<double>(sample.point.distance))
            };
            add_participant(
                atom_index,
                target_evaluation.has_value() ?
                    target_evaluation->jacobian : invalid_jacobian);
            for (const auto & neighbor_sample :
                atom_context.sample_neighbor_list.at(sample_index))
            {
                const auto & neighbor_model{
                    ResolveSecondStageNeighborModel(
                        neighbor_sample,
                        selected_snapshot,
                        unselected_snapshot)
                };
                const auto neighbor_evaluation{
                    detail::EvaluateLocalFittingTransformedResponse(
                        neighbor_model,
                        neighbor_sample.distance)
                };
                const auto jacobian{
                    neighbor_evaluation.has_value() ?
                        neighbor_evaluation->jacobian : invalid_jacobian
                };
                if (neighbor_sample.is_selected)
                {
                    add_participant(neighbor_sample.atom_index, jacobian);
                    continue;
                }
                const auto group_key{
                    context.unselected_atom_list.at(
                        neighbor_sample.atom_index).group_key
                };
                for (std::size_t selected_index = 0;
                    selected_index < context.size();
                    selected_index++)
                {
                    if (data_internal::GetGroupKey(
                            context.at(selected_index).atom) != group_key)
                    {
                        continue;
                    }
                    add_participant(selected_index, jacobian);
                }
            }
            std::vector<detail::LocalFittingCouplingParticipant>
                participant_list;
            participant_list.reserve(jacobian_by_atom_index.size());
            for (const auto & [participant_atom_index, jacobian] :
                jacobian_by_atom_index)
            {
                participant_list.emplace_back(
                    detail::LocalFittingCouplingParticipant{
                        participant_atom_index,
                        jacobian
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
    auto topology{
        weighted_topology.has_value() ?
            std::move(*weighted_topology) : builder.BuildBinary()
    };
    std::vector<detail::LocalFittingCouplingResidueKey>
        residue_key_by_atom_index;
    residue_key_by_atom_index.reserve(context.size());
    for (const auto & atom_context : context)
    {
        residue_key_by_atom_index.emplace_back(
            atom_context.atom->GetChainID(),
            atom_context.atom->GetSequenceID());
    }
    return detail::ApplyLocalFittingCouplingResidueCutoff(
        std::move(topology),
        std::move(residue_key_by_atom_index),
        kLocalFittingCouplingMaximumResidueCount);
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
    SecondStageLocalFittingContext & context,
    const ModelAnalysisView & analysis_view,
    bool & unselected_seed_failure)
{
    unselected_seed_failure = false;
    SecondStageInitialStateBuildResult build_result;
    auto & state{ build_result.state };
    state.resize(context.size());
    std::vector<std::optional<GaussianModel3DWithUncertainty>> group_prior_list(
        context.size());
    std::unordered_map<GroupKey, std::vector<GaussianModel3D>> models_by_group;
    std::vector<GaussianModel3D> global_models;
    global_models.reserve(context.size());

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
        const auto direct_selection{
            detail::SelectSecondStageSeed(
                detail::SecondStageSeedCandidates{
                    result.posterior,
                    group_prior_list.at(i),
                    std::nullopt,
                    std::nullopt
                })
        };
        if (!direct_selection.has_value()) continue;

        models_by_group[group_key].emplace_back(
            direct_selection->model.GetModel());
        global_models.emplace_back(direct_selection->model.GetModel());
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
            detail::SelectSecondStageSeed(
                detail::SecondStageSeedCandidates{
                    result.posterior,
                    group_prior_list.at(i),
                    group_median,
                    global_median
                })
        };
        if (!selection.has_value()) return std::nullopt;

        result.mdpde = selection->model;
        build_result.selection_record_list.emplace_back(
            SecondStageSeedSelectionRecord{
                i,
                selection->source,
                original_model,
                selection->model.GetModel()
            });
    }

    for (std::size_t i = 0;
        i < context.unselected_atom_list.size();
        i++)
    {
        auto & contributor{ context.unselected_atom_list.at(i) };
        std::optional<GaussianModel3DWithUncertainty> group_median;
        const auto group_median_iter{
            median_by_group.find(contributor.group_key)
        };
        if (group_median_iter != median_by_group.end())
        {
            group_median = group_median_iter->second;
        }
        const auto selection{
            detail::SelectSecondStageSeed(
                detail::SecondStageSeedCandidates{
                    std::nullopt,
                    std::nullopt,
                    group_median,
                    global_median
                })
        };
        if (!selection.has_value())
        {
            unselected_seed_failure = true;
            return std::nullopt;
        }

        contributor.initial_seed = selection->model;
        contributor.seed_source = selection->source;
        build_result.unselected_selection_record_list.emplace_back(
            SecondStageSeedSelectionRecord{
                i,
                selection->source,
                GaussianModel3D{},
                selection->model.GetModel()
            });
    }
    return build_result;
}

void LogSecondStageSeedSelections(
    const std::vector<SecondStageSeedSelectionRecord> & selection_record_list,
    const FitOptions & options)
{
    if (options.quiet_mode || selection_record_list.empty()) return;

    constexpr std::array<detail::SecondStageSeedSource, 4> source_list{
        detail::SecondStageSeedSource::GroupPosterior,
        detail::SecondStageSeedSource::GroupPrior,
        detail::SecondStageSeedSource::GroupMedian,
        detail::SecondStageSeedSource::GlobalMedian
    };
    std::array<std::size_t, source_list.size()> source_count{};
    for (const auto & record : selection_record_list)
    {
        source_count.at(static_cast<std::size_t>(record.source))++;
    }

    std::ostringstream summary;
    summary << "Selected second-stage initial seeds = "
        << selection_record_list.size() << ", sources = ";
    for (std::size_t i = 0; i < source_list.size(); i++)
    {
        if (i != 0) summary << ", ";
        summary << GetSecondStageSeedSourceText(source_list.at(i))
            << ":" << source_count.at(i);
    }
    summary << ".";
    Logger::Log(LogLevel::Info, summary.str());

    for (const auto & record : selection_record_list)
    {
        std::ostringstream detail_message;
        detail_message << "Second-stage seed selection: atom index = "
            << record.atom_index
            << ", source = " << GetSecondStageSeedSourceText(record.source)
            << std::scientific << std::setprecision(2)
            << ", original MDPDE A/B/C = "
            << record.original_model.GetAmplitude() << "/"
            << record.original_model.GetWidth() << "/"
            << record.original_model.GetOffset()
            << ", selected A/B/C = "
            << record.selected_model.GetAmplitude() << "/"
            << record.selected_model.GetWidth() << "/"
            << record.selected_model.GetOffset() << ".";
        Logger::Log(LogLevel::Debug, detail_message.str());
    }
}

void LogUnselectedSecondStageSeedSelections(
    const SecondStageLocalFittingContext & context,
    const std::vector<SecondStageSeedSelectionRecord> & selection_record_list,
    const FitOptions & options)
{
    if (options.quiet_mode || selection_record_list.empty()) return;

    std::size_t group_median_count{ 0 };
    std::size_t global_median_count{ 0 };
    for (const auto & record : selection_record_list)
    {
        if (record.source == detail::SecondStageSeedSource::GroupMedian)
        {
            group_median_count++;
        }
        else if (record.source == detail::SecondStageSeedSource::GlobalMedian)
        {
            global_median_count++;
        }
    }

    std::ostringstream summary;
    summary << "Unselected second-stage neighbor seeds = "
        << selection_record_list.size()
        << ", sources = group-median:" << group_median_count
        << ", global-median:" << global_median_count << ".";
    Logger::Log(LogLevel::Info, summary.str());

    for (const auto & record : selection_record_list)
    {
        const auto & contributor{
            context.unselected_atom_list.at(record.atom_index)
        };
        std::ostringstream detail_message;
        detail_message
            << "Unselected second-stage neighbor seed selection: serial ID = "
            << contributor.atom->GetSerialID()
            << ", source = " << GetSecondStageSeedSourceText(record.source)
            << std::scientific << std::setprecision(2)
            << ", seed A/B/C = "
            << record.selected_model.GetAmplitude() << "/"
            << record.selected_model.GetWidth() << "/"
            << record.selected_model.GetOffset() << ".";
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

    const auto & residue_cutoff_summary{ topology.residue_cutoff_summary };
    std::ostringstream residue_cutoff_message;
    residue_cutoff_message
        << "Local-fitting residue cutoff: residues="
        << residue_cutoff_summary.residue_count
        << ", limit=" << kLocalFittingCouplingMaximumResidueCount
        << ", clusters=" << residue_cutoff_summary.cluster_count
        << ", max-residues=" << residue_cutoff_summary.maximum_residue_count
        << ", cutoff-edges=" << residue_cutoff_summary.cut_edge_count << ".";
    Logger::Log(LogLevel::Info, residue_cutoff_message.str());

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

FittedGaussianSnapshot BuildUnselectedContributorSnapshot(
    const SecondStageLocalFittingContext & context,
    const FittedGaussianSnapshot & selected_snapshot)
{
    if (selected_snapshot.size() != context.size())
    {
        throw std::invalid_argument(
            "Second-stage selected contributor snapshot size is inconsistent.");
    }

    std::unordered_map<GroupKey, std::vector<GaussianModel3D>>
        selected_models_by_group;
    selected_models_by_group.reserve(context.size());
    for (std::size_t i = 0; i < context.size(); i++)
    {
        selected_models_by_group[
            data_internal::GetGroupKey(context.at(i).atom)]
            .emplace_back(selected_snapshot.at(i));
    }

    std::unordered_map<GroupKey, GaussianModel3D> median_model_by_group;
    median_model_by_group.reserve(selected_models_by_group.size());
    for (const auto & [group_key, model_list] : selected_models_by_group)
    {
        const auto median_model{
            detail::BuildLocalFittingGaussianParameterMedian(model_list)
        };
        if (median_model.has_value())
        {
            median_model_by_group.emplace(group_key, *median_model);
        }
    }

    FittedGaussianSnapshot snapshot;
    snapshot.reserve(context.unselected_atom_list.size());
    for (const auto & contributor : context.unselected_atom_list)
    {
        const auto median_iter{
            median_model_by_group.find(contributor.group_key)
        };
        if (median_iter != median_model_by_group.end())
        {
            snapshot.emplace_back(median_iter->second);
            continue;
        }
        if (!contributor.initial_seed.has_value())
        {
            throw std::logic_error(
                "Second-stage unselected contributor seed is unavailable.");
        }
        snapshot.emplace_back(contributor.initial_seed->GetModel());
    }
    return snapshot;
}

const GaussianModel3D & ResolveSecondStageNeighborModel(
    const SecondStageNeighborSample & neighbor_sample,
    const FittedGaussianSnapshot & selected_snapshot,
    const FittedGaussianSnapshot & unselected_snapshot)
{
    return neighbor_sample.is_selected ?
        selected_snapshot.at(neighbor_sample.atom_index) :
        unselected_snapshot.at(neighbor_sample.atom_index);
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
    const auto unselected_snapshot{
        BuildUnselectedContributorSnapshot(context, snapshot)
    };
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
        throw std::runtime_error("Joint offset group parameterization is invalid.");
    }

    std::vector<int> active_position_by_atom_index(atom_size, -1);
    std::unordered_map<GroupKey, std::size_t> active_position_by_group_key;
    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto atom_index{ active_index_list.at(i) };
        active_position_by_atom_index.at(atom_index) = static_cast<int>(i);
        active_position_by_group_key.emplace(
            group_key_by_atom_position.at(i),
            i);
    }

    const auto column_count{ parameterization->ParameterCount() };
    std::vector<Eigen::Triplet<double>> triplet_list;
    std::vector<double> response_list;
    Eigen::VectorXd group_column_square_sum{
        Eigen::VectorXd::Zero(column_count)
    };
    std::map<std::pair<Eigen::Index, Eigen::Index>, double>
        group_column_cross_sum_map;
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
        for (std::size_t sample_index = 0;
            sample_index < atom_context.raw_sampling_entries.size();
            sample_index++)
        {
            const auto & sample{ atom_context.raw_sampling_entries.at(sample_index) };
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
                const auto & neighbor_model{
                    ResolveSecondStageNeighborModel(
                        neighbor_sample,
                        snapshot,
                        unselected_snapshot)
                };
                int neighbor_position{ -1 };
                if (neighbor_sample.is_selected)
                {
                    neighbor_position = active_position_by_atom_index.at(
                        neighbor_sample.atom_index);
                }
                else
                {
                    const auto group_key{
                        context.unselected_atom_list.at(
                            neighbor_sample.atom_index).group_key
                    };
                    const auto position_iter{
                        active_position_by_group_key.find(group_key)
                    };
                    if (position_iter != active_position_by_group_key.end())
                    {
                        neighbor_position = static_cast<int>(
                            position_iter->second);
                    }
                }
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
            previous_offset
        };
    }
    auto system{ std::move(build_result.system) };
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
                previous_offset
            };
        }
        return JointOffsetSolveResult{
            status,
            std::move(*atom_offset)
        };
    };
    if (system.response.size() == 0 || system.previous_parameter.size() == 0)
    {
        return JointOffsetSolveResult{
            JointOffsetSolveStatus::EmptySystem,
            previous_offset
        };
    }

    Eigen::VectorXd weight{ Eigen::VectorXd::Ones(system.response.size()) };
    algorithm::WeightedRidgeSolver solver{ system };
    Eigen::VectorXd offset;
    if (!solver.Solve(system, weight, offset))
    {
        return JointOffsetSolveResult{
            JointOffsetSolveStatus::InitialSolveFailed,
            previous_offset
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
            return JointOffsetSolveResult{ JointOffsetSolveStatus::IrlsSolveFailed, previous_offset };
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
            return make_progress_result(JointOffsetSolveStatus::Converged, offset);
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
    const FittedGaussianSnapshot & selected_snapshot,
    const FittedGaussianSnapshot & unselected_snapshot)
{
    const auto & atom_context{ context.at(atom_index) };
    auto response_value{
        static_cast<double>(atom_context.raw_sampling_entries.at(sample_index).response)
    };
    for (const auto & neighbor_sample : atom_context.sample_neighbor_list.at(sample_index))
    {
        response_value -= ResolveSecondStageNeighborModel(
            neighbor_sample,
            selected_snapshot,
            unselected_snapshot).ResponseAtDistance(neighbor_sample.distance);
    }
    return response_value;
}

LocalPotentialSampleList BuildSecondStageAdjustedSamples(
    const SecondStageLocalFittingContext & context,
    std::size_t atom_index,
    const FittedGaussianSnapshot & snapshot)
{
    const auto & atom_context{ context.at(atom_index) };
    const auto unselected_snapshot{
        BuildUnselectedContributorSnapshot(context, snapshot)
    };
    LocalPotentialSampleList adjusted_sampling_entries;
    adjusted_sampling_entries.reserve(atom_context.raw_sampling_entries.size());
    for (std::size_t sample_index = 0;
        sample_index < atom_context.raw_sampling_entries.size();
        sample_index++)
    {
        auto sample{ atom_context.raw_sampling_entries.at(sample_index) };
        sample.response = static_cast<float>(CalculateSecondStageAdjustedResponse(
            context,
            atom_index,
            sample_index,
            snapshot,
            unselected_snapshot));
        adjusted_sampling_entries.emplace_back(sample);
    }
    return adjusted_sampling_entries;
}
struct LocalFittingResidualSample
{
    double adjusted_response{ 0.0 };
    double residual{ 0.0 };
};

std::optional<LocalFittingResidualSample> EvaluateLocalFittingResidualSample(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & state,
    const LocalFittingObjectiveSampleRef & sample_ref,
    const FittedGaussianSnapshot & selected_snapshot,
    const FittedGaussianSnapshot & unselected_snapshot)
{
    const auto & atom_context{ context.at(sample_ref.atom_index) };
    const auto & sample{
        atom_context.raw_sampling_entries.at(sample_ref.sample_index)
    };
    auto adjusted_response{ static_cast<double>(sample.response) };
    for (const auto & neighbor_sample :
        atom_context.sample_neighbor_list.at(sample_ref.sample_index))
    {
        adjusted_response -= ResolveSecondStageNeighborModel(
            neighbor_sample,
            selected_snapshot,
            unselected_snapshot).ResponseAtDistance(neighbor_sample.distance);
    }
    const auto expected_response{
        state.at(sample_ref.atom_index).mdpde.GetModel().ResponseAtDistance(
            static_cast<double>(sample.point.distance))
    };
    const auto residual{ adjusted_response - expected_response };
    if (!std::isfinite(adjusted_response) || !std::isfinite(residual))
    {
        return std::nullopt;
    }
    return LocalFittingResidualSample{ adjusted_response, residual };
}

std::optional<double> BuildFixedLocalFittingObjectiveScale(
    const std::vector<double> & residual_list,
    const std::vector<double> & adjusted_response_list)
{
    if (residual_list.empty() ||
        residual_list.size() != adjusted_response_list.size())
    {
        return std::nullopt;
    }
    const auto scale{
        std::max({
            CalculateMedianAbsoluteDeviationScale(residual_list),
            kLocalFittingObjectiveResidualScaleFloorRatio *
                CalculateMedianAbsoluteDeviationScale(adjusted_response_list),
            kRobustScaleMin
        })
    };
    return numeric_validation::IsFinitePositive(scale) ?
        std::optional<double>{ scale } : std::nullopt;
}

LocalFittingObjectiveDomain BuildLocalFittingObjectiveDomain(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & initial_state,
    const detail::LocalFittingCouplingPartition & partition,
    const FitOptions & options)
{
    LocalFittingObjectiveDomain domain;
    const auto selected_snapshot{ BuildFittedGaussianSnapshot(initial_state) };
    const auto unselected_snapshot{
        BuildUnselectedContributorSnapshot(context, selected_snapshot)
    };
    domain.owner_key_by_atom_index.resize(context.size());
    domain.fit_sample_mask_by_atom.resize(context.size());
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        domain.fit_sample_mask_by_atom.at(atom_index).resize(
            context.at(atom_index).raw_sampling_entries.size(),
            0);
    }
    for (const auto & [key, affected_sample_ref_list] :
        partition.sample_id_list_by_key)
    {
        static_cast<void>(affected_sample_ref_list);
        auto & cluster_domain{ domain.cluster_by_key[key] };
        domain.active_atom_count += key.size();
        std::vector<double> fit_residual_list;
        std::vector<double> fit_response_list;
        std::vector<double> tail_residual_list;
        std::vector<double> tail_response_list;
        for (const auto atom_index : key)
        {
            domain.owner_key_by_atom_index.at(atom_index) = key;
            const auto & raw_sampling_entries{
                context.at(atom_index).raw_sampling_entries
            };
            for (std::size_t sample_index = 0;
                sample_index < raw_sampling_entries.size();
                sample_index++)
            {
                const LocalFittingObjectiveSampleRef sample_ref{
                    atom_index,
                    sample_index
                };
                const auto residual_sample{
                    EvaluateLocalFittingResidualSample(
                        context,
                        initial_state,
                        sample_ref,
                        selected_snapshot,
                        unselected_snapshot)
                };
                const auto distance{
                    static_cast<double>(
                        raw_sampling_entries.at(sample_index).point.distance)
                };
                const auto is_fit_range{
                    distance >= options.distance_min &&
                    distance <= options.distance_max
                };
                domain.fit_sample_mask_by_atom.at(atom_index).at(sample_index) =
                    is_fit_range ? 1 : 0;
                auto & sample_ref_list{
                    is_fit_range ?
                        cluster_domain.fit_sample_ref_list :
                        cluster_domain.tail_sample_ref_list
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
            BuildFixedLocalFittingObjectiveScale(
                fit_residual_list,
                fit_response_list)
        };
        if (!fit_scale.has_value() ||
            fit_residual_list.size() != cluster_domain.fit_sample_ref_list.size())
        {
            continue;
        }
        LocalFittingObjectiveScale scale;
        scale.fit = *fit_scale;
        if (!cluster_domain.tail_sample_ref_list.empty())
        {
            scale.tail = BuildFixedLocalFittingObjectiveScale(
                tail_residual_list,
                tail_response_list);
            if (!scale.tail.has_value() ||
                tail_residual_list.size() !=
                    cluster_domain.tail_sample_ref_list.size())
            {
                continue;
            }
        }
        cluster_domain.scale = scale;
    }
    return domain;
}

void LogLocalFittingObjectiveDomain(
    const LocalFittingObjectiveDomain & domain,
    const FitOptions & options,
    bool is_terminal_reset = false)
{
    if (options.quiet_mode) return;
    std::vector<double> fit_scale_list;
    std::vector<double> tail_scale_list;
    for (const auto & [key, cluster_domain] : domain.cluster_by_key)
    {
        static_cast<void>(key);
        if (!cluster_domain.scale.has_value()) continue;
        fit_scale_list.emplace_back(cluster_domain.scale->fit);
        if (cluster_domain.scale->tail.has_value())
        {
            tail_scale_list.emplace_back(*cluster_domain.scale->tail);
        }
    }
    const auto append_scale_summary = [](
        std::ostringstream & message,
        const std::vector<double> & scale_list)
    {
        if (scale_list.empty())
        {
            message << "unavailable";
            return;
        }
        message
            << array_helper::ComputePercentile(scale_list, 0.5) << "/"
            << array_helper::ComputePercentile(scale_list, 0.99) << "/"
            << *std::max_element(scale_list.begin(), scale_list.end());
    };

    std::ostringstream message;
    message
        << (is_terminal_reset ?
            "Reset second-stage objective domain" :
            "Initialize second-stage objective domain")
        << ": fit/tail/offset weights = "
        << kLocalFittingFitRangeWeight << "/"
        << kLocalFittingTailValidationWeight << "/"
        << kLocalFittingOffsetPlausibilityPenaltyWeight
        << ", clusters = " << domain.cluster_by_key.size()
        << ", active atoms = " << domain.active_atom_count
        << ", unique fit/tail samples = "
        << domain.fit_sample_count << "/"
        << domain.tail_sample_count
        << ", fixed fit scale median/p99/max = ";
    append_scale_summary(message, fit_scale_list);
    message << ", fixed tail scale median/p99/max = ";
    append_scale_summary(message, tail_scale_list);
    message << ".";
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info, message.str());
}

std::optional<detail::LocalFittingObjectiveBreakdown>
EvaluateLocalFittingObjectiveContribution(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & state,
    const LocalFittingClusterKey & changed_key,
    const std::vector<LocalFittingObjectiveSampleRef> & sample_ref_list,
    const LocalFittingObjectiveDomain & domain)
{
    if (domain.active_atom_count == 0) return std::nullopt;
    const auto selected_snapshot{ BuildFittedGaussianSnapshot(state) };
    const auto unselected_snapshot{
        BuildUnselectedContributorSnapshot(context, selected_snapshot)
    };
    detail::LocalFittingObjectiveBreakdown breakdown;
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
        const auto residual_sample{
            EvaluateLocalFittingResidualSample(
                context,
                state,
                sample_ref,
                selected_snapshot,
                unselected_snapshot)
        };
        if (!residual_sample.has_value()) return std::nullopt;
        const auto is_fit_range{
            domain.fit_sample_mask_by_atom.at(sample_ref.atom_index).at(
                sample_ref.sample_index) != 0
        };
        const auto sample_count{
            is_fit_range ?
                owner_iter->second.fit_sample_ref_list.size() :
                owner_iter->second.tail_sample_ref_list.size()
        };
        if (sample_count == 0) return std::nullopt;
        const auto scale{
            is_fit_range ?
                std::optional<double>{ owner_iter->second.scale->fit } :
                owner_iter->second.scale->tail
        };
        if (!scale.has_value()) return std::nullopt;
        const auto loss{
            algorithm::CalculateCauchyLoss(
                residual_sample->residual / *scale,
                kRobustLossCutoffMultiplier)
        };
        const auto coefficient{
            detail::CalculateLocalFittingClusterAtomWeight(
                owner_key.size(),
                domain.active_atom_count) /
            static_cast<double>(sample_count)
        };
        if (is_fit_range)
        {
            breakdown.fit_range_residual_objective +=
                kLocalFittingFitRangeWeight * coefficient * loss;
        }
        else
        {
            breakdown.tail_validation_loss += coefficient * loss;
        }
    }
    for (const auto atom_index : changed_key)
    {
        const auto owner_iter{
            domain.cluster_by_key.find(
                domain.owner_key_by_atom_index.at(atom_index))
        };
        if (owner_iter == domain.cluster_by_key.end() ||
            !owner_iter->second.scale.has_value())
        {
            return std::nullopt;
        }
        const auto & model{ state.at(atom_index).mdpde.GetModel() };
        if (!detail::IsValidSecondStageGaussianModel(model)) return std::nullopt;
        const auto peak_signal{ model.SignalAtDistance(0.0) };
        const auto offset_peak{
            model.GetOffset() * model.OffsetBasisAtDistance(0.0)
        };
        if (!std::isfinite(peak_signal) || !std::isfinite(offset_peak))
        {
            return std::nullopt;
        }
        const auto offset_ratio{
            std::abs(offset_peak) /
            std::max({
                std::abs(peak_signal),
                owner_iter->second.scale->fit,
                kRobustScaleMin
            })
        };
        const auto offset_excess{
            std::max(0.0, offset_ratio - kLocalFittingOffsetPeakRatioMax)
        };
        breakdown.offset_plausibility_penalty +=
            kLocalFittingOffsetPlausibilityPenaltyWeight *
            offset_excess * offset_excess /
            static_cast<double>(domain.active_atom_count);
    }
    return detail::BuildLocalFittingObjectiveBreakdown(
        breakdown.fit_range_residual_objective,
        breakdown.tail_validation_loss,
        breakdown.offset_plausibility_penalty,
        kLocalFittingTailValidationWeight);
}

std::optional<detail::LocalFittingObjectiveBreakdown>
EvaluateLocalFittingAuditObjective(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & state,
    const LocalFittingObjectiveDomain & domain)
{
    detail::LocalFittingObjectiveBreakdown total;
    for (const auto & [key, cluster_domain] : domain.cluster_by_key)
    {
        std::vector<LocalFittingObjectiveSampleRef> owner_sample_ref_list{
            cluster_domain.fit_sample_ref_list
        };
        owner_sample_ref_list.insert(
            owner_sample_ref_list.end(),
            cluster_domain.tail_sample_ref_list.begin(),
            cluster_domain.tail_sample_ref_list.end());
        const auto contribution{
            EvaluateLocalFittingObjectiveContribution(
                context,
                state,
                key,
                owner_sample_ref_list,
                domain)
        };
        if (!contribution.has_value()) return std::nullopt;
        total.fit_range_residual_objective +=
            contribution->fit_range_residual_objective;
        total.tail_validation_loss += contribution->tail_validation_loss;
        total.tail_validation_penalty +=
            contribution->tail_validation_penalty;
        total.offset_plausibility_penalty +=
            contribution->offset_plausibility_penalty;
    }
    total.total_objective =
        total.fit_range_residual_objective +
        total.tail_validation_penalty +
        total.offset_plausibility_penalty;
    return std::isfinite(total.total_objective) ?
        std::optional<detail::LocalFittingObjectiveBreakdown>{ total } :
        std::nullopt;
}

bool IsLocalFittingCombinedObjectiveAcceptable(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & candidate_state,
    const LocalFittingState & previous_state,
    const LocalFittingObjectiveDomain & domain,
    const LocalFittingBestAuditState & audit_state)
{
    const auto candidate_objective{
        EvaluateLocalFittingAuditObjective(context, candidate_state, domain)
    };
    const auto previous_objective{
        EvaluateLocalFittingAuditObjective(context, previous_state, domain)
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
        kLocalFittingObjectiveProgressTolerance);
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
    const LocalFittingObjectiveDomain & domain,
    LocalFittingBestAuditState & audit_state)
{
    const auto candidate_objective{
        EvaluateLocalFittingAuditObjective(context, candidate_state, domain)
    };
    if (!candidate_objective.has_value()) return false;
    if (audit_state.best.has_value() &&
        !detail::IsBetterLocalFittingAuditObjective(
            candidate_objective->total_objective,
            audit_state.best->objective.total_objective,
            kLocalFittingObjectiveStrictTolerance))
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
    const LocalFittingState & initial_state,
    const LocalFittingPolishProvenance & initial_polish_provenance,
    const std::optional<std::size_t> & accepted_iteration,
    const LocalFittingObjectiveDomain & domain)
{
    LocalFittingBestAuditState audit_state;
    static_cast<void>(TryUpdateLocalFittingBestAuditState(
        context,
        initial_state,
        initial_polish_provenance,
        accepted_iteration,
        domain,
        audit_state));
    return audit_state;
}

void ResetLocalFittingBestAuditAfterObjectiveDomainChange(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & validated_state,
    const LocalFittingPolishProvenance & validated_polish_provenance,
    std::size_t accepted_iteration,
    const LocalFittingObjectiveDomain & domain,
    LocalFittingBestAuditState & audit_state)
{
    audit_state = BuildInitialLocalFittingBestAuditState(
        context,
        validated_state,
        validated_polish_provenance,
        accepted_iteration,
        domain);
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
            SummarizeLocalFittingTransformedChanges(assembled_state, previous_state, key)
        };
        if (!IsLocalFittingTransformedPercentileConverged(transformed_change_summary.percentile_stats))
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
    std::unordered_map<GroupKey, std::size_t> local_position_by_group_key;
    for (std::size_t local_position = 0; local_position < key.size(); local_position++)
    {
        const auto atom_index{ key.at(local_position) };
        local_position_by_atom_index.at(atom_index) =
            static_cast<int>(local_position);
        local_position_by_group_key.emplace(
            data_internal::GetGroupKey(context.at(atom_index).atom),
            local_position);
    }
    auto selected_snapshot{ BuildFittedGaussianSnapshot(base_state) };
    for (std::size_t local_position = 0;
        local_position < key.size();
        local_position++)
    {
        selected_snapshot.at(key.at(local_position)) =
            seed_model_list.at(local_position);
    }
    const auto unselected_snapshot{
        BuildUnselectedContributorSnapshot(context, selected_snapshot)
    };

    std::vector<Eigen::Triplet<double>> triplet_list;
    std::vector<double> residual_list;
    residual_list.reserve(sample_ref_list.size());
    for (const auto & sample_ref : sample_ref_list)
    {
        const auto & atom_context{
            context.at(sample_ref.atom_index)
        };
        const auto & sample{
            atom_context.raw_sampling_entries.at(sample_ref.sample_index)
        };
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
        const auto append_unselected_model = [&](
            std::size_t contributor_index,
            double distance) -> bool
        {
            const auto & contributor{
                context.unselected_atom_list.at(contributor_index)
            };
            const auto & model{ unselected_snapshot.at(contributor_index) };
            const auto evaluation{
                detail::EvaluateLocalFittingSharedOffsetResponse(model, distance)
            };
            if (!evaluation.has_value()) return false;
            predicted_response += evaluation->response;

            const auto local_position_iter{
                local_position_by_group_key.find(contributor.group_key)
            };
            if (local_position_iter == local_position_by_group_key.end() ||
                std::abs(evaluation->offset_jacobian) <=
                    std::numeric_limits<double>::epsilon())
            {
                return true;
            }
            triplet_list.emplace_back(
                row_index,
                parameterization.OffsetColumn(local_position_iter->second),
                evaluation->offset_jacobian);
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
            const auto appended{
                neighbor_sample.is_selected ?
                    append_model(
                        neighbor_sample.atom_index,
                        neighbor_sample.distance) :
                    append_unselected_model(
                        neighbor_sample.atom_index,
                        neighbor_sample.distance)
            };
            if (!appended)
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
        for (Eigen::SparseMatrix<double>::InnerIterator iter(system.design_matrix, column_index); iter; ++iter)
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
            static_cast<Eigen::Index>(key.size() * detail::kLocalFittingJointPolishShapeParameterSize)
        };
        if (column_index < offset_column_base)
        {
            const auto local_position{
                static_cast<std::size_t>(column_index) / detail::kLocalFittingJointPolishShapeParameterSize
            };
            parameter_multiplier = ridge_multiplier_list.at(key.at(local_position));
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
    for (std::size_t atom_position = 0; atom_position < key.size(); atom_position++)
    {
        const auto atom_index{ key.at(atom_position) };
        const auto previous{
            detail::EncodeLocalFittingTransformedCoordinates(
                outer_previous_state.at(atom_index).mdpde.GetModel())
        };
        const auto candidate{
            detail::EncodeLocalFittingTransformedCoordinates(candidate_model_list.at(atom_position))
        };
        if (!previous.has_value() || !candidate.has_value())
        {
            return std::nullopt;
        }
        for (std::size_t parameter_index = 0; parameter_index < detail::kTransformedChangeSize; parameter_index++)
        {
            const auto eigen_index{ static_cast<Eigen::Index>(parameter_index) };
            step_norm = std::max(
                step_norm,
                std::abs((*candidate)(eigen_index) - (*previous)(eigen_index)) /
                    kLocalFittingTrustRegionParameterScale.at(parameter_index));
        }
    }
    return std::isfinite(step_norm) ? std::optional<double>{ step_norm } : std::nullopt;
}

struct LocalFittingBaseProposal
{
    LocalFittingState state{};
    double effective_damping{ 0.0 };
    double step_norm{ 0.0 };
};

struct LocalFittingBaseProposalBuildResult
{
    std::optional<LocalFittingBaseProposal> proposal{};
    LocalFittingPreObjectiveFailureReason failure_reason{ LocalFittingPreObjectiveFailureReason::None };
    std::optional<double> attempted_step_norm{};
};

LocalFittingBaseProposalBuildResult
BuildLocalFittingSharedOffsetBaseProposal(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & outer_previous_state,
    const LocalFittingState & raw_state,
    const LocalFittingClusterKey & key,
    double trust_region_radius)
{
    constexpr double trust_region_tolerance{ 1.0e-12 };
    if (key.empty())
    {
        return LocalFittingBaseProposalBuildResult{
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
        group_key_by_atom_position.emplace_back(data_internal::GetGroupKey(context.at(atom_index).atom));
        previous_model_list.emplace_back(outer_previous_state.at(atom_index).mdpde.GetModel());
        raw_model_list.emplace_back(raw_state.at(atom_index).mdpde.GetModel());
    }
    const auto previous_shared_offset_model_list{
        detail::BuildLocalFittingGroupMedianModelList(group_key_by_atom_position, previous_model_list)
    };
    const auto raw_shared_offset_model_list{
        detail::BuildLocalFittingGroupMedianModelList(group_key_by_atom_position, raw_model_list)
    };

    const auto seed_model_list{
        detail::BuildLocalFittingSharedOffsetDampedModelList(
            previous_model_list,
            raw_model_list,
            previous_shared_offset_model_list,
            raw_shared_offset_model_list,
            0.0)
    };
    if (!seed_model_list.has_value())
    {
        return LocalFittingBaseProposalBuildResult{
            std::nullopt,
            LocalFittingPreObjectiveFailureReason::InvalidModel,
            std::nullopt
        };
    }
    const auto seed_step_norm{
        CalculateLocalFittingClusterModelTrustRegionStepNorm(outer_previous_state, key, *seed_model_list)
    };
    if (!seed_step_norm.has_value())
    {
        return LocalFittingBaseProposalBuildResult{
            std::nullopt,
            LocalFittingPreObjectiveFailureReason::InvalidModel,
            std::nullopt
        };
    }
    if (*seed_step_norm > trust_region_radius + trust_region_tolerance)
    {
        return LocalFittingBaseProposalBuildResult{
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
            detail::BuildLocalFittingSharedOffsetDampedModelList(
                previous_model_list,
                raw_model_list,
                previous_shared_offset_model_list,
                raw_shared_offset_model_list,
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
            if (step_norm.has_value() && *step_norm <= trust_region_radius + trust_region_tolerance)
            {
                LocalFittingBaseProposal proposal;
                proposal.state = outer_previous_state;
                proposal.effective_damping = damping;
                proposal.step_norm = *step_norm;
                for (std::size_t atom_position = 0; atom_position < key.size(); atom_position++)
                {
                    const auto atom_index{ key.at(atom_position) };
                    proposal.state.at(atom_index).mdpde =
                        GaussianModel3DWithUncertainty{
                            candidate_model_list->at(atom_position),
                            raw_state.at(atom_index).mdpde
                                .GetStandardDeviationModel()
                        };
                }
                return LocalFittingBaseProposalBuildResult{
                    std::move(proposal),
                    LocalFittingPreObjectiveFailureReason::None,
                    *step_norm
                };
            }
            if (step_norm.has_value()) attempted_step_norm = *step_norm;
        }
        damping *= 0.5;
    }
    return LocalFittingBaseProposalBuildResult{
        std::nullopt,
        LocalFittingPreObjectiveFailureReason::NoCandidateWithinTrustRegion,
        attempted_step_norm
    };
}

bool HasMaterialLocalFittingJointPolishChange(
    const std::vector<GaussianModel3D> & candidate_model_list,
    const std::vector<GaussianModel3D> & seed_model_list)
{
    if (candidate_model_list.size() != seed_model_list.size()) return false;
    for (std::size_t atom_position = 0; atom_position < candidate_model_list.size(); atom_position++)
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
        group_key_by_atom_position.emplace_back(data_internal::GetGroupKey(context.at(atom_index).atom));
        base_model_list.emplace_back(base_state.at(atom_index).mdpde.GetModel());
    }
    const auto parameterization{
        detail::BuildLocalFittingJointPolishParameterization(group_key_by_atom_position, base_model_list)
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
    if (!seed_step_norm.has_value() || *seed_step_norm > trust_region_radius + trust_region_tolerance)
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
            if (step_norm.has_value() && *step_norm <= trust_region_radius + trust_region_tolerance)
            {
                if (!HasMaterialLocalFittingJointPolishChange(*candidate_model_list, *seed_model_list))
                {
                    return std::nullopt;
                }

                LocalFittingJointPolishProposal proposal;
                proposal.state = base_state;
                proposal.effective_damping = damping;
                proposal.step_norm = *step_norm;
                for (std::size_t atom_position = 0; atom_position < key.size(); atom_position++)
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
                        detail::EncodeLocalFittingTransformedCoordinates(base_result.mdpde.GetModel())
                    };
                    const auto candidate_coordinates{
                        detail::EncodeLocalFittingTransformedCoordinates(candidate_model)
                    };
                    if (!base_coordinates.has_value() || !candidate_coordinates.has_value())
                    {
                        return std::nullopt;
                    }
                    if ((base_coordinates->array() != candidate_coordinates->array()).any())
                    {
                        proposal.changed_atom_index_list.emplace_back(atom_index);
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
            auto & result{ candidate_state.at(active_index) };
            result.mdpde = GaussianModel3DWithUncertainty{
                previous_state.at(active_index).mdpde.GetModel(),
                uncertainty_state.at(active_index).mdpde.GetStandardDeviationModel()
            };
            continue;
        }
        const auto candidate_model{
            detail::DecodeLocalFittingTransformedCoordinates(candidate_transformed_estimation)
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

std::optional<LocalFittingState> BuildLocalFittingBacktrackedState(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & previous_state,
    const LocalFittingState & endpoint_state,
    const std::vector<std::size_t> & active_index_list,
    double factor)
{
    if (!std::isfinite(factor) || factor < 0.0 || factor > 1.0)
    {
        throw std::invalid_argument(
            "Local fitting objective backtracking factor must be in [0, 1].");
    }
    std::vector<GroupKey> group_key_by_atom_position;
    std::vector<GaussianModel3D> previous_model_list;
    std::vector<GaussianModel3D> endpoint_model_list;
    group_key_by_atom_position.reserve(active_index_list.size());
    previous_model_list.reserve(active_index_list.size());
    endpoint_model_list.reserve(active_index_list.size());
    for (const auto atom_index : active_index_list)
    {
        group_key_by_atom_position.emplace_back(data_internal::GetGroupKey(context.at(atom_index).atom));
        previous_model_list.emplace_back(previous_state.at(atom_index).mdpde.GetModel());
        endpoint_model_list.emplace_back(endpoint_state.at(atom_index).mdpde.GetModel());
    }
    const auto previous_shared_offset_model_list{
        detail::BuildLocalFittingGroupMedianModelList(group_key_by_atom_position, previous_model_list)
    };
    const auto endpoint_shared_offset_model_list{
        detail::BuildLocalFittingGroupMedianModelList(group_key_by_atom_position, endpoint_model_list)
    };
    const auto candidate_model_list{
        detail::BuildLocalFittingSharedOffsetDampedModelList(
            previous_model_list,
            endpoint_model_list,
            previous_shared_offset_model_list,
            endpoint_shared_offset_model_list,
            factor)
    };
    if (!candidate_model_list.has_value()) return std::nullopt;

    auto candidate_state{ previous_state };
    for (std::size_t atom_position = 0; atom_position < active_index_list.size(); atom_position++)
    {
        const auto atom_index{ active_index_list.at(atom_position) };
        candidate_state.at(atom_index).mdpde =
            GaussianModel3DWithUncertainty{
                candidate_model_list->at(atom_position),
                endpoint_state.at(atom_index).mdpde.GetStandardDeviationModel()
            };
    }
    return candidate_state;
}

double GetLocalFittingMaximumTransformedChange(
    const LocalFittingState & candidate_state,
    const LocalFittingState & previous_state,
    const std::vector<std::size_t> & active_index_list)
{
    const auto summary{
        SummarizeLocalFittingTransformedChanges(
            candidate_state,
            previous_state,
            active_index_list)
    };
    return summary.maximum_list.empty() ? 0.0 :
        *std::max_element(summary.maximum_list.begin(), summary.maximum_list.end());
}

LocalFittingClusterObjectiveState
BuildInitialLocalFittingClusterObjectiveState(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & previous_state,
    const LocalFittingClusterKey & key,
    const std::vector<LocalFittingObjectiveSampleRef> & objective_sample_ref_list,
    const LocalFittingObjectiveDomain & domain)
{
    LocalFittingClusterObjectiveState state;
    state.best_objective = EvaluateLocalFittingObjectiveContribution(
        context, previous_state, key, objective_sample_ref_list, domain);
    return state;
}

void ReconcileLocalFittingClusterObjectiveState(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & previous_state,
    const detail::LocalFittingCouplingPartition & partition,
    const LocalFittingObjectiveDomain & domain,
    LocalFittingClusterObjectiveStateMap & state_by_key)
{
    LocalFittingClusterObjectiveStateMap next_state_by_key;
    for (const auto & [key, objective_sample_ref_list] : partition.sample_id_list_by_key)
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
                objective_sample_ref_list,
                domain));
    }
    state_by_key = std::move(next_state_by_key);
}

bool TryCommitLocalFittingClusterCandidate(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & candidate_state,
    const LocalFittingState & previous_state,
    const LocalFittingClusterKey & key,
    const std::vector<LocalFittingObjectiveSampleRef> & objective_sample_ref_list,
    bool requires_strict_improvement,
    const LocalFittingObjectiveDomain & domain,
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
    const auto domain_iter{ domain.cluster_by_key.find(key) };
    if (domain_iter != domain.cluster_by_key.end())
    {
        diagnostic.fit_sample_count = domain_iter->second.fit_sample_ref_list.size();
        diagnostic.tail_sample_count = domain_iter->second.tail_sample_ref_list.size();
        if (domain_iter->second.scale.has_value())
        {
            diagnostic.fit_scale = domain_iter->second.scale->fit;
            diagnostic.tail_scale = domain_iter->second.scale->tail;
        }
    }
    auto & state{ cluster_objective_state.at(key) };
    diagnostic.candidate_objective =
        EvaluateLocalFittingObjectiveContribution(
            context,
            candidate_state,
            key,
            objective_sample_ref_list,
            domain);
    diagnostic.previous_objective =
        EvaluateLocalFittingObjectiveContribution(
            context,
            previous_state,
            key,
            objective_sample_ref_list,
            domain);
    diagnostic.best_objective = state.best_objective;

    const auto is_objective_deteriorated = [](
        const std::optional<detail::LocalFittingObjectiveBreakdown> & candidate,
        const std::optional<double> & reference)
    {
        if (!reference.has_value()) return false;
        if (!candidate.has_value()) return true;
        return candidate->total_objective > *reference +
            detail::CalculateLocalFittingObjectiveTolerance(
                *reference,
                kLocalFittingObjectiveProgressTolerance);
    };
    if (!diagnostic.candidate_objective.has_value() || !diagnostic.previous_objective.has_value())
    {
        return false;
    }
    diagnostic.rejected_by_previous = is_objective_deteriorated(
        diagnostic.candidate_objective,
        diagnostic.previous_objective->total_objective);
    diagnostic.rejected_by_best = is_objective_deteriorated(
        diagnostic.candidate_objective,
        state.best_objective.has_value() ?
            std::optional<double>{ state.best_objective->total_objective } :
            std::nullopt);
    if (diagnostic.rejected_by_previous || diagnostic.rejected_by_best)
    {
        return false;
    }
    if (requires_strict_improvement &&
        (!diagnostic.candidate_objective.has_value() ||
            !diagnostic.previous_objective.has_value() ||
            !detail::IsBetterLocalFittingAuditObjective(
                diagnostic.candidate_objective->total_objective,
                diagnostic.previous_objective->total_objective,
                kLocalFittingObjectiveStrictTolerance)))
    {
        return false;
    }

    auto is_better_than_best{ !state.best_objective.has_value() };
    if (state.best_objective.has_value())
    {
        const auto candidate_objective_value{
            diagnostic.candidate_objective->total_objective
        };
        const auto best_objective_value{ state.best_objective->total_objective };
        if (detail::IsBetterLocalFittingAuditObjective(
                candidate_objective_value,
                best_objective_value,
                kLocalFittingObjectiveStrictTolerance))
        {
            is_better_than_best = true;
        }
        else if (detail::IsBetterLocalFittingAuditObjective(
                     best_objective_value,
                     candidate_objective_value,
                     kLocalFittingObjectiveStrictTolerance))
        {
            is_better_than_best = false;
        }
        else
        {
            is_better_than_best = maximum_transformed_change < state.best_maximum_transformed_change;
        }
    }
    if (is_better_than_best)
    {
        state.best_objective = diagnostic.candidate_objective;
        state.best_maximum_transformed_change = maximum_transformed_change;
    }
    return true;
}

bool ShouldGrowLocalFittingTrustRegion(const LocalFittingObjectiveAttemptDiagnostic & diagnostic)
{
    return diagnostic.candidate_objective.has_value() &&
        diagnostic.previous_objective.has_value() &&
        diagnostic.trust_region_step_norm >= kLocalFittingTrustRegionBoundaryRatio * diagnostic.trust_region_radius &&
        detail::IsBetterLocalFittingAuditObjective(
            diagnostic.candidate_objective->total_objective,
            diagnostic.previous_objective->total_objective,
            kLocalFittingObjectiveStrictTolerance);
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
    const std::vector<char> & rollback_atom_mask,
    const std::vector<double> & ridge_multiplier_list,
    const std::vector<LocalFittingClusterKey> & unchanged_state_exhausted_key_list,
    const LocalFittingObjectiveDomain & objective_domain,
    LocalFittingClusterObjectiveStateMap & cluster_objective_state,
    const detail::LocalFittingTrustRegionStateSet & trust_region_state)
{
    LocalFittingCandidateSelection selection;
    selection.assembled_state = previous_state;
    selection.assembled_polish_provenance = previous_polish_provenance;
    for (const auto & [key, objective_sample_ref_list] : partition.sample_id_list_by_key)
    {
        const auto is_polish_eligible{
            std::find(
                polish_eligible_key_list.begin(),
                polish_eligible_key_list.end(),
                key) != polish_eligible_key_list.end()
        };
        if (is_polish_eligible) selection.polish_progress.eligible_count++;
        if (std::find(
                unchanged_state_exhausted_key_list.begin(),
                unchanged_state_exhausted_key_list.end(),
                key) != unchanged_state_exhausted_key_list.end())
        {
            if (is_polish_eligible) selection.polish_progress.skipped_count++;
            LocalFittingObjectiveAttemptDiagnostic diagnostic;
            diagnostic.trust_region_radius = trust_region_state.GetRadius(key);
            diagnostic.backtracking_exhausted = true;
            selection.rejected_key_list.emplace_back(key);
            selection.backtracking_exhausted_key_list.emplace_back(key);
            selection.rejected_cluster_diagnostic_list.emplace_back(
                LocalFittingRejectedClusterDiagnostic{
                    key,
                    std::move(diagnostic)
                });
            continue;
        }
        const auto trust_region_radius{ trust_region_state.GetRadius(key) };
        const auto contains_suspicious_atom{
            std::any_of(
                key.begin(),
                key.end(),
                [&](std::size_t atom_index)
                {
                    return rollback_atom_mask.at(atom_index) != 0;
                })
        };
        LocalFittingObjectiveAttemptDiagnostic base_diagnostic;
        base_diagnostic.trust_region_radius = trust_region_radius;
        std::optional<LocalFittingBaseProposal> base_proposal;
        if (!contains_suspicious_atom)
        {
            auto proposal_result{
                BuildLocalFittingSharedOffsetBaseProposal(
                    context,
                    previous_state,
                    raw_state,
                    key,
                    trust_region_radius)
            };
            base_diagnostic.pre_objective_failure_reason = proposal_result.failure_reason;
            if (proposal_result.attempted_step_norm.has_value())
            {
                base_diagnostic.pre_objective_attempted_step_norm = proposal_result.attempted_step_norm;
                base_diagnostic.trust_region_step_norm = *proposal_result.attempted_step_norm;
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
                previous_cluster_estimation_list.emplace_back(previous_transformed_estimation_list.at(atom_index));
                raw_cluster_estimation_list.emplace_back(raw_transformed_estimation_list.at(atom_index));
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
                base_diagnostic.pre_objective_failure_reason = LocalFittingPreObjectiveFailureReason::InvalidModel;
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
        base_diagnostic.effective_damping = base_proposal->effective_damping;
        base_diagnostic.trust_region_step_norm = base_proposal->step_norm;
        base_diagnostic.backtracking_trial_count = 1;
        base_diagnostic.accepted_backtracking_factor = 1.0;
        auto & base_state{ base_proposal->state };
        auto accepted_base_candidate{ TryCommitLocalFittingClusterCandidate(
                context,
                base_state,
                previous_state,
                key,
                objective_sample_ref_list,
                false,
                objective_domain,
                cluster_objective_state,
                base_diagnostic) };
        auto accepted_by_backtracking{ false };
        if (!accepted_base_candidate)
        {
            base_diagnostic.accepted_backtracking_factor.reset();
            const auto endpoint_state{ base_state };
            const auto endpoint_effective_damping{ base_proposal->effective_damping };
            const auto endpoint_step_norm{ base_proposal->step_norm };
            double factor{ 0.5 };
            while (factor >= std::numeric_limits<double>::epsilon())
            {
                auto backtracked_state{
                    BuildLocalFittingBacktrackedState(
                        context,
                        previous_state,
                        endpoint_state,
                        key,
                        factor)
                };
                if (!backtracked_state.has_value())
                {
                    base_diagnostic.is_invalid_model = true;
                    break;
                }
                if (GetLocalFittingMaximumTransformedChange(
                        *backtracked_state,
                        previous_state,
                        key) < kLocalFittingTransformedChangeTolerance)
                {
                    base_diagnostic.backtracking_exhausted = true;
                    break;
                }

                LocalFittingObjectiveAttemptDiagnostic trial_diagnostic;
                trial_diagnostic.effective_damping = endpoint_effective_damping * factor;
                trial_diagnostic.trust_region_radius = trust_region_radius;
                trial_diagnostic.trust_region_step_norm = endpoint_step_norm * factor;
                trial_diagnostic.backtracking_trial_count = base_diagnostic.backtracking_trial_count + 1;
                if (TryCommitLocalFittingClusterCandidate(
                        context,
                        *backtracked_state,
                        previous_state,
                        key,
                        objective_sample_ref_list,
                        false,
                        objective_domain,
                        cluster_objective_state,
                        trial_diagnostic))
                {
                    trial_diagnostic.accepted_backtracking_factor = factor;
                    base_state = std::move(*backtracked_state);
                    base_diagnostic = std::move(trial_diagnostic);
                    accepted_base_candidate = true;
                    accepted_by_backtracking = true;
                    break;
                }
                base_diagnostic = std::move(trial_diagnostic);
                factor *= 0.5;
            }
        }
        if (!accepted_base_candidate)
        {
            if (is_polish_eligible) selection.polish_progress.skipped_count++;
            selection.rejected_key_list.emplace_back(key);
            if (base_diagnostic.backtracking_exhausted)
            {
                selection.backtracking_exhausted_key_list.emplace_back(key);
            }
            selection.rejected_cluster_diagnostic_list.emplace_back(
                LocalFittingRejectedClusterDiagnostic{
                    key,
                    std::move(base_diagnostic)
                });
            continue;
        }
        selection.accepted_key_list.emplace_back(key);
        selection.accepted_cluster_diagnostic_list.emplace_back(
            LocalFittingRejectedClusterDiagnostic{ key, base_diagnostic });
        if (!accepted_by_backtracking &&
            ShouldGrowLocalFittingTrustRegion(base_diagnostic))
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
                (base_transformed->array() == previous_transformed_estimation_list.at(atom_index).array()).all())
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
        polish_diagnostic.effective_damping = polished_candidate->effective_damping;
        polish_diagnostic.trust_region_radius = trust_region_radius;
        polish_diagnostic.trust_region_step_norm = polished_candidate->step_norm;
        if (!TryCommitLocalFittingClusterCandidate(
                context,
                polished_candidate->state,
                base_state,
                key,
                objective_sample_ref_list,
                true,
                objective_domain,
                cluster_objective_state,
                polish_diagnostic))
        {
            selection.polish_progress.rejected_count++;
            continue;
        }
        selection.polish_progress.accepted_count++;
        for (const auto active_index : key)
        {
            selection.assembled_state.at(active_index) = polished_candidate->state.at(active_index);
        }
        for (const auto atom_index : polished_candidate->changed_atom_index_list)
        {
            selection.assembled_polish_provenance.at(atom_index) = 1;
        }
        if (!accepted_by_backtracking &&
            ShouldGrowLocalFittingTrustRegion(polish_diagnostic) &&
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

LocalFittingPolishProvenance BuildBacktrackedLocalFittingPolishProvenance(
    const LocalFittingState & previous_state,
    const LocalFittingPolishProvenance & previous_provenance,
    const LocalFittingPolishProvenance & endpoint_provenance,
    const LocalFittingState & candidate_state,
    const std::vector<std::size_t> & changed_atom_index_list)
{
    auto provenance{ previous_provenance };
    for (const auto atom_index : changed_atom_index_list)
    {
        const auto change{
            detail::CalculateLocalFittingTransformedChange(
                candidate_state.at(atom_index).mdpde.GetModel(),
                previous_state.at(atom_index).mdpde.GetModel())
        };
        const auto has_material_change{
            std::any_of(
                change.value_list.begin(),
                change.value_list.end(),
                [](double value)
                {
                    return std::isfinite(value) &&
                        value >= kLocalFittingTransformedChangeTolerance;
                })
        };
        if (has_material_change)
        {
            provenance.at(atom_index) = endpoint_provenance.at(atom_index);
        }
    }
    return provenance;
}

bool TryBacktrackLocalFittingCombinedCandidate(
    const SecondStageLocalFittingContext & context,
    const detail::LocalFittingCouplingPartition & partition,
    const LocalFittingState & previous_state,
    const LocalFittingPolishProvenance & previous_polish_provenance,
    const LocalFittingObjectiveDomain & objective_domain,
    const LocalFittingBestAuditState & best_audit_state,
    const LocalFittingClusterObjectiveStateMap & committed_objective_state,
    LocalFittingClusterObjectiveStateMap & working_objective_state,
    LocalFittingCandidateSelection & selection)
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

    selection.combined_backtracking_trial_count = 1;
    double factor{ 0.5 };
    while (factor >= std::numeric_limits<double>::epsilon())
    {
        auto candidate_state{
            BuildLocalFittingBacktrackedState(
                context,
                previous_state,
                endpoint_state,
                changed_atom_index_list,
                factor)
        };
        if (!candidate_state.has_value()) return false;
        if (GetLocalFittingMaximumTransformedChange(
                *candidate_state,
                previous_state,
                changed_atom_index_list) < kLocalFittingTransformedChangeTolerance)
        {
            selection.combined_backtracking_exhausted = true;
            return false;
        }

        selection.combined_backtracking_trial_count++;
        auto trial_objective_state{ committed_objective_state };
        auto local_criteria_accepted{ true };
        for (const auto & key : accepted_key_list)
        {
            const auto sample_iter{ partition.sample_id_list_by_key.find(key) };
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
                    *candidate_state,
                    previous_state,
                    key,
                    sample_iter->second,
                    false,
                    objective_domain,
                    trial_objective_state,
                    diagnostic))
            {
                local_criteria_accepted = false;
                break;
            }
        }
        if (local_criteria_accepted &&
            IsLocalFittingCombinedObjectiveAcceptable(
                context,
                *candidate_state,
                previous_state,
                objective_domain,
                best_audit_state))
        {
            selection.assembled_state = std::move(*candidate_state);
            selection.assembled_polish_provenance =
                BuildBacktrackedLocalFittingPolishProvenance(
                    previous_state,
                    previous_polish_provenance,
                    endpoint_provenance,
                    selection.assembled_state,
                    changed_atom_index_list);
            selection.combined_backtracking_factor = factor;
            selection.grow_trust_region_key_list.clear();
            working_objective_state = std::move(trial_objective_state);
            return true;
        }
        factor *= 0.5;
    }
    selection.combined_backtracking_exhausted = true;
    return false;
}

std::optional<LocalAtomRefitResult> FitAtomWithJointOffsetFallback(
    const SecondStageLocalFittingContext & context,
    std::size_t atom_index,
    const LocalGaussianResult & previous_result,
    const FittedGaussianSnapshot & refit_model_snapshot,
    const FitOptions & options)
{
    auto adjusted_sampling_entries{
        BuildSecondStageAdjustedSamples(context, atom_index, refit_model_snapshot)
    };
    const auto & offset_model{ refit_model_snapshot.at(atom_index) };
    const auto & previous_model{ previous_result.mdpde.GetModel() };
    const auto previous_baseline{
        BuildPreviousSuspiciousProfileBaseline(
            adjusted_sampling_entries,
            previous_model,
            options)
    };
    const auto is_post_refit_candidate_acceptable = [&](const GaussianModel3D & model)
    {
        const auto reason{ EvaluateSuspiciousGaussianUpdate(
                adjusted_sampling_entries,
                previous_model,
                model,
                options,
                previous_baseline,
                true) };
        return reason == detail::SuspiciousGaussianReason::None;
    };
    const auto is_offset_only_fallback_acceptable =
        [&](const GaussianModel3D & model)
    {
        const auto reason{ EvaluateSuspiciousGaussianUpdate(
                adjusted_sampling_entries,
                previous_model,
                model,
                options,
                previous_baseline,
                false) };
        return reason == detail::SuspiciousGaussianReason::None;
    };
    try
    {
        auto candidate_result{
            EstimateLocalGaussian(
                adjusted_sampling_entries,
                context.at(atom_index).alpha_r,
                options,
                offset_model)
        };
        const auto candidate_model{ candidate_result.mdpde.GetModel() };
        if (is_post_refit_candidate_acceptable(candidate_model))
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
    if (!is_offset_only_fallback_acceptable(result.mdpde.GetModel()))
    {
        return std::nullopt;
    }
    return LocalAtomRefitResult{ std::move(result), false };
}

void ExpandPostRefitRollbackClusters(
    const std::vector<LocalFittingClusterKey> & cluster_key_list,
    const std::vector<std::size_t> & seed_atom_index_list,
    std::vector<char> & rollback_mask)
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
            rollback_mask.at(atom_index) = 1;
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
    std::map<LocalFittingClusterKey, JointOffsetSolveResult> joint_offset_result_by_key;
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
                current_snapshot.at(atom_index).WithOffset(result.offset(static_cast<Eigen::Index>(i)));
        }
        health_by_key.emplace(key, LocalFittingClusterHealth{ result.status });
    }

    auto iteration_state{ previous_state };
    std::vector<char> rollback_atom_mask(selected_atom_size, 0);
    std::vector<GroupKey> group_key_by_atom_index;
    group_key_by_atom_index.reserve(context.size());
    for (const auto & atom_context : context)
    {
        group_key_by_atom_index.emplace_back(data_internal::GetGroupKey(atom_context.atom));
    }
    for (const auto & [key, result] : joint_offset_result_by_key)
    {
        static_cast<void>(result);
        std::vector<GroupKey> group_key_by_position;
        std::vector<char> suspicious_seed_mask(key.size(), 0);
        group_key_by_position.reserve(key.size());
        for (std::size_t position = 0; position < key.size(); position++)
        {
            const auto atom_index{ key.at(position) };
            group_key_by_position.emplace_back(group_key_by_atom_index.at(atom_index));
            if (detail::EvaluateSuspiciousOffsetUpdate(
                    context.at(atom_index).raw_sampling_entries,
                    previous_state.at(atom_index).mdpde.GetModel(),
                    current_snapshot.at(atom_index),
                    options) != detail::SuspiciousGaussianReason::None)
            {
                suspicious_seed_mask.at(position) = 1;
            }
        }
        const auto cluster_rollback_mask{
            detail::ExpandSuspiciousSharedOffsetGroups(group_key_by_position, suspicious_seed_mask)
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
        current_snapshot.at(atom_index) = previous_state.at(atom_index).mdpde.GetModel();
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
            if (rollback_atom_mask.at(atom_index) != 0) continue;

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

void ApplyLocalFittingState(
    ModelObject & model_object,
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & iteration_state)
{
    const auto fitted_gaussian_snapshot{ BuildFittedGaussianSnapshot(iteration_state) };
    auto analysis{ model_object.EditAnalysis() };
    for (std::size_t i = 0; i < context.size(); i++)
    {
        auto local_editor{
            analysis.EnsureAtomLocalPotential(*context.at(i).atom)
        };
        local_editor.SetGaussianResult(iteration_state.at(i));
        local_editor.SetPeelingSamplingEntries(
            BuildSecondStageAdjustedSamples(context, i, fitted_gaussian_snapshot));
    }
}

LocalFittingOffsetStats SummarizeLocalFittingOffsets(const LocalFittingState & state)
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

    stats.median_absolute_offset = array_helper::ComputeMedian(absolute_offset_list);
    stats.percentile_absolute_offset = array_helper::ComputePercentile(absolute_offset_list, kLocalFittingChangePercentile);
    stats.maximum_absolute_offset = *std::max_element(absolute_offset_list.begin(), absolute_offset_list.end());
    return stats;
}

void AppendLocalFittingOffsetSummary(std::ostringstream & stream, const LocalFittingOffsetStats & stats)
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

void AppendLocalFittingAuditSummary(std::ostringstream & stream, const LocalFittingAuditedState & audited_state)
{
    const auto & objective{ audited_state.objective };
    stream << "; audit best source = ";
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
        << ", fixed audit objective fit/tail-weighted/offset/total = "
        << objective.fit_range_residual_objective << "/"
        << objective.tail_validation_penalty << "/"
        << objective.offset_plausibility_penalty << "/"
        << objective.total_objective
        << ", tail raw/weight = "
        << objective.tail_validation_loss << "/"
        << kLocalFittingTailValidationWeight;
}

void AppendLocalFittingTerminalSummary(std::ostringstream & stream, const LocalFittingTerminalSummary & summary)
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
        << breakdown->fit_range_residual_objective << "/"
        << breakdown->tail_validation_penalty << "/"
        << breakdown->offset_plausibility_penalty << "/"
        << breakdown->total_objective;
}

std::string_view GetLocalFittingPreObjectiveFailureReasonText(
    LocalFittingPreObjectiveFailureReason reason)
{
    switch (reason)
    {
    case LocalFittingPreObjectiveFailureReason::None:
        return "none";
    case LocalFittingPreObjectiveFailureReason::InvalidModel:
        return "invalid-model";
    case LocalFittingPreObjectiveFailureReason::PreviousSharedOffsetProjectionOutsideTrustRegion:
        return "previous-shared-offset-projection-outside-trust-region";
    case LocalFittingPreObjectiveFailureReason::NoCandidateWithinTrustRegion:
        return "no-candidate-within-trust-region";
    }
    return "unknown";
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
            << "Rejected local fitting cluster diagnostics: atoms = "
            << cluster_diagnostic.key.size()
            << ", key first/last = "
            << cluster_diagnostic.key.front() << "/" << cluster_diagnostic.key.back()
            << ", breakdown order = fit/tail-weighted/offset/total";
        Logger::Log(LogLevel::Debug, header.str());

        const auto & diagnostic{ cluster_diagnostic.attempt };
        std::ostringstream message;
        message
            << std::scientific << std::setprecision(2)
            << "  fixed-point effective damping = "
            << diagnostic.effective_damping
            << ", trust radius/step norm = "
            << diagnostic.trust_region_radius << "/";

        if (diagnostic.pre_objective_failure_reason !=
            LocalFittingPreObjectiveFailureReason::None)
        {
            if (diagnostic.pre_objective_attempted_step_norm.has_value())
            {
                message << *diagnostic.pre_objective_attempted_step_norm;
            }
            else
            {
                message << "unavailable";
            }
            message
                << ", status = "
                << GetLocalFittingPreObjectiveFailureReasonText(
                    diagnostic.pre_objective_failure_reason)
                << ", objective = not-evaluated";
            Logger::Log(LogLevel::Debug, message.str());
            continue;
        }

        message << diagnostic.trust_region_step_norm;

        if (diagnostic.is_invalid_model)
        {
            message << ", status = invalid-model";
            Logger::Log(LogLevel::Debug, message.str());
            continue;
        }

        message << ", fit/tail scales = ";
        if (diagnostic.fit_scale.has_value())
        {
            message << *diagnostic.fit_scale;
        }
        else
        {
            message << "unavailable";
        }
        message << "/";
        if (diagnostic.tail_scale.has_value())
        {
            message << *diagnostic.tail_scale;
        }
        else
        {
            message << "empty";
        }
        message
            << ", fit/tail samples = "
            << diagnostic.fit_sample_count << "/"
            << diagnostic.tail_sample_count
            << ", tail raw/weight = ";
        if (diagnostic.candidate_objective.has_value())
        {
            message << diagnostic.candidate_objective->tail_validation_loss;
        }
        else
        {
            message << "unavailable";
        }
        message << "/" << kLocalFittingTailValidationWeight;
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
        message
            << ", backtracking trials/factor/exhausted = "
            << diagnostic.backtracking_trial_count << "/";
        if (diagnostic.accepted_backtracking_factor.has_value())
        {
            message << *diagnostic.accepted_backtracking_factor;
        }
        else
        {
            message << "-";
        }
        message << "/" << (diagnostic.backtracking_exhausted ? "yes" : "no");
        Logger::Log(LogLevel::Debug, message.str());
    }
}

std::string_view GetLocalFittingAllRejectedResolutionText(
    detail::LocalFittingAllRejectedResolution resolution)
{
    switch (resolution)
    {
    case detail::LocalFittingAllRejectedResolution::Retry:
        return "retry";
    case detail::LocalFittingAllRejectedResolution::MaximumIterations:
        return "maximum-iterations";
    case detail::LocalFittingAllRejectedResolution::BacktrackingExhausted:
        return "all-rejected-backtracking-exhausted";
    case detail::LocalFittingAllRejectedResolution::MinimumRadius:
        return "all-rejected-minimum-radius";
    case detail::LocalFittingAllRejectedResolution::NoRetryProgress:
        return "all-rejected-no-retry-progress";
    }
    return "all-rejected-no-retry-progress";
}

void LogLocalFittingAllRejectedResolution(
    const FitOptions & options,
    const detail::LocalFittingRejectedClusterPartition & partition,
    const detail::LocalFittingTrustRegionRadiusUpdate & radius_update,
    detail::LocalFittingAllRejectedResolution resolution)
{
    if (options.quiet_mode || Logger::GetLogLevel() < LogLevel::Debug)
    {
        return;
    }

    Logger::FinishProgressLine();
    std::ostringstream message;
    message
        << "All-rejected local fitting resolution: outcome = "
        << GetLocalFittingAllRejectedResolutionText(resolution)
        << ", exhausted/retryable/radius-changed/radius-saturated = "
        << partition.exhausted_key_list.size() << "/"
        << partition.retryable_key_list.size() << "/"
        << radius_update.changed_key_list.size() << "/"
        << radius_update.saturated_key_list.size() << ".";
    Logger::Log(LogLevel::Debug, message.str());
}

void LogAcceptedLocalFittingBacktrackingDiagnostics(
    const FitOptions & options,
    const LocalFittingCandidateSelection & selection)
{
    if (options.quiet_mode || Logger::GetLogLevel() < LogLevel::Debug)
    {
        return;
    }
    const auto has_local_backtracking{
        std::any_of(
            selection.accepted_cluster_diagnostic_list.begin(),
            selection.accepted_cluster_diagnostic_list.end(),
            [&](const LocalFittingRejectedClusterDiagnostic & diagnostic)
            {
                return diagnostic.attempt.backtracking_trial_count > 1 &&
                    std::find(
                        selection.accepted_key_list.begin(),
                        selection.accepted_key_list.end(),
                        diagnostic.key) != selection.accepted_key_list.end();
            })
    };
    if (!has_local_backtracking &&
        selection.combined_backtracking_trial_count <= 1)
    {
        return;
    }
    Logger::FinishProgressLine();
    for (const auto & cluster_diagnostic :
        selection.accepted_cluster_diagnostic_list)
    {
        const auto & diagnostic{ cluster_diagnostic.attempt };
        if (diagnostic.backtracking_trial_count <= 1 ||
            std::find(
                selection.accepted_key_list.begin(),
                selection.accepted_key_list.end(),
                cluster_diagnostic.key) == selection.accepted_key_list.end())
        {
            continue;
        }
        std::ostringstream message;
        message
            << "Accepted local fitting objective backtracking: atoms = "
            << cluster_diagnostic.key.size()
            << ", key first/last = "
            << cluster_diagnostic.key.front() << "/"
            << cluster_diagnostic.key.back()
            << ", trials/factor = "
            << diagnostic.backtracking_trial_count << "/";
        if (diagnostic.accepted_backtracking_factor.has_value())
        {
            message << *diagnostic.accepted_backtracking_factor;
        }
        else
        {
            message << "-";
        }
        message << ", fixed fit/tail scales = ";
        if (diagnostic.fit_scale.has_value())
        {
            message << *diagnostic.fit_scale;
        }
        else
        {
            message << "unavailable";
        }
        message << "/";
        if (diagnostic.tail_scale.has_value())
        {
            message << *diagnostic.tail_scale;
        }
        else
        {
            message << "empty";
        }
        message << ".";
        Logger::Log(LogLevel::Debug, message.str());
    }
    if (selection.combined_backtracking_trial_count <= 1) return;
    std::ostringstream message;
    message
        << "Combined-objective backtracking: trials/factor/exhausted = "
        << selection.combined_backtracking_trial_count << "/";
    if (selection.combined_backtracking_factor.has_value())
    {
        message << *selection.combined_backtracking_factor;
    }
    else
    {
        message << "-";
    }
    message << "/"
        << (selection.combined_backtracking_exhausted ? "yes" : "no")
        << ".";
    Logger::Log(LogLevel::Debug, message.str());
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
    const auto maximum_iteration_text{ std::to_string(kLocalFittingMaximumIterations) };
    const auto maximum_atom_text{ std::to_string(atom_size) };
    const auto maximum_change_text{
        FormatLocalFittingProgressMaximum(std::numeric_limits<double>::max())
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
        FormatLocalFittingProgressMaximum(progress.accepted_maximum_transformed_change) + "/" +
            FormatLocalFittingProgressMaximum(progress.raw_maximum_transformed_change)
    };
    Logger::ProgressLine(FormatLocalFittingProgressRow(column_widths, cell_list));
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
    detail::LocalFittingFinalStateSource final_state_source,
    const LocalFittingAuditedState * applied_audit_state,
    const LocalFittingTerminalSummary & terminal_summary,
    const LocalFittingOffsetStats & applied_offset_stats)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message << "Reached maximum iteration size";
    AppendLocalFittingTerminalSummary(warning_message, terminal_summary);
    if (final_state_source == detail::LocalFittingFinalStateSource::BestAudit &&
        applied_audit_state != nullptr)
    {
        warning_message << "; applying best validated audit state";
        AppendLocalFittingAuditSummary(warning_message, *applied_audit_state);
    }
    else if (final_state_source == detail::LocalFittingFinalStateSource::LatestValidated)
    {
        warning_message << "; applying latest validated state";
    }
    else
    {
        warning_message << "; no validated state is available";
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
    std::optional<bool> final_uses_polish,
    detail::LocalFittingFinalStateSource final_state_source)
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
    message
        << ", final_state_source="
        << GetLocalFittingFinalStateSourceText(final_state_source)
        << ".";
    Logger::Log(LogLevel::Info, message.str());
}

} // namespace

detail::SuspiciousGaussianReason detail::EvaluateSuspiciousOffsetUpdate(
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
        previous_model,
        candidate_model,
        options,
        previous_baseline,
        false);
}

detail::SuspiciousGaussianReason detail::EvaluateSuspiciousPostRefitUpdate(
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
        previous_model,
        candidate_model,
        options,
        previous_baseline,
        true);
}

std::optional<double> detail::CalculateLocalFittingPeelingRatio(
    const LocalPotentialSampleList & raw_sampling_entries,
    const LocalPotentialSampleList & peeling_sampling_entries,
    bool peeling_applied)
{
    if (!peeling_applied
        || raw_sampling_entries.empty()
        || peeling_sampling_entries.empty())
    {
        return std::nullopt;
    }

    double raw_sum{ 0.0 };
    for (const auto & sample : raw_sampling_entries)
    {
        raw_sum += static_cast<double>(sample.response);
    }
    double peeling_sum{ 0.0 };
    for (const auto & sample : peeling_sampling_entries)
    {
        peeling_sum += static_cast<double>(sample.response);
    }
    if (!std::isfinite(raw_sum) || !std::isfinite(peeling_sum) || raw_sum == 0.0)
    {
        return std::nullopt;
    }

    const auto ratio{ (raw_sum - peeling_sum) / raw_sum };
    return std::isfinite(ratio) ? std::optional<double>{ ratio } : std::nullopt;
}

std::vector<char> detail::ExpandSuspiciousSharedOffsetGroups(
    const std::vector<GroupKey> & group_key_by_position,
    const std::vector<char> & suspicious_seed_mask)
{
    if (group_key_by_position.size() != suspicious_seed_mask.size())
    {
        throw std::invalid_argument("Suspicious shared-offset group input sizes are inconsistent.");
    }

    std::map<GroupKey, bool> has_suspicious_seed_by_group;
    for (std::size_t position = 0; position < group_key_by_position.size(); position++)
    {
        if (suspicious_seed_mask.at(position) == 0) continue;
        has_suspicious_seed_by_group[group_key_by_position.at(position)] = true;
    }

    std::vector<char> rollback_mask(group_key_by_position.size(), 0);
    for (std::size_t position = 0; position < group_key_by_position.size(); position++)
    {
        const auto seed_iter{
            has_suspicious_seed_by_group.find(group_key_by_position.at(position))
        };
        if (seed_iter != has_suspicious_seed_by_group.end())
        {
            rollback_mask.at(position) = 1;
        }
    }
    return rollback_mask;
}

SecondStageLocalFittingDiagnostics RunSecondStageLocalFitting(
    ModelObject & model_object,
    const FitOptions & options)
{
    auto context{ BuildSecondStageLocalFittingContext(model_object, options) };
    auto diagnostics{ BuildSecondStageLocalFittingDiagnostics(context) };
    const auto atom_size{ context.size() };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run 2nd-stage local atom fitting with iterations...");
        Logger::Log(
            LogLevel::Info,
            kApplyLocalFittingBestIteration ?
                "Second-stage best-iteration application: enabled." :
                "Second-stage best-iteration application: disabled.");
    }

    bool unselected_seed_failure{ false };
    auto initial_state_build_result{
        BuildInitialLocalFittingState(
            context,
            model_object.GetAnalysisView(),
            unselected_seed_failure)
    };
    if (!initial_state_build_result.has_value())
    {
        if (!options.quiet_mode)
        {
            Logger::Log(
                LogLevel::Warning,
                unselected_seed_failure ?
                    "Skip 2nd-stage local atom fitting because no valid Gaussian seed "
                    "is available for every unselected neighbor atom." :
                    "Skip 2nd-stage local atom fitting because no valid Gaussian seed "
                    "is available for every selected atom.");
            Logger::Log(
                LogLevel::Info,
                "Second-stage local fitting summary: accepted_iterations=0, "
                "best_iteration=unavailable, stop_reason=" +
                    std::string(unselected_seed_failure ?
                        "no-valid-unselected-neighbor-seed" :
                        "no-valid-seed") +
                    ", "
                "best_audit_objective=unavailable, final_uses_polish=unavailable, "
                "final_state_source=unavailable.");
        }
        return diagnostics;
    }
    LogSecondStageSeedSelections(initial_state_build_result->selection_record_list, options);
    LogUnselectedSecondStageSeedSelections(
        context,
        initial_state_build_result->unselected_selection_record_list,
        options);
    auto previous_state{ std::move(initial_state_build_result->state) };
    LocalFittingPolishProvenance previous_polish_provenance(atom_size, 0);
    const auto coupling_topology{
        BuildLocalFittingCouplingTopology(context, previous_state)
    };
    LogLocalFittingCouplingTopology(coupling_topology, options);
    std::vector<char> terminal_fallback_atom_mask(atom_size, 0);
    const auto initial_active_index_list{
        BuildEligibleLocalFittingActiveIndexList(terminal_fallback_atom_mask)
    };
    const auto initial_cluster_partition{
        detail::BuildLocalFittingCouplingPartition(coupling_topology, initial_active_index_list)
    };
    auto objective_domain{
        BuildLocalFittingObjectiveDomain(
            context,
            previous_state,
            initial_cluster_partition,
            options)
    };
    LogLocalFittingObjectiveDomain(objective_domain, options);
    auto best_audit_state{
        BuildInitialLocalFittingBestAuditState(
            context,
            previous_state,
            previous_polish_provenance,
            std::nullopt,
            objective_domain)
    };

    std::vector<char> rollback_atom_mask(atom_size, 0);
    PersistentTerminalFailureStateMap persistent_terminal_failure_state_by_key;
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
    const auto progress_column_widths{ BuildLocalFittingProgressColumnWidths(atom_size) };
    LogLocalFittingProgressHeader(options, progress_column_widths);

    std::size_t accepted_iteration_count{ 0 };
    std::size_t audit_patience_count{ 0 };
    std::vector<LocalFittingClusterKey> unchanged_state_exhausted_key_list;
    for (std::size_t iter = 0; iter < kLocalFittingMaximumIterations; iter++)
    {
        const auto active_index_list{
            BuildEligibleLocalFittingActiveIndexList(terminal_fallback_atom_mask)
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
            diagnostics.peeling_applied = true;
            LogSecondStageLocalFittingSummary(
                options,
                accepted_iteration_count,
                "terminal-isolation",
                best_audit_state,
                UsesLocalFittingPolish(previous_polish_provenance),
                detail::LocalFittingFinalStateSource::LatestValidated);
            RunGroupPotentialFitting(model_object, options, true);
            return diagnostics;
        }

        const auto cluster_partition{
            detail::BuildLocalFittingCouplingPartition(coupling_topology, active_index_list)
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
            objective_domain,
            cluster_objective_state);
        trust_region_state.Reconcile(cluster_key_list);

        std::vector<double> joint_offset_ridge_multiplier_list(atom_size, 1.0);
        for (std::size_t atom_index = 0; atom_index < atom_size; atom_index++)
        {
            if (rollback_atom_mask.at(atom_index) == 0) continue;
            joint_offset_ridge_multiplier_list.at(atom_index) = kSuspiciousJointOffsetRidgeMultiplier;
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
        rollback_atom_mask = std::move(iteration_result.rollback_atom_mask);
        const auto iteration_suspicious_atom_count{
            static_cast<std::size_t>(std::count_if(
                rollback_atom_mask.begin(),
                rollback_atom_mask.end(),
                [](char is_suspicious)
                {
                    return is_suspicious != 0;
                }))
        };
        const auto has_suspicious_offset_fallback{ iteration_suspicious_atom_count > 0 };

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
                        return rollback_atom_mask.at(atom_index) != 0;
                    })
            };
            if (!contains_suspicious_atom)
            {
                polish_eligible_key_list.emplace_back(key);
            }
        }
        const auto raw_state{ std::move(iteration_result.state) };
        const auto raw_fixed_point_change_summary{
            SummarizeLocalFittingTransformedChanges(raw_state, previous_state, active_index_list)
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
                rollback_atom_mask,
                joint_offset_ridge_multiplier_list,
                unchanged_state_exhausted_key_list,
                objective_domain,
                working_cluster_objective_state,
                trust_region_state)
        };

        const auto needs_combined_objective_guard{
            cluster_partition.boundary_sample_count > 0 && !selection.accepted_key_list.empty()
        };
        const auto combined_changed_key_list{ selection.accepted_key_list };
        auto combined_objective_accepted{
            !needs_combined_objective_guard ||
            IsLocalFittingCombinedObjectiveAcceptable(
                context,
                selection.assembled_state,
                previous_state,
                objective_domain,
                best_audit_state)
        };
        if (!combined_objective_accepted)
        {
            combined_objective_accepted =
                TryBacktrackLocalFittingCombinedCandidate(
                    context,
                    cluster_partition,
                    previous_state,
                    previous_polish_provenance,
                    objective_domain,
                    best_audit_state,
                    cluster_objective_state,
                    working_cluster_objective_state,
                    selection);
        }
        if (!combined_objective_accepted)
        {
            RejectLocalFittingCombinedCandidate(
                previous_state,
                previous_polish_provenance,
                cluster_key_list,
                selection);
            if (selection.combined_backtracking_exhausted)
            {
                for (const auto & key : combined_changed_key_list)
                {
                    if (std::find(
                            selection.backtracking_exhausted_key_list.begin(),
                            selection.backtracking_exhausted_key_list.end(),
                            key) == selection.backtracking_exhausted_key_list.end())
                    {
                        selection.backtracking_exhausted_key_list.emplace_back(key);
                    }
                }
            }
        }
        else
        {
            cluster_objective_state = std::move(working_cluster_objective_state);
        }
        LogAcceptedLocalFittingBacktrackingDiagnostics(options, selection);

        auto assembled_state{ std::move(selection.assembled_state) };
        auto assembled_polish_provenance{ std::move(selection.assembled_polish_provenance) };
        const auto terminal_failure_by_key{
            UpdatePersistentTerminalFailureState(
                selection.accepted_key_list,
                rollback_atom_mask,
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
        auto objective_domain_changed{ false };
        ApplyTerminalFallbackClusters(
            terminal_key_list,
            previous_state,
            previous_polish_provenance,
            terminal_fallback_atom_mask,
            assembled_state,
            assembled_polish_provenance);
        if (!terminal_key_list.empty())
        {
            for (const auto & key : terminal_key_list)
            {
                for (const auto atom_index : key)
                {
                    rollback_atom_mask.at(atom_index) = 0;
                }
            }
            const auto remaining_active_index_list{
                BuildEligibleLocalFittingActiveIndexList(terminal_fallback_atom_mask)
            };
            if (!remaining_active_index_list.empty())
            {
                const auto remaining_partition{
                    detail::BuildLocalFittingCouplingPartition(
                        coupling_topology,
                        remaining_active_index_list)
                };
                objective_domain = BuildLocalFittingObjectiveDomain(
                    context,
                    assembled_state,
                    remaining_partition,
                    options);
                LogLocalFittingObjectiveDomain(
                    objective_domain,
                    options,
                    true);
                cluster_objective_state.clear();
                ReconcileLocalFittingClusterObjectiveState(
                    context,
                    assembled_state,
                    remaining_partition,
                    objective_domain,
                    cluster_objective_state);
                ResetLocalFittingBestAuditAfterObjectiveDomainChange(
                    context,
                    assembled_state,
                    assembled_polish_provenance,
                    accepted_iteration_count + 1,
                    objective_domain,
                    best_audit_state);
                objective_domain_changed = true;
            }
        }

        trust_region_state.Grow(selection.grow_trust_region_key_list);
        const auto rejected_cluster_partition{
            detail::PartitionLocalFittingRejectedClusters(
                selection.rejected_key_list,
                selection.backtracking_exhausted_key_list)
        };
        const auto trust_region_radius_update{
            trust_region_state.Shrink(rejected_cluster_partition.retryable_key_list)
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
            SummarizeLocalFittingProgressMaximum(raw_fixed_point_change_summary.maximum_list)
        };
        if (selection.accepted_key_list.empty())
        {
            const auto all_rejected_resolution{
                detail::ResolveLocalFittingAllRejected(
                    iter + 1 >= kLocalFittingMaximumIterations,
                    rejected_cluster_partition,
                    trust_region_radius_update)
            };
            LogRejectedLocalFittingClusterDiagnostics(options, selection.rejected_cluster_diagnostic_list);
            LogLocalFittingIterationProgress(options, progress_column_widths, progress);
            LogLocalFittingAllRejectedResolution(
                options,
                rejected_cluster_partition,
                trust_region_radius_update,
                all_rejected_resolution);
            if (all_rejected_resolution == detail::LocalFittingAllRejectedResolution::Retry)
            {
                for (const auto & key : rejected_cluster_partition.exhausted_key_list)
                {
                    if (std::find(
                            unchanged_state_exhausted_key_list.begin(),
                            unchanged_state_exhausted_key_list.end(),
                            key) == unchanged_state_exhausted_key_list.end())
                    {
                        unchanged_state_exhausted_key_list.emplace_back(key);
                    }
                }
                continue;
            }

            const auto final_state_selection{
                SelectLocalFittingFinalState(
                    previous_state,
                    previous_polish_provenance,
                    best_audit_state.best)
            };
            ApplyLocalFittingState(
                model_object,
                context,
                *final_state_selection.state);
            diagnostics.peeling_applied = true;
            LogSecondStageLocalFittingSummary(
                options,
                accepted_iteration_count,
                GetLocalFittingAllRejectedResolutionText(all_rejected_resolution),
                best_audit_state,
                UsesLocalFittingPolish(*final_state_selection.polish_provenance),
                final_state_selection.source);
            RunGroupPotentialFitting(model_object, options, true);
            return diagnostics;
        }

        const auto transformed_change_summary{
            SummarizeLocalFittingTransformedChanges(assembled_state, previous_state, active_index_list)
        };

        accepted_iteration_count++;
        const auto improved_best_audit{
            TryUpdateLocalFittingBestAuditState(
                context,
                assembled_state,
                assembled_polish_provenance,
                accepted_iteration_count,
                objective_domain,
                best_audit_state)
        };
        audit_patience_count = objective_domain_changed ? 0 :
            detail::AdvanceLocalFittingAuditPatience(
                audit_patience_count,
                improved_best_audit,
                !selection.rejected_key_list.empty() &&
                    !trust_region_radius_update.changed_key_list.empty());
        progress.accepted_iteration_count = accepted_iteration_count;
        progress.accepted_maximum_transformed_change = SummarizeLocalFittingProgressMaximum(transformed_change_summary.maximum_list);
        LogRejectedLocalFittingClusterDiagnostics(options, selection.rejected_cluster_diagnostic_list);
        LogLocalFittingIterationProgress(options, progress_column_widths, progress);

        if (audit_patience_count >= kLocalFittingAuditPatience)
        {
            const auto final_state_selection{
                SelectLocalFittingFinalState(
                    assembled_state,
                    assembled_polish_provenance,
                    best_audit_state.best)
            };
            ApplyLocalFittingState(
                model_object,
                context,
                *final_state_selection.state);
            diagnostics.peeling_applied = true;
            LogSecondStageLocalFittingSummary(
                options,
                accepted_iteration_count,
                "audit-patience",
                best_audit_state,
                UsesLocalFittingPolish(*final_state_selection.polish_provenance),
                final_state_selection.source);
            RunGroupPotentialFitting(model_object, options, true);
            return diagnostics;
        }

        const auto converged{
            stationarity_ineligible_cluster_count == 0 &&
            !has_suspicious_offset_fallback &&
            selection.rejected_key_list.empty() &&
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
            diagnostics.peeling_applied = true;
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
                UsesLocalFittingPolish(assembled_polish_provenance),
                detail::LocalFittingFinalStateSource::LatestValidated);
            RunGroupPotentialFitting(model_object, options, true);
            return diagnostics;
        }

        if (iter + 1 == kLocalFittingMaximumIterations)
        {
            const auto final_state_selection{
                SelectLocalFittingFinalState(
                    assembled_state,
                    assembled_polish_provenance,
                    best_audit_state.best)
            };
            ApplyLocalFittingState(
                model_object,
                context,
                *final_state_selection.state);
            diagnostics.peeling_applied = true;
            LogLocalFittingMaximumIterations(
                options,
                final_state_selection.source,
                final_state_selection.audit_state,
                terminal_summary,
                SummarizeLocalFittingOffsets(*final_state_selection.state));
            LogSecondStageLocalFittingSummary(
                options,
                accepted_iteration_count,
                "maximum-iterations",
                best_audit_state,
                UsesLocalFittingPolish(*final_state_selection.polish_provenance),
                final_state_selection.source);
            RunGroupPotentialFitting(model_object, options, true);
            return diagnostics;
        }
        unchanged_state_exhausted_key_list.clear();
        previous_state = std::move(assembled_state);
        previous_polish_provenance = std::move(assembled_polish_provenance);
    }
    return diagnostics;
}

} // namespace rhbm_gem::core
