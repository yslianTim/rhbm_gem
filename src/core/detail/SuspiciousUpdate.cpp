#include "core/detail/SuspiciousUpdate.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <utility>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

namespace rhbm_gem::core::detail {

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

enum class SuspiciousProfileAnalysisMode
{
    Candidate,
    PreviousBaseline
};

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

std::optional<StabilizationTerminalDiagnostic>
EvaluateClusterCandidateGuard(
    const SecondStageContext & context,
    const FitOptions & options,
    const SecondStageModelSnapshot & previous_snapshot,
    const ClusterKey & key,
    const FitStateView & candidate_state,
    const SuspiciousBlockActivity & block_activity)
{
    const auto candidate_snapshot{
        BuildSecondStageModelSnapshot(context, candidate_state)
    };
    for (const auto atom_index : key)
    {
        const auto & atom_context{ context.atom_list.at(atom_index) };
        const auto & previous_model{
            previous_snapshot.node.at(atom_index)
        };
        const auto & candidate_model{ candidate_state.GetModel(atom_index) };
        if (!IsValidSecondStageGaussianModel(candidate_model))
            return StabilizationTerminalDiagnostic{ StabilizationTerminalReason::InvalidCandidate };
        if (block_activity.HasActiveOffset(atom_index))
        {
            const auto offset_assessment{
                AssessSuspiciousGaussianUpdate(
                    atom_context.raw_sampling_entries,
                    candidate_model,
                    options,
                    BuildPreviousSuspiciousProfileBaseline(
                        atom_context.raw_sampling_entries,
                        previous_model,
                        options),
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
            BuildSecondStageAdjustedSamples(context, atom_index, previous_snapshot)
        };
        const auto candidate_samples{
            BuildSecondStageAdjustedSamples(context, atom_index, candidate_snapshot)
        };
        const auto shape_assessment{
            AssessSuspiciousGaussianUpdate(
                candidate_samples,
                candidate_model,
                options,
                BuildPreviousSuspiciousProfileBaseline(
                    previous_samples,
                    previous_model,
                    options),
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

std::size_t CountSuspiciousPolishAtoms(
    const SecondStageContext & context,
    const FitOptions & options,
    const std::vector<std::size_t> & atom_index_list,
    const FitStateView & endpoint_state,
    const FitStateView & candidate_state)
{
    const auto endpoint_snapshot{
        BuildSecondStageModelSnapshot(context, endpoint_state)
    };
    const auto candidate_snapshot{
        BuildSecondStageModelSnapshot(context, candidate_state)
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
            BuildSecondStageAdjustedSamples(context, atom_index, endpoint_snapshot)
        };
        const auto candidate_samples{
            BuildSecondStageAdjustedSamples(context, atom_index, candidate_snapshot)
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

} // namespace rhbm_gem::core::detail
