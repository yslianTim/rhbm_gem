#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <ranges>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#include "core/detail/TransformedGaussianModel.hpp"

namespace rhbm_gem::core::detail {

enum class SuspiciousGaussianReason
{
    None,
    InvalidModel,
    NonFiniteResponse,
    OffsetMagnitude,
    CenterSignFlip,
    RadialRebound,
    WidthGrowth,
    AmplitudeOffsetCompensation
};

using SuspiciousUpdateMask = std::vector<char>;

constexpr double kSuspiciousJointOffsetRidgeMultiplier{ 10.0 };

inline std::size_t CountSuspiciousAtoms(const SuspiciousUpdateMask & suspicious_mask)
{
    return static_cast<std::size_t>(
        std::ranges::count_if(suspicious_mask, std::identity{}));
}

inline bool HasSuspiciousAtom(
    const std::vector<std::size_t> & atom_index_list,
    const SuspiciousUpdateMask & suspicious_mask)
{
    return std::any_of(
        atom_index_list.begin(),
        atom_index_list.end(),
        [&](std::size_t atom_index)
        {
            return suspicious_mask.at(atom_index) != 0;
        });
}

inline std::vector<std::size_t> CollectSuspiciousAtomIndices(
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

inline std::vector<double> BuildSuspiciousJointOffsetRidgeMultiplierList(const SuspiciousUpdateMask & suspicious_mask)
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

inline void ClearSuspiciousUpdateMaskForClusters(
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

struct ZeroOffsetProfileDiagnostics
{
    double distance_range{ 0.0 };
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

enum class SuspiciousProfileAnalysisMode
{
    Candidate,
    PreviousBaseline
};

struct SuspiciousUpdateBaseline
{
    GaussianModel3D previous_model{};
    SuspiciousProfileAnalysis previous_analysis{};
};

enum class SuspiciousUpdateMode
{
    OffsetOnly,
    PostRefit
};

inline bool HasSuspiciousCenterSignFlip(
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

inline double CalculateZeroOffsetResponse(
    const LocalPotentialSample & sample,
    const GaussianModel3D & model)
{
    const auto distance{ static_cast<double>(sample.point.distance) };
    const auto model_offset{ model.ResponseAtDistance(distance) - model.SignalAtDistance(distance) };
    return static_cast<double>(sample.response) - model_offset;
}

inline bool IsSameSuspiciousProfileRadius(double lhs, double rhs)
{
    const auto scale{ std::max({ std::abs(lhs), std::abs(rhs), 1.0 }) };
    return std::abs(lhs - rhs) <= kSuspiciousProfileDistanceTolerance * scale;
}

inline SuspiciousProfileAnalysis BuildSuspiciousProfileAnalysis(
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

inline bool HasUsableSuspiciousProfileBaseline(
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

inline bool HasSuspiciousOffsetMagnitude(
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

inline bool HasSuspiciousRadialRebound(
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

inline bool HasSuspiciousWidthGrowth(
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

inline bool HasSuspiciousAmplitudeOffsetCompensation(
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
            kSuspiciousProfileScaleMin
        })
    };
    return signal_delta * offset_delta_response < 0.0 &&
        std::abs(signal_delta) > kSuspiciousCompensationResponseRatio * reference_scale &&
        std::abs(offset_delta_response) > kSuspiciousCompensationResponseRatio * reference_scale;
}

inline SuspiciousGaussianReason EvaluateSuspiciousGaussianUpdate(
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

inline SuspiciousUpdateBaseline BuildPreviousSuspiciousProfileBaseline(
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

inline SuspiciousGaussianReason EvaluateSuspiciousOffsetUpdate(
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

inline SuspiciousUpdateMask ExpandSuspiciousSharedOffsetGroups(
    const std::vector<std::size_t> & group_id_by_position,
    const SuspiciousUpdateMask & suspicious_seed_mask)
{
    if (group_id_by_position.size() != suspicious_seed_mask.size())
    {
        throw std::invalid_argument(
            "Suspicious shared-offset group input sizes are inconsistent.");
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

} // namespace rhbm_gem::core::detail
