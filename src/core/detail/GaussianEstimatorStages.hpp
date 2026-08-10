#pragma once

#include <cstddef>
#include <optional>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>

namespace rhbm_gem::core {

namespace detail {

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

SuspiciousGaussianReason EvaluateSuspiciousOffsetUpdate(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & previous_model,
    const GaussianModel3D & candidate_model,
    const FitOptions & options);

SuspiciousGaussianReason EvaluateSuspiciousPostRefitUpdate(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & previous_model,
    const GaussianModel3D & candidate_model,
    const FitOptions & options);

std::optional<double> CalculateLocalFittingPeelingRatio(
    const LocalPotentialSampleList & raw_sampling_entries,
    const LocalPotentialSampleList & peeling_sampling_entries,
    bool peeling_applied);

std::vector<char> ExpandSuspiciousSharedOffsetGroups(
    const std::vector<GroupKey> & group_key_by_position,
    const std::vector<char> & suspicious_seed_mask);

} // namespace detail

void RunLocalAlphaTraining(
    ModelObject & model_object,
    const FitOptions & options,
    bool use_peeling_sampling_entries);
void RunFixedOffsetLocalFitting(
    ModelObject & model_object,
    const FitOptions & options,
    bool use_peeling_sampling_entries);
void RunGroupPotentialFitting(
    ModelObject & model_object,
    const FitOptions & options,
    bool use_peeling_sampling_entries);

bool RunSecondStageLocalFitting(
    ModelObject & model_object,
    const FitOptions & options);

} // namespace rhbm_gem::core
