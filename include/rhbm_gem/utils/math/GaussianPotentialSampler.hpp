#pragma once

#include <cstddef>
#include <random>

#include <rhbm_gem/utils/math/SamplingTypes.hpp>

struct GaussianModel3D
{
    double amplitude{ 0.0 };
    double width{ 1.0 };
    double intercept{ 0.0 };
};

enum class NeighborhoodRejectedPointPolicy
{
    RemoveRejectedPoints,
    KeepRejectedPointsWithZeroScore
};

struct NeighborhoodSamplingOptions
{
    double radius_min{ 0.0 };
    double radius_max{ 1.0 };
    double neighbor_distance{ 2.0 };
    size_t neighbor_count{ 1 };
    double reject_angle_deg{ 0.0 };
    NeighborhoodRejectedPointPolicy rejected_point_policy{
        NeighborhoodRejectedPointPolicy::RemoveRejectedPoints
    };
};

class GaussianPotentialSampler
{
public:
    GaussianPotentialSampler() = default;
    ~GaussianPotentialSampler() = default;

    LocalPotentialSampleList GenerateRadialSamples(
        size_t sample_count,
        const GaussianModel3D & model,
        double distance_min,
        double distance_max,
        std::mt19937 & generator
    ) const;

    LocalPotentialSampleList GenerateNeighborhoodSamples(
        size_t samples_per_radius,
        const GaussianModel3D & model,
        const NeighborhoodSamplingOptions & options
    ) const;
};
