#include <rhbm_gem/utils/math/GaussianPotentialSampler.hpp>

#include <rhbm_gem/utils/math/GaussianResponseMath.hpp>
#include <rhbm_gem/utils/math/SamplingPointFilter.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
void ValidatePositiveSampleCount(size_t value, const char * name)
{
    if (value == 0)
    {
        throw std::invalid_argument(std::string(name) + " must be positive value");
    }
}

void ValidateDistanceRange(double min_value, double max_value, const char * name)
{
    if (!std::isfinite(min_value) || !std::isfinite(max_value) || min_value < 0.0 ||
        max_value < min_value)
    {
        throw std::invalid_argument(
            std::string(name) + " must be finite and satisfy 0 <= min <= max");
    }
}

void ValidateGaussianModel(const GaussianModel3D & model)
{
    if (!std::isfinite(model.amplitude) || !std::isfinite(model.width) ||
        !std::isfinite(model.intercept))
    {
        throw std::invalid_argument("GaussianModel3D values must be finite.");
    }
    if (model.width <= 0.0)
    {
        throw std::invalid_argument("GaussianModel3D width must be positive.");
    }
}

std::vector<Eigen::VectorXd> BuildNeighborCenterList(
    const NeighborhoodSamplingOptions & options)
{
    constexpr size_t max_neighbor_count{ 4 };
    if (options.neighbor_count > max_neighbor_count)
    {
        throw std::invalid_argument("neighbor_count should be less than or equal to 4");
    }

    std::vector<Eigen::VectorXd> neighbor_center_list;
    neighbor_center_list.reserve(options.neighbor_count);
    if (options.neighbor_count == 0)
    {
        return neighbor_center_list;
    }

    Eigen::VectorXd neighbor_center[max_neighbor_count];
    for (size_t i = 0; i < max_neighbor_count; i++)
    {
        neighbor_center[i] = Eigen::VectorXd::Zero(3);
    }

    if (options.neighbor_count <= 2)
    {
        neighbor_center[0] << 1.0, 0.0, 0.0;
        if (options.neighbor_count == 2)
        {
            neighbor_center[1] << -1.0, 0.0, 0.0;
        }
    }

    if (options.neighbor_count == 3)
    {
        neighbor_center[0] << 1.0, 0.0, 0.0;
        neighbor_center[1] << -0.5, std::sqrt(3) / 2.0, 0.0;
        neighbor_center[2] << -0.5, -std::sqrt(3) / 2.0, 0.0;
    }

    if (options.neighbor_count == 4)
    {
        neighbor_center[0] << 0.0, 0.0, 1.0;
        neighbor_center[1] << 0.0, 2.0 * std::sqrt(2) / 3.0, -1.0 / 3.0;
        neighbor_center[2] << -std::sqrt(6) / 3.0, -std::sqrt(2) / 3.0, -1.0 / 3.0;
        neighbor_center[3] << std::sqrt(6) / 3.0, -std::sqrt(2) / 3.0, -1.0 / 3.0;
    }

    for (size_t i = 0; i < options.neighbor_count; i++)
    {
        neighbor_center[i] *= options.neighbor_distance;
        neighbor_center_list.emplace_back(neighbor_center[i]);
    }

    return neighbor_center_list;
}
} // namespace

LocalPotentialSampleList GaussianPotentialSampler::GenerateRadialSamples(
    size_t sample_count,
    const GaussianModel3D & model,
    double distance_min,
    double distance_max,
    std::mt19937 & generator
) const
{
    ValidatePositiveSampleCount(sample_count, "sample_count");
    ValidateGaussianModel(model);
    ValidateDistanceRange(distance_min, distance_max, "distance range");

    std::uniform_real_distribution<> dist_distance(distance_min, distance_max);
    LocalPotentialSampleList sample_list;
    sample_list.reserve(sample_count);
    for (size_t i = 0; i < sample_count; i++)
    {
        const auto distance{ dist_distance(generator) };
        const auto response{
            model.amplitude *
                rhbm_gem::GaussianResponseMath::GetGaussianResponseAtDistance(
                    distance,
                    model.width) +
            model.intercept
        };
        sample_list.emplace_back(LocalPotentialSample{
            static_cast<float>(distance),
            static_cast<float>(response),
            1.0
        });
    }

    return sample_list;
}

LocalPotentialSampleList GaussianPotentialSampler::GenerateNeighborhoodSamples(
    size_t samples_per_radius,
    const GaussianModel3D & model,
    const NeighborhoodSamplingOptions & options
) const
{
    ValidatePositiveSampleCount(samples_per_radius, "samples_per_radius");
    ValidateGaussianModel(model);
    ValidateDistanceRange(options.radius_min, options.radius_max, "radius range");
    if (!std::isfinite(options.neighbor_distance) || options.neighbor_distance <= 0.0)
    {
        throw std::invalid_argument("neighbor_distance must be finite and positive.");
    }

    auto neighbor_center_list{ BuildNeighborCenterList(options) };
    const Eigen::VectorXd atom_center{ Eigen::VectorXd::Zero(3) };

    SphereSampler sampler;
    sampler.SetSamplingProfile(
        SphereSamplingProfile::FibonacciDeterministic(
            SphereDistanceRange{ options.radius_min, options.radius_max },
            0.1,
            static_cast<unsigned int>(samples_per_radius),
            false
        )
    );
    const auto sampling_points{ sampler.GenerateSamplingPoints({ 0.0f, 0.0f, 0.0f }) };
    const auto sampling_scores{
        rhbm_gem::BuildSamplingPointScoreList(
            sampling_points,
            neighbor_center_list,
            options.reject_angle_deg
        )
    };

    LocalPotentialSampleList sample_list;
    if (options.rejected_point_policy ==
        NeighborhoodRejectedPointPolicy::KeepRejectedPointsWithZeroScore)
    {
        sample_list.reserve(sampling_points.size());
    }
    else
    {
        sample_list.reserve(static_cast<size_t>(std::count_if(
            sampling_scores.begin(),
            sampling_scores.end(),
            [](float score)
            {
                return score > 0.0f;
            })));
    }

    for (size_t i = 0; i < sampling_points.size(); i++)
    {
        const auto & sampling_point{ sampling_points.at(i) };
        const auto sampling_score{ sampling_scores.at(i) };
        if (
            options.rejected_point_policy ==
                NeighborhoodRejectedPointPolicy::RemoveRejectedPoints &&
            sampling_score <= 0.0f)
        {
            continue;
        }

        Eigen::VectorXd point{ Eigen::VectorXd::Zero(3) };
        point(0) = sampling_point.position[0];
        point(1) = sampling_point.position[1];
        point(2) = sampling_point.position[2];

        const auto response{
            model.amplitude *
                rhbm_gem::GaussianResponseMath::GetGaussianResponseAtPointWithNeighborhood(
                    point,
                    atom_center,
                    neighbor_center_list,
                    model.width
                ) +
            model.intercept
        };
        sample_list.emplace_back(LocalPotentialSample{
            sampling_point.distance,
            static_cast<float>(response),
            sampling_score,
            sampling_point.position
        });
    }

    return sample_list;
}
