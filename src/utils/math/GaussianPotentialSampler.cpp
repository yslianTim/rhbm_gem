#include <rhbm_gem/utils/math/GaussianPotentialSampler.hpp>

#include <rhbm_gem/utils/math/GaussianResponseMath.hpp>
#include <rhbm_gem/utils/math/LocalPotentialSampleScoring.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace
{
std::vector<Eigen::VectorXd> BuildNeighborCenterList(
    const NeighborhoodSamplingOptions & options)
{
    constexpr size_t max_neighbor_count{ 4 };
    rhbm_gem::NumericValidation::RequireAtMost(
        options.neighbor_count,
        max_neighbor_count,
        "neighbor_count");

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
    rhbm_gem::NumericValidation::RequirePositive(sample_count, "sample_count");
    rhbm_gem::NumericValidation::RequireFinite(model.amplitude, "GaussianModel3D amplitude");
    rhbm_gem::NumericValidation::RequireFinitePositive(model.width, "GaussianModel3D width");
    rhbm_gem::NumericValidation::RequireFinite(model.intercept, "GaussianModel3D intercept");
    rhbm_gem::NumericValidation::RequireFiniteNonNegativeRange(
        distance_min,
        distance_max,
        "distance range");

    std::uniform_real_distribution<> dist_distance(distance_min, distance_max);
    const auto sampling_scores{
        rhbm_gem::BuildDefaultLocalPotentialSampleScoreList(sample_count)
    };
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
            sampling_scores.at(i)
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
    rhbm_gem::NumericValidation::RequirePositive(samples_per_radius, "samples_per_radius");
    rhbm_gem::NumericValidation::RequireFinite(model.amplitude, "GaussianModel3D amplitude");
    rhbm_gem::NumericValidation::RequireFinitePositive(model.width, "GaussianModel3D width");
    rhbm_gem::NumericValidation::RequireFinite(model.intercept, "GaussianModel3D intercept");
    rhbm_gem::NumericValidation::RequireFiniteNonNegativeRange(
        options.radius_min,
        options.radius_max,
        "radius range");
    rhbm_gem::NumericValidation::RequireFinitePositive(
        options.neighbor_distance,
        "neighbor_distance");

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
        rhbm_gem::BuildLocalPotentialSampleScoreList(
            sampling_points, neighbor_center_list, options.reject_angle_deg)
        //rhbm_gem::BuildLocalPotentialCleanSampleScoreList(
        //    sampling_points, neighbor_center_list)
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
