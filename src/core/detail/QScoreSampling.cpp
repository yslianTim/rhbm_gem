#include "QScoreSampling.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>
#include <algorithm>

namespace rhbm_gem::core::detail {
namespace {
constexpr int kMaximumRadialPointAttempts{ 50 };
constexpr std::size_t kCandidatePointIncrement{ 2 };
constexpr double kNeighborExclusionRadiusRatio{ 0.9 };
constexpr double kNeighborSearchRadiusRatio{ 2.0 };
} // namespace

bool IsInsideNeighborExclusion(
    const std::array<double, 3> & sample_position,
    const std::array<double, 3> & neighbor_position,
    double exclusion_radius_squared)
{
    numeric_validation::RequireAllFinite(sample_position, "sample position");
    numeric_validation::RequireAllFinite(neighbor_position, "neighbor position");
    numeric_validation::RequireFiniteNonNegative(exclusion_radius_squared, "exclusion radius squared");
    double distance_square{ 0.0 };
    for (std::size_t axis = 0; axis < sample_position.size(); ++axis)
    {
        const auto difference{ sample_position.at(axis) - neighbor_position.at(axis) };
        distance_square += difference * difference;
    }
    return distance_square < exclusion_radius_squared;
}

SamplingPointList SampleQScoreShell(
    const AtomObject & atom,
    const ModelObject & model,
    double radius,
    int num_points)
{
    numeric_validation::RequireFinitePositive(radius, "GetRadialPointsForQScore radius");
    numeric_validation::RequireAtLeast(num_points, 2, "GetRadialPointsForQScore num_points");

    auto neighbor_atoms{
        model.FindNeighborAtoms(atom, kNeighborSearchRadiusRatio * radius)
    };
    neighbor_atoms.erase(
        std::remove_if(
            neighbor_atoms.begin(),
            neighbor_atoms.end(),
            [](const AtomObject * neighbor_atom)
            {
                return neighbor_atom->GetElement() == Element::HYDROGEN;
            }),
        neighbor_atoms.end());
    const auto exclusion_radius{
        kNeighborExclusionRadiusRatio * radius
    };
    const auto exclusion_radius_square{
        exclusion_radius * exclusion_radius
    };
    const auto target_point_count{
        static_cast<std::size_t>(num_points)
    };

    for (int attempt = 0; attempt < kMaximumRadialPointAttempts; ++attempt)
    {
        const auto candidate_point_count{
            target_point_count + static_cast<std::size_t>(attempt) * kCandidatePointIncrement
        };
        auto accepted_points{
            sphere_sampler::GenerateSpiralShell(atom.GetPositionRef(), radius, candidate_point_count)
        };
        std::erase_if(accepted_points, [&](const SamplingPoint & point)
        {
            return std::any_of(neighbor_atoms.begin(), neighbor_atoms.end(), [&](const AtomObject * neighbor)
            {
                return IsInsideNeighborExclusion(
                    point.position, neighbor->GetPositionRef(), exclusion_radius_square);
            });
        });

        if (accepted_points.size() >= target_point_count)
        {
            return accepted_points;
        }
    }

    return {};
}

} // namespace rhbm_gem::core::detail
