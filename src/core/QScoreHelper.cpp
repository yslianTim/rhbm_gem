#include <rhbm_gem/core/QScoreHelper.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

namespace rhbm_gem::core {
namespace {

constexpr int kMaximumRadialPointAttempts{ 50 };
constexpr std::size_t kCandidatePointIncrement{ 2 };
constexpr double kNeighborExclusionRadiusRatio{ 0.9 };
constexpr double kNeighborSearchRadiusRatio{ 2.0 };
constexpr double kSpiralSpacing{ 3.6 };

using RadialPoint = std::array<double, 3>;

std::vector<RadialPoint> GenerateSpiralSpherePoints(
    const std::array<float, 3> & center,
    double radius,
    std::size_t point_count)
{
    std::vector<RadialPoint> points;
    points.reserve(point_count);

    double previous_theta{ 0.0 };
    for (std::size_t i = 0; i < point_count; ++i)
    {
        const auto height{
            -1.0 + 2.0 * static_cast<double>(i) /
                static_cast<double>(point_count - 1)
        };
        const auto phi{ std::acos(height) };
        double theta{ 0.0 };
        if (i > 0 && i + 1 < point_count)
        {
            theta = std::fmod(
                previous_theta +
                    kSpiralSpacing /
                    std::sqrt(static_cast<double>(point_count) * (1.0 - height * height)),
                Constants::two_pi);
        }
        previous_theta = theta;

        const auto sin_phi{ std::sin(phi) };
        points.emplace_back(RadialPoint{
            static_cast<double>(center.at(0)) + radius * sin_phi * std::cos(theta),
            static_cast<double>(center.at(1)) + radius * sin_phi * std::sin(theta),
            static_cast<double>(center.at(2)) + radius * std::cos(phi)
        });
    }
    return points;
}

} // namespace

std::tuple<float, float> GetReferenceGaussianParameters(const MapObject & map_object)
{
    const auto reference_high{
        std::min(
            map_object.GetMapValueMean() + 10.0f * map_object.GetMapValueSD(),
            map_object.GetMapValueMax())
    };
    const auto offset{
        std::max(
            map_object.GetMapValueMean() - map_object.GetMapValueSD(),
            map_object.GetMapValueMin())
    };
    return std::make_tuple(reference_high - offset, offset);
}

SamplingPointList GetRadialPointsForQScore(
    const AtomObject & atom,
    const ModelObject & model,
    double radius,
    int num_points)
{
    numeric_validation::RequireFinitePositive(
        radius,
        "GetRadialPointsForQScore radius");
    numeric_validation::RequireAtLeast(
        num_points,
        2,
        "GetRadialPointsForQScore num_points");

    auto neighbor_atoms{
        model.FindNeighborAtoms(
            atom,
            kNeighborSearchRadiusRatio * radius)
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
            target_point_count +
                static_cast<std::size_t>(attempt) * kCandidatePointIncrement
        };
        const auto candidate_points{
            GenerateSpiralSpherePoints(
                atom.GetPositionRef(),
                radius,
                candidate_point_count)
        };

        SamplingPointList accepted_points;
        accepted_points.reserve(candidate_points.size());
        for (const auto & candidate_point : candidate_points)
        {
            bool is_rejected{ false };
            for (const auto * neighbor_atom : neighbor_atoms)
            {
                const auto & neighbor_position{
                    neighbor_atom->GetPositionRef()
                };
                double distance_square{ 0.0 };
                for (std::size_t axis = 0; axis < candidate_point.size(); ++axis)
                {
                    const auto difference{
                        candidate_point.at(axis) -
                            static_cast<double>(neighbor_position.at(axis))
                    };
                    distance_square += difference * difference;
                }
                if (distance_square < exclusion_radius_square)
                {
                    is_rejected = true;
                    break;
                }
            }
            if (!is_rejected)
            {
                accepted_points.emplace_back(SamplingPoint{
                    static_cast<float>(radius),
                    {
                        static_cast<float>(candidate_point.at(0)),
                        static_cast<float>(candidate_point.at(1)),
                        static_cast<float>(candidate_point.at(2))
                    },
                    true
                });
            }
        }

        if (accepted_points.size() >= target_point_count)
        {
            return accepted_points;
        }
    }

    return {};
}

} // namespace rhbm_gem::core
