#include "QScoreSampling.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace rhbm_gem::core::detail {
namespace {
constexpr int kMaximumRadialPointAttempts{ 50 };
constexpr std::size_t kCandidatePointIncrement{ 2 };
constexpr double kNeighborExclusionRadiusRatio{ 0.9 };
constexpr double kNeighborSearchRadiusRatio{ 2.0 };
constexpr double kSpiralSpacing{ 3.6 };
}
std::vector<RadialPoint> GenerateSpiralSpherePoints(
    const std::array<double, 3> & center,
    double radius,
    std::size_t point_count)
{
    if (point_count < 2) throw std::invalid_argument("Spiral sphere requires at least two points.");
    std::vector<RadialPoint> points;
    points.reserve(point_count);

    double previous_theta{ 0.0 };
    for (std::size_t i = 0; i < point_count; ++i)
    {
        const auto height{
            -1.0 + 2.0 * static_cast<double>(i) / static_cast<double>(point_count - 1)
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
            center.at(0) + radius * sin_phi * std::cos(theta),
            center.at(1) + radius * sin_phi * std::sin(theta),
            center.at(2) + radius * std::cos(phi)
        });
    }
    return points;
}

QScoreShellSamples SampleQScoreShell(
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

    QScoreShellSamples result;
    result.diagnostics.radius = radius;
    result.diagnostics.requested_count = target_point_count;
    for (int attempt = 0; attempt < kMaximumRadialPointAttempts; ++attempt)
    {
        const auto candidate_point_count{
            target_point_count + static_cast<std::size_t>(attempt) * kCandidatePointIncrement
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
                        candidate_point.at(axis) - neighbor_position.at(axis)
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
                    radius,
                    candidate_point,
                    true
                });
            }
        }

        result.diagnostics.attempt_count = static_cast<std::size_t>(attempt + 1);
        result.diagnostics.candidate_count = candidate_point_count;
        result.diagnostics.accepted_count = accepted_points.size();
        if (accepted_points.size() >= target_point_count)
        {
            result.diagnostics.successful = true;
            result.points = std::move(accepted_points);
            return result;
        }
    }

    return result;
}

QScoreSeedSamples BuildQScoreSeedSamples(const AtomObject & atom, const ModelObject & model)
{
    QScoreSeedSamples result;
    result.points.push_back(SamplingPoint{0.0, atom.GetPositionRef(), true});
    for (double radius = 0.1; radius < 2.0 + 0.01; radius += 0.1)
    {
        auto shell = SampleQScoreShell(atom, model, radius, 8);
        result.shells.push_back(shell.diagnostics);
        result.points.insert(result.points.end(), shell.points.begin(), shell.points.end());
    }
    return result;
}
} // namespace rhbm_gem::core::detail
