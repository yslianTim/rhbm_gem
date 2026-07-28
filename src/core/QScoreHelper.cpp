#include <rhbm_gem/core/QScoreHelper.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <unordered_map>
#include <utility>
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
constexpr double kMapBoundaryMarginRatio{ 0.5 };
constexpr double kMaximumRadiusTolerance{ 0.01 };

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
            static_cast<double>(center.at(0)) + radius * sin_phi * std::cos(theta),
            static_cast<double>(center.at(1)) + radius * sin_phi * std::sin(theta),
            static_cast<double>(center.at(2)) + radius * std::cos(phi)
        });
    }
    return points;
}

double InterpolateMapValueTrilinear(const MapObject & map, const std::array<float, 3> & position)
{
    const auto grid_size{ map.GetGridSize() };
    const auto grid_spacing{ map.GetGridSpacing() };
    const auto origin{ map.GetOrigin() };

    std::array<int, 3> lower_indices;
    std::array<int, 3> upper_indices;
    std::array<double, 3> fractions;
    for (std::size_t axis = 0; axis < position.size(); ++axis)
    {
        const auto lower_boundary{
            static_cast<double>(origin.at(axis)) -
                kMapBoundaryMarginRatio * static_cast<double>(grid_spacing.at(axis))
        };
        const auto upper_boundary{
            static_cast<double>(origin.at(axis)) +
                (static_cast<double>(grid_size.at(axis)) - kMapBoundaryMarginRatio) *
                    static_cast<double>(grid_spacing.at(axis))
        };
        const auto coordinate{ static_cast<double>(position.at(axis)) };
        if (!std::isfinite(coordinate) ||
            coordinate < lower_boundary ||
            coordinate > upper_boundary)
        {
            throw std::out_of_range(
                "CalculateQScoreForAtom sampling position is outside the map boundary.");
        }

        const auto grid_coordinate{
            std::clamp(
                (coordinate - static_cast<double>(origin.at(axis))) /
                    static_cast<double>(grid_spacing.at(axis)),
                0.0,
                static_cast<double>(grid_size.at(axis) - 1))
        };
        lower_indices.at(axis) = static_cast<int>(std::floor(grid_coordinate));
        upper_indices.at(axis) = std::min(lower_indices.at(axis) + 1, grid_size.at(axis) - 1);
        fractions.at(axis) = grid_coordinate - static_cast<double>(lower_indices.at(axis));
    }

    const auto interpolate{
        [](double lower_value, double upper_value, double fraction)
        {
            return lower_value + fraction * (upper_value - lower_value);
        }
    };

    std::array<std::array<double, 2>, 2> x_interpolated_values{};
    for (std::size_t y_offset = 0; y_offset < 2; ++y_offset)
    {
        const auto y_index{
            y_offset == 0 ? lower_indices.at(1) : upper_indices.at(1)
        };
        for (std::size_t z_offset = 0; z_offset < 2; ++z_offset)
        {
            const auto z_index{
                z_offset == 0 ? lower_indices.at(2) : upper_indices.at(2)
            };
            const auto lower_value{
                static_cast<double>(map.GetMapValue(
                    lower_indices.at(0), y_index, z_index))
            };
            const auto upper_value{
                static_cast<double>(map.GetMapValue(
                    upper_indices.at(0), y_index, z_index))
            };
            x_interpolated_values.at(y_offset).at(z_offset) = interpolate(
                lower_value,
                upper_value,
                fractions.at(0));
        }
    }

    std::array<double, 2> y_interpolated_values{};
    for (std::size_t z_offset = 0; z_offset < 2; ++z_offset)
    {
        y_interpolated_values.at(z_offset) = interpolate(
            x_interpolated_values.at(0).at(z_offset),
            x_interpolated_values.at(1).at(z_offset),
            fractions.at(1));
    }
    return interpolate(
        y_interpolated_values.at(0),
        y_interpolated_values.at(1),
        fractions.at(2));
}

double CalculateMeanSubtractedCorrelation(
    const std::vector<double> & map_values,
    const std::vector<double> & reference_values)
{
    if (map_values.empty())
    {
        return 0.0;
    }

    double map_value_sum{ 0.0 };
    double reference_value_sum{ 0.0 };
    for (std::size_t i = 0; i < map_values.size(); ++i)
    {
        map_value_sum += map_values.at(i);
        reference_value_sum += reference_values.at(i);
    }
    const auto sample_count{ static_cast<double>(map_values.size()) };
    const auto map_value_mean{ map_value_sum / sample_count };
    const auto reference_value_mean{ reference_value_sum / sample_count };

    double numerator{ 0.0 };
    double map_value_norm_square{ 0.0 };
    double reference_value_norm_square{ 0.0 };
    for (std::size_t i = 0; i < map_values.size(); ++i)
    {
        const auto map_value_difference{
            map_values.at(i) - map_value_mean
        };
        const auto reference_value_difference{
            reference_values.at(i) - reference_value_mean
        };
        numerator += map_value_difference * reference_value_difference;
        map_value_norm_square += map_value_difference * map_value_difference;
        reference_value_norm_square += reference_value_difference * reference_value_difference;
    }

    const auto denominator{
        std::sqrt(map_value_norm_square) * std::sqrt(reference_value_norm_square)
    };
    if (denominator == 0.0 ||
        !std::isfinite(denominator) ||
        !std::isfinite(numerator))
    {
        return 0.0;
    }

    const auto q_score{ numerator / denominator };
    return std::isfinite(q_score) ? q_score : 0.0;
}

} // namespace

std::tuple<double, double> GetReferenceGaussianParameters(const MapObject & map_object)
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
    return std::make_tuple(static_cast<double>(reference_high - offset), static_cast<double>(offset));
}

SamplingPointList GetRadialPointsForQScore(
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
                        candidate_point.at(axis) - static_cast<double>(neighbor_position.at(axis))
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

double CalculateQScoreForAtom(
    const AtomObject & atom,
    const MapObject & map,
    const ModelObject & model,
    double height,
    double offset,
    double sigma,
    int num_points)
{
    numeric_validation::RequireFinitePositive(sigma, "CalculateQScoreForAtom sigma");
    numeric_validation::RequireAtLeast(num_points, 2, "CalculateQScoreForAtom num_points");
    
    const auto center_map_value{
        InterpolateMapValueTrilinear(map, atom.GetPositionRef())
    };
    const auto center_reference_value{ height + offset };

    std::vector<double> map_values(static_cast<std::size_t>(num_points), center_map_value);
    std::vector<double> reference_values(static_cast<std::size_t>(num_points), center_reference_value);

    for (double radius = 0.1; radius < 2.0 + kMaximumRadiusTolerance; radius += 0.1)
    {
        const auto radial_points{
            GetRadialPointsForQScore(atom, model, radius, num_points)
        };
        if (radial_points.empty()) continue;
        const auto reference_value{
            height * std::exp(-0.5 * radius * radius / (sigma * sigma)) + offset
        };
        for (const auto & radial_point : radial_points)
        {
            map_values.emplace_back(InterpolateMapValueTrilinear(map, radial_point.position));
            reference_values.emplace_back(reference_value);
        }
    }

    return CalculateMeanSubtractedCorrelation(map_values, reference_values);
}

double CalculateQScoreForAtom(
    const LocalPotentialSampleList & sampling_entries,
    double height,
    double offset,
    double sigma)
{
    if (sampling_entries.empty()) return 0.0;
    if (std::isfinite(sigma) == false || sigma <= 0.0) return 0.0;

    std::vector<double> map_values;
    std::vector<double> reference_values;
    map_values.reserve(sampling_entries.size());
    reference_values.reserve(sampling_entries.size());
    for (const auto & sample : sampling_entries)
    {
        const auto radius{ static_cast<double>(sample.point.distance) };
        const auto reference_value{
            height * std::exp(-0.5 * radius * radius / (sigma * sigma)) + offset
        };
        map_values.emplace_back(sample.response);
        reference_values.emplace_back(reference_value);
    }

    return CalculateMeanSubtractedCorrelation(map_values, reference_values);
}

namespace {

double CalculateAverageQScoresImpl(
    const MapObject & map,
    const ModelObject & model,
    std::unordered_map<int, double> * q_scores_by_serial_id)
{
    const auto [height, offset]{ GetReferenceGaussianParameters(map) };
    std::unordered_map<int, double> computed_q_scores;
    if (q_scores_by_serial_id != nullptr)
    {
        computed_q_scores.reserve(model.GetNumberOfAtom());
    }

    double q_score_sum{ 0.0 };
    std::size_t atom_count{ 0 };
    for (const auto & atom : model.GetAtomList())
    {
        if (atom->GetElement() == Element::HYDROGEN) continue;
        const auto q_score{ CalculateQScoreForAtom(*atom, map, model, height, offset) };
        q_score_sum += q_score;
        if (q_scores_by_serial_id != nullptr)
        {
            computed_q_scores[atom->GetSerialID()] = q_score;
        }
        atom_count++;
    }

    if (q_scores_by_serial_id != nullptr)
    {
        *q_scores_by_serial_id = std::move(computed_q_scores);
    }
    return atom_count == 0 ? 0.0 : q_score_sum / static_cast<double>(atom_count);
}

} // namespace

double CalculateAverageQScores(const MapObject & map, const ModelObject & model)
{
    return CalculateAverageQScoresImpl(map, model, nullptr);
}

double CalculateAverageQScores(
    const MapObject & map,
    const ModelObject & model,
    std::unordered_map<int, double> & q_scores_by_serial_id)
{
    return CalculateAverageQScoresImpl(map, model, &q_scores_by_serial_id);
}

} // namespace rhbm_gem::core
