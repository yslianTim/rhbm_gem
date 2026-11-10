#include <rhbm_gem/core/QScoreHelper.hpp>
#include "core/sampling/MapInterpolation.hpp"

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
#include <rhbm_gem/utils/math/NumericValidation.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>

namespace rhbm_gem::core {
namespace {

constexpr double kMaximumRadiusTolerance{ 0.01 };

constexpr int kMaximumRadialPointAttempts{ 50 };
constexpr std::size_t kCandidatePointIncrement{ 2 };
constexpr double kNeighborExclusionRadiusRatio{ 0.9 };
constexpr double kNeighborSearchRadiusRatio{ 2.0 };

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

double InterpolateQScoreResponse(const MapObject & map, const std::array<double, 3> & position)
{
    try
    {
        return detail::InterpolateTrilinear(map, position);
    }
    catch (const std::out_of_range &)
    {
        throw std::out_of_range(
            "CalculateQScoreForAtom sampling position is outside the map boundary.");
    }
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
            map_object.GetMapValueMean() + 10.0 * map_object.GetMapValueSD(),
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
    numeric_validation::RequireFinitePositive(radius, "GetRadialPointsForQScore radius");
    numeric_validation::RequireAtLeast(num_points, 2, "GetRadialPointsForQScore num_points");

    auto neighbor_atoms{
        model.FindNeighborAtoms(atom, kNeighborSearchRadiusRatio * radius)
    };
    std::erase_if(neighbor_atoms, [](const AtomObject * neighbor)
    {
        return neighbor->GetElement() == Element::HYDROGEN;
    });
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
        InterpolateQScoreResponse(map, atom.GetPositionRef())
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
            map_values.emplace_back(InterpolateQScoreResponse(map, radial_point.position));
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
        const auto radius{ sample.point.distance };
        const auto reference_value{
            height * std::exp(-0.5 * radius * radius / (sigma * sigma)) + offset
        };
        map_values.emplace_back(sample.response);
        reference_values.emplace_back(reference_value);
    }

    return CalculateMeanSubtractedCorrelation(map_values, reference_values);
}

double CalculateAverageQScores(
    const MapObject & map,
    const ModelObject & model,
    double reference_height,
    double reference_offset,
    std::unordered_map<int, double> & q_scores_by_serial_id)
{
    std::unordered_map<int, double> computed_q_scores;
    computed_q_scores.reserve(model.GetNumberOfAtom());

    double q_score_sum{ 0.0 };
    std::size_t atom_count{ 0 };
    for (const auto & atom : model.GetAtomList())
    {
        if (atom->GetElement() == Element::HYDROGEN) continue;
        const auto q_score{
            CalculateQScoreForAtom(*atom, map, model, reference_height, reference_offset)
        };
        q_score_sum += q_score;
        computed_q_scores[atom->GetSerialID()] = q_score;
        atom_count++;
    }

    q_scores_by_serial_id = std::move(computed_q_scores);
    return atom_count == 0 ? 0.0 : q_score_sum / static_cast<double>(atom_count);
}

} // namespace rhbm_gem::core
