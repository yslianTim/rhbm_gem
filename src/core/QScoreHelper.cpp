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

LocalPotentialSampleList PrepareQScoreSamples(
    const AtomObject & atom, const MapObject & map, const ModelObject & model, int num_points)
{
    numeric_validation::RequireAtLeast(num_points, 2, "CalculateQScoreForAtom num_points");
    const auto center_map_value{InterpolateQScoreResponse(map, atom.GetPositionRef())};
    LocalPotentialSampleList samples(static_cast<std::size_t>(num_points),
        LocalPotentialSample{center_map_value, {0.0, atom.GetPositionRef(), true}});
    for (double radius = 0.1; radius < 2.0 + kMaximumRadiusTolerance; radius += 0.1)
    {
        const auto radial_points{GetRadialPointsForQScore(atom, model, radius, num_points)};
        for (const auto & point : radial_points)
            samples.push_back({InterpolateQScoreResponse(map, point.position), point});
    }
    return samples;
}

std::vector<double> StandardQScoreReference(
    const LocalPotentialSampleList & samples, double height, double offset, double sigma)
{
    std::vector<double> values;
    values.reserve(samples.size());
    double previous_radius=0.0;
    double value=height+offset;
    for (const auto & sample : samples)
    {
        const auto radius=sample.point.distance;
        if (radius!=previous_radius)
        {
            value=height * std::exp(-0.5 * radius * radius / (sigma * sigma)) + offset;
            previous_radius=radius;
        }
        values.push_back(value);
    }
    return values;
}

QScoreCorrelation CalculateQScoreCorrelation(
    const std::vector<double> & observations, const std::vector<double> & reference)
{
    if (observations.empty() || observations.size()!=reference.size())
        return {std::nullopt, "sample-unavailable"};
    double sum_y{}, sum_m{};
    for (std::size_t i=0;i<observations.size();++i)
    {
        if (!std::isfinite(observations[i]) || !std::isfinite(reference[i]))
            return {std::nullopt, "nonfinite-response"};
        sum_y+=observations[i]; sum_m+=reference[i];
    }
    if (std::all_of(observations.begin(),observations.end(),[&](double y){return y==observations.front();}))
        return {std::nullopt, "observation-zero-variance"};
    if (std::all_of(reference.begin(),reference.end(),[&](double m){return m==reference.front();}))
        return {std::nullopt, "reference-zero-variance"};
    const double mean_y=sum_y/static_cast<double>(observations.size());
    const double mean_m=sum_m/static_cast<double>(reference.size());
    double variance_y{},variance_m{};
    for (std::size_t i=0;i<observations.size();++i)
    {
        const double dy=observations[i]-mean_y,dm=reference[i]-mean_m;
        variance_y+=dy*dy; variance_m+=dm*dm;
    }
    if (variance_y==0) return {std::nullopt,"observation-zero-variance"};
    if (variance_m==0) return {std::nullopt,"reference-zero-variance"};
    const double denominator=std::sqrt(variance_y)*std::sqrt(variance_m);
    if (!std::isfinite(denominator) || denominator==0)
        return {std::nullopt,"nonfinite-correlation"};
    // Standard retains its original arithmetic and legacy missing-data contract.
    const double score=CalculateMeanSubtractedCorrelation(observations,reference);
    if (!std::isfinite(score) || std::abs(score)>1.0+1.0e-12)
        return {std::nullopt,"nonfinite-correlation"};
    return {score,{}};
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
    const auto samples=PrepareQScoreSamples(atom,map,model,num_points);
    std::vector<double> map_values;
    map_values.reserve(samples.size());
    for (const auto & sample:samples) map_values.push_back(sample.response);
    return CalculateMeanSubtractedCorrelation(map_values,StandardQScoreReference(samples,height,offset,sigma));
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
