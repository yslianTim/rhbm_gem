#include <rhbm_gem/utils/domain/SampleFilter.hpp>

#include <cmath>
#include <map>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

namespace rhbm_gem::sample_filter {
namespace {

bool IsOwnedByNeighbor(
    const std::array<double, 3> & sample_position,
    const std::array<double, 3> & local_position,
    const std::vector<std::array<double, 3>> & valid_neighbor_list)
{
    bool is_bad_point{ false };
    const auto to_local_distance{
        array_helper::ComputeNorm(sample_position, local_position)
    };
    for (const auto & neighbor_position : valid_neighbor_list)
    {
        const auto to_neighbor_distance{
            array_helper::ComputeNorm(sample_position, neighbor_position)
        };
        if (to_local_distance > to_neighbor_distance) is_bad_point = true;
    }
    return is_bad_point;
}

bool IsInsideNeighborCone(
    const std::array<double, 3> & sample_position,
    const std::array<double, 3> & local_position,
    const std::vector<std::array<double, 3>> & valid_neighbor_list,
    double angle)
{
    const auto cos_threshold{ std::cos(angle * Constants::pi / 180.0) };
    const auto sampling_unit_vector{
        array_helper::ComputeVector(local_position, sample_position, true)
    };

    for (const auto & neighbor_position : valid_neighbor_list)
    {
        const auto neighbor_unit_vector{
            array_helper::ComputeVector(local_position, neighbor_position, true)
        };
        const auto cos_theta{
            array_helper::ComputeDotProduct(sampling_unit_vector, neighbor_unit_vector)
        };
        if (cos_theta > cos_threshold) return true;
    }
    return false;
}

} // namespace

void FilterSamplingPointList(
    SamplingPointList & sample_point_list,
    const std::array<double, 3> & local_position,
    const std::vector<std::array<double, 3>> & reject_position_list,
    double angle)
{
    numeric_validation::RequireFiniteInclusiveRange(angle, 0.0, 180.0, "angle");

    for (auto & sample_point : sample_point_list)
    {
        sample_point.is_selected = true;
    }

    if (reject_position_list.empty())
    {
        return;
    }

    std::vector<std::array<double, 3>> valid_neighbor_list;
    valid_neighbor_list.reserve(reject_position_list.size());
    for (const auto & reject_position : reject_position_list)
    {
        numeric_validation::RequireAllFinite(reject_position, "reject positions");
        if (!numeric_validation::IsNonEqual(local_position, reject_position)) continue;
        valid_neighbor_list.emplace_back(reject_position);
    }

    for (auto & sample_point : sample_point_list)
    {
        const auto & sample_position{ sample_point.position };
        if (IsOwnedByNeighbor(sample_position, local_position, valid_neighbor_list) ||
            IsInsideNeighborCone(sample_position, local_position, valid_neighbor_list, angle))
        {
            sample_point.is_selected = false;
        }
    }
}

LocalPotentialSampleList BuildMedianResponseSampleEntriesByRadius(
    const LocalPotentialSampleList & sample_entries)
{
    std::map<double, std::vector<double>> response_by_radius;
    for (const auto & sample : sample_entries)
    {
        response_by_radius[sample.point.distance].emplace_back(sample.response);
    }

    LocalPotentialSampleList median_sample_entries;
    median_sample_entries.reserve(response_by_radius.size());
    for (const auto & [radius, response_list] : response_by_radius)
    {
        median_sample_entries.emplace_back(
            LocalPotentialSample{
                array_helper::ComputeMedian(response_list),
                SamplingPoint{ radius }
            }
        );
    }
    return median_sample_entries;
}

} // namespace rhbm_gem::sample_filter
