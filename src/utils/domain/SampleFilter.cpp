#include <rhbm_gem/utils/domain/SampleFilter.hpp>

#include <cmath>
#include <cstddef>
#include <map>
#include <utility>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

namespace rhbm_gem::sample_filter {
namespace {

struct ResponseRankedLocalPotentialSample
{
    LocalPotentialSample sample{};
    std::size_t sequence{ 0 };

    bool operator<(const ResponseRankedLocalPotentialSample & other) const
    {
        if (sample.response != other.sample.response)
        {
            return sample.response < other.sample.response;
        }
        return sequence < other.sequence;
    }
};

bool IsOwnedByNeighbor(
    const std::array<float, 3> & sample_position,
    const std::array<float, 3> & local_position,
    const std::vector<std::array<float, 3>> & valid_neighbor_list)
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
    const std::array<float, 3> & sample_position,
    const std::array<float, 3> & local_position,
    const std::vector<std::array<float, 3>> & valid_neighbor_list,
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

LocalPotentialSampleList KeepLowestResponseRatioByDistance(
    LocalPotentialSampleList sample_list,
    double retained_ratio)
{
    using ResponseRankedLocalPotentialSampleList = std::vector<ResponseRankedLocalPotentialSample>;

    std::map<float, ResponseRankedLocalPotentialSampleList> samples_by_distance;
    for (auto & sample : sample_list)
    {
        auto & distance_samples{ samples_by_distance[sample.point.distance] };
        distance_samples.emplace_back(
            ResponseRankedLocalPotentialSample{ std::move(sample), distance_samples.size() });
    }

    LocalPotentialSampleList retained_samples;
    for (auto & distance_entry : samples_by_distance)
    {
        auto retained_distance_samples{
            array_helper::ComputeSmallestProportionValues(distance_entry.second, retained_ratio)
        };
        for (auto & retained_sample : retained_distance_samples)
        {
            retained_samples.emplace_back(std::move(retained_sample.sample));
        }
    }

    return retained_samples;
}

} // namespace

void FilterSamplingPointList(
    SamplingPointList & sample_point_list,
    const std::array<float, 3> & local_position,
    const std::vector<std::array<float, 3>> & reject_position_list,
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

    std::vector<std::array<float, 3>> valid_neighbor_list;
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

LocalPotentialSampleList FilterLocalPotentialSampleList(LocalPotentialSampleList sample_list)
{
    constexpr double kRetainedRatio{ 0.5 };
    return KeepLowestResponseRatioByDistance(std::move(sample_list), kRetainedRatio);
}

} // namespace rhbm_gem::sample_filter
