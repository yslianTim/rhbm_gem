#include <rhbm_gem/utils/hrl/LocalPotentialSeries.hpp>

#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include <cmath>
#include <vector>

namespace rhbm_gem::local_potential_series
{

std::tuple<float, float> ComputeDistanceRange(
    const LocalPotentialSampleList & sampling_entries,
    double margin_rate)
{
    std::vector<float> distance_array;
    distance_array.reserve(sampling_entries.size());
    for (const auto & sample : sampling_entries)
    {
        distance_array.emplace_back(sample.point.distance);
    }
    return array_helper::ComputeScalingRangeTuple(
        distance_array, static_cast<float>(margin_rate));
}

std::tuple<float, float> ComputeResponseRange(
    const LocalPotentialSampleList & sampling_entries,
    double margin_rate)
{
    std::vector<float> map_value_array;
    map_value_array.reserve(sampling_entries.size());
    for (const auto & sample : sampling_entries)
    {
        map_value_array.emplace_back(sample.response);
    }
    return array_helper::ComputeScalingRangeTuple(
        map_value_array, static_cast<float>(margin_rate));
}

SeriesPointList BuildBinnedDistanceResponseSeries(
    const LocalPotentialSampleList & sampling_entries,
    int bin_size,
    double x_min,
    double x_max)
{
    const auto bin_spacing{ (x_max - x_min) / static_cast<double>(bin_size) };
    std::vector<std::vector<float>> bin_map(static_cast<size_t>(bin_size));
    for (const auto & sample : sampling_entries)
    {
        const auto distance{ sample.point.distance };
        if (distance < x_min || distance >= x_max)
        {
            continue;
        }
        const auto shifted_distance{ static_cast<double>(distance) - x_min };
        const auto bin_index{ static_cast<int>(std::floor(shifted_distance / bin_spacing)) };
        if (bin_index < 0 || bin_index >= bin_size)
        {
            continue;
        }
        bin_map.at(static_cast<size_t>(bin_index)).emplace_back(sample.response);
    }

    SeriesPointList binned_distance_response_series;
    binned_distance_response_series.reserve(static_cast<size_t>(bin_size));
    for (int i = 0; i < bin_size; i++)
    {
        const auto x_value{ static_cast<float>(x_min + (i + 0.5) * bin_spacing) };
        const auto & bin_values{ bin_map.at(static_cast<size_t>(i)) };
        const auto y_value{ bin_values.empty() ? 0.0f : array_helper::ComputeMedian(bin_values) };
        binned_distance_response_series.emplace_back(SeriesPoint{ { x_value }, y_value });
    }
    return binned_distance_response_series;
}

double ComputeMapValueNearCenter(const LocalPotentialSampleList & sampling_entries)
{
    if (sampling_entries.empty()) return 0.0;
    int count{ 0 };
    double sum{ 0.0 };
    for (const auto & sample : sampling_entries)
    {
        if (sample.point.distance > 0.05f) continue;
        sum += sample.response;
        count++;
    }
    if (count == 0) return 0.0;
    return sum/static_cast<double>(count);
}

} // namespace rhbm_gem::local_potential_series
