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

double ComputeQScore(
    const LocalPotentialSampleList & sampling_entries,
    const LocalGaussianResult & gaussian_result,
    QScoreReference reference)
{
    if (sampling_entries.empty())
    {
        return 0.0;
    }

    GaussianModel3D reference_model;
    bool has_reference_model{ false };
    if (reference == QScoreReference::Fixed)
    {
        constexpr double peak_intensity{ 0.05 };
        constexpr double width{ 0.6 };
        constexpr double offset{ -0.005 };
        const auto unit_peak_response{
            GaussianModel3D{ 1.0, width, 0.0 }.ResponseAtDistance(0.0)
        };
        if (std::isfinite(unit_peak_response) && unit_peak_response != 0.0)
        {
            reference_model = GaussianModel3D{
                peak_intensity / unit_peak_response,
                width,
                offset
            };
            has_reference_model = true;
        }
    }
    else if (reference == QScoreReference::MDPDE)
    {
        const auto & estimate{ gaussian_result.mdpde.GetModel() };
        reference_model = estimate.WithOffset(0.0);
        has_reference_model = true;
    }

    if (!has_reference_model ||
        std::isfinite(reference_model.GetAmplitude()) == false ||
        std::isfinite(reference_model.GetWidth()) == false ||
        std::isfinite(reference_model.GetOffset()) == false ||
        reference_model.GetWidth() <= 0.0)
    {
        return 0.0;
    }

    std::vector<float> map_value_list;
    std::vector<float> reference_value_list;
    map_value_list.reserve(sampling_entries.size());
    reference_value_list.reserve(sampling_entries.size());
    for (const auto & sample : sampling_entries)
    {
        map_value_list.emplace_back(sample.response);
        reference_value_list.emplace_back(
            static_cast<float>(
                reference_model.ResponseAtDistance(
                    static_cast<double>(sample.point.distance))));
    }
    if (map_value_list.empty())
    {
        return 0.0;
    }

    float map_value_mean{ array_helper::ComputeMean(map_value_list.data(), map_value_list.size()) };
    float reference_value_mean{
        array_helper::ComputeMean(reference_value_list.data(), reference_value_list.size())
    };

    auto numerator{ 0.0 };
    auto map_value_norm_square{ 0.0 };
    auto reference_value_norm_square{ 0.0 };
    for (size_t i = 0; i < map_value_list.size(); i++)
    {
        auto map_value_diff{ map_value_list.at(i) - map_value_mean };
        auto reference_value_diff{ reference_value_list.at(i) - reference_value_mean };
        numerator += map_value_diff * reference_value_diff;
        map_value_norm_square += map_value_diff * map_value_diff;
        reference_value_norm_square += reference_value_diff * reference_value_diff;
    }
    const auto denominator{
        std::sqrt(map_value_norm_square) * std::sqrt(reference_value_norm_square)
    };
    if (denominator == 0.0 ||
        std::isfinite(denominator) == false ||
        std::isfinite(numerator) == false)
    {
        return 0.0;
    }

    const auto q_score{ numerator/denominator };
    return std::isfinite(q_score) ? q_score : 0.0;
}

} // namespace rhbm_gem::local_potential_series
