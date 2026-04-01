#pragma once

#include <cmath>
#include <map>
#include <tuple>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/utils/math/GausLinearTransformHelper.hpp>

namespace rhbm_gem::series_ops {

inline std::vector<std::tuple<float, float>> BuildLinearModelDistanceAndMapValueList(
    const LocalPotentialEntry & entry)
{
    Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
    const auto & data_array{ entry.GetDistanceAndMapValueList() };
    std::vector<std::tuple<float, float>> linear_model_distance_and_map_value_list;
    linear_model_distance_and_map_value_list.reserve(data_array.size());
    for (const auto & [distance, map_value] : data_array)
    {
        if (map_value <= 0.0f) continue;
        const auto data_vector{
            GausLinearTransformHelper::BuildLinearModelDataVector(
                distance, map_value, model_par_init)
        };
        linear_model_distance_and_map_value_list.emplace_back(
            std::make_tuple(data_vector(1), data_vector(2)));
    }
    return linear_model_distance_and_map_value_list;
}

inline std::vector<std::tuple<float, float>> BuildBinnedDistanceAndMapValueList(
    const LocalPotentialEntry & entry,
    int bin_size = 15,
    double x_min = 0.0,
    double x_max = 1.5)
{
    const auto bin_spacing{ (x_max - x_min) / static_cast<double>(bin_size) };
    std::map<int, std::vector<float>> bin_map;
    for (const auto & [distance, map_value] : entry.GetDistanceAndMapValueList())
    {
        const auto bin_index{ static_cast<int>(std::floor(distance / bin_spacing)) };
        bin_map[bin_index].emplace_back(map_value);
    }

    std::vector<std::tuple<float, float>> binned_distance_and_map_value_list;
    binned_distance_and_map_value_list.reserve(static_cast<size_t>(bin_size));
    for (int i = 0; i < bin_size; i++)
    {
        const auto x_value{ static_cast<float>(x_min + (i + 0.5) * bin_spacing) };
        const auto y_value{ (bin_map.find(i) == bin_map.end()) ?
            0.0f : ArrayStats<float>::ComputeMedian(bin_map.at(i))
        };
        binned_distance_and_map_value_list.emplace_back(std::make_tuple(x_value, y_value));
    }
    return binned_distance_and_map_value_list;
}

inline std::tuple<float, float> ComputeDistanceRange(
    const LocalPotentialEntry & entry,
    double margin_rate = 0.0)
{
    const auto & data_array{ entry.GetDistanceAndMapValueList() };
    std::vector<float> distance_array;
    distance_array.reserve(data_array.size());
    for (const auto & [distance, map_value] : data_array)
    {
        (void)map_value;
        distance_array.emplace_back(distance);
    }
    return ArrayStats<float>::ComputeScalingRangeTuple(
        distance_array, static_cast<float>(margin_rate));
}

inline std::tuple<float, float> ComputeMapValueRange(
    const LocalPotentialEntry & entry,
    double margin_rate = 0.0)
{
    const auto & data_array{ entry.GetDistanceAndMapValueList() };
    std::vector<float> map_value_array;
    map_value_array.reserve(data_array.size());
    for (const auto & [distance, map_value] : data_array)
    {
        (void)distance;
        map_value_array.emplace_back(map_value);
    }
    return ArrayStats<float>::ComputeScalingRangeTuple(
        map_value_array, static_cast<float>(margin_rate));
}

inline std::tuple<float, float> ComputeBinnedMapValueRange(
    const LocalPotentialEntry & entry,
    int bin_size = 15,
    double x_min = 0.0,
    double x_max = 1.5,
    double margin_rate = 0.0)
{
    const auto data_array{ BuildBinnedDistanceAndMapValueList(entry, bin_size, x_min, x_max) };
    std::vector<float> map_value_array;
    map_value_array.reserve(data_array.size());
    for (const auto & [distance, map_value] : data_array)
    {
        (void)distance;
        map_value_array.emplace_back(map_value);
    }
    return ArrayStats<float>::ComputeScalingRangeTuple(
        map_value_array, static_cast<float>(margin_rate));
}

} // namespace rhbm_gem::series_ops
