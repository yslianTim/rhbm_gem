#pragma once

#include <cmath>
#include <map>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/utils/math/GausLinearTransformHelper.hpp>

namespace rhbm_gem {

class LocalPotentialView
{
    const LocalPotentialEntry & m_local_entry;

public:
    explicit LocalPotentialView(const AtomObject & atom_object) :
        m_local_entry(*RequireLocalEntry(atom_object.GetLocalPotentialEntry(), "Atom local entry")) {}

    explicit LocalPotentialView(const BondObject & bond_object) :
        m_local_entry(*RequireLocalEntry(bond_object.GetLocalPotentialEntry(), "Bond local entry")) {}

    explicit LocalPotentialView(const LocalPotentialEntry & local_entry) :
        m_local_entry(local_entry) {}

    const LocalPotentialEntry & GetLocalEntry() const { return m_local_entry; }

    bool IsOutlier(const std::string & class_key) const
    {
        return m_local_entry.GetOutlierTag(class_key);
    }

    std::vector<std::tuple<float, float>> GetLinearModelDistanceAndMapValueList() const
    {
        Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
        const auto & data_array{ GetDistanceAndMapValueList() };
        std::vector<std::tuple<float, float>> linear_model_distance_and_map_value_list;
        linear_model_distance_and_map_value_list.reserve(data_array.size());
        for (auto & [distance, map_value] : data_array)
        {
            if (map_value <= 0.0f) continue;
            auto data_vector{
                GausLinearTransformHelper::BuildLinearModelDataVector(
                    distance, map_value, model_par_init)
            };
            linear_model_distance_and_map_value_list.emplace_back(
                std::make_tuple(data_vector(1), data_vector(2)));
        }
        return linear_model_distance_and_map_value_list;
    }

    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList() const
    {
        return m_local_entry.GetDistanceAndMapValueList();
    }

    std::vector<std::tuple<float, float>> GetBinnedDistanceAndMapValueList(
        int bin_size=15,
        double x_min=0.0,
        double x_max=1.5) const
    {
        auto bin_spacing{ (x_max - x_min) / static_cast<double>(bin_size) };
        std::map<int, std::vector<float>> bin_map;
        for (auto & [distance, map_value] : GetDistanceAndMapValueList())
        {
            auto bin_index{ static_cast<int>(std::floor(distance / bin_spacing)) };
            bin_map[bin_index].emplace_back(map_value);
        }

        std::vector<std::tuple<float, float>> binned_distance_and_map_value_list;
        binned_distance_and_map_value_list.reserve(static_cast<size_t>(bin_size));
        for (int i = 0; i < bin_size; i++)
        {
            auto x_value{ static_cast<float>(x_min + (i + 0.5) * bin_spacing) };
            auto y_value{ (bin_map.find(i) == bin_map.end()) ?
                0.0f : ArrayStats<float>::ComputeMedian(bin_map.at(i))
            };
            binned_distance_and_map_value_list.emplace_back(std::make_tuple(x_value, y_value));
        }
        return binned_distance_and_map_value_list;
    }

    std::tuple<float, float> GetDistanceRange(double margin_rate=0.0) const
    {
        const auto & data_array{ GetDistanceAndMapValueList() };
        std::vector<float> distance_array;
        distance_array.reserve(static_cast<size_t>(data_array.size()));
        for (auto & [distance, map_value] : data_array)
        {
            (void)map_value;
            distance_array.emplace_back(distance);
        }
        return ArrayStats<float>::ComputeScalingRangeTuple(
            distance_array, static_cast<float>(margin_rate));
    }

    std::tuple<float, float> GetMapValueRange(double margin_rate=0.0) const
    {
        const auto & data_array{ GetDistanceAndMapValueList() };
        std::vector<float> map_value_array;
        map_value_array.reserve(data_array.size());
        for (auto & [distance, map_value] : data_array)
        {
            (void)distance;
            map_value_array.emplace_back(map_value);
        }
        return ArrayStats<float>::ComputeScalingRangeTuple(
            map_value_array, static_cast<float>(margin_rate));
    }

    std::tuple<float, float> GetBinnedMapValueRange(
        int bin_size=15,
        double x_min=0.0,
        double x_max=1.5,
        double margin_rate=0.0) const
    {
        auto data_array{ GetBinnedDistanceAndMapValueList(bin_size, x_min, x_max) };
        std::vector<float> map_value_array;
        map_value_array.reserve(data_array.size());
        for (auto & [distance, map_value] : data_array)
        {
            (void)distance;
            map_value_array.emplace_back(map_value);
        }
        return ArrayStats<float>::ComputeScalingRangeTuple(
            map_value_array, static_cast<float>(margin_rate));
    }

    const GaussianEstimate & GetEstimateMDPDE() const { return m_local_entry.GetEstimateMDPDE(); }
    const GaussianPosterior & GetPosterior(const std::string & key) const { return m_local_entry.GetPosterior(key); }
    double GetAlphaR() const { return m_local_entry.GetAlphaR(); }

private:
    static const LocalPotentialEntry * RequireLocalEntry(
        const LocalPotentialEntry * entry,
        const char * context)
    {
        if (entry == nullptr)
        {
            throw std::runtime_error(std::string(context) + " is not available.");
        }
        return entry;
    }
};

} // namespace rhbm_gem
