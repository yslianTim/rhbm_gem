#pragma once

#include "IPlotDataBuilder.hpp"

#include <rhbm_gem/data/object/PotentialEntryQuery.hpp>

namespace rhbm_gem {

class QueryPlotDataBuilder final : public IPlotDataBuilder
{
public:
    explicit QueryPlotDataBuilder(const PotentialEntryQuery & query) :
        m_query{ &query }
    {
    }

    std::vector<std::tuple<float, float>> GetLinearModelDistanceAndMapValueList() const override
    {
        return m_query->GetLinearModelDistanceAndMapValueList();
    }

    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList() const override
    {
        return m_query->GetDistanceAndMapValueList();
    }

    std::vector<std::tuple<float, float>> GetBinnedDistanceAndMapValueList(
        int bin_size,
        double x_min,
        double x_max) const override
    {
        return m_query->GetBinnedDistanceAndMapValueList(bin_size, x_min, x_max);
    }

private:
    const PotentialEntryQuery * m_query;
};

} // namespace rhbm_gem
