#pragma once

#include <tuple>
#include <vector>

namespace rhbm_gem {

class IPlotDataBuilder
{
public:
    virtual ~IPlotDataBuilder() = default;

    virtual std::vector<std::tuple<float, float>> GetLinearModelDistanceAndMapValueList() const = 0;
    virtual const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList() const = 0;
    virtual std::vector<std::tuple<float, float>> GetBinnedDistanceAndMapValueList(
        int bin_size,
        double x_min,
        double x_max) const = 0;
};

} // namespace rhbm_gem
