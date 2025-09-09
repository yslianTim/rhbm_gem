#pragma once

#include <array>
#include <tuple>
#include <vector>

class SamplerBase
{
public:
    virtual ~SamplerBase() = default;

    virtual void Print(void) const = 0;
    virtual std::vector<std::tuple<float, std::array<float, 3>>> GenerateSamplingPoints(
        const std::array<float, 3> & position,
        const std::array<float, 3> & axis_vector) const = 0;
};
