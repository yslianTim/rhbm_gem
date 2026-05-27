#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>

#include <rhbm_gem/utils/math/NumericValidation.hpp>

namespace rhbm_gem::data_internal {

inline size_t CountVoxelCount(const std::array<int, 3> & array_size)
{
    try
    {
        numeric_validation::RequireAllPositive(array_size, "map dimensions");
    }
    catch (const std::invalid_argument &)
    {
        throw std::runtime_error("Map dimensions must be positive.");
    }
    return static_cast<size_t>(array_size[0]) *
           static_cast<size_t>(array_size[1]) *
           static_cast<size_t>(array_size[2]);
}

std::unique_ptr<float[]> ReorderToCanonicalXYZ(
    std::unique_ptr<float[]> raw_data,
    const std::array<int, 3> & array_size,
    const std::array<int, 3> & axis_order);

} // namespace rhbm_gem::data_internal
