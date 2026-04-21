#pragma once

#include <array>
#include <cstddef>

#include <rhbm_gem/utils/math/NumericValidation.hpp>

namespace rhbm_gem::map_io {

inline size_t CountVoxelCount(const std::array<int, 3> & array_size)
{
    NumericValidation::RequireAllPositive(array_size, "map dimensions");
    return static_cast<size_t>(array_size[0]) *
           static_cast<size_t>(array_size[1]) *
           static_cast<size_t>(array_size[2]);
}

} // namespace rhbm_gem::map_io
