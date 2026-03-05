#pragma once

#include <array>
#include <cstddef>
#include <memory>

namespace rhbm_gem::map_io {

bool IsCanonicalAxisOrder(const std::array<int, 3> & axis_order);

std::unique_ptr<float[]> ReorderToCanonicalXYZ(
    const float * raw_data,
    const std::array<int, 3> & array_size,
    const std::array<int, 3> & axis_order);

size_t CountVoxelCount(const std::array<int, 3> & array_size);

} // namespace rhbm_gem::map_io
