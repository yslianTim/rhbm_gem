#pragma once

#include <array>
#include <memory>

namespace rhbm_gem::map_io {

std::unique_ptr<float[]> ReorderToCanonicalXYZ(
    std::unique_ptr<float[]> raw_data,
    const std::array<int, 3>& array_size,
    const std::array<int, 3>& axis_order);

} // namespace rhbm_gem::map_io
