#include "internal/file/MapAxisOrderHelper.hpp"

#include <cstddef>
#include <stdexcept>

namespace rhbm_gem::map_io {

namespace {

size_t CountVoxelCount(const std::array<int, 3>& array_size) {
    if (array_size[0] <= 0 || array_size[1] <= 0 || array_size[2] <= 0) {
        throw std::runtime_error("Map dimensions must be positive.");
    }
    return static_cast<size_t>(array_size[0]) *
           static_cast<size_t>(array_size[1]) *
           static_cast<size_t>(array_size[2]);
}

bool IsCanonicalAxisOrder(const std::array<int, 3>& axis_order) {
    return axis_order[0] == 1 && axis_order[1] == 2 && axis_order[2] == 3;
}

} // namespace

std::unique_ptr<float[]> ReorderToCanonicalXYZ(
    std::unique_ptr<float[]> raw_data,
    const std::array<int, 3>& array_size,
    const std::array<int, 3>& axis_order) {
    if (raw_data == nullptr) {
        throw std::runtime_error("Cannot reorder null map data pointer.");
    }

    const auto voxel_count{CountVoxelCount(array_size)};
    if (IsCanonicalAxisOrder(axis_order)) {
        return raw_data;
    }

    auto reordered_array{std::make_unique<float[]>(voxel_count)};
    std::array<int, 3> axis_to_index{-1, -1, -1};
    for (size_t i = 0; i < axis_order.size(); i++) {
        const int axis_value{axis_order[i]};
        if (axis_value < 1 || axis_value > 3) {
            throw std::runtime_error("Invalid axis mapping in map header.");
        }
        axis_to_index[static_cast<size_t>(axis_value - 1)] = static_cast<int>(i);
    }
    if (axis_to_index[0] < 0 || axis_to_index[1] < 0 || axis_to_index[2] < 0) {
        throw std::runtime_error("Incomplete axis mapping in map header.");
    }

    std::array<size_t, 3> dimensions{
        static_cast<size_t>(array_size[static_cast<size_t>(axis_to_index[0])]),
        static_cast<size_t>(array_size[static_cast<size_t>(axis_to_index[1])]),
        static_cast<size_t>(array_size[static_cast<size_t>(axis_to_index[2])])};
    std::array<size_t, 3> source_stride{};
    size_t stride_acc{1};
    for (size_t i = 0; i < axis_order.size(); ++i) {
        source_stride[static_cast<size_t>(axis_order[i] - 1)] = stride_acc;
        stride_acc *= static_cast<size_t>(array_size[i]);
    }

    for (size_t z = 0; z < dimensions[2]; z++) {
        const size_t src_off_z{z * source_stride[2]};
        const size_t dst_off_z{z * dimensions[0] * dimensions[1]};
        for (size_t y = 0; y < dimensions[1]; y++) {
            const size_t src_off_y{src_off_z + y * source_stride[1]};
            const size_t dst_off_y{dst_off_z + y * dimensions[0]};
            for (size_t x = 0; x < dimensions[0]; x++) {
                reordered_array[dst_off_y + x] = raw_data[src_off_y + x * source_stride[0]];
            }
        }
    }

    return reordered_array;
}

} // namespace rhbm_gem::map_io
