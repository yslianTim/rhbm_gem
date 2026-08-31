#include "MapHelper.hpp"

#include <algorithm>
#include <cstddef>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>

namespace rhbm_gem::data_internal {

namespace {

constexpr int kFloat32Mode = 2;

bool IsCanonicalAxisOrder(const std::array<int, 3> & axis_order)
{
    return axis_order[0] == 1 && axis_order[1] == 2 && axis_order[2] == 3;
}

std::array<int, 3> BuildAxisToIndex(const std::array<int, 3> & axis_order)
{
    std::array<int, 3> axis_to_index{ -1, -1, -1 };
    for (size_t i = 0; i < axis_order.size(); ++i)
    {
        const auto axis_value{ axis_order[i] };
        if (axis_value < 1 || axis_value > 3)
        {
            throw std::runtime_error("Invalid axis mapping in map header.");
        }
        auto & index{ axis_to_index[static_cast<size_t>(axis_value - 1)] };
        if (index >= 0)
        {
            throw std::runtime_error("Duplicate axis mapping in map header.");
        }
        index = static_cast<int>(i);
    }
    return axis_to_index;
}

void RequireFloat32Mode(int mode)
{
    if (mode != kFloat32Mode)
    {
        throw std::runtime_error(
            "Only float32 map mode (2) is supported; received mode "
            + std::to_string(mode) + ".");
    }
}

} // namespace

std::unique_ptr<double[]> ReorderToCanonicalXYZ(
    std::unique_ptr<double[]> raw_data,
    const std::array<int, 3> & array_size,
    const std::array<int, 3> & axis_order)
{
    if (raw_data == nullptr)
    {
        throw std::runtime_error("Cannot reorder null map data pointer.");
    }

    const auto voxel_count{ CountVoxelCount(array_size) };
    if (IsCanonicalAxisOrder(axis_order))
    {
        return raw_data;
    }

    auto reordered_array{ std::make_unique<double[]>(voxel_count) };
    const auto axis_to_index{ BuildAxisToIndex(axis_order) };

    std::array<size_t, 3> dimensions{
        static_cast<size_t>(array_size[static_cast<size_t>(axis_to_index[0])]),
        static_cast<size_t>(array_size[static_cast<size_t>(axis_to_index[1])]),
        static_cast<size_t>(array_size[static_cast<size_t>(axis_to_index[2])])
    };
    std::array<size_t, 3> source_stride{};
    size_t stride_acc{ 1 };
    for (size_t i = 0; i < axis_order.size(); ++i)
    {
        source_stride[static_cast<size_t>(axis_order[i] - 1)] = stride_acc;
        stride_acc *= static_cast<size_t>(array_size[i]);
    }

    for (size_t z = 0; z < dimensions[2]; z++)
    {
        const size_t src_off_z{ z * source_stride[2] };
        const size_t dst_off_z{ z * dimensions[0] * dimensions[1] };
        for (size_t y = 0; y < dimensions[1]; y++)
        {
            const size_t src_off_y{ src_off_z + y * source_stride[1] };
            const size_t dst_off_y{ dst_off_z + y * dimensions[0] };
            for (size_t x = 0; x < dimensions[0]; x++)
            {
                reordered_array[dst_off_y + x] = raw_data[src_off_y + x * source_stride[0]];
            }
        }
    }

    return reordered_array;
}

std::unique_ptr<double[]> ReadFloat32VoxelData(
    std::istream & stream,
    std::streamoff data_offset,
    const std::array<int, 3> & array_size,
    const std::array<int, 3> & axis_order,
    int mode)
{
    RequireFloat32Mode(mode);
    const auto voxel_count{ CountVoxelCount(array_size) };
    stream.seekg(data_offset, std::ios::beg);
    if (!stream)
    {
        throw std::runtime_error("Failed to seek to map voxel data.");
    }

    auto float32_data{ std::make_unique<float[]>(voxel_count) };
    stream.read(
        reinterpret_cast<char *>(float32_data.get()),
        static_cast<std::streamsize>(voxel_count * sizeof(float)));
    if (!stream)
    {
        throw std::runtime_error("Failed to read float32 map voxel data.");
    }
    auto raw_data{ std::make_unique<double[]>(voxel_count) };
    for (size_t i = 0; i < voxel_count; ++i)
    {
        raw_data[i] = static_cast<double>(float32_data[i]);
    }
    return ReorderToCanonicalXYZ(std::move(raw_data), array_size, axis_order);
}

void WriteFloat32VoxelData(
    std::ostream & stream,
    std::streamoff data_offset,
    const double * data,
    size_t data_size,
    const std::array<int, 3> & array_size,
    int mode)
{
    RequireFloat32Mode(mode);
    if (data == nullptr)
    {
        throw std::runtime_error("Cannot write a null map voxel array.");
    }
    const auto voxel_count{ CountVoxelCount(array_size) };
    if (data_size != voxel_count)
    {
        throw std::runtime_error("Map voxel count does not match the header dimensions.");
    }

    stream.seekp(data_offset, std::ios::beg);
    if (!stream)
    {
        throw std::runtime_error("Failed to seek to map voxel data.");
    }
    auto float32_data{ std::make_unique<float[]>(voxel_count) };
    for (size_t i = 0; i < voxel_count; ++i)
    {
        float32_data[i] = static_cast<float>(data[i]);
    }
    stream.write(
        reinterpret_cast<const char *>(float32_data.get()),
        static_cast<std::streamsize>(voxel_count * sizeof(float)));
    if (!stream)
    {
        throw std::runtime_error("Failed to write float32 map voxel data.");
    }
}

void ReorderHeaderAxesToCanonical(
    int (&array_size)[3],
    int (&location_index)[3],
    int (&grid_size)[3],
    float (&cell_dimensions)[3],
    float (&cell_angles)[3],
    int (&axis_order)[3],
    float * origin)
{
    const std::array<int, 3> source_axis_order{
        axis_order[0], axis_order[1], axis_order[2] };
    const auto axis_to_index{ BuildAxisToIndex(source_axis_order) };
    if (IsCanonicalAxisOrder(source_axis_order))
    {
        return;
    }

    const std::array<int, 3> source_array_size{
        array_size[0], array_size[1], array_size[2] };
    const std::array<int, 3> source_location_index{
        location_index[0], location_index[1], location_index[2] };
    const std::array<int, 3> source_grid_size{
        grid_size[0], grid_size[1], grid_size[2] };
    const std::array<float, 3> source_cell_dimensions{
        cell_dimensions[0], cell_dimensions[1], cell_dimensions[2] };
    const std::array<float, 3> source_cell_angles{
        cell_angles[0], cell_angles[1], cell_angles[2] };
    std::array<float, 3> source_origin{};
    if (origin != nullptr)
    {
        std::copy_n(origin, 3, source_origin.begin());
    }

    for (size_t i = 0; i < axis_to_index.size(); ++i)
    {
        const auto source_index{ static_cast<size_t>(axis_to_index[i]) };
        array_size[i] = source_array_size[source_index];
        location_index[i] = source_location_index[source_index];
        grid_size[i] = source_grid_size[source_index];
        cell_dimensions[i] = source_cell_dimensions[source_index];
        cell_angles[i] = source_cell_angles[source_index];
        if (origin != nullptr)
        {
            origin[i] = source_origin[source_index];
        }
        axis_order[i] = static_cast<int>(i + 1);
    }
}

} // namespace rhbm_gem::data_internal
