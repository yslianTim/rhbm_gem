#pragma once

#include <array>
#include <cstddef>
#include <iosfwd>
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

std::unique_ptr<double[]> ReorderToCanonicalXYZ(
    std::unique_ptr<double[]> raw_data,
    const std::array<int, 3> & array_size,
    const std::array<int, 3> & axis_order);

std::unique_ptr<double[]> ReadFloat32VoxelData(
    std::istream & stream,
    std::streamoff data_offset,
    const std::array<int, 3> & array_size,
    const std::array<int, 3> & axis_order,
    int mode);

void WriteFloat32VoxelData(
    std::ostream & stream,
    std::streamoff data_offset,
    const double * data,
    size_t data_size,
    const std::array<int, 3> & array_size,
    int mode);

void ReorderHeaderAxesToCanonical(
    int (&array_size)[3],
    int (&location_index)[3],
    int (&grid_size)[3],
    float (&cell_dimensions)[3],
    float (&cell_angles)[3],
    int (&axis_order)[3],
    float * origin = nullptr);

} // namespace rhbm_gem::data_internal
