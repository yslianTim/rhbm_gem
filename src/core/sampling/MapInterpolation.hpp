#pragma once
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <algorithm>
#include <array>
#include <cstddef>
#include <cmath>
#include <stdexcept>

namespace rhbm_gem::core::detail {

struct TricubicStencil
{
    std::array<int, 3> base, dimensions;
    std::array<double, 3> local;
    std::array<int, 3> Node(int i, int j, int k) const
    {
        return {std::clamp(base[0] + i, 0, dimensions[0] - 1),
            std::clamp(base[1] + j, 0, dimensions[1] - 1),
            std::clamp(base[2] + k, 0, dimensions[2] - 1)};
    }

};
inline TricubicStencil MakeTricubicStencil(const MapObject & map, const std::array<double, 3> & position)
{
    for (const auto coordinate : position)
        if (!std::isfinite(coordinate))
            throw std::out_of_range("Map interpolation position is nonfinite.");
    TricubicStencil stencil{map.GetIndexFromPosition(position), map.GetGridSize(), {}};
    const auto origin = map.GetOrigin(), spacing = map.GetGridSpacing();
    for (std::size_t i = 0; i < 3; ++i)
        stencil.local[i] = (position[i] - origin[i] - static_cast<double>(stencil.base[i]) * spacing[i]) / spacing[i];
    return stencil;
}
inline constexpr auto CubicInterpolate = [](double p0, double p1, double p2, double p3, double t)
    {
        double a0{ p1 };
        double a1{ 0.5 * (p2 - p0) };
        double a2{ 0.5 * (-p3 + 4.0 * p2 - 5.0 * p1 + 2.0 * p0) };
        double a3{ 0.5 * (p3 - 3.0 * p2 + 3.0 * p1 - p0) };
        return a3 * t * t * t + a2 * t * t + a1 * t + a0;
    };


template<class Read>
double InterpolateTricubic(const TricubicStencil & stencil, Read read)
{
    const auto & local = stencil.local;
    std::array<std::array<std::array<double, 4>, 4>, 4> values;
    for (int i = -1; i < 3; ++i)
        for (int j = -1; j < 3; ++j)
            for (int k = -1; k < 3; ++k)
                values[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(j + 1)][static_cast<std::size_t>(k + 1)] = read(stencil.Node(i, j, k));
    std::array<std::array<double, 4>, 4> temp_y;
    for (size_t j = 0; j < 4; j++)
    {
        for (size_t k = 0; k < 4; k++)
        {
            temp_y[j][k] = CubicInterpolate(
                values[0][j][k], values[1][j][k], values[2][j][k], values[3][j][k], local.at(0));
        }
    }

    std::array<double, 4> temp_z;
    for (size_t k = 0; k < 4; k++)
    {
        temp_z[k] = CubicInterpolate(
            temp_y[0][k], temp_y[1][k], temp_y[2][k], temp_y[3][k], local.at(1));
    }

    return CubicInterpolate(temp_z[0], temp_z[1], temp_z[2], temp_z[3], local.at(2));
}

inline double InterpolateTrilinear(const MapObject & map, const std::array<double, 3> & position)
{
    const auto grid_size{ map.GetGridSize() };
    const auto grid_spacing{ map.GetGridSpacing() };
    const auto origin{ map.GetOrigin() };

    std::array<int, 3> lower_indices;
    std::array<int, 3> upper_indices;
    std::array<double, 3> fractions;
    for (std::size_t axis = 0; axis < position.size(); ++axis)
    {
        const auto lower_boundary{
            origin.at(axis) - 0.5 * grid_spacing.at(axis)
        };
        const auto upper_boundary{
            origin.at(axis) +
                (static_cast<double>(grid_size.at(axis)) - 0.5) *
                    grid_spacing.at(axis)
        };
        const auto coordinate{ position.at(axis) };
        if (!std::isfinite(coordinate) ||
            coordinate < lower_boundary ||
            coordinate > upper_boundary)
        {
            throw std::out_of_range(
                "Map interpolation position is outside the map boundary.");
        }

        const auto grid_coordinate{
            std::clamp(
                (coordinate - origin.at(axis)) / grid_spacing.at(axis),
                0.0,
                static_cast<double>(grid_size.at(axis) - 1))
        };
        lower_indices.at(axis) = static_cast<int>(std::floor(grid_coordinate));
        upper_indices.at(axis) = std::min(lower_indices.at(axis) + 1, grid_size.at(axis) - 1);
        fractions.at(axis) = grid_coordinate - static_cast<double>(lower_indices.at(axis));
    }

    const auto interpolate{
        [](double lower_value, double upper_value, double fraction)
        {
            return lower_value + fraction * (upper_value - lower_value);
        }
    };

    std::array<std::array<double, 2>, 2> x_interpolated_values{};
    for (std::size_t y_offset = 0; y_offset < 2; ++y_offset)
    {
        const auto y_index{
            y_offset == 0 ? lower_indices.at(1) : upper_indices.at(1)
        };
        for (std::size_t z_offset = 0; z_offset < 2; ++z_offset)
        {
            const auto z_index{
                z_offset == 0 ? lower_indices.at(2) : upper_indices.at(2)
            };
            const auto lower_value{
                map.GetMapValue(lower_indices.at(0), y_index, z_index)
            };
            const auto upper_value{
                map.GetMapValue(upper_indices.at(0), y_index, z_index)
            };
            x_interpolated_values.at(y_offset).at(z_offset) = interpolate(
                lower_value,
                upper_value,
                fractions.at(0));
        }
    }

    std::array<double, 2> y_interpolated_values{};
    for (std::size_t z_offset = 0; z_offset < 2; ++z_offset)
    {
        y_interpolated_values.at(z_offset) = interpolate(
            x_interpolated_values.at(0).at(z_offset),
            x_interpolated_values.at(1).at(z_offset),
            fractions.at(1));
    }
    return interpolate(
        y_interpolated_values.at(0),
        y_interpolated_values.at(1),
        fractions.at(2));
}

inline LocalPotentialSampleList BuildLocalPotentialSampleList(
    const MapObject & map_object,
    const SamplingPointList & sample_point_list)
{
    LocalPotentialSampleList sampling_data_list;
    sampling_data_list.reserve(sample_point_list.size());
    for (const auto & sampling_point : sample_point_list)
    {
        const auto map_value{InterpolateTricubic(MakeTricubicStencil(map_object, sampling_point.position),
            [&](const auto & node) { return map_object.GetMapValue(node[0], node[1], node[2]); })};
        sampling_data_list.emplace_back(LocalPotentialSample{
            map_value,
            sampling_point
        });
    }
    return sampling_data_list;
}
} // namespace rhbm_gem::core::detail
