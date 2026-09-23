#pragma once
#include <rhbm_gem/data/object/MapObject.hpp>
#include <algorithm>
#include <array>
#include <map>

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
    std::size_t Index(const std::array<int, 3> & p) const
    {
        return static_cast<std::size_t>(p[0]) + static_cast<std::size_t>(dimensions[0]) *
            (static_cast<std::size_t>(p[1]) + static_cast<std::size_t>(dimensions[1]) * static_cast<std::size_t>(p[2]));
    }
};
inline TricubicStencil MakeTricubicStencil(const MapObject & map, const std::array<double, 3> & position)
{
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

// Merge clamped duplicate nodes before testing whether a node affects a sample.
inline std::map<std::size_t, double> TricubicWeights(const TricubicStencil & stencil)
{
    std::array<std::array<double, 4>, 3> weights;
    for (std::size_t axis = 0; axis < 3; ++axis)
        for (std::size_t node = 0; node < 4; ++node)
            weights[axis][node] = CubicInterpolate(node == 0, node == 1, node == 2, node == 3, stencil.local[axis]);
    std::map<std::size_t, double> result;
    for (int i = -1; i < 3; ++i) for (int j = -1; j < 3; ++j) for (int k = -1; k < 3; ++k)
        result[stencil.Index(stencil.Node(i, j, k))] += weights[0][static_cast<std::size_t>(i + 1)] * weights[1][static_cast<std::size_t>(j + 1)] * weights[2][static_cast<std::size_t>(k + 1)];
    return result;
}
} // namespace rhbm_gem::core::detail
