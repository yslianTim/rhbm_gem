#pragma once

#include <array>
#include <cmath>

namespace rhbm_gem::core::simulation {
// Versioned arithmetic for simulation and immutable experimental domains.
inline std::array<double, 3> GridPosition(const std::array<int, 3> & index,
    const std::array<double, 3> & spacing, const std::array<double, 3> & origin)
{
    return { std::fma(static_cast<double>(index[0]), spacing[0], origin[0]),
        std::fma(static_cast<double>(index[1]), spacing[1], origin[1]),
        std::fma(static_cast<double>(index[2]), spacing[2], origin[2]) };
}

inline double SupportSquare(const std::array<double, 3> & a, const std::array<double, 3> & b)
{
    const double x{ a[0] - b[0] }, y{ a[1] - b[1] }, z{ a[2] - b[2] };
    return std::fma(z, z, std::fma(y, y, x * x));
}
} // namespace rhbm_gem::core::simulation
