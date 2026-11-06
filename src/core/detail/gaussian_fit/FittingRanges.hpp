#pragma once

namespace rhbm_gem::core::detail {

// Radial ranges in angstroms. Signal and tail may overlap or leave a gap.
inline constexpr double kSignalDistanceMax{ 1.0 };
inline constexpr double kTailDistanceMin{ 1.2 };
inline constexpr double kTailDistanceMax{ 2.0 };

static_assert(kSignalDistanceMax >= 0.0);
static_assert(kTailDistanceMin >= 0.0 && kTailDistanceMin <= kTailDistanceMax);

constexpr bool IsSignalDistance(double distance)
{
    return distance >= 0.0 && distance <= kSignalDistanceMax;
}

constexpr bool IsTailDistance(double distance)
{
    return distance >= kTailDistanceMin && distance <= kTailDistanceMax;
}

} // namespace rhbm_gem::core::detail
