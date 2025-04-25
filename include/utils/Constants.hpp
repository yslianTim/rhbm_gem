#pragma once

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

namespace Constants
{
    inline constexpr double pi{ M_PI };
    inline constexpr double two_pi{ 2.0 * pi };
}