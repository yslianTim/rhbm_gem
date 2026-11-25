#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

namespace rhbm_gem::core::joint_component {

struct WidthBounds
{
    std::optional<double> upper, lower;

    static double RoundoffTolerance(double eta)
    {
        return 64.0*std::numeric_limits<double>::epsilon()*std::max(1.0,std::abs(eta));
    }

    bool IsValid() const
    {
        return upper && std::isfinite(*upper) &&
            (!lower || (std::isfinite(*lower) && *lower < *upper));
    }

    bool Contains(double eta) const
    {
        return upper && eta <= *upper && (!lower || eta >= *lower);
    }

    bool ContainsWithRoundoff(double eta) const
    {
        return upper && *upper-eta>=-RoundoffTolerance(*upper) &&
            (!lower || eta-*lower>=-RoundoffTolerance(*lower));
    }
};

} // namespace rhbm_gem::core::joint_component
