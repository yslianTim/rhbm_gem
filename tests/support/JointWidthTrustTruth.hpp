#pragma once

#include <cmath>
#include <stdexcept>

namespace joint_test::width_trust {

struct Truth
{
    double amplitude{},width{},charge{};
};

inline Truth ToNormalizedDomain(Truth truth,double divisor)
{
    if(!(divisor>0.0) || !std::isfinite(divisor))
        throw std::invalid_argument("Invalid map normalization divisor.");
    truth.amplitude/=divisor;
    truth.charge/=divisor;
    return truth;
}

inline Truth ToOriginalDomain(Truth truth,double divisor)
{
    if(!(divisor>0.0) || !std::isfinite(divisor))
        throw std::invalid_argument("Invalid map normalization divisor.");
    truth.amplitude*=divisor;
    truth.charge*=divisor;
    return truth;
}

} // namespace joint_test::width_trust
