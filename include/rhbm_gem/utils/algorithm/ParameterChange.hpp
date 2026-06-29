#pragma once

#include <algorithm>
#include <vector>

namespace rhbm_gem::algorithm {

struct ParameterChange
{
    std::vector<double> value_list{};
};

inline double GetMaximumParameterChange(const ParameterChange & change)
{
    if (change.value_list.empty()) return 0.0;
    return *std::max_element(change.value_list.begin(), change.value_list.end());
}

} // namespace rhbm_gem::algorithm
