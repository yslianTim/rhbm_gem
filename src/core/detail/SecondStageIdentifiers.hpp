#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace rhbm_gem::core::detail {

using ClusterKey = std::vector<std::size_t>;
using ResidueKey = std::pair<std::string, int>;

struct SampleRef
{
    std::size_t atom_index{ 0 };
    std::size_t sample_index{ 0 };
};

} // namespace rhbm_gem::core::detail
