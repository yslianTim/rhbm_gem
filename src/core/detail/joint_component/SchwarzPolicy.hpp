#pragma once

#include <cstddef>

namespace rhbm_gem::core::joint_component {

struct SchwarzPolicy
{
    std::size_t core_atoms{128},overlap_hops{1},max_block_atoms{512};
    std::size_t storage_bytes{512ULL*1024*1024},scratch_bytes{256ULL*1024*1024};
};

} // namespace rhbm_gem::core::joint_component
