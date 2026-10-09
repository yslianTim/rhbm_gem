#pragma once

#include <cstddef>

namespace rhbm_gem::core::joint_component {

enum class FixedNeighborBlockOrder {Forward,Reverse};

// The production policy contains only the qualified route controls. Numerical
// tolerances and diagnostic switches remain implementation details.
struct FixedNeighborSearchPolicy
{
    std::size_t core_atoms{12};
    std::size_t maximum_sweeps{30};
    FixedNeighborBlockOrder order{FixedNeighborBlockOrder::Forward};
};

} // namespace rhbm_gem::core::joint_component
