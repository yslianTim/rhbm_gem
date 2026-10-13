#pragma once

#include <cstddef>
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <functional>
#endif

namespace rhbm_gem::core::joint_component {

enum class FixedNeighborBlockOrder {Forward,Reverse};
struct FixedNeighborPolicy;

// The production policy contains only the qualified policy controls. Numerical
// tolerances and diagnostic switches remain implementation details.
struct FixedNeighborSearchPolicy
{
    std::size_t core_atoms{12};
    std::size_t maximum_sweeps{200};
    FixedNeighborBlockOrder order{FixedNeighborBlockOrder::Forward};
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    std::function<void(FixedNeighborPolicy &)> configure_search_for_testing;
#endif
};

} // namespace rhbm_gem::core::joint_component
