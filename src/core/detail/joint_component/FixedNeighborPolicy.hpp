#pragma once

#include <cstddef>
#include <functional>
#include <string>

namespace rhbm_gem::core::joint_component {

enum class FixedNeighborBlockOrder {Forward,Reverse};
struct FixedNeighborPolicy;

// The search controls are the qualified production policy. The component
// identity and census callback are supplied internally for failure-only logs.
struct FixedNeighborSearchPolicy
{
    std::size_t core_atoms{12};
    std::size_t maximum_sweeps{200};
    FixedNeighborBlockOrder order{FixedNeighborBlockOrder::Forward};
    std::size_t component_index{};
    std::string component_id;
    std::function<void(std::size_t &,std::size_t &,std::size_t &)> halo_census;
    std::function<void(FixedNeighborPolicy &)> configure_search_for_testing;
};

} // namespace rhbm_gem::core::joint_component
