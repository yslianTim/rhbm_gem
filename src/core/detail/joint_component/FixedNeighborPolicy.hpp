#pragma once

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include "ConvergenceCriteria.hpp"

namespace rhbm_gem::core::joint_component {

enum class FixedNeighborBlockOrder {Forward,Reverse};
struct FixedNeighborPolicy;
struct BlockCoordinateState;
struct FixedNeighborBlockSweep;
struct FixedNeighborSearchResult;

// The search controls are the qualified production policy. The component
// identity and census callback are supplied internally for failure-only logs.
struct FixedNeighborSearchPolicy
{
    std::size_t core_atoms{12};
    std::size_t maximum_sweeps{200};
    JointConvergenceCriteria convergence_criteria{PracticalJointConvergenceCriteria};
    bool collect_telemetry{};
    bool collect_diagnostics{};
    FixedNeighborBlockOrder order{FixedNeighborBlockOrder::Forward};
    std::size_t component_index{};
    std::string component_id;
    std::function<void(std::size_t &,std::size_t &)> halo_census;
    std::function<bool(int)> claim_parameter_stability_block;
    std::function<void(std::size_t,std::size_t,const BlockCoordinateState &,
        const FixedNeighborBlockSweep &)> state_observer;
    std::function<void(std::string_view,const FixedNeighborSearchResult &)> search_result_observer;
    std::function<void(FixedNeighborPolicy &)> configure_search_for_testing;
    std::optional<double> eta_upper_bound;
};

} // namespace rhbm_gem::core::joint_component
