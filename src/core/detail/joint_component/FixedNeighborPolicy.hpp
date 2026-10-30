#pragma once

#include <cstddef>
#include <optional>
#include <string>
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <functional>
#endif
#include "ConvergenceCriteria.hpp"

namespace rhbm_gem::core::joint_component {

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
struct FixedNeighborPolicy;
struct FixedNeighborSearchResult;
#endif

inline constexpr std::size_t ProductionFixedNeighborCoreAtoms=12;
inline constexpr std::size_t ProductionFixedNeighborNominalSweeps=200;
inline constexpr std::size_t ProductionFixedNeighborMaximumSweeps=300;

// The production search accepts its calibrated width bound.
// Component identity is supplied internally for failure-only logs.
struct FixedNeighborSearchPolicy
{
    std::size_t component_index{};
    std::string component_id;
    std::optional<double> eta_upper_bound;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    // Reuse production orchestration/assembly with a scoped test search policy.
    std::function<void(FixedNeighborPolicy &)> configure_for_testing;
    std::function<void(const FixedNeighborSearchResult &)> observe_search_for_testing;
#endif
};

} // namespace rhbm_gem::core::joint_component
