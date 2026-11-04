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
struct BcdSearchPolicy;
struct BcdSearchResult;
#endif

inline constexpr std::size_t ProductionBcdCoreAtoms=12;
inline constexpr std::size_t ProductionBcdNominalSweeps=200;
inline constexpr std::size_t ProductionBcdMaximumSweeps=300;

// The production search accepts its calibrated width bound.
// Component identity is supplied internally for failure-only logs.
struct BcdExecutionPolicy
{
    std::size_t component_index{};
    std::string component_id;
    std::optional<double> eta_upper_bound;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    // Reuse production orchestration/assembly with a scoped test search policy.
    std::function<void(BcdSearchPolicy &)> configure_for_testing;
    std::function<void(const BcdSearchResult &)> observe_search_for_testing;
#endif
};

} // namespace rhbm_gem::core::joint_component
