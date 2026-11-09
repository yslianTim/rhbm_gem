#pragma once

#include <cstddef>
#include <optional>
#include <string>
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <functional>
#endif

namespace rhbm_gem::core::joint_component {

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
struct BcdSearchPolicy;
#endif

inline constexpr std::size_t ProductionBcdCoreAtoms{ 12 };
inline constexpr std::size_t ProductionBcdNominalSweeps{ 200 };
inline constexpr std::size_t ProductionBcdMaximumSweeps{ 300 };

// Component execution supplies the calibrated width bound and failure-log identity.
// Block search settings remain in BcdSearchPolicy.
struct BcdExecutionPolicy
{
    std::size_t component_index{};
    std::string component_id;
    std::optional<double> eta_upper_bound;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    // Reuse production orchestration/assembly with a scoped test search policy.
    std::function<void(BcdSearchPolicy &)> configure_for_testing;
#endif
};

} // namespace rhbm_gem::core::joint_component
