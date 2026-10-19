#pragma once

#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"

#include <functional>
#include <string>

namespace joint_test::width_trust {

enum class Policy {Baseline,ConditioningAware,ConditioningAndCancellationAware};

const char * Name(Policy policy);

std::function<std::string(std::size_t,std::size_t,
    const rhbm_gem::core::joint_component::FixedNeighborProfileTrial &,
    const rhbm_gem::core::joint_component::FixedNeighborProfileTrial &)> Make(Policy policy);

} // namespace joint_test::width_trust
