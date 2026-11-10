#pragma once

#include "core/estimator/joint_component/GaussSeidelBcd.hpp"

namespace joint_component_test {
rhbm_gem::core::joint_component::BcdSearchResult SearchBcdFixture(
    const rhbm_gem::core::JointProblem &,rhbm_gem::core::joint_component::VectorRef,
    const rhbm_gem::core::joint_component::BcdSearchPolicy &,
    std::vector<rhbm_gem::core::joint_component::BcdSweepRecord> * sweeps=nullptr);
}
