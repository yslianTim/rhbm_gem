#pragma once
#include "core/detail/joint_component/JointComponentEstimator.hpp"
#include <string>
namespace joint_component_test {
// Deterministic lattice support and observations for numerical regression tests.
rhbm_gem::core::JointProblemInput SyntheticJointWorkload(const std::string & topology,int atoms);
}
