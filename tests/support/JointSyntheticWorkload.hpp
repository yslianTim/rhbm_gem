#pragma once
#include "core/detail/joint_component/JointComponentEstimator.hpp"
#include <string>
namespace second_stage_test {
// Frozen support, no finite map box. Versioned benchmark data, not production.
rhbm_gem::core::JointProblemInput SyntheticJointWorkload(const std::string & topology,int atoms);
std::string SyntheticJointWorkloadHash(const rhbm_gem::core::JointProblemInput &);
}
