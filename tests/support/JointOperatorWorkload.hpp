#pragma once
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <string>
namespace second_stage_test {
// Frozen support, no finite map box. Versioned benchmark data, not production.
rhbm_gem::core::JointProblemInput OperatorWorkload(const std::string & topology,int atoms);
std::string OperatorWorkloadHash(const rhbm_gem::core::JointProblemInput &);
}
