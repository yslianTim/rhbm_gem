#pragma once
#include "support/JointComponentChecks.hpp"
#include <filesystem>
namespace second_stage_test::matched::joint_abc {
void RunJointRuntime(const std::string & model,const std::string & map,const std::string & output);
void RunPhysicalJointRuntime(const std::string & output,bool inputs_only=false);
}
