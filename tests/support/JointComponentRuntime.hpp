#pragma once
#include "support/JointABCComponents.hpp"
#include <filesystem>
namespace second_stage_test::matched::joint_abc {
void WriteEndpointAudit(const Domain &,const Eigen::VectorXd &,const boost::json::object &,
    const EvaluationContext &,const std::filesystem::path &);
void WriteLocalAudit(const ComponentView &,const Eigen::VectorXd &,const boost::json::object &,
    const EvaluationContext &,const std::filesystem::path &);
void RunJointRuntime(const std::string & model,const std::string & map,const std::string & output);
void RunPhysicalJointRuntime(const std::string & output,bool inputs_only=false);
}
