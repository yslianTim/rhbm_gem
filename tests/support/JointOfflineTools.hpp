#pragma once
#include "support/JointComponentChecks.hpp"
#include <filesystem>
namespace second_stage_test::matched::joint_abc {
void ComponentLocalBundleRerun(const std::string &,const std::string &);
void WriteLocalAudit(const ComponentView &,const Eigen::VectorXd &,const boost::json::object &,
    const EvaluationContext &,const std::filesystem::path &);
}
