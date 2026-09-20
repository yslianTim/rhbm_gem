#pragma once
#include "support/JointTestNumerics.hpp"
namespace second_stage_test::matched::certification {
boost::json::object AssessWithDerivativeAudit(const joint_abc::Domain &,const Eigen::VectorXd &,
    const Eigen::VectorXd &,const joint_abc::EvaluationContext &,const Eigen::VectorXd * beta=nullptr);
}
