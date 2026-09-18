#pragma once
#include "support/JointABCProfile.hpp"

namespace second_stage_test::matched::certification {
boost::json::object PrecisionAudit(const joint_abc::Domain &, const Eigen::VectorXd &,
    const joint_abc::Evaluation &, const Eigen::MatrixXd & directions);
boost::json::object BoundaryAudit(const joint_abc::Domain &, const Eigen::VectorXd &,
    const Eigen::VectorXd & eta, const Eigen::VectorXd & beta);
} // namespace second_stage_test::matched::certification
