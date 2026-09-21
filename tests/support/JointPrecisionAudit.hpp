#pragma once
#include "support/JointTestNumerics.hpp"

namespace second_stage_test::matched::certification {
boost::json::object PrecisionNormalization(const Eigen::VectorXd & parent_observations);
void ResetPrecisionCache();
boost::json::object PrecisionCacheCosts();
boost::json::object PrecisionAudit(const joint_abc::Domain &, const Eigen::VectorXd &,
    const joint_abc::Evaluation &, const Eigen::MatrixXd & directions,
    const joint_abc::EvaluationContext * = nullptr, bool block_reference=false);
boost::json::object BoundaryAudit(const joint_abc::Domain &, const Eigen::VectorXd &,
    const Eigen::VectorXd & eta, const Eigen::VectorXd & beta,
    const joint_abc::EvaluationContext * = nullptr);
// Independently reprofile unchanged observations at two widths, at 50/100 digits.
boost::json::object PrecisionProfileChange(const joint_abc::Domain &, const Eigen::VectorXd &,
    const Eigen::VectorXd & base_eta, const Eigen::VectorXd & eta);
} // namespace second_stage_test::matched::certification
