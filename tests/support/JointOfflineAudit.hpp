#pragma once
#include "support/JointTestNumerics.hpp"
#include <filesystem>

namespace second_stage_test::matched::certification {
void Audit(const joint_abc::Domain &, const Eigen::VectorXd &, const boost::json::object & fit,
    const std::filesystem::path & output, const joint_abc::EvaluationContext * = nullptr);
} // namespace second_stage_test::matched::certification
