#pragma once
#include "support/JointABCProfile.hpp"
#include <filesystem>

namespace second_stage_test::matched::certification {
joint_abc::Domain Snapshot(const unique_grid::Grid &, const std::vector<Atom> &,
    const std::filesystem::path & output);
void Audit(const joint_abc::Domain &, const Eigen::VectorXd &, const boost::json::object & fit,
    const std::filesystem::path & output, const joint_abc::EvaluationContext * = nullptr);
void AuditDirectory(const std::filesystem::path & directory);
} // namespace second_stage_test::matched::certification
