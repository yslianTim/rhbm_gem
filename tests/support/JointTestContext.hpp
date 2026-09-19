#pragma once
#include "core/detail/joint_component/Numerics.hpp"
#include <Eigen/Core>
#include <boost/json.hpp>
#include <memory>
#include <string>
#include <vector>

namespace second_stage_test::matched::joint_abc {
namespace runtime=rhbm_gem::core::joint_component;
using runtime::RankPolicy;
using runtime::LinearPolicy;
using runtime::AuditPlan;
using runtime::EvaluationContext;
// This is the fixture audit registry. Numerical kernels never inspect case names.
AuditPlan RegisteredAudit(Eigen::Index atoms, const std::string & dataset = "",
    const std::string & case_name = "");
EvaluationContext MakeContext(const Eigen::VectorXd &, Eigen::Index atoms,
    const std::string & snapshot_hash = "", const AuditPlan * = nullptr);
boost::json::object ContextEvidence(const EvaluationContext &);
} // namespace second_stage_test::matched::joint_abc
