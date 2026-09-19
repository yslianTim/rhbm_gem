#pragma once
#include "support/JointABCCertification.hpp"

namespace second_stage_test::matched::certification {
struct LocalAuditState
{
    joint_abc::EvaluationContext context;
    boost::json::object fit, scope;
};
// Only the component's actual trusted state and immutable parent-derived context
// are needed. Inherited global directions are deliberately discarded.
LocalAuditState PrepareLocalAudit(const joint_abc::Domain &, const Eigen::VectorXd &,
    const boost::json::object &, const joint_abc::EvaluationContext &);
}
