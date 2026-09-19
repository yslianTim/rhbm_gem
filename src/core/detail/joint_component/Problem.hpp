#pragma once
#include "Numerics.hpp"
#include <rhbm_gem/core/JointComponentEstimator.hpp>
namespace rhbm_gem::core::joint_component {
struct ProblemData
{
    JointProblemInput input;
    Domain domain{0,{}};
    Vector y;
    EvaluationContext context;
    ComponentPartition partition;
};
}
namespace rhbm_gem::core {
struct JointProblemAccess
{
    static const joint_component::ProblemData & Get(const JointProblem & p) {return *p.m_data;}
};
}
