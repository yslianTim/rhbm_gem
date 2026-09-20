#pragma once
#include "Numerics.hpp"
#include <rhbm_gem/core/JointComponentEstimator.hpp>
namespace rhbm_gem::core::joint_component {
struct ProblemData
{
    std::shared_ptr<const JointProblemInput> input;
    Domain domain{0,{}};
    VectorMap y;
    EvaluationContext context;
    ComponentPartition partition;
    explicit ProblemData(std::shared_ptr<const JointProblemInput> snapshot):input(std::move(snapshot)),domain(input),
        y(input->observations.data(),static_cast<Eigen::Index>(input->observations.size())) {}
};
}
namespace rhbm_gem::core {
struct JointProblemAccess
{
    static const joint_component::ProblemData & Get(const JointProblem & p) {return *p.m_data;}
};
}
