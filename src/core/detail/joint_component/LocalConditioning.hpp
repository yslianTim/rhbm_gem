#pragma once
#include "Numerics.hpp"
#include "Problem.hpp"

namespace rhbm_gem::core::joint_component {
struct LocalConditioningAssessment
{
    JointBlockConditioning evidence;
    MemoryEstimate memory_estimate;
    double seconds{};
};
LocalConditioningAssessment AssessFixedNeighborLocalConditioning(
    const Domain &,VectorRef,const Vector &,const Vector &,const EvaluationContext &,
    std::string,std::size_t=0);
}
