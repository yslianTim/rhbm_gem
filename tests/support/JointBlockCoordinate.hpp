#pragma once
#include "core/detail/joint_component/StructuralPartition.hpp"
#include "core/detail/joint_component/Numerics.hpp"

namespace second_stage_test::joint_block {
namespace n=rhbm_gem::core::joint_component;
struct BlockCoordinateState
{
    n::Vector observations,eta,beta,prediction,residual;
    double scale{},objective{};
};
struct ConditionalBlockProblem
{
    n::Indices rows;
    n::Vector effective_response,old_core_prediction;
    double scale{};
};

double ParentScale(n::VectorRef observations);
BlockCoordinateState Initialize(const rhbm_gem::core::JointProblemInput &,const rhbm_gem::JointParameterLayout &,
    n::VectorRef observations,n::Vector eta,n::Vector beta,double parent_scale);
ConditionalBlockProblem BuildConditionalProblem(const rhbm_gem::core::JointProblemInput &,
    const n::StructuralCore &,const BlockCoordinateState &);
double ConditionalObjective(const rhbm_gem::core::JointProblemInput &,const n::StructuralCore &,
    const ConditionalBlockProblem &,n::VectorRef eta,n::VectorRef beta);
void ReplaceBlock(const rhbm_gem::core::JointProblemInput &,const n::StructuralCore &,
    n::VectorRef eta,n::VectorRef beta,BlockCoordinateState &);
BlockCoordinateState Replay(const rhbm_gem::core::JointProblemInput &,const rhbm_gem::JointParameterLayout &,
    n::VectorRef observations,n::VectorRef eta,n::VectorRef beta,double parent_scale);
}
