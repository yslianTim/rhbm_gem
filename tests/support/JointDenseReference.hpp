#pragma once
#include "core/detail/joint_component/Numerics.hpp"
namespace second_stage_test::matched::joint_abc {
namespace runtime=rhbm_gem::core::joint_component;
struct DenseDifferential {runtime::Matrix projected,jacobian; bool valid{}; std::string reason;};
DenseDifferential DenseDifferentiate(const runtime::Evaluation &,double,const runtime::EvaluationContext *,double=-1);
DenseDifferential MaterializeDerivative(const runtime::Evaluation &,double,const runtime::EvaluationContext *,double=-1,Eigen::Index=8192);
runtime::Vector DenseLocalCorrection(const runtime::Evaluation &,const DenseDifferential &,const runtime::EvaluationContext &,double=-1);
}
