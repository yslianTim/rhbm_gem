#pragma once
#include "Problem.hpp"
#include "CompactSvd.hpp"
namespace rhbm_gem::core::joint_component {
// Orthogonal reduction of the complete, column-scaled observation Jacobian.
// Selector rows annihilating its null space are equivalent to nuisance projection.
struct TargetGeometry
{
    CompactSvdResult svd;
    Vector scales;
    Evaluation endpoint;
    bool valid{},boundary{},threshold_sensitive{};
    std::string reason;
};
MemoryEstimate TargetIdentifiabilityMemoryEstimate(const JointProblemInput &,const JointParameterLayout &,std::size_t);
std::string TargetRankContextIdentity(const JointProblemInput &,const JointParameterLayout &,const JointState &,std::size_t);
JointTargetEvidence DeriveTargetIdentifiability(const JointProblem &,const JointParameterLayout &,
    const std::optional<JointState> &,const std::vector<JointCheck> &,const std::vector<JointRankEvidence> &,
    JointEvidenceScope,std::size_t);
TargetGeometry BuildTargetGeometry(const Domain &,VectorRef,const JointState &,const EvaluationContext &);
JointTargetEvidence AssessTargetState(const JointProblem &,const JointParameterLayout &,
    const std::optional<JointState> &,const std::vector<JointCheck> &,JointEvidenceScope,std::size_t);
JointCheckStatus TargetConvergence(const std::optional<JointTargetEvidence> &,bool);
void AddTargetEvidence(JointFitResult &);
}
