#pragma once
#include "Numerics.hpp"
#include <rhbm_gem/core/JointComponentEstimator.hpp>
namespace rhbm_gem::core::joint_component {
JointInitialization InitializeContributors(MapObject &,ModelObject &,const JointProblem &);
struct ProblemData
{
    std::shared_ptr<const JointProblemInput> input;
    Domain domain{0,{}};
    VectorMap y;
    EvaluationContext context;
    ComponentPartition partition;
    JointParameterLayout layout;
    explicit ProblemData(std::shared_ptr<const JointProblemInput> snapshot):input(std::move(snapshot)),domain(input),
        y(input->observations.data(),static_cast<Eigen::Index>(input->observations.size())) {}
};
JointParameterLayout BuildParameterLayout(const JointProblemInput &);
JointFitResult FitObservableComponents(const JointProblem &,const std::vector<double> &,const SearchPolicy & = {});
JointFitResult FitWithSearchPolicy(const JointProblem &,const std::vector<double> &,const SearchPolicy &);
std::vector<JointRankEvidence> AssessmentRanks(const Assessment &,JointEvidenceScope);
Domain ProfileDomain(const Domain &,const JointParameterLayout &);
EvaluationContext ProfileContext(const EvaluationContext &,const JointParameterLayout &,Eigen::Index);
}
namespace rhbm_gem::core {
struct JointProblemAccess
{
    static const joint_component::ProblemData & Get(const JointProblem & p) {return *p.m_data;}
};
}
