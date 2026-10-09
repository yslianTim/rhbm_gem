#pragma once
#include "Numerics.hpp"
#include "JointProgress.hpp"
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
struct PreparedComponent
{
    std::shared_ptr<const JointProblemInput> input;
    JointParameterLayout parent_layout,local_layout;
    ComponentView view;
    Domain domain{0,{}};
    Indices local_to_parent_atoms,local_to_parent_rows,parent_to_local_atoms,parent_to_local_rows;
    Identities atom_ids,row_ids;
};
PreparedComponent PrepareComponent(const JointParameterLayout &,const ComponentView &);
ComponentResult SolveFixedNeighborComponentView(
    const ComponentView &,const JointParameterLayout &,const Vector &,const EvaluationContext &,
    const FixedNeighborSearchPolicy &,const JointProgressObserver & = {},
    const JointProgressComponent * = nullptr);
JointFitResult FitObservableComponents(const JointProblem &,const std::vector<double> &,const FixedNeighborSearchPolicy & = {},
    const JointProgressObserver & = {});
JointFitResult FitFixedNeighborComponents(const JointProblem &,const std::vector<double> &,const FixedNeighborSearchPolicy &,
    const JointProgressObserver & = {});
JointFitResult FitFixedNeighborComponentsImpl(const JointProblem &,const std::vector<double> &,const FixedNeighborSearchPolicy &,
    const JointProgressObserver & = {});
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
