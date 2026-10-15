#pragma once
#include "Numerics.hpp"
#include "JointProgress.hpp"
#include "FixedNeighborBlockCoordinate.hpp"
#include "ResourceWork.hpp"
#include <rhbm_gem/core/JointComponentEstimator.hpp>
namespace rhbm_gem::core::joint_component {
JointInitialization InitializeContributors(MapObject &,ModelObject &,const JointProblem &);
struct ProblemData
{
    std::shared_ptr<const JointProblemInput> input;
    std::shared_ptr<const std::vector<double>> optimization_observations;
    std::vector<JointFixedHaloParameter> fixed_halos;
    std::vector<std::size_t> target_to_contributor,contributor_to_target;
    Domain domain{0,{}};
    VectorMap y;
    Eigen::VectorXd fixed_prediction;
    EvaluationContext context;
    ComponentPartition partition;
    JointParameterLayout layout;
    MemoryEstimate construction_memory;
    explicit ProblemData(std::shared_ptr<const JointProblemInput> snapshot,
        std::vector<JointFixedHaloParameter> fixed={},
        std::shared_ptr<const std::vector<double>> observations={})
        :input(std::move(snapshot)),optimization_observations(observations ? std::move(observations) :
            std::shared_ptr<const std::vector<double>>(input,&input->observations)),
        fixed_halos(std::move(fixed)),domain(input),
        y(optimization_observations->data(),static_cast<Eigen::Index>(optimization_observations->size())),
        fixed_prediction(Eigen::VectorXd::Zero(static_cast<Eigen::Index>(input->observations.size()))) {}
};
struct HaloSupportDiagnostic
{
    std::size_t atom_index{};
    std::string atom_id,selection_role{"halo"};
    std::size_t original_observation_support_count{},profiled_observation_support_count{},
        distinct_squared_distance_count{};
    bool singleton{},exact_common_response{},insufficient_full_abc_dimension{};
};
JointParameterLayout BuildParameterLayout(const JointProblemInput &);
std::vector<HaloSupportDiagnostic> DiagnoseHaloSupport(const JointProblemInput &,
    const JointParameterLayout &,const std::vector<Eigen::Index> * component_atoms=nullptr);
struct ComponentProblemView
{
    JointParameterLayout local_layout;
    Domain domain{0,{}};
    Eigen::Index parent_row_count{};
    Indices local_to_parent_atoms,local_to_parent_rows;
    Identities atom_ids,row_ids;
};
ComponentProblemView MakeComponentProblemView(const JointParameterLayout &,const ComponentView &);
FixedNeighborComponentResult SolveFixedNeighborComponentView(
    const ComponentView &,const JointParameterLayout &,const Vector &,const EvaluationContext &,
    const FixedNeighborSearchPolicy &,const JointProgressObserver & = {},
    const JointProgressComponent * = nullptr);
struct FixedNeighborComponentBatch
{
    std::vector<JointComponentResult> components;
    // The single-component FullABC path can reuse its certified assessment
    // during assembly. Other search and certification state is released as
    // each component finishes.
    std::optional<Assessment> single_component_assessment;
    JointFitCosts costs;
};
FixedNeighborComponentBatch OrchestrateFixedNeighborComponents(
    const JointProblem &,const std::vector<double> &,const FixedNeighborSearchPolicy &,
    const JointProgressObserver &);
JointFitResult FitFixedNeighborComponents(const JointProblem &,const std::vector<double> &,const FixedNeighborSearchPolicy &,
    const JointProgressObserver & = {});
JointFitResult FitFixedNeighborComponentsImpl(const JointProblem &,const std::vector<double> &,const FixedNeighborSearchPolicy &,
    const JointProgressObserver & = {});
}
namespace rhbm_gem::core {
struct JointProblemAccess
{
    static const joint_component::ProblemData & Get(const JointProblem & p) {return *p.m_data;}
};
}
