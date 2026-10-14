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
    Domain domain{0,{}};
    VectorMap y;
    EvaluationContext context;
    ComponentPartition partition;
    JointParameterLayout layout;
    MemoryEstimate construction_memory;
    explicit ProblemData(std::shared_ptr<const JointProblemInput> snapshot):input(std::move(snapshot)),domain(input),
        y(input->observations.data(),static_cast<Eigen::Index>(input->observations.size())) {}
};
struct HaloSupportDiagnostic
{
    std::size_t atom_index{};
    std::string atom_id,selection_role{"halo"};
    std::size_t original_observation_support_count{},profiled_observation_support_count{},
        distinct_squared_distance_count{},ac_structural_rank_upper_bound{};
    bool singleton{},exact_common_response{},insufficient_full_abc_dimension{};
    std::optional<bool> gaussian_column_nonzero_norm,charge_column_nonzero_norm;
    std::optional<std::size_t> fixed_width_ac_numerical_rank;
    std::optional<double> column_normalized_condition_indicator;
};
JointParameterLayout BuildParameterLayout(const JointProblemInput &);
std::vector<HaloSupportDiagnostic> DiagnoseHaloSupport(const JointProblemInput &,
    const JointParameterLayout &,const std::vector<double> & widths,
    const std::vector<Eigen::Index> * component_atoms=nullptr,Eigen::Index rank_context_rows=-1);
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
JointFitResult FitObservableComponents(const JointProblem &,const std::vector<double> &,const FixedNeighborSearchPolicy & = {},
    const JointProgressObserver & = {});
JointFitResult FitFixedNeighborComponents(const JointProblem &,const std::vector<double> &,const FixedNeighborSearchPolicy &,
    const JointProgressObserver & = {});
JointFitResult FitFixedNeighborComponentsImpl(const JointProblem &,const std::vector<double> &,const FixedNeighborSearchPolicy &,
    const JointProgressObserver & = {});
std::vector<JointRankEvidence> AssessmentRanks(const Assessment &,JointEvidenceScope,
    const JointProblemInput * = nullptr,const JointParameterLayout * = nullptr,std::size_t rank_context_rows=0);
std::vector<JointCheck> AnalyticNuisanceEvidence(JointEvidenceScope);
Domain ProfileDomain(const Domain &,const JointParameterLayout &);
EvaluationContext ProfileContext(const EvaluationContext &,const JointParameterLayout &,Eigen::Index);
}
namespace rhbm_gem::core {
struct JointProblemAccess
{
    static const joint_component::ProblemData & Get(const JointProblem & p) {return *p.m_data;}
};
}
