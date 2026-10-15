#include <gtest/gtest.h>
#include <boost/json.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/data/io/JointAnalysisFileIO.hpp>
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/TargetEvidence.hpp"
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/FittingWorkset.hpp"
#include "core/detail/gaussian_fit/PreparedLocalGaussianFit.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include "data/io/detail/JointResultJson.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include "core/detail/PostFitPeeling.hpp"
#include "core/detail/MapSampler.hpp"
#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <fstream>

namespace {
namespace core=rhbm_gem::core;
namespace n=core::joint_component;
constexpr double missing_width=std::numeric_limits<double>::quiet_NaN();
n::Domain FullDomain(const n::Domain & domain,const rhbm_gem::JointParameterLayout & layout)
{
    EXPECT_EQ(layout.full_atoms.size(),domain.atoms.size());
    EXPECT_EQ(layout.informative_rows.size(),static_cast<std::size_t>(domain.rows));
    return domain;
}
n::EvaluationContext FullContext(const n::EvaluationContext & context,
    const rhbm_gem::JointParameterLayout &,Eigen::Index)
{return context;}
core::JointProblemInput Input(bool noise=false)
{
    core::JointProblemInput input;
    input.atom_ids={"1","2","3"}; input.support.resize(3);
    input.selection_domain.emplace(); input.selection_domain->target_indices={0};
    for(std::size_t row=0;row<24;++row)
    {
        const double square=std::pow(.05+.1*static_cast<double>(row),2);
        const auto basis=n::EvaluateKernel(square,.6,2.5);
        input.row_ids.push_back(std::to_string(row));
        input.observations.push_back(2*basis.gaussian+.1*basis.charge+(noise ? 1e-5*std::sin(static_cast<double>(row)) : 0));
        input.support[0].push_back({row,square});
    }
    input.support[1]={{23,6.25}}; input.support[2]={{23,6.24}};
    input.observations[23]+=.7;
    return input;
}
struct OriginalTargetGeometryOracle
{
    n::Vector scales,singular_values;
    n::Matrix right_vectors;
    double threshold{};
    Eigen::Index rank{};
    bool valid{},threshold_sensitive{};
    std::string reason;
};
OriginalTargetGeometryOracle OriginalGeometryOracle(const n::Domain & domain,n::VectorRef y,
    const rhbm_gem::JointState & state,const n::EvaluationContext & context)
{
    OriginalTargetGeometryOracle out;
    const auto columns=static_cast<Eigen::Index>(state.b.size());
    if(!columns || state.ac.size()!=static_cast<std::size_t>(2*columns) ||
        state.log_b.size()!=static_cast<std::size_t>(columns))
    {out.reason="missing-full-parameter-state"; return out;}
    const n::VectorMap eta(state.log_b.data(),columns),beta(state.ac.data(),2*columns);
    for(Eigen::Index atom=0;atom<columns;++atom)
        if(beta(2*atom)==0) {out.reason="active-amplitude-boundary"; return out;}
    const auto endpoint=n::EvaluateState(domain,y,eta,beta,context);
    if(!endpoint.valid || !endpoint.certificate.kkt_passed) {out.reason="invalid-saved-state"; return out;}
    n::Matrix jacobian=n::Matrix::Zero(y.size(),3*columns);
    for(Eigen::Index atom=0;atom<columns;++atom)
        for(int k=0;k<2;++k)
        {
            for(n::Sparse::InnerIterator value(endpoint.x,2*atom+k);value;++value)
                jacobian(value.row(),3*atom+k)=value.value();
            for(n::Sparse::InnerIterator value(endpoint.derivative,2*atom+k);value;++value)
                jacobian(value.row(),3*atom+2)+=value.value()*beta(2*atom+k);
        }
    out.scales=jacobian.colwise().norm();
    if(!out.scales.allFinite()) {out.reason="nonfinite-jacobian"; return out;}
    for(Eigen::Index col=0;col<out.scales.size();++col)
    {
        if(out.scales(col)==0) out.scales(col)=1;
        jacobian.col(col)/=out.scales(col);
    }
    Eigen::JacobiSVD<n::Matrix> svd(jacobian,Eigen::ComputeFullV);
    if(svd.info()!=Eigen::Success || !svd.singularValues().allFinite())
    {out.reason="jacobian-factorization-failed"; return out;}
    svd.setThreshold(context.rank.Relative(jacobian.cols()));
    out.singular_values=svd.singularValues(); out.right_vectors=svd.matrixV();
    out.threshold=svd.threshold()*out.singular_values(0); out.rank=svd.rank();
    out.threshold_sensitive=out.threshold>0 &&
        ((out.singular_values.array()>=.5*out.threshold) && (out.singular_values.array()<=2*out.threshold)).any();
    out.valid=true; return out;
}
}
TEST(JointObservableProfileTest, FixedLayoutAndSharedRows)
{
    const core::JointProblem problem(Input()); const auto & layout=problem.ParameterLayout();
    EXPECT_EQ(layout.full_atoms,(std::vector<std::size_t>{0,1,2}));
    EXPECT_TRUE(layout.groups.empty()); EXPECT_EQ(layout.informative_rows.size(),24);
    auto selected=Input(); selected.selection_domain->target_indices={0,1};
    EXPECT_EQ(core::JointProblem(selected).ParameterLayout().full_atoms,(std::vector<std::size_t>{0,1,2}));
    auto legacy=Input(); legacy.selection_domain.reset();
    EXPECT_TRUE(core::JointProblem(legacy).ParameterLayout().groups.empty());
    EXPECT_EQ(core::JointProblem(legacy).ParameterLayout().full_atoms,(std::vector<std::size_t>{0,1,2}));
    auto multiple=Input(); multiple.support[2].push_back({22,6.0});
    EXPECT_EQ(core::JointProblem(multiple).ParameterLayout().full_atoms,(std::vector<std::size_t>{0,1,2}));
    auto separate=Input(); separate.support[2][0].row=22;
    EXPECT_TRUE(core::JointProblem(separate).ParameterLayout().groups.empty());
}
TEST(JointObservableProfileTest, HaloSupportCensusSeparatesSupportDistanceAndRankEvidence)
{
    core::JointProblemInput input;
    input.atom_ids={"target","singleton","equal-two","distinct-two","equal-three"};
    input.row_ids={"r0","r1","r2","r3","r4"}; input.observations.assign(5,0.0);
    input.support.resize(5);
    input.support[0]={{0,.1},{1,.2},{2,.3},{3,.4},{4,.5}};
    input.support[1]={{4,.6}};
    input.support[2]={{1,.7},{2,.7}};
    input.support[3]={{0,.2},{1,1.2}};
    input.support[4]={{2,.9},{3,.9},{4,.9}};
    input.selection_domain.emplace(); input.selection_domain->target_indices={0};

    const auto layout=n::BuildParameterLayout(input);
    const auto census=n::DiagnoseHaloSupport(input,layout,{.6,.6,.6,.6,.6});
    ASSERT_EQ(census.size(),4u);
    EXPECT_EQ(census[0].selection_role,"halo");
    EXPECT_EQ(census[0].atom_id,"singleton");
    EXPECT_TRUE(census[0].singleton);
    EXPECT_EQ(census[0].original_observation_support_count,1u);
    EXPECT_EQ(census[0].profiled_observation_support_count,1u);
    EXPECT_EQ(census[0].distinct_squared_distance_count,1u);
    EXPECT_TRUE(census[0].exact_common_response);
    EXPECT_TRUE(census[0].insufficient_full_abc_dimension);
    ASSERT_TRUE(census[0].fixed_width_ac_numerical_rank);
    EXPECT_EQ(*census[0].fixed_width_ac_numerical_rank,1u);

    EXPECT_EQ(census[1].original_observation_support_count,2u);
    EXPECT_EQ(census[1].distinct_squared_distance_count,1u);
    EXPECT_TRUE(census[1].exact_common_response);
    EXPECT_TRUE(census[1].insufficient_full_abc_dimension);
    ASSERT_TRUE(census[1].fixed_width_ac_numerical_rank);
    EXPECT_EQ(*census[1].fixed_width_ac_numerical_rank,1u);

    EXPECT_EQ(census[2].original_observation_support_count,2u);
    EXPECT_EQ(census[2].distinct_squared_distance_count,2u);
    EXPECT_FALSE(census[2].exact_common_response);
    ASSERT_TRUE(census[2].fixed_width_ac_numerical_rank);
    EXPECT_EQ(*census[2].fixed_width_ac_numerical_rank,2u);
    EXPECT_TRUE(census[2].gaussian_column_nonzero_norm.value_or(false));
    EXPECT_TRUE(census[2].charge_column_nonzero_norm.value_or(false));

    EXPECT_EQ(census[3].original_observation_support_count,3u);
    EXPECT_EQ(census[3].distinct_squared_distance_count,1u);
    EXPECT_FALSE(census[3].insufficient_full_abc_dimension);
    ASSERT_TRUE(census[3].fixed_width_ac_numerical_rank);
    EXPECT_EQ(*census[3].fixed_width_ac_numerical_rank,1u);
}
TEST(JointObservableProfileTest, IndividuallyFullRankHalosCanHaveDeficientJointDesign)
{
    core::JointProblemInput input;
    input.atom_ids={"target","halo-a","halo-b"}; input.observations.assign(3,0.0);
    input.row_ids={"r0","r1","r2"}; input.support.resize(3);
    input.support[0]={{0,.1}};
    input.support[1]={{0,.2},{1,.7},{2,1.3}};
    input.support[2]=input.support[1];
    input.selection_domain.emplace(); input.selection_domain->target_indices={0};
    const auto layout=n::BuildParameterLayout(input);
    const auto census=n::DiagnoseHaloSupport(input,layout,{.6,.6,.6});
    ASSERT_EQ(census.size(),2u);
    ASSERT_TRUE(census[0].fixed_width_ac_numerical_rank);
    ASSERT_TRUE(census[1].fixed_width_ac_numerical_rank);
    EXPECT_EQ(*census[0].fixed_width_ac_numerical_rank,2u);
    EXPECT_EQ(*census[1].fixed_width_ac_numerical_rank,2u);

    n::Matrix joint_design(3,4);
    for(Eigen::Index row=0;row<3;++row)
    {
        const auto basis=n::EvaluateKernel(input.support[1][static_cast<std::size_t>(row)].squared_distance,.6,2.5);
        joint_design.row(row)<<basis.gaussian,basis.charge,basis.gaussian,basis.charge;
    }
    const auto spectrum=n::ComputeSpectrum(joint_design,{3,4,0},4,true);
    ASSERT_TRUE(spectrum.available);
    EXPECT_EQ(spectrum.rank,2);
}
TEST(JointObservableProfileTest, FullAbcLayoutKeepsAllContributorRows)
{
    const core::JointProblem problem(Input(true)); const auto & data=core::JointProblemAccess::Get(problem);
    const auto & layout=problem.ParameterLayout();
    EXPECT_EQ(layout.full_atoms,(std::vector<std::size_t>{0,1,2}));
    EXPECT_TRUE(layout.groups.empty());
    EXPECT_EQ(layout.informative_rows.size(),problem.Input().observations.size());
    EXPECT_EQ(data.domain.atoms.size(),problem.Input().atom_ids.size());
    EXPECT_EQ(data.domain.rows,static_cast<Eigen::Index>(problem.Input().observations.size()));
}
TEST(JointObservableProfileTest, ReconstructionAndUnavailableHaloParameters)
{
    auto input=Input(); input.atom_ids.resize(1); input.support.resize(1); input.observations[23]-=.7;
    const core::JointProblem problem(std::move(input)); const auto fit=core::FitJointComponents(problem,{.65});
    ASSERT_TRUE(fit.assembled_state); ASSERT_TRUE(fit.prediction);
    EXPECT_EQ(fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(fit.components[0].atoms.size(),1); EXPECT_EQ(fit.components[0].state->b.size(),1);
    EXPECT_FALSE(fit.components[0].layout);
    ASSERT_TRUE(fit.components[0].candidate_diagnostic);
    EXPECT_TRUE(fit.components[0].candidate_diagnostic->candidate_available);
    EXPECT_EQ(fit.components[0].candidate_diagnostic->candidate_eta,fit.components[0].state->log_b);
    EXPECT_EQ(fit.components[0].candidate_diagnostic->candidate_ac,fit.components[0].state->ac);
    auto invalid_state=*fit.assembled_state; invalid_state.ac[0]+=.01;
    const auto invalid=n::AssessTargetState(problem,problem.ParameterLayout(),invalid_state,fit.evidence,fit.ranks,
        rhbm_gem::JointEvidenceScope::AssembledGlobal,problem.Input().observations.size());
    EXPECT_EQ(invalid.atoms[0].reason,"missing-rank-evidence");
    auto inconsistent_state=*fit.assembled_state; inconsistent_state.b[0]*=1.01;
    const auto inconsistent=n::AssessTargetState(problem,problem.ParameterLayout(),inconsistent_state,fit.evidence,fit.ranks,
        rhbm_gem::JointEvidenceScope::AssembledGlobal,problem.Input().observations.size());
    EXPECT_EQ(inconsistent.atoms[0].reason,"inconsistent-saved-widths");
    EXPECT_EQ(fit.components[0].candidate_diagnostic->certification_status,rhbm_gem::JointCheckStatus::Passed);
    const auto & problem_data=core::JointProblemAccess::Get(problem);
    const auto candidate_domain=problem_data.domain;
    const auto candidate_context=problem_data.context;
    n::Vector candidate_y=problem_data.y;
    n::Vector candidate_eta(static_cast<Eigen::Index>(fit.components[0].candidate_diagnostic->candidate_eta.size()));
    n::Vector candidate_ac(static_cast<Eigen::Index>(fit.components[0].candidate_diagnostic->candidate_ac.size()));
    for(std::size_t k=0;k<fit.components[0].candidate_diagnostic->candidate_eta.size();++k)
        candidate_eta(static_cast<Eigen::Index>(k))=fit.components[0].candidate_diagnostic->candidate_eta[k];
    for(std::size_t k=0;k<fit.components[0].candidate_diagnostic->candidate_ac.size();++k)
        candidate_ac(static_cast<Eigen::Index>(k))=fit.components[0].candidate_diagnostic->candidate_ac[k];
    const auto reconstructed=n::EvaluateStateStreaming(candidate_domain,candidate_y,candidate_eta,candidate_ac,candidate_context);
    ASSERT_TRUE(reconstructed.valid);
    ASSERT_TRUE(fit.components[0].candidate_diagnostic->candidate_objective);
    EXPECT_NEAR(*fit.components[0].candidate_diagnostic->candidate_objective,
        reconstructed.certificate.objective/(fit.observation_scale*fit.observation_scale),1e-14);
    EXPECT_EQ(boost::json::parse(rhbm_gem::joint_result_io::Encode(core::CaptureJointAnalysisResult(fit)))
        .at("components").as_array()[0].at("candidate_diagnostic").at("candidate_status").as_string(),"certified");
    EXPECT_NEAR(fit.assembled_state->b[0],.6,1e-10);
    EXPECT_TRUE(fit.assembled_state->nuisance_amplitudes.empty());
    for(std::size_t r=0;r<24;++r) EXPECT_NEAR(fit.prediction->at(r),problem.Input().observations[r],1e-10);
    const auto saved=core::CaptureJointAnalysisResult(fit);
    const auto stages=rhbm_gem::data_internal::BuildJointStageEstimates(saved,"test");
    EXPECT_TRUE(stages.at(1).point);
    const auto decoded=rhbm_gem::joint_result_io::Decode(rhbm_gem::joint_result_io::Encode(saved));
    EXPECT_FALSE(decoded.layout);
    EXPECT_EQ(decoded.parameterization_contract,"full-abc-v1");
    EXPECT_TRUE(decoded.assembled_state->nuisance_amplitudes.empty());
    ASSERT_TRUE(decoded.components[0].candidate_diagnostic);
    EXPECT_EQ(decoded.components[0].candidate_diagnostic->candidate_eta,
        saved.components[0].candidate_diagnostic->candidate_eta);
    EXPECT_EQ(decoded.components[0].candidate_diagnostic->certification_status,
        saved.components[0].candidate_diagnostic->certification_status);
    auto malformed=boost::json::parse(rhbm_gem::joint_result_io::Encode(saved));
    malformed.at("target_evidence").as_object()["column_scales"]=boost::json::array{1.};
    EXPECT_THROW(rhbm_gem::joint_result_io::Decode(boost::json::serialize(malformed)),std::invalid_argument);
}
TEST(JointObservableProfileTest, SharedOrchestrationReturnsCertifiedStatesBeforeAssembly)
{
    auto input=Input(); input.atom_ids.resize(1); input.support.resize(1); input.observations[23]-=.7;
    const core::JointProblem problem(std::move(input));
    const std::vector<double> initial{.65};
    std::vector<n::JointProgressEvent> events;
    const auto batch=n::OrchestrateFixedNeighborComponents(problem,initial,{},
        [&](const auto & event) {events.push_back(event);});
    ASSERT_EQ(batch.components.size(),1);
    ASSERT_TRUE(batch.components.front().state);
    EXPECT_FALSE(batch.components.front().layout);
    EXPECT_TRUE(batch.single_component_assessment);
    EXPECT_EQ(events.back().phase,n::JointProgressPhase::ComponentCompleted);
    EXPECT_TRUE(std::none_of(events.begin(),events.end(),[](const auto & event) {
        return event.phase==n::JointProgressPhase::AssemblyStarted ||
            event.phase==n::JointProgressPhase::AssemblyCompleted;
    }));

    const auto fit=n::FitFixedNeighborComponents(problem,initial,{});
    ASSERT_TRUE(fit.assembled_state); ASSERT_TRUE(fit.prediction); ASSERT_TRUE(fit.objective);
    ASSERT_TRUE(fit.components.front().state);
    EXPECT_EQ(fit.components.front().state->ac,batch.components.front().state->ac);
    EXPECT_EQ(fit.components.front().state->log_b,batch.components.front().state->log_b);
    EXPECT_TRUE(fit.assembled_state->nuisance_amplitudes.empty());
    n::Vector residual(static_cast<Eigen::Index>(problem.Input().observations.size()));
    for(std::size_t row=0;row<problem.Input().observations.size();++row)
    {
        EXPECT_NEAR(fit.prediction->at(row),problem.Input().observations[row],1e-10);
        residual(static_cast<Eigen::Index>(row))=fit.prediction->at(row)-problem.Input().observations[row];
    }
    EXPECT_NEAR(*fit.objective,.5*residual.squaredNorm()/(fit.observation_scale*fit.observation_scale),1e-15);
}
TEST(JointObservableProfileTest, SweepBudgetRetainsCandidateWithoutCompletingSearch)
{
    auto input=Input(); input.atom_ids.resize(1); input.support.resize(1); input.observations[23]-=.7;
    const core::JointProblem problem(std::move(input));
    n::FixedNeighborSearchPolicy policy; policy.maximum_sweeps=1;
    const auto fit=n::FitFixedNeighborComponents(problem,{.65},policy);
    ASSERT_EQ(fit.components.size(),1u);
    const auto & component=fit.components.front();
    ASSERT_TRUE(component.candidate_diagnostic);
    EXPECT_TRUE(component.candidate_diagnostic->candidate_available);
    EXPECT_FALSE(component.search_completed);
    EXPECT_EQ(component.stop_reason,"block-sweep-budget");
    EXPECT_EQ(component.candidate_diagnostic->candidate_eta.size(),problem.ParameterLayout().full_atoms.size());
    EXPECT_NE(component.candidate_diagnostic->certification_status,rhbm_gem::JointCheckStatus::NotRun);
}
TEST(JointObservableProfileTest, RemainingUnobservedTargetStillFails)
{
    auto input=Input(); input.support[0]={{23,1}};
    const core::JointProblem problem(input);
    const auto fit=core::FitJointComponents(problem,{.6,.7,.8});
    EXPECT_FALSE(fit.assembled_state); EXPECT_NE(fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
}
TEST(JointObservableProfileTest, RemainingDependentTargetIsNotReducedAgain)
{
    auto input=Input(); input.support[0]={{0,1},{1,1},{23,1}};
    const core::JointProblem problem(input);
    EXPECT_EQ(problem.ParameterLayout().full_atoms,(std::vector<std::size_t>{0,1,2}));
    const auto fit=core::FitJointComponents(problem,{.6,.7,.8});
    ASSERT_EQ(fit.components.size(),1); EXPECT_FALSE(fit.components[0].state);
    EXPECT_FALSE(fit.assembled_state);
    EXPECT_NE(fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
}
TEST(JointObservableProfileTest, RemainingTwoRowInteriorWidthIsUnidentified)
{
    auto input=Input(); input.atom_ids.resize(1); input.support.resize(1);
    input.atom_ids.push_back("4"); input.support.push_back({{0,.4},{1,1.2}});
    input.selection_domain->target_indices.push_back(1);
    for(const auto & support:input.support[1])
    {
        const auto basis=n::EvaluateKernel(support.squared_distance,.65,2.5);
        input.observations[support.row]+=basis.gaussian+.2*basis.charge;
    }
    const core::JointProblem problem(input); const auto & data=core::JointProblemAccess::Get(problem);
    const auto & layout=problem.ParameterLayout();
    EXPECT_EQ(layout.full_atoms,(std::vector<std::size_t>{0,1}));
    const auto context=FullContext(data.context,layout,data.y.size());
    const n::Vector eta=Eigen::Vector2d(.6,.65).array().log();
    const auto assessment=n::AssessProfile(FullDomain(data.domain,layout),data.y,eta,context);
    ASSERT_TRUE(assessment.primary.valid); ASSERT_TRUE(assessment.independent_replay.passed);
    EXPECT_GT(assessment.primary.beta(0),0); EXPECT_GT(assessment.primary.beta(2),0);
    EXPECT_FALSE(assessment.identified); EXPECT_EQ(assessment.failure,"width-unidentified");
}
TEST(JointObservableProfileTest, FrozenMedianDonorsAndOriginalFailures)
{
    rhbm_gem::JointInitialization initialization; initialization.b={1,3,missing_width,-1,missing_width}; initialization.atoms.resize(5);
    for(std::size_t k=0;k<2;++k)
    {auto & atom=initialization.atoms[k]; atom.original_b=initialization.b[k]; atom.reason="valid-width"; atom.seed_source="fitted";}
    initialization.atoms[2].reason="initialization-exception: test";
    initialization.atoms[3].reason="invalid-width"; initialization.atoms[3].original_b=-1;
    initialization.atoms[4].reason="invalid-width";
    core::detail::ApplyJointSeedFallback(initialization);
    EXPECT_TRUE(initialization.valid); EXPECT_EQ(initialization.b[2],2); EXPECT_EQ(initialization.b[3],2);
    EXPECT_EQ(initialization.atoms[2].reason,"initialization-exception: test");
    EXPECT_FALSE(initialization.atoms[2].original_b); EXPECT_EQ(initialization.atoms[3].original_b,-1);
    EXPECT_EQ(initialization.atoms[2].seed_source,"median-fallback"); EXPECT_EQ(initialization.atoms[2].donor_count,2);
    EXPECT_EQ(initialization.b[4],2);
    EXPECT_EQ(initialization.atoms[4].seed_source,"median-fallback");
    EXPECT_EQ(initialization.atoms[4].donor_count,2u);
    rhbm_gem::JointInitialization empty; empty.b={missing_width}; empty.atoms.resize(1); empty.atoms[0].reason="invalid-width";
    core::detail::ApplyJointSeedFallback(empty); EXPECT_FALSE(empty.valid); EXPECT_TRUE(std::isnan(empty.b[0]));
}
TEST(JointObservableProfileTest, LegacyV3DoesNotInferReductionAndMalformedLayoutRejected)
{
    auto input=Input(); input.atom_ids.resize(1); input.support.resize(1); input.observations[23]-=.7;
    const auto saved=core::CaptureJointAnalysisResult(core::FitJointComponents(core::JointProblem(input),{.6}));
    auto json=boost::json::parse(rhbm_gem::joint_result_io::Encode(saved)).as_object();
    json["schema_version"]=3; json.erase("layout"); json.erase("parameterization_contract");
    auto & solver=json.at("metadata").as_object().at("solver").as_object();
    solver["contract_version"]="joint-solver-provenance-v3";
    solver["fixed_neighbor_policy_version"]="fixed-neighbor-production-v2";
    const auto old=rhbm_gem::joint_result_io::Decode(boost::json::serialize(json));
    EXPECT_FALSE(old.layout); EXPECT_EQ(old.parameterization_contract,"full-abc-v1");
    auto legacy_input=Input();
    const core::JointProblem legacy_problem(legacy_input);
    const auto current=core::CaptureJointAnalysisResult(core::FitJointComponents(legacy_problem,{.6,.7,.8}));
    auto legacy_json=boost::json::parse(rhbm_gem::joint_result_io::Encode(current)).as_object();
    boost::json::array full_atoms{0};
    boost::json::array informative_rows; for(std::size_t row=0;row<23;++row) informative_rows.push_back(row);
    boost::json::array groups{boost::json::object{{"row",23},{"atoms",boost::json::array{1,2}}}};
    const boost::json::object legacy_layout{{"full_atoms",full_atoms},{"informative_rows",informative_rows},{"groups",groups}};
    legacy_json["schema_version"]=10;
    legacy_json["parameterization_contract"]="singleton-halo-profile-v1";
    legacy_json["layout"]=legacy_layout;
    legacy_json["assembled_state"]=nullptr;
    legacy_json["target_evidence"]=nullptr;
    legacy_json["ranks"]=boost::json::array{};
    auto & component=legacy_json.at("components").as_array()[0].as_object();
    component["layout"]=legacy_layout;
    component["state"]=nullptr;
    component["candidate_diagnostic"]=nullptr;
    component["target_evidence"]=nullptr;
    component["ranks"]=boost::json::array{};
    auto & metadata=legacy_json.at("metadata").as_object();
    auto & legacy_solver=metadata.at("solver").as_object();
    legacy_solver["contract_version"]="joint-solver-provenance-v5";
    legacy_solver["fixed_neighbor_policy_version"]="fixed-neighbor-production-v4";
    const auto legacy=rhbm_gem::joint_result_io::Decode(boost::json::serialize(legacy_json));
    ASSERT_TRUE(legacy.layout);
    EXPECT_EQ(legacy.parameterization_contract,"singleton-halo-profile-v1");
    EXPECT_EQ(legacy.layout->groups.size(),1u);
    EXPECT_EQ(legacy.layout->groups.front().atoms,(std::vector<std::size_t>{1,2}));
    EXPECT_EQ(rhbm_gem::joint_result_io::Decode(rhbm_gem::joint_result_io::Encode(legacy)).layout->groups.front().atoms,
        (std::vector<std::size_t>{1,2}));
    auto bad=legacy_json;
    bad.at("layout").as_object().at("groups").as_array()[0].as_object()["atoms"]=boost::json::array{0,1,2};
    EXPECT_THROW(rhbm_gem::joint_result_io::Decode(boost::json::serialize(bad)),std::invalid_argument);
}
TEST(JointObservableProfileTest, IndependentAnalyticComponentSurvivesInvalidTargetSeed)
{
    auto input=Input(); input.row_ids.push_back("24"); input.observations.push_back(5); input.support[2]={{24,2}};
    const auto fit=core::FitJointComponents(core::JointProblem(input),{.6,.6,.6});
    ASSERT_EQ(fit.components.size(),2); EXPECT_FALSE(fit.assembled_state);
    EXPECT_FALSE(fit.components[0].state); EXPECT_FALSE(fit.components[1].state);
    EXPECT_NE(fit.components[1].stop_reason,"analytic-nuisance-only");
    EXPECT_TRUE(fit.components[0].layout == std::nullopt);
    EXPECT_TRUE(fit.components[1].layout == std::nullopt);
}
TEST(JointObservableProfileTest, FullAbcPersistenceExportPeelingAndInitialization)
{
    using namespace rhbm_gem;
    auto input=Input(); input.atom_ids.resize(1); input.support.resize(1); input.observations[23]-=.7;
    const core::JointProblem problem(std::move(input));
    const auto result=core::CaptureJointAnalysisResult(core::FitJointComponents(problem,{.6}));
    std::vector<std::unique_ptr<AtomObject>> atoms;
    for(int id=1;id<=3;++id)
    {auto atom=std::make_unique<AtomObject>(); atom->SetSerialID(id); atom->SetElement(Element::CARBON); atom->SetPosition(0,0,0); atoms.push_back(std::move(atom));}
    ModelObject model(std::move(atoms)); model.SelectAtoms([](const auto & a){return a.GetSerialID()==1;});
    auto values=std::make_unique<double[]>(24*4*4);
    for(std::size_t r=0;r<24;++r) values[r]=problem.Input().observations[r];
    MapObject map({24,4,4},{1,1,1},{0,0,0},std::move(values));
    auto editor=model.EditAnalysis(); editor.ApplyJointResult(result,"observable-test");
    const SamplingPointList points{{23,{23,0,0},false}};
    const auto raw=core::detail::BuildLocalPotentialSampleList(map,points);
    editor.SetAtomLocalRawSamplingEntries(*model.FindAtomPtr(1),raw);
    const auto peeled=core::detail::BuildPostFitPeelingSamples(map,model,problem,std::nullopt,&result);
    ASSERT_TRUE(peeled.at(1).samples[0].response);
    editor.SetJointResult(result);
    const command_test::ScopedTempDir directory{"observable_roundtrip"};
    DataRepository repository(directory.path()/"saved.sqlite"); repository.SaveModel(model,"test");
    const auto loaded=repository.LoadModel("test");
    ASSERT_TRUE(loaded->GetAnalysisView().GetJointResult());
    EXPECT_EQ(joint_result_io::Encode(*loaded->GetAnalysisView().GetJointResult()),joint_result_io::Encode(result));
    EXPECT_THROW(AtomLocalPotentialView::For(*loaded->FindAtomPtr(2)).GetStageEstimate(FittingStage::Second),std::runtime_error);
    WriteJointAnalysisResult(result,directory.path()/"saved.json",directory.path()/"saved.csv");
    EXPECT_FALSE(std::filesystem::exists(directory.path()/"saved.contributions.csv"));
    std::ifstream csv(directory.path()/"saved.csv"); std::string line;
    std::getline(csv,line); std::getline(csv,line);
    EXPECT_TRUE(line.starts_with("\"1\",\""+result.components[0].id+"\","));
    const auto workset=core::detail::MakeJointFittingWorkset(model,problem);
    std::vector<int> visited;
    core::detail::FirstStageObserverForTesting()=[&](int id,std::string_view){visited.push_back(id);};
    core::FitOptions options; options.quiet_mode=true; options.thread_size=1;
    const auto first_stage_atoms=core::detail::CollectFirstStageAtoms(workset);
    core::RunPotentialSamplingWorkflow(map,model,first_stage_atoms,options.sampling_method,1);
    const auto initialization=core::detail::RunJointFirstStageInitializationFromPreparedSamples(model,workset,options);
    core::detail::FirstStageObserverForTesting()={};
    EXPECT_EQ(visited,(std::vector<int>{1}));
    EXPECT_EQ(initialization.atoms.size(),1u);
}

TEST(JointObservableProfileTest, PostJointLocalSamplesSeparateFixedChargeAndFilterLogResponses)
{
    using namespace rhbm_gem;
    auto input = Input();
    input.atom_ids.resize(1);
    input.support.resize(1);
    input.observations[23] -= .7;
    input.support[0].front().squared_distance = 0.0;
    input.support[0].back().squared_distance = 6.25;
    const core::JointProblem problem(std::move(input));
    const auto base = core::CaptureJointAnalysisResult(core::FitJointComponents(problem, {.6}));
    ASSERT_TRUE(base.components.front().state);
    std::vector<std::unique_ptr<AtomObject>> atoms;
    auto atom = std::make_unique<AtomObject>();
    atom->SetSerialID(1);
    atom->SetElement(Element::CARBON);
    atom->SetPosition(0, 0, 0);
    atoms.push_back(std::move(atom));
    ModelObject model(std::move(atoms));
    model.SelectAtoms([](const auto &) { return true; });
    auto editor = model.EditAnalysis();

    constexpr double joint_width = .85;
    constexpr double local_width = .42;
    constexpr double local_amplitude = 2.4;
    constexpr double negative_response = -.1;
    const std::array<double, 3> charges{0.0, .35, -.35};
    for (const double charge : charges)
    {
        auto result = base;
        result.components.front().state->b[0] = joint_width;
        result.components.front().state->log_b[0] = std::log(joint_width);
        result.components.front().state->ac[0] = 1.7;
        result.components.front().state->ac[1] = charge;
        editor.ApplyJointResult(result, "charge-correction-test");

        const auto center_joint = n::EvaluateKernel(0.0, joint_width, 2.5);
        const auto boundary_joint = n::EvaluateKernel(6.25, joint_width, 2.5);
        const auto center_local = n::EvaluateKernel(0.0, local_width, 2.5);
        const auto boundary_local = n::EvaluateKernel(6.25, local_width, 2.5);
        auto values = std::make_unique<double[]>(24 * 4 * 4);
        std::fill_n(values.get(), 24 * 4 * 4, 0.0);
        values[0] = local_amplitude * center_local.gaussian + charge * center_joint.charge;
        values[23] = local_amplitude * boundary_local.gaussian + charge * boundary_joint.charge;
        MapObject map({24, 4, 4}, {1, 1, 1}, {0, 0, 0}, std::move(values));

        const SamplingPointList points{
            {0.0, {0, 0, 0}, true},
            {2.5, {23, 0, 0}, true},
            {0.0, {0, 0, 0}, true}};
        auto raw = core::detail::BuildLocalPotentialSampleList(map, points);
        raw[2].response = charge * center_joint.charge + negative_response;
        editor.SetAtomLocalRawSamplingEntries(*model.FindAtomPtr(1), raw);
        const auto peeling = core::detail::BuildPostFitPeelingSamples(
            map, model, problem, std::nullopt, &result).at(1);
        const auto prepared = core::detail::BuildPostJointLocalSamples(
            map, model, problem, 1, peeling, result);
        ASSERT_EQ(prepared.samples.size(), 3u);
        EXPECT_DOUBLE_EQ(prepared.fixed_charge.width, joint_width);
        EXPECT_NE(prepared.fixed_charge.width, local_width);
        for (std::size_t i = 0; i < prepared.samples.size(); ++i)
        {
            const auto & sample = prepared.samples[i];
            EXPECT_EQ(sample.raw_sample_index, i);
            EXPECT_DOUBLE_EQ(sample.distance, points[i].distance);
            ASSERT_TRUE(sample.peeled_response);
            ASSERT_TRUE(sample.fixed_charge_response);
            ASSERT_TRUE(sample.gaussian_response);
            EXPECT_DOUBLE_EQ(*sample.gaussian_response,
                *sample.peeled_response - *sample.fixed_charge_response);
        }
        EXPECT_NEAR(*prepared.samples[0].fixed_charge_response,
            charge * center_joint.charge, 1e-12);
        EXPECT_NEAR(*prepared.samples[1].fixed_charge_response,
            charge * boundary_joint.charge, 1e-12);
        EXPECT_NEAR(*prepared.samples[0].gaussian_response,
            local_amplitude * center_local.gaussian, 1e-12);
        EXPECT_NEAR(*prepared.samples[1].gaussian_response,
            local_amplitude * boundary_local.gaussian, 1e-12);
        EXPECT_DOUBLE_EQ(*prepared.samples[2].gaussian_response, negative_response);

        const auto fitting_samples = prepared.BuildGaussianSamples(raw);
        ASSERT_EQ(fitting_samples.size(), 3u);
        EXPECT_DOUBLE_EQ(fitting_samples.back().response, negative_response);
        std::vector<double> responses;
        for (const auto & sample : fitting_samples) responses.push_back(sample.response);
        const core::detail::PreparedLocalGaussianDesign design{fitting_samples, 0.0, 2.5};
        const auto dataset = design.BuildDataset(responses, GaussianModel3D{0.0, 1.0, 0.0});
        EXPECT_EQ(dataset.y.size(), 2);
        EXPECT_TRUE(dataset.y.allFinite());
    }
}

namespace {
core::JointProblemInput TwoRowHalo(bool selected=false,bool noise=false)
{
    auto input=Input(noise);
    input.atom_ids.push_back("4"); input.support.push_back({{0,.4},{1,1.2}});
    if(selected) input.selection_domain->target_indices.push_back(3);
    for(const auto & p:input.support[3])
    {
        const auto basis=n::EvaluateKernel(p.squared_distance,.65,2.5);
        input.observations[p.row]+=basis.gaussian+.2*basis.charge;
    }
    return input;
}
}
TEST(JointObservableProfileTest, HaloNonuniquenessDoesNotCertifyItsParameters)
{
    const core::JointProblem problem(TwoRowHalo());
    const auto fit=core::FitJointComponents(problem,{.6,.7,.8,.65});
    ASSERT_FALSE(fit.assembled_state);
    ASSERT_TRUE(fit.components[0].candidate_diagnostic);
    EXPECT_TRUE(fit.components[0].candidate_diagnostic->candidate_available);
    EXPECT_NE(fit.components[0].candidate_diagnostic->endpoint_verification_status,
        rhbm_gem::JointCheckStatus::Passed);
    EXPECT_NE(fit.components[0].candidate_diagnostic->numerical_qualification_status,
        rhbm_gem::JointCheckStatus::Passed);
    EXPECT_NE(fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(fit.TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::Unavailable);
    const auto saved=core::CaptureJointAnalysisResult(fit);
    const auto stages=rhbm_gem::data_internal::BuildJointStageEstimates(saved,"target-test");
    EXPECT_FALSE(stages.at(1).point); EXPECT_FALSE(stages.at(4).point);
    const auto text=rhbm_gem::joint_result_io::Encode(saved);
    auto value=boost::json::parse(text);
    EXPECT_EQ(value.at("schema_version").as_int64(),11);
    auto downgraded=value;
    downgraded.as_object()["schema_version"]=7;
    EXPECT_THROW(rhbm_gem::joint_result_io::Decode(boost::json::serialize(downgraded)),std::invalid_argument);
    const auto restored=rhbm_gem::joint_result_io::Decode(text);
    EXPECT_EQ(restored.target_runtime_convergence,saved.target_runtime_convergence);
    ASSERT_TRUE(restored.components[0].candidate_diagnostic);
    EXPECT_NE(restored.components[0].candidate_diagnostic->endpoint_verification_status,
        rhbm_gem::JointCheckStatus::Passed);
    EXPECT_FALSE(restored.components[0].state);
    value.as_object()["schema_version"]=4;
    value.as_object().at("selection_domain").as_object()["contributor_policy"]="all-non-hydrogen";
    value.as_object().at("selection_domain").as_object()["halo_eligibility"]=nullptr;
    auto & solver=value.as_object().at("metadata").as_object().at("solver").as_object();
    solver["contract_version"]="joint-solver-provenance-v4";
    solver["fixed_neighbor_policy_version"]="fixed-neighbor-production-v3";
    for(auto & component:value.as_object().at("components").as_array())
        component.as_object().erase("candidate_diagnostic");
    const auto old=rhbm_gem::joint_result_io::Decode(boost::json::serialize(value));
    EXPECT_FALSE(old.target_evidence); EXPECT_FALSE(old.components[0].target_evidence);
    EXPECT_FALSE(old.components[0].candidate_diagnostic);
    EXPECT_EQ(old.target_runtime_convergence,rhbm_gem::JointCheckStatus::NotRun);
}
TEST(JointObservableProfileTest, TwoRowGeometryDoesNotPromoteAnUnqualifiedTarget)
{
    const core::JointProblem problem(TwoRowHalo(true));
    const auto fit=core::FitJointComponents(problem,{.6,.7,.8,.65});
    EXPECT_FALSE(fit.assembled_state);
    ASSERT_TRUE(fit.components[0].candidate_diagnostic);
    const auto & candidate=*fit.components[0].candidate_diagnostic;
    EXPECT_TRUE(candidate.candidate_available);
    EXPECT_EQ(fit.TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::Unavailable);
    rhbm_gem::JointState candidate_state; candidate_state.ac=candidate.candidate_ac;
    candidate_state.log_b=candidate.candidate_eta; candidate_state.b.resize(candidate.candidate_eta.size());
    for(std::size_t k=0;k<candidate.candidate_eta.size();++k)
        candidate_state.b[k]=std::exp(candidate.candidate_eta[k]);
    const auto evidence=n::AssessTargetState(problem,problem.ParameterLayout(),candidate_state,
        fit.components[0].evidence,fit.components[0].ranks,rhbm_gem::JointEvidenceScope::ComponentLocal,
        fit.components[0].rows.size());
    const auto target=std::find_if(evidence.checks.begin(),evidence.checks.end(),[](const auto & check) {
        return check.name=="target-identifiability";
    });
    ASSERT_NE(target,evidence.checks.end());
    EXPECT_EQ(target->status,rhbm_gem::JointCheckStatus::Unavailable);
    EXPECT_EQ(evidence.atoms[3].status,rhbm_gem::JointCheckStatus::Unavailable);
    EXPECT_EQ(evidence.atoms[3].status,rhbm_gem::JointCheckStatus::Unavailable);
}
TEST(JointObservableProfileTest, UnqualifiedHaloRepresentativesDoNotBecomePublicEndpoints)
{
    const core::JointProblem problem(TwoRowHalo(false,true));
    const auto first=core::FitJointComponents(problem,{.6,.7,.8,.65});
    const auto second=core::FitJointComponents(problem,{.6,.7,.8,.7});
    ASSERT_EQ(first.components.size(),1);
    ASSERT_EQ(second.components.size(),1);
    ASSERT_TRUE(first.components[0].candidate_diagnostic);
    ASSERT_TRUE(second.components[0].candidate_diagnostic);
    const auto & a=*first.components[0].candidate_diagnostic;
    const auto & b=*second.components[0].candidate_diagnostic;
    ASSERT_TRUE(a.candidate_available && b.candidate_available);
    ASSERT_GE(a.candidate_ac.size(),2); ASSERT_GE(b.candidate_ac.size(),2);
    ASSERT_FALSE(first.assembled_state); ASSERT_FALSE(second.assembled_state);
    EXPECT_EQ(first.TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::Unavailable);
    EXPECT_EQ(second.TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::Unavailable);
    EXPECT_NEAR(a.candidate_ac[0],b.candidate_ac[0],1e-10);
    EXPECT_NEAR(a.candidate_ac[1],b.candidate_ac[1],1e-10);
    EXPECT_NEAR(std::exp(a.candidate_eta[0]),std::exp(b.candidate_eta[0]),1e-10);
}
TEST(JointObservableProfileTest, ActiveBoundaryDoesNotAcquireTargetCertificate)
{
    auto input=Input(); input.atom_ids.resize(1); input.support.resize(1); input.observations[23]-=.7;
    const core::JointProblem problem(std::move(input));
    auto fit=core::FitJointComponents(problem,{.6});
    ASSERT_TRUE(fit.assembled_state);
    auto state=*fit.assembled_state; state.ac[0]=0;
    const auto evidence=n::AssessTargetState(problem,problem.ParameterLayout(),state,fit.evidence,fit.ranks,
        rhbm_gem::JointEvidenceScope::AssembledGlobal,problem.Input().observations.size());
    EXPECT_EQ(n::TargetConvergence(evidence,true),rhbm_gem::JointCheckStatus::Unavailable);
}
TEST(JointObservableProfileTest, TargetContractKeepsCorrectionRequiredOnlyForHistoricalEvidence)
{
    using Status=rhbm_gem::JointCheckStatus;
    using Scope=rhbm_gem::JointEvidenceScope;
    rhbm_gem::JointTargetEvidence evidence;
    evidence.contract="target-estimability-v1";
    evidence.checks={
        {"conditional-ac-optimality",Status::Passed,Scope::ComponentLocal,{},{},{}},
        {"independent-replay",Status::Passed,Scope::ComponentLocal,{},{},{}},
        {"width-stationarity",Status::Passed,Scope::ComponentLocal,{},{},{}},
        {"design-identifiability",Status::Passed,Scope::ComponentLocal,{},{},{}},
        {"local-correction",Status::Failed,Scope::ComponentLocal,{},{},{}},
        {"numerical-identifiability",Status::Passed,Scope::ComponentLocal,{},{},{}},
        {"target-identifiability",Status::Passed,Scope::ComponentLocal,{},{},{}},
        {"target-local-correction",Status::Failed,Scope::ComponentLocal,{},{},{}}};
    EXPECT_EQ(n::TargetConvergence(evidence,true),Status::Failed);
    evidence.contract="target-estimability-v2";
    EXPECT_EQ(n::TargetConvergence(evidence,true),Status::Passed);
    evidence.contract="target-estimability-v3";
    EXPECT_EQ(n::TargetConvergence(evidence,true),Status::Passed);
    evidence.contract="target-estimability-v4";
    EXPECT_EQ(n::TargetConvergence(evidence,true),Status::Passed);
}
TEST(JointObservableProfileTest, TargetEvidenceRequiresSelectionMetadata)
{
    auto input=Input(); input.atom_ids.resize(1); input.support.resize(1); input.selection_domain.reset();
    const auto fit=core::FitJointComponents(core::JointProblem(input),{.6});
    EXPECT_EQ(fit.TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::NotRun);
    EXPECT_FALSE(fit.target_evidence);
}

TEST(JointObservableProfileTest, FullAbcPointEvidenceMatchesFreshAssessment)
{
    auto input=Input(); input.atom_ids.resize(1); input.support.resize(1); input.observations[23]-=.7;
    const core::JointProblem problem(std::move(input));
    n::AssessmentWorkForTesting()={};
    const auto fit=core::FitJointComponents(problem,{.6});
    ASSERT_EQ(fit.components.size(),1);
    ASSERT_TRUE(fit.components[0].target_evidence);
    ASSERT_TRUE(fit.target_evidence);
    EXPECT_EQ(fit.target_evidence->contract,rhbm_gem::JointTargetPointValidationContractVersion);
    const auto & component=fit.components.front();
    static const std::vector<rhbm_gem::JointBlockParameterStability> no_stability;
    const auto & stability=component.parameter_stability ? component.parameter_stability->blocks : no_stability;
    const auto fresh=n::AssessTargetPointState(problem,problem.ParameterLayout(),fit.assembled_state,
        fit.RuntimeConvergence(),component.local_conditioning,stability,
        rhbm_gem::JointEvidenceScope::AssembledGlobal,problem.Input().observations.size());
    ASSERT_EQ(fit.target_evidence->checks.size(),fresh.checks.size());
    for(std::size_t i=0;i<fresh.checks.size();++i)
    {
        EXPECT_EQ(fit.target_evidence->checks[i].name,fresh.checks[i].name);
        EXPECT_EQ(fit.target_evidence->checks[i].status,fresh.checks[i].status);
        EXPECT_EQ(fit.target_evidence->checks[i].scope,fresh.checks[i].scope);
        EXPECT_EQ(fit.target_evidence->checks[i].reason,fresh.checks[i].reason);
    }
    ASSERT_EQ(fit.target_evidence->point_atoms.size(),fresh.point_atoms.size());
    for(std::size_t i=0;i<fresh.point_atoms.size();++i)
    {
        EXPECT_EQ(fit.target_evidence->point_atoms[i].atom,fresh.point_atoms[i].atom);
        EXPECT_EQ(fit.target_evidence->point_atoms[i].parameter_validity,fresh.point_atoms[i].parameter_validity);
        EXPECT_EQ(fit.target_evidence->point_atoms[i].local_conditioning,fresh.point_atoms[i].local_conditioning);
        EXPECT_EQ(fit.target_evidence->point_atoms[i].parameter_stability,fresh.point_atoms[i].parameter_stability);
        EXPECT_EQ(fit.target_evidence->point_atoms[i].point_quality,fresh.point_atoms[i].point_quality);
    }
    EXPECT_TRUE(fit.target_evidence->atoms.empty());
    EXPECT_TRUE(fit.target_evidence->ranks.empty());
    EXPECT_EQ(problem.ObservationScale(),fit.observation_scale);
}

TEST(JointObservableProfileTest, TargetGeometryNumericalOracleCapturesOldContract)
{
    using Status=rhbm_gem::JointCheckStatus;
    auto input=Input(true); input.atom_ids.resize(1); input.support.resize(1); input.observations[23]-=.7;
    const core::JointProblem problem(std::move(input));
    const auto fit=core::FitJointComponents(problem,{.65});
    ASSERT_EQ(fit.components.size(),1u);
    ASSERT_TRUE(fit.components.front().state);
    ASSERT_TRUE(fit.components.front().target_evidence);
    const auto & component=fit.components.front();
    const auto & layout=problem.ParameterLayout();
    const auto & data=core::JointProblemAccess::Get(problem);
    const auto domain=FullDomain(data.domain,layout);
    const auto context=FullContext(data.context,layout,data.y.size());
    const n::Vector y=n::SelectValues(data.y,{layout.informative_rows.begin(),layout.informative_rows.end()});
    const n::Vector eta=Eigen::Map<const n::Vector>(component.state->log_b.data(),
        static_cast<Eigen::Index>(component.state->log_b.size()));
    const n::Vector beta=Eigen::Map<const n::Vector>(component.state->ac.data(),
        static_cast<Eigen::Index>(component.state->ac.size()));
    const auto endpoint=n::EvaluateState(domain,y,eta,beta,context);
    ASSERT_TRUE(endpoint.valid);
    n::Matrix jacobian= n::Matrix::Zero(y.size(),3*static_cast<Eigen::Index>(layout.full_atoms.size()));
    for(Eigen::Index atom=0;atom<static_cast<Eigen::Index>(layout.full_atoms.size());++atom)
    {
        jacobian.middleCols(3*atom,2)=endpoint.x.middleCols(2*atom,2);
        jacobian.col(3*atom+2)=endpoint.derivative.col(2*atom)*beta(2*atom)+
            endpoint.derivative.col(2*atom+1)*beta(2*atom+1);
    }
    const n::Vector scales=jacobian.colwise().norm();
    ASSERT_TRUE((scales.array()>0).all());
    const n::Matrix scaled=jacobian*scales.cwiseInverse().asDiagonal();
    const Eigen::JacobiSVD<n::Matrix> oracle(scaled,Eigen::ComputeThinU|Eigen::ComputeFullV);
    const auto old_geometry=OriginalGeometryOracle(domain,y,*component.state,context);
    ASSERT_TRUE(old_geometry.valid)<<old_geometry.reason;
    EXPECT_FALSE(old_geometry.threshold_sensitive);
    EXPECT_EQ(old_geometry.rank,3);
    EXPECT_GT(old_geometry.threshold,0);
    EXPECT_EQ(old_geometry.scales.size(),3);
    EXPECT_TRUE(old_geometry.scales.isApprox(scales,1e-12));
    ASSERT_EQ(old_geometry.singular_values.size(),oracle.singularValues().size());
    EXPECT_TRUE(old_geometry.singular_values.isApprox(oracle.singularValues(),1e-11));
    EXPECT_TRUE(component.ranks.empty());
    EXPECT_TRUE(component.target_evidence->ranks.empty());
    EXPECT_TRUE(component.target_evidence->atoms.empty());
    ASSERT_EQ(component.target_evidence->point_atoms.size(),1u);
    EXPECT_EQ(component.target_evidence->point_atoms.front().parameter_validity,Status::Passed);
    EXPECT_EQ(component.target_evidence->point_atoms.front().atom,0u);
    const auto old_nullspace=old_geometry.right_vectors.rightCols(old_geometry.right_vectors.cols()-old_geometry.rank);
    EXPECT_LE(old_nullspace.topRows(3).norm(),1e-10);
    EXPECT_EQ(component.TargetRuntimeConvergence(),Status::Passed);
    auto boundary=*component.state;
    boundary.ac[0]=0;
    const auto boundary_geometry=OriginalGeometryOracle(domain,y,boundary,context);
    EXPECT_FALSE(boundary_geometry.valid);
    EXPECT_EQ(boundary_geometry.reason,"active-amplitude-boundary");
}

TEST(JointObservableProfileTest, TargetGeometryOracleMapsFullAbcHaloCoordinates)
{
    auto input=Input();
    input.atom_ids.resize(1); input.support.resize(1); input.observations[23]-=.7;
    input.atom_ids.push_back("4"); input.support.emplace_back();
    for(std::size_t row=0;row<23;++row)
    {
        const double square=std::pow(.1+.09*static_cast<double>(row),2);
        input.support.back().push_back({row,square});
        const auto basis=n::EvaluateKernel(square,.9,2.5);
        input.observations[row]+=.7*basis.gaussian+.15*basis.charge;
    }
    const core::JointProblem problem(std::move(input));
    EXPECT_EQ(problem.ParameterLayout().full_atoms,(std::vector<std::size_t>{0,1}));
    const auto fit=core::FitJointComponents(problem,{.6,.9});
    ASSERT_TRUE(fit.assembled_state);
    ASSERT_TRUE(fit.target_evidence);
    ASSERT_EQ(fit.target_evidence->point_atoms.size(),2u);
    const auto target=std::find_if(fit.target_evidence->point_atoms.begin(),fit.target_evidence->point_atoms.end(),[](const auto & atom) {
        return atom.atom==0;
    });
    const auto halo=std::find_if(fit.target_evidence->point_atoms.begin(),fit.target_evidence->point_atoms.end(),[](const auto & atom) {
        return atom.atom==1;
    });
    ASSERT_NE(target,fit.target_evidence->point_atoms.end()); ASSERT_NE(halo,fit.target_evidence->point_atoms.end());
    EXPECT_EQ(target->parameter_validity,rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(halo->parameter_validity,rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(target->component_convergence,rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(halo->component_convergence,rhbm_gem::JointCheckStatus::Passed);
    EXPECT_TRUE(std::binary_search(problem.Input().selection_domain->target_indices.begin(),
        problem.Input().selection_domain->target_indices.end(),0u));
    EXPECT_FALSE(std::binary_search(problem.Input().selection_domain->target_indices.begin(),
        problem.Input().selection_domain->target_indices.end(),1u));
    const auto & data=core::JointProblemAccess::Get(problem);
    const auto & layout=problem.ParameterLayout();
    const auto domain=FullDomain(data.domain,layout);
    const auto context=FullContext(data.context,layout,data.y.size());
    const auto y=n::SelectValues(data.y,n::Indices{layout.informative_rows.begin(),layout.informative_rows.end()});
    const auto old_geometry=OriginalGeometryOracle(domain,y,*fit.assembled_state,context);
    ASSERT_TRUE(old_geometry.valid)<<old_geometry.reason;
    EXPECT_FALSE(old_geometry.threshold_sensitive);
    EXPECT_EQ(old_geometry.rank,6);
    ASSERT_TRUE(fit.target_evidence);
    EXPECT_TRUE(fit.target_evidence->ranks.empty());
    EXPECT_EQ(target->parameter_validity,rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(halo->parameter_validity,rhbm_gem::JointCheckStatus::Passed);
    const auto old_nullspace=old_geometry.right_vectors.rightCols(old_geometry.right_vectors.cols()-old_geometry.rank);
    EXPECT_EQ(old_nullspace.cols(),0);
}

TEST(JointObservableProfileTest, PoorLocalQualityDoesNotRemoveRetainedHaloPoint)
{
    auto input=Input(); input.atom_ids.resize(1); input.support.resize(1); input.observations[23]-=.7;
    input.atom_ids.push_back("4"); input.support.emplace_back();
    for(std::size_t row=0;row<23;++row)
    {
        const double square=std::pow(.1+.09*static_cast<double>(row),2);
        input.support.back().push_back({row,square});
        const auto basis=n::EvaluateKernel(square,.9,2.5);
        input.observations[row]+=.7*basis.gaussian+.15*basis.charge;
    }
    const core::JointProblem problem(std::move(input));
    auto saved=core::CaptureJointAnalysisResult(core::FitJointComponents(problem,{.6,.9}));
    ASSERT_EQ(saved.runtime_convergence,rhbm_gem::JointCheckStatus::Passed);
    auto & points=saved.components.front().target_evidence->point_atoms;
    const auto halo=std::find_if(points.begin(),points.end(),[](const auto & point){return point.atom==1;});
    ASSERT_NE(halo,points.end());
    halo->local_conditioning="Warning";
    halo->parameter_stability="Sensitive";
    halo->point_quality="Warning";
    halo->reason="conditioning-or-stability-warning";
    EXPECT_EQ(rhbm_gem::core::joint_component::TargetConvergence(
        saved.components.front().target_evidence,true),rhbm_gem::JointCheckStatus::Passed);
    const auto stages=rhbm_gem::data_internal::BuildJointStageEstimates(saved,"point-quality-test");
    ASSERT_TRUE(stages.at(4).point);
    EXPECT_EQ(stages.at(4).source.role,rhbm_gem::FittingRole::Halo);
    EXPECT_EQ(stages.at(4).convergence,rhbm_gem::JointCheckStatus::Passed);
}

TEST(JointObservableProfileTest, DerivedTargetEvidenceRequiresMatchingQualifiedRankEvidence)
{
    using Status=rhbm_gem::JointCheckStatus;
    using Scope=rhbm_gem::JointEvidenceScope;
    auto input=Input(true); input.atom_ids.resize(1); input.support.resize(1); input.observations[23]-=.7;
    const core::JointProblem problem(std::move(input));
    const auto fit=core::FitJointComponents(problem,{.65});
    ASSERT_TRUE(fit.assembled_state);
    const auto & data=core::JointProblemAccess::Get(problem);
    const n::Vector eta=Eigen::Map<const n::Vector>(fit.assembled_state->log_b.data(),
        static_cast<Eigen::Index>(fit.assembled_state->log_b.size()));
    const n::Vector beta=Eigen::Map<const n::Vector>(fit.assembled_state->ac.data(),
        static_cast<Eigen::Index>(fit.assembled_state->ac.size()));
    const auto legacy_assessment=n::AssessProfile(data.domain,data.y,eta,data.context,&beta);
    auto legacy_checks=n::AssessmentEvidence(legacy_assessment,Scope::AssembledGlobal);
    for(const auto & check:fit.evidence)
        if(check.name=="assembly-mapping" || check.name=="assembled-consistency") legacy_checks.push_back(check);
    const auto legacy_ranks=n::AssessmentRanks(legacy_assessment,Scope::AssembledGlobal,data.input.get(),
        &problem.ParameterLayout(),problem.Input().observations.size());
    const auto evidence=n::DeriveTargetIdentifiability(problem,problem.ParameterLayout(),fit.assembled_state,
        legacy_checks,legacy_ranks,Scope::AssembledGlobal,problem.Input().observations.size());
    ASSERT_EQ(evidence.contract,"target-estimability-v4");
    ASSERT_TRUE(evidence.column_scales.empty());
    ASSERT_EQ(evidence.ranks.size(),2u);
    EXPECT_EQ(evidence.ranks[0].name,"design");
    EXPECT_EQ(evidence.ranks[1].name,"projected-width");
    ASSERT_EQ(evidence.atoms.size(),1u);
    EXPECT_EQ(evidence.atoms[0].status,Status::Passed);
    EXPECT_EQ(evidence.atoms[0].reason,"derived-from-component-identifiability");
    EXPECT_FALSE(evidence.atoms[0].null_space_leakage);
    EXPECT_EQ(n::TargetConvergence(evidence,true),Status::Passed);

    const auto derive_with=[&](const std::optional<rhbm_gem::JointState> & state,
        const std::vector<rhbm_gem::JointCheck> & checks,const std::vector<rhbm_gem::JointRankEvidence> & ranks,
        const rhbm_gem::JointParameterLayout & layout,std::size_t rows) {
        return n::DeriveTargetIdentifiability(problem,layout,state,checks,ranks,Scope::AssembledGlobal,rows);
    };
    const auto derive=[&](const std::optional<rhbm_gem::JointState> & state,
        const std::vector<rhbm_gem::JointCheck> & checks,const std::vector<rhbm_gem::JointRankEvidence> & ranks) {
        return derive_with(state,checks,ranks,problem.ParameterLayout(),problem.Input().observations.size());
    };
    const auto status_reason=[](const auto & value) {
        const auto check=std::find_if(value.checks.begin(),value.checks.end(),[](const auto & item) {
            return item.name=="target-identifiability";
        });
        return check==value.checks.end() ? std::string{} : check->reason;
    };
    EXPECT_EQ(status_reason(derive({},legacy_checks,legacy_ranks)),"missing-component-state");
    auto no_design=legacy_ranks;
    no_design.erase(std::remove_if(no_design.begin(),no_design.end(),[](const auto & rank) {
        return rank.name=="design";
    }),no_design.end());
    EXPECT_EQ(status_reason(derive(fit.assembled_state,legacy_checks,no_design)),"missing-rank-evidence");
    auto no_width=legacy_ranks;
    no_width.erase(std::remove_if(no_width.begin(),no_width.end(),[](const auto & rank) {
        return rank.name=="projected-width";
    }),no_width.end());
    EXPECT_EQ(status_reason(derive(fit.assembled_state,legacy_checks,no_width)),"missing-rank-evidence");
    auto deficient=legacy_ranks;
    std::find_if(deficient.begin(),deficient.end(),[](const auto & rank) {
        return rank.name=="projected-width";
    })->rank=0;
    EXPECT_EQ(status_reason(derive(fit.assembled_state,legacy_checks,deficient)),"component-rank-deficient");
    auto near_threshold=legacy_ranks;
    auto design=std::find_if(near_threshold.begin(),near_threshold.end(),[](const auto & rank) {
        return rank.name=="design";
    });
    ASSERT_NE(design,near_threshold.end()); ASSERT_GE(design->singular_values.size(),2u);
    design->singular_values.back()=std::nextafter(design->threshold,std::numeric_limits<double>::infinity());
    EXPECT_EQ(status_reason(derive(fit.assembled_state,legacy_checks,near_threshold)),"rank-threshold-sensitive");
    auto changed_state=*fit.assembled_state; changed_state.ac[0]+=.01;
    EXPECT_EQ(status_reason(derive(changed_state,legacy_checks,legacy_ranks)),"rank-context-mismatch");
    auto wrong_layout=problem.ParameterLayout(); wrong_layout.informative_rows.pop_back();
    EXPECT_EQ(status_reason(derive_with(fit.assembled_state,legacy_checks,legacy_ranks,wrong_layout,
        problem.Input().observations.size())),"wrong-parameter-layout");
    auto nonfinite=*fit.assembled_state; nonfinite.log_b[0]=std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(status_reason(derive(nonfinite,legacy_checks,legacy_ranks)),"nonfinite-parameter-state");
    auto boundary=*fit.assembled_state; boundary.ac[0]=0;
    EXPECT_EQ(status_reason(derive(boundary,legacy_checks,legacy_ranks)),"active-amplitude-boundary");
    auto missing_check=legacy_checks;
    missing_check.erase(std::remove_if(missing_check.begin(),missing_check.end(),[](const auto & check) {
        return check.name=="independent-replay";
    }),missing_check.end());
    EXPECT_EQ(status_reason(derive(fit.assembled_state,missing_check,legacy_ranks)),"missing-component-qualification");
    EXPECT_EQ(status_reason(derive_with(fit.assembled_state,legacy_checks,legacy_ranks,
        problem.ParameterLayout(),problem.Input().observations.size()-1)),"observation-domain-mismatch");

    auto other_input=problem.Input(); other_input.observations[0]+=1e-6;
    const core::JointProblem other_problem(std::move(other_input));
    const auto other=n::DeriveTargetIdentifiability(other_problem,other_problem.ParameterLayout(),fit.assembled_state,
        legacy_checks,legacy_ranks,Scope::AssembledGlobal,other_problem.Input().observations.size());
    EXPECT_EQ(status_reason(other),"rank-context-mismatch");
}

TEST(JointObservableProfileTest, DerivedTargetEvidencePreservesNotRunWithoutSelectionMetadata)
{
    auto input=Input(); input.selection_domain.reset(); input.atom_ids.resize(1); input.support.resize(1);
    const core::JointProblem problem(std::move(input));
    const auto fit=core::FitJointComponents(problem,{.6});
    EXPECT_FALSE(fit.target_evidence);
    EXPECT_EQ(fit.TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::NotRun);
}

TEST(JointObservableProfileTest, TargetMemoryShapeSeparatesTargetsFromHaloColumns)
{
    const auto input=TwoRowHalo(false);
    const core::JointProblem problem(input);
    const auto estimate=n::TargetIdentifiabilityMemoryEstimate(input,problem.ParameterLayout(),
        problem.ParameterLayout().informative_rows.size());
    EXPECT_EQ(estimate.shape.component_atoms,4u);
    EXPECT_EQ(estimate.shape.component_target_atoms,1u);
    EXPECT_EQ(estimate.shape.component_halo_atoms,3u);
    EXPECT_EQ(estimate.shape.component_design_columns,8u);
    EXPECT_EQ(estimate.shape.component_width_columns,4u);
    EXPECT_EQ(estimate.stage,n::MemoryStage::TargetIdentifiability);
    EXPECT_TRUE(estimate.attributed_stages.empty());
    EXPECT_EQ(estimate.dense_factor_bytes,0u);
    EXPECT_EQ(estimate.right_singular_vector_bytes,0u);
}

TEST(JointObservableProfileTest, LargeDerivedTargetEvidenceHasLinearMemoryAndNoGeometryStage)
{
    constexpr std::size_t atoms=37'406;
    core::JointProblemInput input; input.observations.assign(8,0.0); input.row_ids.reserve(8);
    for(std::size_t row=0;row<input.observations.size();++row) input.row_ids.push_back("row-"+std::to_string(row));
    input.atom_ids.reserve(atoms); input.support.resize(atoms);
    for(std::size_t atom=0;atom<atoms;++atom)
    {
        input.atom_ids.push_back("atom-"+std::to_string(atom)); input.support[atom].push_back({0,0.0});
    }
    input.selection_domain.emplace(); input.selection_domain->target_indices.reserve(atoms);
    for(std::size_t atom=0;atom<atoms;++atom) input.selection_domain->target_indices.push_back(atom);
    const core::JointProblem problem(std::move(input));
    const auto & data=core::JointProblemAccess::Get(problem);
    rhbm_gem::JointState state; state.ac.assign(2*atoms,1.0); state.b.assign(atoms,.55);
    state.log_b.assign(atoms,std::log(.55));
    state.width_gradient.assign(atoms,0.0);
    const auto memory=n::TargetIdentifiabilityMemoryEstimate(problem.Input(),problem.ParameterLayout(),8);
    EXPECT_EQ(memory.status,n::MemoryStatus::Pass);
    EXPECT_EQ(memory.shape.component_target_atoms,atoms);
    EXPECT_LT(memory.estimated_peak_bytes,50'000'000u);
    EXPECT_EQ(memory.dense_factor_bytes,0u);
    n::AssessmentWorkForTesting()={};
    const auto evidence=n::AssessTargetState(problem,problem.ParameterLayout(),state,{},{},
        rhbm_gem::JointEvidenceScope::AssembledGlobal,data.input->observations.size());
    const auto identifiability=std::find_if(evidence.checks.begin(),evidence.checks.end(),[](const auto & check) {
        return check.name=="target-identifiability";
    });
    ASSERT_NE(identifiability,evidence.checks.end());
    EXPECT_EQ(identifiability->status,rhbm_gem::JointCheckStatus::Unavailable);
    EXPECT_EQ(identifiability->reason,"missing-rank-evidence");
    EXPECT_EQ(evidence.atoms.size(),atoms);
    EXPECT_TRUE(std::all_of(evidence.atoms.begin(),evidence.atoms.end(),[](const auto & atom) {
        return atom.status==rhbm_gem::JointCheckStatus::Unavailable;
    }));
    EXPECT_TRUE(evidence.ranks.empty());
    EXPECT_EQ(n::TargetConvergence(evidence,true),rhbm_gem::JointCheckStatus::Unavailable);
    const auto & stages=n::AssessmentWorkForTesting().stages;
    const auto preflight=std::find_if(stages.begin(),stages.end(),[](const auto & stage) {
        return stage.name=="derived-target-memory-preflight";
    });
    ASSERT_NE(preflight,stages.end());
    EXPECT_EQ(preflight->memory_status,n::MemoryStatus::Pass);
    EXPECT_EQ(preflight->calls,1); EXPECT_EQ(preflight->completed_calls,1);
    EXPECT_LT(preflight->estimated_peak_bytes,preflight->budget_bytes);
    EXPECT_TRUE(preflight->attributed_memory_stages.empty());
    const auto derived=std::find_if(stages.begin(),stages.end(),[](const auto & stage) {
        return stage.name=="derived-target-identifiability";
    });
    ASSERT_NE(derived,stages.end());
    EXPECT_EQ(derived->calls,1); EXPECT_EQ(derived->completed_calls,1);
    EXPECT_EQ(std::find_if(stages.begin(),stages.end(),[](const auto & stage) {
        return stage.name=="target-local-correction";
    }),stages.end());
}

TEST(JointObservableProfileTest, OldGeometryOracleAndDerivedEvidenceKeepSeparateThresholdContracts)
{
    auto input=Input(); input.atom_ids.resize(1); input.support.resize(1); input.observations[23]-=.7;
    const core::JointProblem problem(std::move(input));
    const auto fit=core::FitJointComponents(problem,{.6});
    ASSERT_TRUE(fit.assembled_state);
    const auto & data=core::JointProblemAccess::Get(problem); const auto & layout=problem.ParameterLayout();
    auto context=FullContext(data.context,layout,24);
    const auto domain=FullDomain(data.domain,layout);
    const auto old_geometry=OriginalGeometryOracle(domain,data.y,*fit.assembled_state,context);
    ASSERT_TRUE(old_geometry.valid);
    auto near_context=context;
    near_context.rank.rows=static_cast<Eigen::Index>(old_geometry.singular_values.tail(1)(0)/
        (std::numeric_limits<double>::epsilon()*old_geometry.singular_values(0)));
    const auto old_near=OriginalGeometryOracle(domain,data.y,*fit.assembled_state,near_context);
    ASSERT_TRUE(old_near.valid);
    EXPECT_TRUE(old_near.threshold_sensitive);
    ASSERT_TRUE(fit.target_evidence);
    EXPECT_EQ(fit.TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(fit.target_evidence->point_atoms.front().component_convergence,rhbm_gem::JointCheckStatus::Passed);

    const n::Vector eta=Eigen::Map<const n::Vector>(fit.assembled_state->log_b.data(),
        static_cast<Eigen::Index>(fit.assembled_state->log_b.size()));
    const n::Vector beta=Eigen::Map<const n::Vector>(fit.assembled_state->ac.data(),
        static_cast<Eigen::Index>(fit.assembled_state->ac.size()));
    const auto legacy_assessment=n::AssessProfile(data.domain,data.y,eta,data.context,&beta);
    auto checks=n::AssessmentEvidence(legacy_assessment,rhbm_gem::JointEvidenceScope::AssembledGlobal);
    for(const auto & check:fit.evidence)
        if(check.name=="assembly-mapping" || check.name=="assembled-consistency") checks.push_back(check);
    auto near_ranks=n::AssessmentRanks(legacy_assessment,rhbm_gem::JointEvidenceScope::AssembledGlobal,
        data.input.get(),&layout,problem.Input().observations.size());
    const auto design=std::find_if(near_ranks.begin(),near_ranks.end(),[](const auto & rank) {
        return rank.name=="design";
    });
    ASSERT_NE(design,near_ranks.end()); ASSERT_GE(design->singular_values.size(),2u);
    design->singular_values.back()=std::nextafter(design->threshold,std::numeric_limits<double>::infinity());
    const auto evidence=n::AssessTargetState(problem,layout,fit.assembled_state,checks,near_ranks,
        rhbm_gem::JointEvidenceScope::AssembledGlobal,problem.Input().observations.size());
    EXPECT_EQ(n::TargetConvergence(evidence,true),rhbm_gem::JointCheckStatus::Unavailable);
    EXPECT_EQ(evidence.atoms[0].reason,"rank-threshold-sensitive");
}

TEST(JointObservableProfileTest, TargetEvidenceUsesEveryComponentAndTheSavedState)
{
    core::JointProblemInput input;
    input.atom_ids={"1","2"}; input.support.resize(2);
    input.selection_domain.emplace(); input.selection_domain->target_indices={0,1};
    for(int atom=0;atom<2;++atom) for(int k=0;k<40;++k)
    {
        const double square=.003*k*k;
        const auto basis=n::EvaluateKernel(square,.5,2.5);
        input.support[static_cast<std::size_t>(atom)].push_back({input.observations.size(),square});
        input.row_ids.push_back(std::to_string(input.observations.size()));
        input.observations.push_back(2*basis.gaussian+.2*basis.charge);
    }
    const core::JointProblem problem(std::move(input));
    const auto fit=core::FitJointComponents(problem,{.5,.5});
    ASSERT_EQ(fit.components.size(),2);
    ASSERT_TRUE(fit.assembled_state);
    ASSERT_TRUE(fit.objective);
    EXPECT_EQ(fit.TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    ASSERT_TRUE(fit.target_evidence);
    for(const auto & component:fit.components)
    {
        ASSERT_TRUE(component.state);
        ASSERT_TRUE(component.candidate_diagnostic);
        EXPECT_TRUE(component.candidate_diagnostic->candidate_available);
        EXPECT_EQ(component.candidate_diagnostic->endpoint_verification_status,rhbm_gem::JointCheckStatus::Passed);
        EXPECT_EQ(component.candidate_diagnostic->numerical_qualification_status,rhbm_gem::JointCheckStatus::Passed);
    }
    std::vector<rhbm_gem::JointBlockConditioning> conditioning;
    std::vector<rhbm_gem::JointBlockParameterStability> stability;
    for(const auto & component:fit.components)
    {
        conditioning.insert(conditioning.end(),component.local_conditioning.begin(),component.local_conditioning.end());
        if(component.parameter_stability)
            stability.insert(stability.end(),component.parameter_stability->blocks.begin(),component.parameter_stability->blocks.end());
    }
    const auto fresh_target=n::AssessTargetPointState(problem,problem.ParameterLayout(),fit.assembled_state,
        fit.RuntimeConvergence(),conditioning,stability,rhbm_gem::JointEvidenceScope::AssembledGlobal,
        problem.Input().observations.size());
    ASSERT_EQ(fit.target_evidence->checks.size(),fresh_target.checks.size());
    for(std::size_t k=0;k<fresh_target.checks.size();++k)
    {
        EXPECT_EQ(fit.target_evidence->checks[k].name,fresh_target.checks[k].name);
        EXPECT_EQ(fit.target_evidence->checks[k].status,fresh_target.checks[k].status);
        EXPECT_EQ(fit.target_evidence->checks[k].scope,fresh_target.checks[k].scope);
        EXPECT_EQ(fit.target_evidence->checks[k].reason,fresh_target.checks[k].reason);
    }
    ASSERT_EQ(fit.target_evidence->point_atoms.size(),fresh_target.point_atoms.size());
    for(std::size_t k=0;k<fresh_target.point_atoms.size();++k)
    {
        EXPECT_EQ(fit.target_evidence->point_atoms[k].atom,fresh_target.point_atoms[k].atom);
        EXPECT_EQ(fit.target_evidence->point_atoms[k].parameter_validity,fresh_target.point_atoms[k].parameter_validity);
        EXPECT_EQ(fit.target_evidence->point_atoms[k].point_quality,fresh_target.point_atoms[k].point_quality);
    }
    const auto saved=core::CaptureJointAnalysisResult(fit);
    auto json=boost::json::parse(rhbm_gem::joint_result_io::Encode(saved));
    EXPECT_EQ(json.at("schema_version").as_int64(),11);
    json.at("target_evidence").as_object()["column_scales"]=boost::json::array{1.};
    EXPECT_THROW(rhbm_gem::joint_result_io::Decode(boost::json::serialize(json)),std::invalid_argument);
    const auto restored=rhbm_gem::joint_result_io::Decode(rhbm_gem::joint_result_io::Encode(saved));
    EXPECT_TRUE(restored.assembled_state);
    ASSERT_EQ(restored.components.size(),2);
    for(std::size_t k=0;k<2;++k)
    {
        ASSERT_TRUE(restored.components[k].candidate_diagnostic);
        EXPECT_EQ(restored.components[k].candidate_diagnostic->endpoint_verification_status,
            fit.components[k].candidate_diagnostic->endpoint_verification_status);
        EXPECT_TRUE(restored.components[k].state);
    }
    ASSERT_TRUE(restored.target_evidence);
    ASSERT_EQ(restored.target_evidence->point_atoms.size(),fit.target_evidence->point_atoms.size());
    EXPECT_EQ(restored.target_evidence->point_atoms[0].point_quality,fit.target_evidence->point_atoms[0].point_quality);
}

TEST(JointObservableProfileTest, UnobservedTargetHasNoTargetCertificate)
{
    core::JointProblemInput input; input.atom_ids={"1"}; input.support.resize(1);
    input.selection_domain.emplace(); input.selection_domain->target_indices={0};
    const auto fit=core::FitJointComponents(core::JointProblem(input),{.6});
    EXPECT_EQ(fit.TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::Unavailable);
    EXPECT_EQ(fit.components.at(0).TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::Unavailable);
}
