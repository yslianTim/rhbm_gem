#include <gtest/gtest.h>
#include <boost/json.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/data/io/JointAnalysisFileIO.hpp>
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/TargetEvidence.hpp"
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/FittingWorkset.hpp"
#include "core/detail/JointUncertainty.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include "data/io/detail/JointResultJson.hpp"
#include <algorithm>
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
}
TEST(JointObservableProfileTest, FixedLayoutAndSharedRows)
{
    const core::JointProblem problem(Input()); const auto & layout=problem.ParameterLayout();
    EXPECT_EQ(layout.full_atoms,(std::vector<std::size_t>{0}));
    ASSERT_EQ(layout.groups.size(),1); EXPECT_EQ(layout.groups[0].atoms,(std::vector<std::size_t>{1,2}));
    EXPECT_EQ(layout.groups[0].row,23); EXPECT_EQ(layout.informative_rows.size(),23);
    auto selected=Input(); selected.selection_domain->target_indices={0,1};
    EXPECT_EQ(core::JointProblem(selected).ParameterLayout().full_atoms,(std::vector<std::size_t>{0,1}));
    auto legacy=Input(); legacy.selection_domain.reset();
    EXPECT_TRUE(core::JointProblem(legacy).ParameterLayout().groups.empty());
    auto multiple=Input(); multiple.support[2].push_back({22,6.0});
    EXPECT_EQ(core::JointProblem(multiple).ParameterLayout().full_atoms,(std::vector<std::size_t>{0,2}));
    auto separate=Input(); separate.support[2][0].row=22;
    EXPECT_EQ(core::JointProblem(separate).ParameterLayout().groups.size(),2);
}
TEST(JointObservableProfileTest, ExplicitNuisanceAndNonzeroResidualDerivative)
{
    const core::JointProblem problem(Input(true)); const auto & data=core::JointProblemAccess::Get(problem);
    const auto & layout=problem.ParameterLayout();
    const auto domain=n::ProfileDomain(data.domain,layout);
    EXPECT_EQ(domain.atoms.Storage(),data.domain.atoms.Storage());
    ASSERT_EQ(domain.atoms[0].size(),23);
    const auto context=n::ProfileContext(data.context,layout,data.y.size());
    EXPECT_EQ(context.scale,data.context.scale); EXPECT_EQ(context.rank.rows,24);
    const n::Vector y=data.y.head(23),eta=n::Vector::Constant(1,std::log(.7));
    const auto evaluation=n::EvaluateProfile(domain,y,eta,false,&context); ASSERT_TRUE(evaluation.valid);
    n::Matrix design=n::Matrix::Zero(24,3);
    for(std::size_t r=0;r<24;++r)
    {
        const auto basis=n::EvaluateKernel(problem.Input().support[0][r].squared_distance,.7,2.5);
        design(static_cast<Eigen::Index>(r),0)=basis.gaussian; design(static_cast<Eigen::Index>(r),1)=basis.charge;
    }
    design(23,2)=1;
    const n::Vector explicit_solution=design.colPivHouseholderQr().solve(data.y);
    EXPECT_GT(explicit_solution(0),0);
    EXPECT_NEAR((explicit_solution.head(2)-evaluation.beta).norm(),0,1e-11);
    const n::Vector residual=design*explicit_solution-data.y;
    EXPECT_NEAR(residual.squaredNorm(),evaluation.residual.squaredNorm(),1e-12);
    const auto differential=n::PrepareDerivative(evaluation,context.scale,&context);
    ASSERT_TRUE(differential.valid); n::Matrix projected,jacobian;
    differential.Rows(0,23,projected,jacobian);
    const double step=1e-5;
    const auto plus=n::EvaluateProfile(domain,y,eta.array()+step,false,&context);
    const auto minus=n::EvaluateProfile(domain,y,eta.array()-step,false,&context);
    ASSERT_TRUE(plus.valid && minus.valid);
    EXPECT_GT(evaluation.residual.norm(),1e-4);
    EXPECT_LT(((plus.residual-minus.residual)/(2*step*context.scale)-jacobian.col(0)).norm(),1e-8);
    const double gradient=(plus.residual.squaredNorm()-minus.residual.squaredNorm())/(4*step*context.scale*context.scale);
    EXPECT_NEAR(gradient,evaluation.gradient(0),1e-9);
}
TEST(JointObservableProfileTest, ReconstructionAndUnavailableHaloParameters)
{
    const core::JointProblem problem(Input()); const auto fit=core::FitJointComponents(problem,{.65,missing_width,missing_width});
    ASSERT_TRUE(fit.assembled_state); ASSERT_TRUE(fit.prediction);
    EXPECT_EQ(fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(fit.components[0].atoms.size(),3); EXPECT_EQ(fit.components[0].state->b.size(),1);
    ASSERT_TRUE(fit.components[0].candidate_diagnostic);
    EXPECT_TRUE(fit.components[0].candidate_diagnostic->candidate_available);
    EXPECT_EQ(fit.components[0].candidate_diagnostic->candidate_eta,fit.components[0].state->log_b);
    EXPECT_EQ(fit.components[0].candidate_diagnostic->candidate_ac,fit.components[0].state->ac);
    auto invalid_state=*fit.assembled_state; invalid_state.ac[0]+=.01;
    const auto invalid=n::AssessTargetState(problem,problem.ParameterLayout(),invalid_state,fit.evidence,
        rhbm_gem::JointEvidenceScope::AssembledGlobal,problem.Input().observations.size());
    EXPECT_EQ(invalid.atoms[0].reason,"invalid-saved-state");
    auto inconsistent_state=*fit.assembled_state; inconsistent_state.b[0]*=1.01;
    const auto inconsistent=n::AssessTargetState(problem,problem.ParameterLayout(),inconsistent_state,fit.evidence,
        rhbm_gem::JointEvidenceScope::AssembledGlobal,problem.Input().observations.size());
    EXPECT_EQ(inconsistent.atoms[0].reason,"inconsistent-saved-widths");
    EXPECT_EQ(fit.components[0].candidate_diagnostic->certification_status,rhbm_gem::JointCheckStatus::Passed);
    const auto & problem_data=core::JointProblemAccess::Get(problem);
    const auto & component_layout=*fit.components[0].layout;
    const auto candidate_domain=n::ProfileDomain(problem_data.domain,component_layout);
    const auto candidate_context=n::ProfileContext(problem_data.context,component_layout,problem_data.domain.rows);
    n::Vector candidate_y(static_cast<Eigen::Index>(component_layout.informative_rows.size()));
    for(std::size_t k=0;k<component_layout.informative_rows.size();++k)
        candidate_y(static_cast<Eigen::Index>(k))=problem_data.input->observations[component_layout.informative_rows[k]];
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
    EXPECT_NEAR(fit.assembled_state->nuisance_amplitudes.at(0),.7,1e-10);
    for(std::size_t r=0;r<24;++r) EXPECT_NEAR(fit.prediction->at(r),problem.Input().observations[r],1e-10);
    const auto saved=core::CaptureJointAnalysisResult(fit);
    const auto stages=rhbm_gem::data_internal::BuildJointStageEstimates(saved,"test");
    EXPECT_TRUE(stages.at(1).point); EXPECT_FALSE(stages.at(2).point); EXPECT_FALSE(stages.at(3).point);
    EXPECT_EQ(stages.at(2).reason,"observable-contribution-only");
    const auto decoded=rhbm_gem::joint_result_io::Decode(rhbm_gem::joint_result_io::Encode(saved));
    EXPECT_EQ(decoded.assembled_state->nuisance_amplitudes,saved.assembled_state->nuisance_amplitudes);
    EXPECT_EQ(decoded.layout->groups[0].atoms,saved.layout->groups[0].atoms);
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
    const core::JointProblem problem(Input());
    const std::vector<double> initial{.65,missing_width,missing_width};
    std::vector<n::JointProgressEvent> events;
    const auto batch=n::OrchestrateFixedNeighborComponents(problem,initial,{},
        [&](const auto & event) {events.push_back(event);});
    ASSERT_EQ(batch.components.size(),1);
    ASSERT_TRUE(batch.components.front().state);
    ASSERT_TRUE(batch.components.front().layout);
    EXPECT_FALSE(batch.single_component_assessment);
    EXPECT_EQ(events.back().phase,n::JointProgressPhase::ComponentCompleted);
    EXPECT_TRUE(std::none_of(events.begin(),events.end(),[](const auto & event) {
        return event.phase==n::JointProgressPhase::AssemblyStarted ||
            event.phase==n::JointProgressPhase::AssemblyCompleted;
    }));

    const auto fit=n::FitObservableComponents(problem,initial);
    ASSERT_TRUE(fit.assembled_state); ASSERT_TRUE(fit.prediction); ASSERT_TRUE(fit.objective);
    ASSERT_TRUE(fit.components.front().state);
    EXPECT_EQ(fit.components.front().state->ac,batch.components.front().state->ac);
    EXPECT_EQ(fit.components.front().state->log_b,batch.components.front().state->log_b);
    ASSERT_EQ(fit.assembled_state->nuisance_amplitudes.size(),1);
    EXPECT_NEAR(fit.assembled_state->nuisance_amplitudes.front(),.7,1e-10);
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
    const core::JointProblem problem(Input());
    n::FixedNeighborSearchPolicy policy; policy.maximum_sweeps=1;
    const auto fit=n::FitFixedNeighborComponents(problem,{.65,missing_width,missing_width},policy);
    ASSERT_EQ(fit.components.size(),1u);
    const auto & component=fit.components.front();
    ASSERT_TRUE(component.candidate_diagnostic);
    EXPECT_TRUE(component.candidate_diagnostic->candidate_available);
    EXPECT_FALSE(component.search_completed);
    EXPECT_EQ(component.stop_reason,"block-sweep-budget");
    EXPECT_EQ(component.candidate_diagnostic->candidate_eta.size(),component.layout->full_atoms.size());
    EXPECT_NE(component.candidate_diagnostic->certification_status,rhbm_gem::JointCheckStatus::NotRun);
}
TEST(JointObservableProfileTest, RemainingUnobservedTargetStillFails)
{
    auto input=Input(); input.support[0]={{23,1}};
    const core::JointProblem problem(input);
    const auto fit=core::FitJointComponents(problem,{.6,missing_width,missing_width});
    EXPECT_FALSE(fit.assembled_state); EXPECT_NE(fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
}
TEST(JointObservableProfileTest, RemainingDependentTargetIsNotReducedAgain)
{
    auto input=Input(); input.support[0]={{0,1},{1,1},{23,1}};
    const core::JointProblem problem(input);
    EXPECT_EQ(problem.ParameterLayout().full_atoms,(std::vector<std::size_t>{0}));
    const auto fit=core::FitJointComponents(problem,{.6,missing_width,missing_width});
    ASSERT_EQ(fit.components.size(),1); EXPECT_FALSE(fit.components[0].state);
    EXPECT_FALSE(fit.assembled_state);
    EXPECT_NE(fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
}
TEST(JointObservableProfileTest, RemainingTwoRowInteriorWidthIsUnidentified)
{
    auto input=Input(); input.atom_ids.push_back("4"); input.support.push_back({{0,.4},{1,1.2}});
    input.selection_domain->target_indices.push_back(3);
    for(const auto & support:input.support[3])
    {
        const auto basis=n::EvaluateKernel(support.squared_distance,.65,2.5);
        input.observations[support.row]+=basis.gaussian+.2*basis.charge;
    }
    const core::JointProblem problem(input); const auto & data=core::JointProblemAccess::Get(problem);
    const auto & layout=problem.ParameterLayout();
    EXPECT_EQ(layout.full_atoms,(std::vector<std::size_t>{0,3}));
    const auto context=n::ProfileContext(data.context,layout,data.y.size());
    const n::Vector eta=Eigen::Vector2d(.6,.65).array().log();
    const auto assessment=n::AssessProfile(n::ProfileDomain(data.domain,layout),data.y.head(23),eta,context);
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
    initialization.atoms[4].reason="not-required-observable-contribution";
    core::detail::ApplyJointSeedFallback(initialization);
    EXPECT_TRUE(initialization.valid); EXPECT_EQ(initialization.b[2],2); EXPECT_EQ(initialization.b[3],2);
    EXPECT_EQ(initialization.atoms[2].reason,"initialization-exception: test");
    EXPECT_FALSE(initialization.atoms[2].original_b); EXPECT_EQ(initialization.atoms[3].original_b,-1);
    EXPECT_EQ(initialization.atoms[2].seed_source,"median-fallback"); EXPECT_EQ(initialization.atoms[2].donor_count,2);
    EXPECT_TRUE(std::isnan(initialization.b[4]));
    rhbm_gem::JointInitialization empty; empty.b={missing_width}; empty.atoms.resize(1); empty.atoms[0].reason="invalid-width";
    core::detail::ApplyJointSeedFallback(empty); EXPECT_FALSE(empty.valid); EXPECT_TRUE(std::isnan(empty.b[0]));
}
TEST(JointObservableProfileTest, ProfileCovarianceMatchesAugmentedJacobian)
{
    const core::JointProblem problem(Input(true)); const auto fit=core::FitJointComponents(problem,{.65,missing_width,missing_width});
    ASSERT_TRUE(fit.assembled_state); ASSERT_EQ(fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    const auto result=core::CaptureJointAnalysisResult(fit);
    const auto uncertainties=core::detail::ComputeJointUncertainty(problem,result);
    const auto & uncertainty=uncertainties.at(1); ASSERT_TRUE(uncertainty.covariance) << uncertainty.reason;
    EXPECT_EQ(uncertainty.degrees_of_freedom,20);
    const auto & state=*fit.assembled_state; n::Matrix jacobian=n::Matrix::Zero(24,4);
    for(std::size_t r=0;r<24;++r)
    {
        const auto basis=n::EvaluateKernel(problem.Input().support[0][r].squared_distance,state.b[0],2.5);
        jacobian(static_cast<Eigen::Index>(r),0)=basis.gaussian; jacobian(static_cast<Eigen::Index>(r),1)=basis.charge;
        jacobian(static_cast<Eigen::Index>(r),2)=state.ac[0]*basis.gaussian_log_width+state.ac[1]*basis.charge_log_width;
    }
    jacobian(23,3)=1;
    const Eigen::JacobiSVD<n::Matrix> svd(jacobian,Eigen::ComputeThinV);
    const n::Matrix factor=svd.matrixV()*svd.singularValues().cwiseInverse().asDiagonal();
    const n::Matrix covariance=*uncertainty.residual_variance*factor*factor.transpose();
    EXPECT_LT(((*uncertainty.covariance-covariance.topLeftCorner<3,3>()).norm()/covariance.topLeftCorner<3,3>().norm()),1e-10);
    // A saved convergence flag cannot bypass the regular interior covariance gate.
    auto boundary=result; boundary.components[0].state->ac[0]=0;
    const auto unavailable=core::detail::ComputeJointUncertainty(problem,boundary);
    EXPECT_FALSE(unavailable.at(1).covariance);
    EXPECT_EQ(unavailable.at(1).reason,"active-amplitude-boundary");
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
    const auto reduced=core::CaptureJointAnalysisResult(core::FitJointComponents(core::JointProblem(Input()),{.6,missing_width,missing_width}));
    auto bad=boost::json::parse(rhbm_gem::joint_result_io::Encode(reduced)).as_object();
    bad.at("layout").as_object().at("groups").as_array()[0].as_object()["atoms"]=boost::json::array{0,1,2};
    EXPECT_THROW(rhbm_gem::joint_result_io::Decode(boost::json::serialize(bad)),std::invalid_argument);
}
TEST(JointObservableProfileTest, IndependentAnalyticComponentSurvivesInvalidTargetSeed)
{
    auto input=Input(); input.row_ids.push_back("24"); input.observations.push_back(5); input.support[2]={{24,2}};
    const auto fit=core::FitJointComponents(core::JointProblem(input),{missing_width,missing_width,missing_width});
    ASSERT_EQ(fit.components.size(),2); EXPECT_FALSE(fit.assembled_state);
    EXPECT_FALSE(fit.components[0].state); ASSERT_TRUE(fit.components[1].state);
    EXPECT_EQ(fit.components[1].state->nuisance_amplitudes,(std::vector<double>{5}));
    EXPECT_TRUE(fit.available_row_mask[24]); EXPECT_FALSE(fit.available_row_mask[23]);
}
TEST(JointObservableProfileTest, PersistenceExportPeelingAndSkippedInitialization)
{
    using namespace rhbm_gem;
    const core::JointProblem problem(Input());
    const auto result=core::CaptureJointAnalysisResult(core::FitJointComponents(problem,{.6,missing_width,missing_width}));
    std::vector<std::unique_ptr<AtomObject>> atoms;
    for(int id=1;id<=3;++id)
    {auto atom=std::make_unique<AtomObject>(); atom->SetSerialID(id); atom->SetElement(Element::CARBON); atom->SetPosition(0,0,0); atoms.push_back(std::move(atom));}
    ModelObject model(std::move(atoms)); model.SelectAtoms([](const auto & a){return a.GetSerialID()==1;});
    auto values=std::make_unique<double[]>(24*4*4);
    for(std::size_t r=0;r<24;++r) values[r]=problem.Input().observations[r];
    MapObject map({24,4,4},{1,1,1},{0,0,0},std::move(values));
    auto editor=model.EditAnalysis(); editor.ApplyJointResult(result,"observable-test");
    const SamplingPointList points{{23,{23,0,0},false},{1,{1,1,1},false}};
    const auto raw=core::detail::BuildLocalPotentialSampleList(map,points);
    editor.SetAtomLocalRawSamplingEntries(*model.FindAtomPtr(1),raw);
    const auto peeled=core::detail::BuildPostFitPeelingSamples(map,model,problem,std::nullopt,&result);
    ASSERT_TRUE(peeled.at(1).samples[0].response);
    EXPECT_NEAR(*peeled.at(1).samples[0].response,raw[0].response-.7,1e-10);
    EXPECT_FALSE(peeled.at(1).samples[1].response);
    EXPECT_EQ(peeled.at(1).samples[1].reason,"outside-joint-domain");
    editor.SetJointResult(result);
    const command_test::ScopedTempDir directory{"observable_roundtrip"};
    DataRepository repository(directory.path()/"saved.sqlite"); repository.SaveModel(model,"test");
    const auto loaded=repository.LoadModel("test");
    ASSERT_TRUE(loaded->GetAnalysisView().GetJointResult());
    EXPECT_EQ(joint_result_io::Encode(*loaded->GetAnalysisView().GetJointResult()),joint_result_io::Encode(result));
    EXPECT_FALSE(AtomLocalPotentialView::For(*loaded->FindAtomPtr(2)).GetStageEstimate(FittingStage::Second).point);
    WriteJointAnalysisResult(result,directory.path()/"saved.json",directory.path()/"saved.csv");
    std::ifstream groups(directory.path()/"saved.contributions.csv");
    const std::string group_text(std::istreambuf_iterator<char>(groups),{});
    EXPECT_EQ(std::count(group_text.begin(),group_text.end(),'\n'),2); EXPECT_NE(group_text.find("2;3"),std::string::npos);
    const auto last=group_text.find_last_of(','),before=group_text.find_last_of(',',last-1);
    EXPECT_DOUBLE_EQ(std::stod(group_text.substr(before+1,last-before-1)),result.assembled_state->nuisance_amplitudes[0]);
    std::ifstream csv(directory.path()/"saved.csv"); std::string line;
    std::getline(csv,line); std::getline(csv,line); // Header and FullABC target.
    for(const auto id:{"2","3"})
    {
        std::getline(csv,line);
        EXPECT_TRUE(line.starts_with("\""+std::string(id)+"\",\""+result.components[0].id+"\",,,,0,"));
    }
    const auto workset=core::detail::MakeJointFittingWorkset(model,problem);
    std::vector<int> visited;
    core::detail::FirstStageObserverForTesting()=[&](int id,std::string_view){visited.push_back(id);};
    core::FitOptions options; options.quiet_mode=true; options.thread_size=1;
    const auto first_stage_atoms=core::detail::CollectFirstStageAtoms(workset);
    core::RunPotentialSamplingWorkflow(map,model,first_stage_atoms,options.sampling_method,1);
    const auto initialization=core::detail::RunJointFirstStageInitializationFromPreparedSamples(model,workset,options);
    core::detail::FirstStageObserverForTesting()={};
    EXPECT_TRUE(std::all_of(visited.begin(),visited.end(),[](int id){return id==1;}));
    EXPECT_EQ(initialization.atoms[1].reason,"not-required-observable-contribution");
    EXPECT_EQ(initialization.atoms[2].reason,"not-required-observable-contribution");
    EXPECT_FALSE(initialization.atoms[1].original_b); EXPECT_FALSE(initialization.atoms[1].used_b);
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
    const auto fit=core::FitJointComponents(problem,{.6,missing_width,missing_width,.65});
    ASSERT_FALSE(fit.assembled_state);
    ASSERT_TRUE(fit.components[0].candidate_diagnostic);
    EXPECT_TRUE(fit.components[0].candidate_diagnostic->candidate_available);
    EXPECT_EQ(fit.components[0].candidate_diagnostic->endpoint_verification_status,
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
    EXPECT_EQ(value.at("schema_version").as_int64(),8);
    auto downgraded=value;
    downgraded.as_object()["schema_version"]=7;
    EXPECT_THROW(rhbm_gem::joint_result_io::Decode(boost::json::serialize(downgraded)),std::invalid_argument);
    const auto restored=rhbm_gem::joint_result_io::Decode(text);
    EXPECT_EQ(restored.target_runtime_convergence,saved.target_runtime_convergence);
    ASSERT_TRUE(restored.components[0].candidate_diagnostic);
    EXPECT_EQ(restored.components[0].candidate_diagnostic->endpoint_verification_status,
        rhbm_gem::JointCheckStatus::Passed);
    EXPECT_FALSE(restored.components[0].state);
    value.as_object()["schema_version"]=4;
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
    const auto fit=core::FitJointComponents(problem,{.6,missing_width,missing_width,.65});
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
        fit.components[0].evidence,rhbm_gem::JointEvidenceScope::ComponentLocal,fit.components[0].rows.size());
    const auto target=std::find_if(evidence.checks.begin(),evidence.checks.end(),[](const auto & check) {
        return check.name=="target-identifiability";
    });
    ASSERT_NE(target,evidence.checks.end());
    EXPECT_EQ(target->status,rhbm_gem::JointCheckStatus::Failed);
    EXPECT_EQ(evidence.atoms[3].status,rhbm_gem::JointCheckStatus::Failed);
    EXPECT_EQ(evidence.atoms[3].reason,"nonunique-parameter-directions");
}
TEST(JointObservableProfileTest, UnqualifiedHaloRepresentativesDoNotProduceTargetCovariance)
{
    const core::JointProblem problem(TwoRowHalo(false,true));
    const auto first=core::FitJointComponents(problem,{.6,missing_width,missing_width,.65});
    const auto second=core::FitJointComponents(problem,{.6,missing_width,missing_width,.7});
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
    const auto ca=core::detail::ComputeJointUncertainty(problem,core::CaptureJointAnalysisResult(first));
    const auto cb=core::detail::ComputeJointUncertainty(problem,core::CaptureJointAnalysisResult(second));
    EXPECT_FALSE(ca.at(1).covariance); EXPECT_FALSE(cb.at(1).covariance);
}
TEST(JointObservableProfileTest, ActiveBoundaryDoesNotAcquireTargetCertificate)
{
    const core::JointProblem problem(Input());
    auto fit=core::FitJointComponents(problem,{.6,missing_width,missing_width});
    ASSERT_TRUE(fit.assembled_state);
    auto state=*fit.assembled_state; state.ac[0]=0;
    const auto evidence=n::AssessTargetState(problem,problem.ParameterLayout(),state,fit.evidence,
        rhbm_gem::JointEvidenceScope::AssembledGlobal,problem.Input().observations.size());
    EXPECT_EQ(n::TargetConvergence(evidence,true),rhbm_gem::JointCheckStatus::Unavailable);
    auto saved=core::CaptureJointAnalysisResult(fit);
    saved.components[0].state->ac[0]=0;
    const auto uncertainty=core::detail::ComputeJointUncertainty(problem,saved);
    EXPECT_FALSE(uncertainty.at(1).covariance);
    EXPECT_EQ(uncertainty.at(1).reason,"active-amplitude-boundary");
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
}
TEST(JointObservableProfileTest, TargetEvidenceRequiresSelectionMetadata)
{
    auto input=Input(); input.atom_ids.resize(1); input.support.resize(1); input.selection_domain.reset();
    const auto fit=core::FitJointComponents(core::JointProblem(input),{.6});
    EXPECT_EQ(fit.TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::NotRun);
    EXPECT_FALSE(fit.target_evidence);
}

TEST(JointObservableProfileTest, ExactSingleComponentTargetEvidenceIsReused)
{
    const core::JointProblem problem(Input());
    n::AssessmentWorkForTesting()={};
    const auto fit=core::FitJointComponents(problem,{.6,missing_width,missing_width});
    ASSERT_EQ(fit.components.size(),1);
    ASSERT_TRUE(fit.components[0].target_evidence);
    ASSERT_TRUE(fit.target_evidence);
    EXPECT_EQ(fit.target_evidence->contract,"target-estimability-v2");
    auto geometry_calls=[&] {
        const auto & stages=n::AssessmentWorkForTesting().stages;
        const auto found=std::find_if(stages.begin(),stages.end(),[](const auto & stage) {
            return stage.name=="target-geometry";
        });
        return found==stages.end() ? 0 : found->calls;
    };
    EXPECT_EQ(geometry_calls(),1);
    for(const auto & check:fit.target_evidence->checks)
        EXPECT_NE(check.name,"target-local-correction");
    for(const auto & rank:fit.target_evidence->ranks)
        EXPECT_NE(rank.name,"undamped-full-profile-step");

    const auto & data=core::JointProblemAccess::Get(problem);
    const auto fresh=n::AssessTargetState(problem,problem.ParameterLayout(),fit.assembled_state,fit.evidence,
        rhbm_gem::JointEvidenceScope::AssembledGlobal,problem.Input().observations.size());
    EXPECT_EQ(geometry_calls(),2);
    ASSERT_EQ(fit.target_evidence->checks.size(),fresh.checks.size());
    for(std::size_t i=0;i<fresh.checks.size();++i)
    {
        EXPECT_EQ(fit.target_evidence->checks[i].name,fresh.checks[i].name);
        EXPECT_EQ(fit.target_evidence->checks[i].status,fresh.checks[i].status);
        EXPECT_EQ(fit.target_evidence->checks[i].scope,fresh.checks[i].scope);
        EXPECT_EQ(fit.target_evidence->checks[i].reason,fresh.checks[i].reason);
    }
    ASSERT_EQ(fit.target_evidence->atoms.size(),fresh.atoms.size());
    for(std::size_t i=0;i<fresh.atoms.size();++i)
    {
        EXPECT_EQ(fit.target_evidence->atoms[i].atom,fresh.atoms[i].atom);
        EXPECT_EQ(fit.target_evidence->atoms[i].status,fresh.atoms[i].status);
        EXPECT_EQ(fit.target_evidence->atoms[i].null_space_leakage,fresh.atoms[i].null_space_leakage);
    }
    EXPECT_EQ(data.context.scale,fit.observation_scale);
}

TEST(JointObservableProfileTest, TargetGeometryNumericalOracleCapturesOldContract)
{
    using Status=rhbm_gem::JointCheckStatus;
    const core::JointProblem problem(Input(true));
    const auto fit=core::FitJointComponents(problem,{.6,missing_width,missing_width});
    ASSERT_EQ(fit.components.size(),1u);
    ASSERT_TRUE(fit.components.front().state);
    ASSERT_TRUE(fit.components.front().target_evidence);
    const auto & component=fit.components.front();
    const auto & layout=*component.layout;
    const auto & data=core::JointProblemAccess::Get(problem);
    const auto domain=n::ProfileDomain(data.domain,layout);
    const auto context=n::ProfileContext(data.context,layout,data.y.size());
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
    const auto geometry=n::BuildTargetGeometry(domain,y,*component.state,context);
    ASSERT_TRUE(geometry.valid)<<geometry.reason;
    EXPECT_FALSE(geometry.threshold_sensitive);
    EXPECT_EQ(geometry.svd.rank,3);
    EXPECT_GT(geometry.svd.threshold,0);
    EXPECT_EQ(geometry.scales.size(),3);
    EXPECT_TRUE(geometry.scales.isApprox(scales,1e-12));
    ASSERT_EQ(geometry.svd.singular_values.size(),oracle.singularValues().size());
    EXPECT_TRUE(geometry.svd.singular_values.isApprox(oracle.singularValues(),1e-11));
    const auto design=std::find_if(component.ranks.begin(),component.ranks.end(),[](const auto & rank) {
        return rank.name=="design";
    });
    const auto widths=std::find_if(component.ranks.begin(),component.ranks.end(),[](const auto & rank) {
        return rank.name=="projected-width";
    });
    ASSERT_NE(design,component.ranks.end()); ASSERT_NE(widths,component.ranks.end());
    EXPECT_EQ(design->rank,2u); EXPECT_EQ(widths->rank,1u);
    ASSERT_EQ(component.target_evidence->ranks.size(),1u);
    EXPECT_EQ(component.target_evidence->ranks.front().name,"column-scaled-observation-jacobian");
    EXPECT_EQ(component.target_evidence->ranks.front().rank,3u);
    ASSERT_EQ(component.target_evidence->atoms.size(),1u);
    EXPECT_EQ(component.target_evidence->atoms.front().status,Status::Passed);
    ASSERT_TRUE(component.target_evidence->atoms.front().null_space_leakage);
    EXPECT_LE(*component.target_evidence->atoms.front().null_space_leakage,1e-10);
    EXPECT_EQ(component.TargetRuntimeConvergence(),Status::Passed);
    const auto uncertainty=core::detail::ComputeJointUncertainty(problem,core::CaptureJointAnalysisResult(fit));
    EXPECT_EQ(uncertainty.at(1).status,rhbm_gem::EvidenceStatus::Available);
    EXPECT_TRUE(uncertainty.at(1).covariance);

    auto boundary=*component.state;
    boundary.ac[0]=0;
    const auto boundary_geometry=n::BuildTargetGeometry(domain,y,boundary,context);
    EXPECT_FALSE(boundary_geometry.valid);
    EXPECT_EQ(boundary_geometry.reason,"active-amplitude-boundary");
}

TEST(JointObservableProfileTest, TargetGeometryOracleMapsFullAbcHaloCoordinates)
{
    auto input=Input();
    input.atom_ids.push_back("4"); input.support.emplace_back();
    for(std::size_t row=0;row<23;++row)
    {
        const double square=std::pow(.1+.09*static_cast<double>(row),2);
        input.support.back().push_back({row,square});
        const auto basis=n::EvaluateKernel(square,.9,2.5);
        input.observations[row]+=.7*basis.gaussian+.15*basis.charge;
    }
    const core::JointProblem problem(std::move(input));
    EXPECT_EQ(problem.ParameterLayout().full_atoms,(std::vector<std::size_t>{0,3}));
    const auto fit=core::FitJointComponents(problem,{.6,missing_width,missing_width,.9});
    ASSERT_TRUE(fit.assembled_state);
    ASSERT_TRUE(fit.target_evidence);
    ASSERT_EQ(fit.target_evidence->atoms.size(),4u);
    const auto target=std::find_if(fit.target_evidence->atoms.begin(),fit.target_evidence->atoms.end(),[](const auto & atom) {
        return atom.atom==0;
    });
    const auto halo=std::find_if(fit.target_evidence->atoms.begin(),fit.target_evidence->atoms.end(),[](const auto & atom) {
        return atom.atom==3;
    });
    ASSERT_NE(target,fit.target_evidence->atoms.end()); ASSERT_NE(halo,fit.target_evidence->atoms.end());
    EXPECT_EQ(target->status,rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(halo->status,rhbm_gem::JointCheckStatus::Passed);
    EXPECT_TRUE(std::binary_search(problem.Input().selection_domain->target_indices.begin(),
        problem.Input().selection_domain->target_indices.end(),0u));
    EXPECT_FALSE(std::binary_search(problem.Input().selection_domain->target_indices.begin(),
        problem.Input().selection_domain->target_indices.end(),3u));
}

TEST(JointObservableProfileTest, DerivedTargetEvidenceRequiresMatchingQualifiedRankEvidence)
{
    using Status=rhbm_gem::JointCheckStatus;
    using Scope=rhbm_gem::JointEvidenceScope;
    const core::JointProblem problem(Input(true));
    const auto fit=core::FitJointComponents(problem,{.6,missing_width,missing_width});
    ASSERT_TRUE(fit.assembled_state);
    const auto evidence=n::DeriveTargetIdentifiability(problem,problem.ParameterLayout(),fit.assembled_state,
        fit.evidence,fit.ranks,Scope::AssembledGlobal,problem.Input().observations.size());
    ASSERT_EQ(evidence.contract,"target-estimability-v3");
    ASSERT_TRUE(evidence.column_scales.empty());
    ASSERT_EQ(evidence.ranks.size(),2u);
    EXPECT_EQ(evidence.ranks[0].name,"design");
    EXPECT_EQ(evidence.ranks[1].name,"projected-width");
    ASSERT_EQ(evidence.atoms.size(),3u);
    EXPECT_EQ(evidence.atoms[0].status,Status::Passed);
    EXPECT_EQ(evidence.atoms[0].reason,"derived-from-component-identifiability");
    EXPECT_FALSE(evidence.atoms[0].null_space_leakage);
    EXPECT_EQ(evidence.atoms[1].status,Status::Unavailable);
    EXPECT_EQ(evidence.atoms[1].reason,"observable-contribution-only");
    EXPECT_EQ(evidence.atoms[2].reason,"observable-contribution-only");
    EXPECT_EQ(n::TargetConvergence(evidence,true),Status::Passed);

    const auto derive_with=[&](const std::optional<rhbm_gem::JointState> & state,
        const std::vector<rhbm_gem::JointCheck> & checks,const std::vector<rhbm_gem::JointRankEvidence> & ranks,
        const core::JointParameterLayout & layout,std::size_t rows) {
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
    EXPECT_EQ(status_reason(derive({},fit.evidence,fit.ranks)),"missing-component-state");
    auto no_design=fit.ranks;
    no_design.erase(std::remove_if(no_design.begin(),no_design.end(),[](const auto & rank) {
        return rank.name=="design";
    }),no_design.end());
    EXPECT_EQ(status_reason(derive(fit.assembled_state,fit.evidence,no_design)),"missing-rank-evidence");
    auto no_width=fit.ranks;
    no_width.erase(std::remove_if(no_width.begin(),no_width.end(),[](const auto & rank) {
        return rank.name=="projected-width";
    }),no_width.end());
    EXPECT_EQ(status_reason(derive(fit.assembled_state,fit.evidence,no_width)),"missing-rank-evidence");
    auto deficient=fit.ranks;
    std::find_if(deficient.begin(),deficient.end(),[](const auto & rank) {
        return rank.name=="projected-width";
    })->rank=0;
    EXPECT_EQ(status_reason(derive(fit.assembled_state,fit.evidence,deficient)),"component-rank-deficient");
    auto near_threshold=fit.ranks;
    auto design=std::find_if(near_threshold.begin(),near_threshold.end(),[](const auto & rank) {
        return rank.name=="design";
    });
    ASSERT_NE(design,near_threshold.end()); ASSERT_GE(design->singular_values.size(),2u);
    design->singular_values.back()=std::nextafter(design->threshold,std::numeric_limits<double>::infinity());
    EXPECT_EQ(status_reason(derive(fit.assembled_state,fit.evidence,near_threshold)),"rank-threshold-sensitive");
    auto changed_state=*fit.assembled_state; changed_state.ac[0]+=.01;
    EXPECT_EQ(status_reason(derive(changed_state,fit.evidence,fit.ranks)),"rank-context-mismatch");
    auto wrong_layout=problem.ParameterLayout(); wrong_layout.informative_rows.pop_back();
    EXPECT_EQ(status_reason(derive_with(fit.assembled_state,fit.evidence,fit.ranks,wrong_layout,
        problem.Input().observations.size())),"wrong-parameter-layout");
    auto nonfinite=*fit.assembled_state; nonfinite.log_b[0]=std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(status_reason(derive(nonfinite,fit.evidence,fit.ranks)),"nonfinite-parameter-state");
    auto boundary=*fit.assembled_state; boundary.ac[0]=0;
    EXPECT_EQ(status_reason(derive(boundary,fit.evidence,fit.ranks)),"active-amplitude-boundary");
    auto missing_check=fit.evidence;
    missing_check.erase(std::remove_if(missing_check.begin(),missing_check.end(),[](const auto & check) {
        return check.name=="independent-replay";
    }),missing_check.end());
    EXPECT_EQ(status_reason(derive(fit.assembled_state,missing_check,fit.ranks)),"missing-component-qualification");
    EXPECT_EQ(status_reason(derive_with(fit.assembled_state,fit.evidence,fit.ranks,
        problem.ParameterLayout(),problem.Input().observations.size()-1)),"observation-domain-mismatch");

    auto other_input=Input(true); other_input.observations[0]+=1e-6;
    const core::JointProblem other_problem(std::move(other_input));
    const auto other=n::DeriveTargetIdentifiability(other_problem,other_problem.ParameterLayout(),fit.assembled_state,
        fit.evidence,fit.ranks,Scope::AssembledGlobal,other_problem.Input().observations.size());
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
    EXPECT_EQ(estimate.shape.component_atoms,2u);
    EXPECT_EQ(estimate.shape.component_target_atoms,1u);
    EXPECT_EQ(estimate.shape.component_halo_atoms,1u);
    EXPECT_EQ(estimate.shape.component_design_columns,4u);
    EXPECT_EQ(estimate.shape.component_width_columns,2u);
    EXPECT_EQ(estimate.shape.component_target_columns,3u);
    EXPECT_EQ(estimate.shape.component_jacobian_columns,6u);
    ASSERT_EQ(estimate.attributed_stages.size(),1u);
    EXPECT_EQ(estimate.attributed_stages[0].factor_rows,6u);
    EXPECT_EQ(estimate.attributed_stages[0].stage,n::MemoryStage::TargetGeometry);
}

TEST(JointObservableProfileTest, LargeTargetIdentifiabilityPreflightReturnsUnavailableEvidence)
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
    n::AssessmentWorkForTesting()={};
    const auto evidence=n::AssessTargetState(problem,problem.ParameterLayout(),state,{},
        rhbm_gem::JointEvidenceScope::AssembledGlobal,data.input->observations.size());
    const auto identifiability=std::find_if(evidence.checks.begin(),evidence.checks.end(),[](const auto & check) {
        return check.name=="target-identifiability";
    });
    ASSERT_NE(identifiability,evidence.checks.end());
    EXPECT_EQ(identifiability->status,rhbm_gem::JointCheckStatus::Unavailable);
    EXPECT_EQ(identifiability->reason,"resource-limit:target-identifiability");
    EXPECT_EQ(evidence.atoms.size(),atoms);
    EXPECT_TRUE(std::all_of(evidence.atoms.begin(),evidence.atoms.end(),[](const auto & atom) {
        return atom.status==rhbm_gem::JointCheckStatus::Unavailable;
    }));
    EXPECT_TRUE(evidence.ranks.empty());
    EXPECT_EQ(n::TargetConvergence(evidence,true),rhbm_gem::JointCheckStatus::Unavailable);
    const auto & stages=n::AssessmentWorkForTesting().stages;
    const auto preflight=std::find_if(stages.begin(),stages.end(),[](const auto & stage) {
        return stage.name=="target-memory-preflight";
    });
    ASSERT_NE(preflight,stages.end());
    EXPECT_EQ(preflight->memory_status,n::MemoryStatus::Reject);
    EXPECT_EQ(preflight->calls,1); EXPECT_EQ(preflight->completed_calls,1);
    EXPECT_GT(preflight->estimated_peak_bytes,preflight->budget_bytes);
    ASSERT_EQ(preflight->attributed_memory_stages.size(),1u);
    EXPECT_EQ(preflight->attributed_memory_stages[0].stage,n::MemoryStage::TargetGeometry);
    EXPECT_EQ(std::find_if(stages.begin(),stages.end(),[](const auto & stage) {
        return stage.name=="target-geometry";
    }),stages.end());
    EXPECT_EQ(std::find_if(stages.begin(),stages.end(),[](const auto & stage) {
        return stage.name=="target-local-correction";
    }),stages.end());
}

TEST(JointObservableProfileTest, TargetEvidenceIsUnavailableAtTheRankThreshold)
{
    const core::JointProblem problem(Input());
    const auto fit=core::FitJointComponents(problem,{.6,missing_width,missing_width});
    ASSERT_TRUE(fit.assembled_state);
    const auto & data=core::JointProblemAccess::Get(problem); const auto & layout=problem.ParameterLayout();
    auto context=n::ProfileContext(data.context,layout,24);
    const auto geometry=n::BuildTargetGeometry(n::ProfileDomain(data.domain,layout),data.y.head(23),*fit.assembled_state,context);
    ASSERT_TRUE(geometry.valid);
    const auto rows=static_cast<std::size_t>(geometry.svd.singular_values.tail(1)(0)/
        (std::numeric_limits<double>::epsilon()*geometry.svd.singular_values(0)));
    const auto evidence=n::AssessTargetState(problem,layout,fit.assembled_state,fit.evidence,rhbm_gem::JointEvidenceScope::AssembledGlobal,rows);
    EXPECT_EQ(n::TargetConvergence(evidence,true),rhbm_gem::JointCheckStatus::Unavailable);
    EXPECT_EQ(evidence.atoms[0].reason,"rank-threshold-sensitive");
}

TEST(JointObservableProfileTest, TargetEvidenceUsesEveryComponentAndTheSavedState)
{
    auto input=TwoRowHalo(false,true);
    input.support[3].push_back({2,.8});
    input.observations[2]+=n::EvaluateKernel(.8,.65,2.5).gaussian+
        .2*n::EvaluateKernel(.8,.65,2.5).charge;
    const auto copy=input;
    for(std::size_t a=0;a<copy.atom_ids.size();++a)
    {
        input.atom_ids.push_back(std::to_string(a+5)); input.support.emplace_back();
        for(const auto & p:copy.support[a]) input.support.back().push_back({p.row+24,p.squared_distance});
    }
    input.selection_domain->target_indices.push_back(4);
    for(std::size_t r=0;r<24;++r)
    {input.row_ids.push_back(std::to_string(r+24)); input.observations.push_back(copy.observations[r]);}
    const core::JointProblem problem(input);
    const auto fit=core::FitJointComponents(problem,{.6,missing_width,missing_width,.65,.6,missing_width,missing_width,.65});
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
    const auto fresh_target=n::AssessTargetState(problem,problem.ParameterLayout(),fit.assembled_state,fit.evidence,
        rhbm_gem::JointEvidenceScope::AssembledGlobal,problem.Input().observations.size());
    ASSERT_EQ(fit.target_evidence->checks.size(),fresh_target.checks.size());
    for(std::size_t k=0;k<fresh_target.checks.size();++k)
    {
        EXPECT_EQ(fit.target_evidence->checks[k].name,fresh_target.checks[k].name);
        EXPECT_EQ(fit.target_evidence->checks[k].status,fresh_target.checks[k].status);
        EXPECT_EQ(fit.target_evidence->checks[k].scope,fresh_target.checks[k].scope);
        EXPECT_EQ(fit.target_evidence->checks[k].reason,fresh_target.checks[k].reason);
    }
    ASSERT_EQ(fit.target_evidence->atoms.size(),fresh_target.atoms.size());
    for(std::size_t k=0;k<fresh_target.atoms.size();++k)
    {
        EXPECT_EQ(fit.target_evidence->atoms[k].atom,fresh_target.atoms[k].atom);
        EXPECT_EQ(fit.target_evidence->atoms[k].status,fresh_target.atoms[k].status);
        EXPECT_EQ(fit.target_evidence->atoms[k].null_space_leakage,fresh_target.atoms[k].null_space_leakage);
    }
    const auto saved=core::CaptureJointAnalysisResult(fit);
    const auto covariance=core::detail::ComputeJointUncertainty(problem,saved);
    EXPECT_TRUE(covariance.at(1).covariance); EXPECT_TRUE(covariance.at(5).covariance);
    auto json=boost::json::parse(rhbm_gem::joint_result_io::Encode(saved));
    EXPECT_EQ(json.at("schema_version").as_int64(),8);
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
}

TEST(JointObservableProfileTest, UnobservedTargetHasNoTargetCertificate)
{
    core::JointProblemInput input; input.atom_ids={"1"}; input.support.resize(1);
    input.selection_domain.emplace(); input.selection_domain->target_indices={0};
    const auto fit=core::FitJointComponents(core::JointProblem(input),{.6});
    EXPECT_EQ(fit.TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::Unavailable);
    EXPECT_EQ(fit.components.at(0).TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::Unavailable);
}
