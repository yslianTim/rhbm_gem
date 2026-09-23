#include <gtest/gtest.h>
#include <boost/json.hpp>
#include <rhbm_gem/data/io/JointAnalysisFileIO.hpp>
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/JointUncertainty.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include "data/io/detail/JointResultJson.hpp"
#include <limits>
#include "core/detail/PostFitPeeling.hpp"
#include "support/ForwardModelExperiment.hpp"
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
    ASSERT_TRUE(assessment.primary.valid); ASSERT_TRUE(assessment.reference.valid);
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
    const auto raw=second_stage_test::SampleExperimentPoints(map,points);
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
    const auto initialization=core::detail::RunFirstStage(model,workset,options,core::detail::FirstStageMode::SampleContributorsIsolated,&map);
    core::detail::FirstStageObserverForTesting()={};
    EXPECT_TRUE(std::all_of(visited.begin(),visited.end(),[](int id){return id==1;}));
    EXPECT_EQ(initialization.atoms[1].reason,"not-required-observable-contribution");
    EXPECT_EQ(initialization.atoms[2].reason,"not-required-observable-contribution");
    EXPECT_FALSE(initialization.atoms[1].original_b); EXPECT_FALSE(initialization.atoms[1].used_b);
}
