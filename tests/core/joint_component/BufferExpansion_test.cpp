#include <gtest/gtest.h>

#include "core/estimator/joint_component/BufferExpansion.hpp"
#include "core/command/detail/SimulationGeometry.hpp"
#include "core/estimator/InitialStage.hpp"
#include "core/estimator/PotentialFittingWorkflow.hpp"
#include "core/estimator/joint_component/JointComponentEstimator.hpp"
#include "core/estimator/joint_component/Numerics.hpp"
#include "core/estimator/joint_component/PostFitPeeling.hpp"
#include "core/estimator/joint_component/ResourceWork.hpp"
#include "core/estimator/joint_component/WidthCalibration.hpp"
#include "data/io/detail/JointResultJson.hpp"
#include "data/io/detail/StageResultJson.hpp"
#include "support/CommandTestHelpers.hpp"

#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <boost/json.hpp>
#include <limits>

#ifdef RHBM_GEM_ENABLE_UMAP
#include "core/classifier/LocalFittingFeatures.hpp"
#endif

namespace {
using namespace rhbm_gem;
using namespace rhbm_gem::core;
using namespace rhbm_gem::core::joint_component;

TEST(BufferExpansionTest, ScoresCompleteBlockersWithFullCoreDenominators)
{
    JointProblemInput input;
    input.observations.resize(10);
    input.support={
        {{0,0.0},{1,.25},{2,.25},{3,.25},{4,.25},{5,.25},{6,2.},{7,2.},{8,2.},{9,2.}},
        {},
        {{1,.1},{2,.2},{6,6.25}},
        {{3,.2},{5,.3},{7,6.25},{9,3.0}},
        {{5,.5},{8,.5}}
    };
    // Buffer 2 has two core incidences at 1/6 each. Buffer 3 has one core
    // incidence at 1/6 and one at 1/1. Row 5 has two blockers.
    const std::array<std::size_t,2> outputs{0,1};
    const auto choice=SelectBuffer(input,outputs);
    EXPECT_EQ(choice.buffer,2u);
    EXPECT_EQ(choice.recovered_rows,3u);
    input.support[1]={{1,.25001},{3,.25}};
    const auto shared=SelectBuffer(input,outputs);
    EXPECT_EQ(shared.buffer,3u);
    EXPECT_EQ(shared.recovered_rows,3u);
    EXPECT_EQ(shared.retained_rows,(std::vector<bool>{true,false,false,true,true,false,false,true,false,true}));

    input.observations.resize(11);
    input.support[0]={{0,0.0},{1,.25},{2,.25},{4,.25},{5,2.},{6,2.},{7,2.},{8,2.},{9,2.},{10,.25}};
    input.support[1]={{0,0.0},{3,.25}};
    // The candidates recover 2/5 and 1/2. A blocked-only denominator would tie.
    EXPECT_EQ(SelectBuffer(input,outputs).buffer,3u);
}

TEST(BufferExpansionTest, UsesCanonicalTiesAndScreensDistinctRadialSupport)
{
    JointProblemInput input;
    input.observations.resize(9);
    const double close=1.0+4*std::numeric_limits<double>::epsilon();
    input.support={
        {{0,.25},{1,.25},{2,.25},{3,.3},{4,.4},{5,.5},{6,.6},{7,.7},{8,.8}},
        {{0,1.0},{3,close},{4,1.0}},
        {{1,1.0},{5,1.2},{5,1.2}},
        {{2,1.0},{6,1.5},{7,2.0}}
    };
    const std::array<std::size_t,1> outputs{0};
    EXPECT_EQ(SelectBuffer(input,outputs).buffer,3u);
    input.support[2].push_back({8,1.1});
    EXPECT_EQ(SelectBuffer(input,outputs).buffer,2u);
    input.support[1][1].squared_distance=1.0+16*std::numeric_limits<double>::epsilon();
    EXPECT_EQ(SelectBuffer(input,outputs).buffer,1u);
    input.support[0][0].squared_distance=input.support[0][1].squared_distance=
        input.support[0][2].squared_distance=.251;
    const auto no_core=SelectBuffer(input,outputs);
    EXPECT_FALSE(no_core.buffer);
    EXPECT_EQ(no_core.recovered_rows,0u);
    EXPECT_EQ(std::count(no_core.retained_rows.begin(),no_core.retained_rows.end(),true),0);
}

TEST(BufferExpansionTest, ScratchEstimateUsesCheckedArithmeticAndBuildBudget)
{
    const auto scratch=EstimateBufferSelectionScratchBytes(100,10);
    ASSERT_TRUE(scratch);
    EXPECT_GT(*scratch,100*sizeof(std::size_t));
    const auto maximum=std::numeric_limits<std::size_t>::max();
    EXPECT_FALSE(EstimateBufferSelectionScratchBytes(maximum,1));
    EXPECT_FALSE(EstimateBufferSelectionScratchBytes(1,maximum));
    MemoryShape shape;
    shape.problem_scope=true;
    shape.resident_bytes=JointMemoryPlanningBudgetBytes;
    ASSERT_TRUE(CheckedAddBytes(shape.resident_bytes,*scratch,shape.resident_bytes));
    EXPECT_EQ(PlanMemory(shape,MemoryStage::ProblemConstruction).status,MemoryStatus::Reject);
    EXPECT_EQ(JointMemoryPlanningBudgetBytes,3'200'000'000u);
    EXPECT_EQ(JointMemoryHardLimitBytes,4'000'000'000u);
    shape.resident_bytes=JointMemoryHardLimitBytes+1;
    const auto hard_limit=PlanMemory(shape,MemoryStage::ProblemConstruction,maximum);
    EXPECT_EQ(hard_limit.budget_bytes,JointMemoryHardLimitBytes);
    EXPECT_EQ(hard_limit.status,MemoryStatus::Reject);
}

std::unique_ptr<AtomObject> Carbon(int id,std::array<double,3> position)
{
    auto atom=std::make_unique<AtomObject>();
    atom->SetSerialID(id);
    atom->SetElement(Element::CARBON);
    atom->SetPosition(position);
    atom->SetChainID("A");
    atom->SetComponentID("ALA");
    atom->SetComponentKey(static_cast<ComponentKey>(id));
    atom->SetAtomKey(static_cast<AtomKey>(7));
    atom->SetAtomID("CA");
    return atom;
}

TEST(BufferExpansionTest, RestoresEmptyBaselineOnFixedNativeDomain)
{
    std::vector<std::unique_ptr<AtomObject>> atoms;
    atoms.push_back(Carbon(40,{2,0,0}));
    atoms.push_back(Carbon(12,{2.25,0,0}));
    atoms.push_back(Carbon(3,{2,0,0}));
    atoms.back()->SetElement(Element::HYDROGEN);
    ModelObject model(std::move(atoms));
    model.SelectAtoms([](const auto & atom){return atom.GetSerialID()==40;});
    const MapObject map({9,1,1},{.5,.5,.5},{0,0,0});
    const auto problem=BuildJointProblem(map,model);
    const auto & original=problem.OriginalInput();
    EXPECT_EQ(original.atom_ids,(std::vector<std::string>{"12","40"}));
    EXPECT_EQ(original.selection_domain->output_indices,(std::vector<std::size_t>{1}));
    EXPECT_EQ(original.selection_domain->buffer_indices,(std::vector<std::size_t>{0}));
    EXPECT_EQ(problem.FittingToOriginalContributors(),(std::vector<std::size_t>{0,1}));
    EXPECT_EQ(problem.Input().row_ids,original.row_ids);
    EXPECT_EQ(problem.FittingToOriginalRows(),(std::vector<std::size_t>{0,1,2,3,4,5,6,7,8}));
    EXPECT_EQ(original.selection_domain->buffer_expansion->recovered_rows,9u);
    ASSERT_EQ(model.GetSelectedAtoms().size(),1u);
    EXPECT_EQ(model.GetSelectedAtoms().front()->GetSerialID(),40);
}

TEST(BufferExpansionTest, RejectsInvalidContributorIdentity)
{
    std::vector<std::unique_ptr<AtomObject>> atoms;
    atoms.push_back(Carbon(40,{2,0,0}));
    atoms.push_back(Carbon(0,{2.25,0,0}));
    ModelObject model(std::move(atoms));
    model.SelectAtoms([](const auto & atom){return atom.GetSerialID()==40;});
    const MapObject map({9,1,1},{.5,.5,.5},{0,0,0});
    EXPECT_THROW(BuildJointProblem(map,model),std::invalid_argument);
}

TEST(BufferExpansionTest, WidthFallbackAndCalibrationUseOutputDonors)
{
    JointInitialization initialization;
    for(std::size_t i=0;i<21;++i)
    {
        const bool buffer=i==3;
        const auto width=buffer ? 12.0 : .5;
        auto & atom=initialization.atoms.emplace_back();
        atom.role=buffer ? FittingRole::BufferTarget : FittingRole::OutputTarget;
        atom.reason="valid-width";
        atom.original_tau=width;
        atom.mdpde[1]=width;
        atom.seed_source="fitted";
        atom.native_status=static_cast<int>(RHBMEstimationStatus::SUCCESS);
        atom.initial_stage_fit_qualified=true;
        initialization.tau.push_back(width);
    }
    initialization.tau[0]=0;
    initialization.atoms[0].reason="invalid-width";
    rhbm_gem::core::detail::ApplyJointSeedFallback(initialization);
    EXPECT_EQ(initialization.atoms[0].donor_count,19u);
    EXPECT_EQ(initialization.tau[0],.5);
    EXPECT_EQ(CalibrateGlobalWidth(initialization).qualified_donor_count,19u);

    initialization.atoms[0]=initialization.atoms[1];
    EXPECT_EQ(CalibrateGlobalWidth(initialization).qualified_donor_coverage,1.0);
    const auto calibration=CalibrateGlobalWidth(initialization);
    EXPECT_EQ(calibration.median_width,.5);
    EXPECT_EQ(calibration.target_count,20u);
    EXPECT_EQ(MakeBoundedInitialization(initialization,calibration,.125).tau_values.at(3),2.0);

    initialization.tau[3]=0;
    initialization.atoms[3].reason="invalid-width";
    EXPECT_EQ(CalibrateGlobalWidth(initialization).unavailable_reason,"invalid-target-initial-width");
    rhbm_gem::core::detail::ApplyJointSeedFallback(initialization);
    EXPECT_EQ(initialization.atoms[3].donor_count,20u);
    EXPECT_EQ(initialization.atoms[3].seed_source,"median-fallback");
    EXPECT_EQ(initialization.tau[3],.5);
    for(auto & atom:initialization.atoms) atom.reason="invalid-width";
    initialization.tau[3]=0;
    rhbm_gem::core::detail::ApplyJointSeedFallback(initialization);
    EXPECT_FALSE(initialization.valid);
}

TEST(BufferExpansionTest, CertifiedBufferSurvivesWorkflowAndSQLiteRoundTrip)
{
    std::vector<std::unique_ptr<AtomObject>> atoms;
    for(int i=0;i<20;++i)
        atoms.push_back(Carbon(i==0 ? 10 : 30+3*i,{7.0*i,0,0}));
    atoms.push_back(Carbon(20,{2,0,0}));
    atoms.push_back(Carbon(15,{3.75,1.5,0}));
    constexpr std::array<int,3> dimensions{281,17,17};
    constexpr std::array<double,3> spacing{.5,.5,.5},origin{-3,-4,-4};
    const auto row_count=static_cast<std::size_t>(dimensions[0])*dimensions[1]*dimensions[2];
    MapObject map(dimensions,spacing,origin);
    auto values=std::make_unique<double[]>(row_count);
    for(std::size_t row=0;row<row_count;++row)
    {
        const auto position=map.GetGridPosition(row);
        for(std::size_t i=0;i<atoms.size();++i)
        {
            const double square=simulation::SupportSquare(position,atoms[i]->GetPosition());
            if(square>6.25) continue;
            const auto index=static_cast<double>(i);
            const auto basis=EvaluateKernel(square,.5+.01*index,2.5);
            values[row]+=(2+.02*index)*basis.gaussian+.1*basis.charge;
        }
    }
    map.SetMapValueArray(std::move(values));
    ModelObject model(std::move(atoms));
    model.SelectAtoms([](const auto & atom){return atom.GetSerialID()!=15 && atom.GetSerialID()!=20;});

    const auto problem=BuildJointProblem(map,model);
    const auto & original=problem.OriginalInput();
    ASSERT_EQ(original.selection_domain->buffer_indices,(std::vector<std::size_t>{2}));
    EXPECT_EQ(problem.Input().selection_domain->buffer_indices,(std::vector<std::size_t>{1}));
    ASSERT_EQ(original.atom_ids.at(2),"20");
    EXPECT_EQ(problem.FittingToOriginalContributors().at(1),2u);
    EXPECT_EQ(original.selection_domain->buffer_expansion->external_halos,1u);
    EXPECT_GT(original.selection_domain->buffer_expansion->recovered_rows,0u);
    EXPECT_LT(problem.Input().row_ids.size(),original.row_ids.size());
    for(const auto & entry:original.support.at(1))
        EXPECT_FALSE(std::binary_search(problem.FittingToOriginalRows().begin(),problem.FittingToOriginalRows().end(),entry.row));

    FitOptions options;
    options.quiet_mode=true;
    options.thread_size=1;
    rhbm_gem::core::detail::JointPostprocessingMetrics metrics;
    rhbm_gem::core::detail::RunPotentialFittingWorkflowForSelection(map,model,options,&metrics);
    const auto & result=model.GetAnalysisView().GetJointResult();
    ASSERT_TRUE(result);
    ASSERT_TRUE(IsQualifiedJointAnalysisResult(*result));
    EXPECT_EQ(result->initialization.atoms.at(1).role,FittingRole::BufferTarget);
    EXPECT_GT(result->initialization.atoms.at(1).sample_count,0u);
    EXPECT_EQ(result->initialization.atoms.at(1).seed_source,"fitted");
    ASSERT_TRUE(result->assembled_state);
    const auto & state=*result->assembled_state;
    EXPECT_NEAR(state.ac.at(2),2.4,1e-4);
    EXPECT_NEAR(state.tau.at(1),.7,1e-4);
    EXPECT_NEAR(state.ac.at(3),.1,1e-4);
    EXPECT_TRUE(std::any_of(result->components.begin(),result->components.end(),[](const auto & component) {
        return component.atoms==std::vector<std::size_t>{0,2};
    }));
    EXPECT_TRUE(std::any_of(result->evidence.begin(),result->evidence.end(),[](const auto & check) {
        return check.name=="global-width-calibration" && check.reason.ends_with("qualified-donors=20/20");
    }));
    for(const auto * name:{"independent-replay","endpoint-verification","width-stationarity"})
        EXPECT_TRUE(std::any_of(result->evidence.begin(),result->evidence.end(),[name](const auto & check) {
            return check.name==name && check.status==JointCheckStatus::Passed;
        }));
    ASSERT_TRUE(result->contributor_reconstruction);
    ASSERT_EQ(result->contributor_reconstruction->replacements.size(),1u);
    const auto & halo=result->contributor_reconstruction->replacements.front();
    EXPECT_EQ(halo.atom_id,"15");
    EXPECT_NE(std::find(halo.donor_atom_ids.begin(),halo.donor_atom_ids.end(),"20"),halo.donor_atom_ids.end());

    const command_test::ScopedTempDir dir{"joint_buffer_expansion"};
    const DataRepository repository{dir.path()/"buffer.sqlite"};
    repository.SaveModel(model,"fit");
    const auto loaded=repository.LoadModel("fit");
    ASSERT_EQ(joint_result_io::Encode(*loaded->GetAnalysisView().GetJointResult()),joint_result_io::Encode(*result));
    auto stage=boost::json::parse(stage_result_io::Encode(model));
    auto loaded_stage=boost::json::parse(stage_result_io::Encode(*loaded));
    for(const auto & [array,identity]:{std::pair{"atoms","id"},std::pair{"groups","key"}})
    {
        auto & expected=stage.at(array).as_array();
        auto & actual=loaded_stage.at(array).as_array();
        const auto by_identity=[identity](const auto & a,const auto & b) {
            return boost::json::value_to<std::uint64_t>(a.at(identity))<
                boost::json::value_to<std::uint64_t>(b.at(identity));
        };
        std::sort(expected.begin(),expected.end(),by_identity);
        std::sort(actual.begin(),actual.end(),by_identity);
        EXPECT_TRUE(actual==expected) << array;
    }
    ASSERT_EQ(loaded->GetSelectedAtoms().size(),20u);
    const auto buffer=AtomLocalPotentialView::For(*loaded->FindAtomPtr(20));
    EXPECT_EQ(buffer.GetInitialStageEstimate().source.role,FittingRole::BufferTarget);
    const auto & point=buffer.GetJointStageEstimate().point;
    ASSERT_TRUE(point);
    EXPECT_EQ(buffer.GetJointStageEstimate().source.role,FittingRole::BufferTarget);
    EXPECT_EQ(point->GetAmplitude(),state.ac.at(2));
    EXPECT_EQ(point->GetWidth(),state.tau.at(1));
    EXPECT_EQ(point->GetOffset(),state.ac.at(3));
    EXPECT_FALSE(buffer.GetPostFitPeeling());
    EXPECT_FALSE(buffer.GetPostJointLocalResult());
    EXPECT_FALSE(buffer.GetPostJointQScore());
    EXPECT_FALSE(buffer.GetGroupMemberResult());
    EXPECT_EQ(metrics.qscore_targets,20u);

    const auto target=AtomLocalPotentialView::For(*loaded->FindAtomPtr(10));
    ASSERT_TRUE(target.GetPostFitPeeling());
    EXPECT_TRUE(target.GetPostJointLocalResult());
    const auto & reloaded_snapshot=*loaded->GetAnalysisView().GetJointResult();
    const rhbm_gem::core::detail::JointPostprocessingContext context(map,*loaded,problem,reloaded_snapshot,
        nullptr,nullptr,rhbm_gem::core::detail::BuildJointPeelingContributorUniverse(*loaded,problem,reloaded_snapshot));
    const auto rebuilt=context.BuildPostFitPeelingBatch(original.selection_domain->output_indices);
    ASSERT_EQ(rebuilt.size(),20u);
    ASSERT_EQ(rebuilt.front().first,10);
    const auto & samples=target.GetPostFitPeeling()->samples;
    const auto center=std::find_if(samples.begin(),samples.end(),[](const auto & sample) {
        return sample.observation.squared_distance==0;
    });
    ASSERT_NE(center,samples.end());
    ASSERT_TRUE(center->response);
    const auto buffer_basis=EvaluateKernel(4.0,point->GetWidth(),2.5);
    EXPECT_NEAR(*center->response,center->observation.map_value-
        (point->GetAmplitude()*buffer_basis.gaussian+point->GetOffset()*buffer_basis.charge),1e-14);
    const auto sample_index=static_cast<std::size_t>(center-samples.begin());
    EXPECT_EQ(rebuilt.front().second.samples.at(sample_index).response,center->response);

    ModelObject direct_model(*loaded);
    const auto direct=EstimateJointComponents(map,direct_model);
    ASSERT_TRUE(direct.prediction);
    const auto & fitted_rows=problem.FittingToOriginalRows();
    const auto fitted_center=std::lower_bound(fitted_rows.begin(),fitted_rows.end(),center->observation.original_row_index);
    ASSERT_NE(fitted_center,fitted_rows.end());
    ASSERT_EQ(*fitted_center,center->observation.original_row_index);
    EXPECT_NEAR(direct.prediction->at(static_cast<std::size_t>(fitted_center-fitted_rows.begin())),
        center->observation.map_value,1e-10);

#ifdef RHBM_GEM_ENABLE_UMAP
    const auto features=classifier::BuildLocalFittingFeatureRows(*loaded);
    ASSERT_EQ(features.size(),20u);
    for(std::size_t i=0;i<features.size();++i)
        EXPECT_EQ(features[i].serial_id,i==0 ? 10 : 30+3*static_cast<int>(i));
#endif

    const auto encoded=boost::json::parse(joint_result_io::Encode(*result));
    auto invalid=encoded;
    invalid.at("selection_domain").as_object()["buffer_indices"]=boost::json::array{0};
    EXPECT_THROW(joint_result_io::Decode(boost::json::serialize(invalid)),std::invalid_argument);
    invalid.at("selection_domain").as_object()["buffer_indices"]=boost::json::array{2,2};
    EXPECT_THROW(joint_result_io::Decode(boost::json::serialize(invalid)),std::invalid_argument);
    invalid=encoded;
    invalid.at("fitting_to_original_contributors").as_array()[1]=1;
    EXPECT_THROW(joint_result_io::Decode(boost::json::serialize(invalid)),std::invalid_argument);
    invalid=encoded;
    invalid.at("assembled_state").at("ac").as_array()[2]=-1;
    EXPECT_THROW(joint_result_io::Decode(boost::json::serialize(invalid)),std::invalid_argument);
    invalid=encoded;
    for(auto & component:invalid.at("components").as_array())
        if(std::any_of(component.at("atoms").as_array().begin(),component.at("atoms").as_array().end(),
                [](const auto & atom){return atom==2;}))
            component.as_object()["runtime_convergence"]="failed";
    EXPECT_THROW(joint_result_io::Decode(boost::json::serialize(invalid)),std::invalid_argument);

    auto & buffer_stage=loaded_stage.at("atoms").as_array().at(2).at("joint").at("source").as_object();
    buffer_stage["role"]=99;
    EXPECT_THROW(stage_result_io::Decode(*loaded,boost::json::serialize(loaded_stage)),std::invalid_argument);
    buffer_stage["role"]=static_cast<int>(FittingRole::OutputTarget);
    EXPECT_THROW(stage_result_io::Decode(*loaded,boost::json::serialize(loaded_stage)),std::invalid_argument);
}
}
