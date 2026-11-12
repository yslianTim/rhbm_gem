#include <gtest/gtest.h>
#include "support/InitialStageSynthetic.hpp"
#include "support/CommandTestHelpers.hpp"
#include "core/estimator/InitialStage.hpp"
#include "core/estimator/PotentialFittingWorkflow.hpp"
#include "core/estimator/joint_component/Problem.hpp"
#include "core/estimator/joint_component/WidthCalibration.hpp"
#include "data/io/detail/JointResultJson.hpp"
#include "data/io/detail/StageResultJson.hpp"
#include "core/sampling/MapSampler.hpp"
#include <rhbm_gem/core/QScoreHelper.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <boost/json.hpp>
#include <bit>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>

namespace {
using namespace rhbm_gem;
using namespace rhbm_gem::core;
namespace d=rhbm_gem::core::detail;
namespace n=rhbm_gem::core::joint_component;
namespace j=boost::json;

j::object Baseline()
{
    const auto path=std::filesystem::path(__FILE__).parent_path().parent_path().parent_path()/
        "fixtures/joint_component/initial_stage.json";
    std::ifstream stream(path);
    if(!stream.good()) throw std::runtime_error("Missing InitialStage regression fixture");
    j::parse_options options; options.numbers=j::number_precision::precise;
    return j::parse(std::string(std::istreambuf_iterator<char>(stream),{}),{},options).as_object();
}

j::array Checks(const std::vector<JointCheck> & checks)
{
    j::array out;
    for(const auto & c:checks)
    {
        j::object entry{{"name",c.name},{"status",static_cast<int>(c.status)},{"reason",c.reason}};
        out.push_back(entry);
    }
    return out;
}

// FNV-1a over IEEE-754 words, low byte first, independent of host endianness.
// Separate fingerprints make exact raw geometry/mask/response regressions visible.
j::object SampleFingerprint(const LocalPotentialSampleList & samples)
{
    std::array<std::uint64_t,4> hashes; hashes.fill(14695981039346656037ull);
    const auto add=[&](std::size_t field,std::uint64_t bits) {
        for(int byte=0;byte<8;++byte) { hashes[field]^=(bits&255); hashes[field]*=1099511628211ull; bits>>=8; }
    };
    std::size_t selected=0;
    for(const auto & s:samples)
    {
        for(double v:s.point.position) add(0,std::bit_cast<std::uint64_t>(v));
        add(1,std::bit_cast<std::uint64_t>(s.point.distance));
        add(2,s.point.is_selected); add(3,std::bit_cast<std::uint64_t>(s.response));
        selected+=s.point.is_selected;
    }
    return {{"positions",std::to_string(hashes[0])},{"distances",std::to_string(hashes[1])},
        {"selection",std::to_string(hashes[2])},{"responses",std::to_string(hashes[3])},
        {"raw",samples.size()},{"selected",selected}};
}

TEST(InitialStageTest, ExactProductionBaselineAndCertifiedEndpoint)
{
    const auto baseline=Baseline();
    auto model=joint_component_test::InitialStageSyntheticModel(5.5);
    auto map=joint_component_test::InitialStageSyntheticMap(*model,.4);
    const auto problem=BuildPureHaloSupportExclusionProblem(*map,*model);
    const auto workset=d::BuildInitialStageWorkset(*model,problem);
    model->EditAnalysis().InitializeFromSelection();
    FitOptions options; options.quiet_mode=true;
    const auto init=d::RunInitialStage(*map,*model,workset,options);
    ASSERT_TRUE(init.valid); ASSERT_EQ(init.atoms.size(),24u);
    for(std::size_t i=0;i<init.atoms.size();++i)
    {
        SCOPED_TRACE(i);
        const auto & atom=init.atoms[i];
        const auto & expected=baseline.at("atoms").as_array()[i].as_object();
        EXPECT_EQ(atom.id,j::value_to<std::string>(expected.at("id")));
        EXPECT_EQ(atom.alpha,j::value_to<double>(expected.at("alpha")));
        for(std::size_t k=0;k<3;++k)
        {
            EXPECT_EQ(atom.ols[k],j::value_to<double>(expected.at("ols").as_array()[k]));
            EXPECT_EQ(atom.mdpde[k],j::value_to<double>(expected.at("mdpde").as_array()[k]));
        }
        EXPECT_TRUE(atom.target); EXPECT_TRUE(atom.initial_stage_fit_qualified);
        EXPECT_EQ(atom.initial_stage_qualification,RHBMSolveQualification::NativeSuccess);
        EXPECT_EQ(atom.reason,"valid-width"); EXPECT_EQ(atom.seed_source,"fitted");
        EXPECT_EQ(atom.original_tau,atom.mdpde[1]); EXPECT_EQ(atom.used_tau,atom.mdpde[1]);
        EXPECT_EQ(init.tau[i],atom.mdpde[1]); EXPECT_EQ(atom.sample_count,200u);
    }
    const auto calibration=n::CalibrateGlobalWidth(init,"baseline");
    const auto & expected_calibration=baseline.at("calibration").as_object();
    EXPECT_EQ(calibration.status,"Available");
    EXPECT_EQ(calibration.qualified_donor_count,24u); EXPECT_EQ(calibration.native_success_count,24u);
    EXPECT_EQ(calibration.qualified_donor_coverage,1.); EXPECT_EQ(calibration.fallback_count,0u);
    EXPECT_EQ(calibration.refined_success_count,0u); EXPECT_EQ(calibration.unqualified_count,0u);
    EXPECT_EQ(calibration.initial_violations4,0u);
    for(const auto & pair:{std::pair{"bound",calibration.tau_max4},
        {"p50",calibration.median_width},{"p90",calibration.p90},{"p95",calibration.p95},{"p99",calibration.p99}})
    { ASSERT_TRUE(pair.second); EXPECT_EQ(*pair.second,j::value_to<double>(expected_calibration.at(pair.first))); }
    const auto bounded=n::MakeBoundedInitialization(init,calibration);
    ASSERT_TRUE(bounded.available); EXPECT_EQ(bounded.tau_values,init.tau); EXPECT_EQ(bounded.projection_count,0u);
    const auto fit=FitJointComponentsInitialized(problem,init);
    ASSERT_TRUE(fit.assembled_state); ASSERT_TRUE(fit.objective);
    const auto & expected=baseline.at("joint").as_object();
    EXPECT_EQ(fit.objective,j::value_to<double>(expected.at("objective")));
    EXPECT_EQ(fit.assembled_state->tau,j::value_to<std::vector<double>>(expected.at("tau")));
    EXPECT_EQ(fit.assembled_state->ac,j::value_to<std::vector<double>>(expected.at("ac")));
    EXPECT_EQ(fit.assembled_state->eta_gradient,j::value_to<std::vector<double>>(expected.at("eta_gradient")));
    EXPECT_EQ(fit.RuntimeConvergence(),JointCheckStatus::Passed);
    EXPECT_EQ(fit.TargetRuntimeConvergence(),JointCheckStatus::Passed);
    EXPECT_EQ(Checks(fit.evidence),expected.at("checks"));
    ASSERT_EQ(fit.components.size(),expected.at("components").as_array().size());
    for(std::size_t i=0;i<fit.components.size();++i)
    {
        const auto & component=fit.components[i];
        const auto & reference=expected.at("components").as_array()[i].as_object();
        EXPECT_EQ(component.id,j::value_to<std::string>(reference.at("id")));
        EXPECT_EQ(component.search_sweeps,j::value_to<std::size_t>(reference.at("sweeps")));
        EXPECT_EQ(component.stop_reason,j::value_to<std::string>(reference.at("stop")));
        EXPECT_EQ(Checks(component.evidence),expected.at("component_checks"));
    }
    // The original Debug correlation differs from Release by one ULP.
#ifdef NDEBUG
    const auto qscore_key="qscore";
#else
    const auto qscore_key="qscore_debug";
#endif
    EXPECT_EQ(CalculateQScoreForAtom(*workset.contributors[0],*map,*model,1.,0.),
        j::value_to<double>(baseline.at(qscore_key)));
    d::JointPostprocessingMetrics metrics;
    d::CompletePotentialFittingWorkflow(*map,*model,problem,CaptureJointAnalysisResult(fit),options,metrics);
    EXPECT_EQ(metrics.available_peeling_sample_count,metrics.grid_peeling_sample_count);
    j::array downstream;
    for(const auto * atom:workset.contributors)
    {
        const auto view=AtomLocalPotentialView::For(*atom);
        ASSERT_TRUE(view.GetPostFitPeeling()); ASSERT_TRUE(view.GetPostJointLocalResult());
        EXPECT_EQ(view.GetPostJointLocalResult()->status,EvidenceStatus::Available);
        const auto & local=*view.GetPostJointLocalResult();
        const auto ols=local.fit.ols.GetModel().ToVector(),mdpde=local.fit.mdpde.GetModel().ToVector();
        LocalPotentialSampleList peeled;
        for(const auto & sample:view.GetPostFitPeeling()->samples)
        {
            ASSERT_TRUE(sample.response);
            peeled.push_back({*sample.response,{sample.observation.squared_distance,sample.observation.position,true}});
        }
        downstream.push_back(j::object{{"peeled",SampleFingerprint(peeled)},
            {"ols",j::array{ols[0],ols[1],ols[2]}},{"mdpde",j::array{mdpde[0],mdpde[1],mdpde[2]}},
            {"retained",local.effective_sample_count},{"group",view.GetGroupMemberResult().has_value()}});
    }
    EXPECT_EQ(downstream,baseline.at("downstream"));
    const auto snapshot=model->GetAnalysisView().GetJointResult(); ASSERT_TRUE(snapshot);
    const auto encoded=joint_result_io::Encode(*snapshot);
    EXPECT_EQ(joint_result_io::Encode(joint_result_io::Decode(encoded)),encoded);
}

TEST(InitialStageTest, ExactRawSamplingAcrossWorkersAndCrowdedSelection)
{
    j::array captured;
    for(double separation:{5.5,1.4})
    {
        auto model=joint_component_test::InitialStageSyntheticModel(separation);
        auto map=joint_component_test::InitialStageSyntheticMap(*model,.4);
        model->EditAnalysis().InitializeFromSelection();
        const auto workset=d::BuildInitialStageWorkset(*model,BuildPureHaloSupportExclusionProblem(*map,*model));
        std::vector<LocalPotentialSampleList> serial;
        RunPotentialSamplingWorkflow(*map,*model,workset.contributors,1,{},true);
        for(const auto * atom:workset.contributors)
        {
            serial.push_back(AtomLocalPotentialView::For(*atom).GetSamplingEntriesView());
            captured.push_back(SampleFingerprint(serial.back()));
        }
        FitOptions options; options.quiet_mode=true; options.joint_initialization_thread_size=4;
        const auto init=d::RunInitialStage(*map,*model,workset,options);
        const auto calibration=n::CalibrateGlobalWidth(init,"sampling-regression");
        EXPECT_EQ(calibration.qualified_donor_count,24u); EXPECT_EQ(calibration.fallback_count,0u);
        if(separation==1.4)
        {
            double square=0.;
            for(std::size_t i=0;i<init.atoms.size();++i)
            {
                const double error=init.atoms[i].mdpde[1]-(.35+.005*static_cast<double>(i));
                square+=error*error;
            }
            EXPECT_LT(std::sqrt(square/24.),.03); // Original crowded F0 truth contract.
        }
        for(std::size_t i=0;i<serial.size();++i)
        {
            const auto & raw=AtomLocalPotentialView::For(*workset.contributors[i]).GetSamplingEntriesView();
            ASSERT_EQ(raw.size(),serial[i].size());
            for(std::size_t k=0;k<raw.size();++k)
            {
                EXPECT_EQ(raw[k].point.position,serial[i][k].point.position);
                EXPECT_EQ(raw[k].point.distance,serial[i][k].point.distance);
                EXPECT_EQ(raw[k].point.is_selected,serial[i][k].point.is_selected);
                EXPECT_EQ(raw[k].response,serial[i][k].response);
            }
        }
    }
    EXPECT_EQ(captured,Baseline().at("sampling_fingerprints"));
}

TEST(InitialStageTest, MedianFallbackPreservesRawEvidenceAndExcludesCalibrationDonors)
{
    JointInitialization init;
    for(int i=0;i<24;++i)
    {
        JointInitializationAtom atom; atom.id=std::to_string(i); atom.target=true;
        atom.mdpde={1.,1.+.01*i,0.}; atom.ols={1.,.5,0.};
        atom.original_tau=atom.mdpde[1]; atom.reason="valid-width"; atom.seed_source="fitted";
        atom.native_status=static_cast<int>(RHBMEstimationStatus::SUCCESS);
        atom.initial_stage_fit_qualified=true; atom.initial_stage_qualification=RHBMSolveQualification::NativeSuccess;
        init.tau.push_back(atom.mdpde[1]); init.atoms.push_back(atom);
    }
    init.tau[0]=0.; init.atoms[0].original_tau=0.; init.atoms[0].mdpde[1]=0.; init.atoms[0].reason="invalid-width";
    // Seed fallback and calibration deliberately have different donor contracts.
    init.atoms[1].initial_stage_fit_qualified=false;
    const auto before=init.atoms;
    d::ApplyJointSeedFallback(init);
    EXPECT_TRUE(init.valid); EXPECT_EQ(init.tau[0],1.12); EXPECT_EQ(init.atoms[0].donor_count,23u);
    EXPECT_EQ(init.atoms[0].seed_source,"median-fallback"); EXPECT_EQ(init.atoms[0].used_tau,1.12);
    EXPECT_EQ(init.atoms[0].original_tau,0.); EXPECT_EQ(init.atoms[0].mdpde,before[0].mdpde);
    EXPECT_EQ(init.atoms[0].ols,before[0].ols);
    const auto calibration=n::CalibrateGlobalWidth(init,"fallback");
    EXPECT_EQ(calibration.qualified_donor_count,22u); EXPECT_EQ(calibration.fallback_count,1u);
    EXPECT_EQ(calibration.unavailable_reason,"qualified-initial-stage-target-coverage-below-0.95");
    init.atoms.back().target=false; init.tau.back()=0.;
    d::ApplyJointSeedFallback(init);
    EXPECT_TRUE(init.valid); EXPECT_EQ(init.tau.back(),0.);
    for(auto & atom:init.atoms) atom.reason="invalid-width";
    std::fill(init.tau.begin(),init.tau.end(),0.);
    d::ApplyJointSeedFallback(init); EXPECT_FALSE(init.valid);
}
TEST(InitialStageTest, WidthTruthRecoveryAndQualificationAtTwoVoxelSpacings)
{
    for(double spacing:{.2,.4})
    {
        auto model=joint_component_test::InitialStageSyntheticModel(5.5);
        auto map=joint_component_test::InitialStageSyntheticMap(*model,spacing);
        const auto problem=BuildPureHaloSupportExclusionProblem(*map,*model);
        model->EditAnalysis().InitializeFromSelection();
        FitOptions options; options.quiet_mode=true;
        const auto workset=d::BuildInitialStageWorkset(*model,problem);
        const auto init=d::RunInitialStage(*map,*model,workset,options);
        const auto calibration=n::CalibrateGlobalWidth(init,"synthetic-truth-regression");
        EXPECT_EQ(calibration.status,"Available"); EXPECT_EQ(calibration.qualified_donor_count,24u);
        EXPECT_EQ(calibration.qualified_donor_coverage,1.); EXPECT_EQ(calibration.fallback_count,0u);
        double square=0.;
        for(std::size_t i=0;i<init.atoms.size();++i)
        {
            const double error=init.atoms[i].mdpde[1]-(.35+.005*static_cast<double>(i));
            square+=error*error;
        }
        // Truth-based bounds cover interpolation error, not just solver self-consistency.
        EXPECT_LT(std::sqrt(square/24.),spacing==.2 ? .002 : .02);
    }
}

TEST(InitialStageTest, RawSamplingPersistenceAndGridInvalidationIndependence)
{
    auto model=joint_component_test::InitialStageSyntheticModel(5.5);
    auto map=joint_component_test::InitialStageSyntheticMap(*model,.4);
    FitOptions options; options.quiet_mode=true;
    d::RunPotentialFittingWorkflowForSelection(*map,*model,options);
    ASSERT_TRUE(model->GetAnalysisView().GetJointResult());
    const auto original=AtomLocalPotentialView::For(*model->FindAtomPtr(1));
    const LocalPotentialSampleList raw=original.GetSamplingEntriesView();
    ASSERT_EQ(raw.size(),200u); ASSERT_TRUE(original.GetPostFitPeeling());
    ASSERT_TRUE(original.GetPostJointLocalResult());
    ASSERT_TRUE(original.GetGroupMemberResult());
    const auto peeling_count=original.GetPostFitPeeling()->samples.size();
    const auto retained=original.GetPostJointLocalResult()->effective_sample_count;
    command_test::ScopedTempDir directory("joint-initial-stage");
    DataRepository repository(directory.path()/"samples.sqlite");
    repository.SaveModel(*model,"F0"); auto restored=repository.LoadModel("F0");
    EXPECT_EQ(stage_result_io::Encode(*restored),stage_result_io::Encode(*model));
    const auto view=AtomLocalPotentialView::For(*restored->FindAtomPtr(1));
    const LocalPotentialSampleList saved=view.GetSamplingEntriesView();
    ASSERT_EQ(saved.size(),raw.size());
    for (std::size_t i=0;i<raw.size();++i)
    {
        EXPECT_EQ(saved[i].point.position,raw[i].point.position);
        EXPECT_DOUBLE_EQ(saved[i].point.distance,raw[i].point.distance);
        EXPECT_DOUBLE_EQ(saved[i].response,raw[i].response);
        EXPECT_EQ(saved[i].point.is_selected,raw[i].point.is_selected);
    }
    restored->EditAnalysis().SetAtomLocalSamplingEntries(*restored->FindAtomPtr(1),
        LocalPotentialSampleList{saved.front()});
    ASSERT_TRUE(view.GetPostFitPeeling()); ASSERT_TRUE(view.GetPostJointLocalResult());
    EXPECT_EQ(view.GetPostFitPeeling()->samples.size(),peeling_count);
    EXPECT_EQ(view.GetPostJointLocalResult()->effective_sample_count,retained);
    EXPECT_EQ(view.GetSamplingEntriesView().size(),1u);
}


} // namespace
