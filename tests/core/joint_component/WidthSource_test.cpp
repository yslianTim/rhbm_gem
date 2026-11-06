#include <gtest/gtest.h>
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/FittingWorkset.hpp"
#include "core/detail/PotentialFittingWorkflow.hpp"
#include "core/detail/joint_component/WidthCalibration.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/FirstStageSamplingFixture.hpp"
#include "support/ResearchWidthInitialization.hpp"
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include "data/io/detail/JointResultJson.hpp"
#include <boost/json.hpp>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>

namespace {
using namespace rhbm_gem;
using namespace rhbm_gem::core;
namespace d=rhbm_gem::core::detail;
namespace n=rhbm_gem::core::joint_component;
namespace j=boost::json;

JointInitialization SourceFixture(std::size_t count=24)
{
    JointInitialization out; out.valid=true;
    for(std::size_t i=0;i<count;++i)
    {
        JointInitializationAtom atom; atom.id=std::to_string(i); atom.target=true;
        atom.ols={1.,.5+.01*static_cast<double>(i),0.};
        atom.mdpde={1.,1.+.01*static_cast<double>(i),0.};
        atom.original_tau=atom.used_tau=atom.mdpde[1]; atom.reason="valid-width"; atom.seed_source="fitted";
        atom.native_status=static_cast<int>(RHBMEstimationStatus::SUCCESS);
        atom.first_stage_fit_qualified=true; atom.first_stage_qualification=RHBMSolveQualification::NativeSuccess;
        out.tau.push_back(atom.mdpde[1]); out.atoms.push_back(atom);
    }
    return out;
}

TEST(WidthSourceTest, SelectsBothSourcesDeterministicallyWithoutChangingEvidenceOrRoles)
{
    auto input=SourceFixture(); input.atoms.back().target=false;
    const auto before=input;
    const auto m=width_research::SelectWidths(input,width_research::WidthSource::MDPDE);
    const auto o=width_research::SelectWidths(input,width_research::WidthSource::OLS);
    EXPECT_EQ(m.tau,input.tau);
    EXPECT_EQ(o.tau,width_research::SelectWidths(input,width_research::WidthSource::OLS).tau);
    EXPECT_TRUE(m.available); EXPECT_TRUE(o.available);
    for(std::size_t i=0;i<input.atoms.size();++i)
    {
        EXPECT_EQ(o.tau[i],before.atoms[i].ols[1]);
        EXPECT_EQ(input.atoms[i].mdpde,before.atoms[i].mdpde);
        EXPECT_EQ(input.atoms[i].ols,before.atoms[i].ols);
        EXPECT_EQ(input.atoms[i].original_tau,before.atoms[i].original_tau);
        EXPECT_EQ(input.atoms[i].used_tau,before.atoms[i].used_tau);
        EXPECT_EQ(input.atoms[i].seed_source,before.atoms[i].seed_source);
        EXPECT_EQ(o.atoms[i].target,before.atoms[i].target);
    }
    EXPECT_FALSE(o.atoms.back().calibration_donor);
    EXPECT_EQ(o.atoms.back().donor_exclusion,"halo");
    const auto bridge=width_research::ProjectionInput(input,o);
    EXPECT_EQ(bridge.atoms[0].mdpde,input.atoms[0].mdpde);
    EXPECT_EQ(bridge.atoms[0].original_tau,input.atoms[0].ols[1]);
}

TEST(WidthSourceTest, UsesOnlySameSourceMedianForInvalidWidthsAndExcludesFallbackDonors)
{
    auto input=SourceFixture();
    input.atoms[0].ols[1]=std::numeric_limits<double>::quiet_NaN();
    input.atoms[1].ols[1]=0.; input.atoms[2].ols[1]=-1.;
    input.atoms[3].ols[0]=std::numeric_limits<double>::infinity();
    const auto o=width_research::SelectWidths(input,width_research::WidthSource::OLS);
    ASSERT_TRUE(o.fallback_median); EXPECT_EQ(o.fallback_donor_count,21u);
    EXPECT_EQ(o.fallback_count,3u); EXPECT_TRUE(o.available);
    for(std::size_t i=0;i<3;++i)
    { EXPECT_EQ(o.tau[i],*o.fallback_median); EXPECT_FALSE(o.atoms[i].calibration_donor); }
    EXPECT_FALSE(o.atoms[3].parameters_valid); EXPECT_TRUE(o.atoms[3].width_valid);
    input.atoms[0].mdpde[1]=0.; input.atoms[0].original_tau=0.;
    input.atoms[0].seed_source="median-fallback"; input.tau[0]=1.12;
    const auto m=width_research::SelectWidths(input,width_research::WidthSource::MDPDE);
    EXPECT_EQ(m.tau,input.tau); EXPECT_FALSE(m.atoms[0].width_valid);
    EXPECT_FALSE(m.atoms[0].calibration_donor); EXPECT_EQ(m.fallback_count,1u);
}

TEST(WidthSourceTest, ReportsUnavailableWithoutDonorsAndDoesNotFallbackHalos)
{
    auto input=SourceFixture();
    for(auto & atom:input.atoms) atom.ols[1]=0.;
    const auto o=width_research::SelectWidths(input,width_research::WidthSource::OLS);
    EXPECT_FALSE(o.available); EXPECT_FALSE(o.fallback_median); EXPECT_EQ(o.fallback_count,0u);
    input=SourceFixture(); input.atoms.back().target=false; input.atoms.back().ols[1]=0.;
    const auto halo=width_research::SelectWidths(input,width_research::WidthSource::OLS);
    EXPECT_TRUE(halo.available); EXPECT_EQ(halo.tau.back(),0.); EXPECT_EQ(halo.fallback_count,0u);
    input.atoms[0].first_stage_fit_qualified=false;
    const auto gate=width_research::SelectWidths(input,width_research::WidthSource::OLS);
    EXPECT_TRUE(gate.atoms[0].parameters_valid); EXPECT_FALSE(gate.atoms[0].calibration_donor);
    EXPECT_EQ(gate.atoms[0].donor_exclusion,"unqualified-first-stage");
}

TEST(WidthSourceTest, KeepsSeedAndCalibrationOrthogonalAndProjectsOnlyWorkingValues)
{
    const auto input=SourceFixture();
    const auto m=width_research::SelectWidths(input,width_research::WidthSource::MDPDE);
    const auto o=width_research::SelectWidths(input,width_research::WidthSource::OLS);
    const auto cm=width_research::CalibrateWidths(input,m),co=width_research::CalibrateWidths(input,o);
    ASSERT_EQ(cm.status,"Available"); ASSERT_EQ(co.status,"Available");
    EXPECT_EQ(cm.tau_max4,n::CalibrateGlobalWidth(input,"production").tau_max4);
    EXPECT_EQ(cm.p95,width_research::CalibrateWidths(input,m,true).p95);
    EXPECT_EQ(co.p95,width_research::CalibrateWidths(input,o,true).p95);
    const auto mm=n::MakeBoundedInitialization(width_research::ProjectionInput(input,m),cm);
    const auto om=n::MakeBoundedInitialization(width_research::ProjectionInput(input,o),cm);
    const auto mo=n::MakeBoundedInitialization(width_research::ProjectionInput(input,m),co);
    const auto oo=n::MakeBoundedInitialization(width_research::ProjectionInput(input,o),co);
    EXPECT_EQ(mm.tau_max,om.tau_max); EXPECT_EQ(mo.tau_max,oo.tau_max);
    EXPECT_EQ(mm.tau_values,mo.tau_values); EXPECT_EQ(om.tau_values,oo.tau_values);
    auto extreme=o; extreme.tau[0]=10.; extreme.atoms[0].used=10.;
    const auto projected=n::MakeBoundedInitialization(width_research::ProjectionInput(input,extreme),cm);
    EXPECT_EQ(projected.projection_count,1u); EXPECT_EQ(projected.tau_values[0],*cm.tau_max4);
    EXPECT_EQ(extreme.tau[0],10.); EXPECT_EQ(input.atoms[0].ols[1],.5);
}

TEST(WidthSourceTest, RetainsSourceSpecificDonorsMinimumCoverageAndMatchedCohort)
{
    auto input=SourceFixture(); input.atoms[0].ols[1]=0.;
    auto o=width_research::SelectWidths(input,width_research::WidthSource::OLS);
    auto m=width_research::SelectWidths(input,width_research::WidthSource::MDPDE);
    const auto co=width_research::CalibrateWidths(input,o);
    EXPECT_EQ(co.qualified_donor_count,23u); EXPECT_EQ(co.fallback_count,1u);
    EXPECT_EQ(co.qualified_donor_coverage,23./24.); EXPECT_EQ(co.status,"Available");
    EXPECT_EQ(width_research::CalibrateWidths(input,m).qualified_donor_count,24u);
    EXPECT_EQ(width_research::CalibrateWidths(input,m,true).qualified_donor_count,23u);
    input.atoms[1].ols[1]=0.; o=width_research::SelectWidths(input,width_research::WidthSource::OLS);
    EXPECT_EQ(width_research::CalibrateWidths(input,o).unavailable_reason,
        "qualified-first-stage-target-coverage-below-0.95");
    input=SourceFixture(19); o=width_research::SelectWidths(input,width_research::WidthSource::OLS);
    EXPECT_EQ(width_research::CalibrateWidths(input,o).unavailable_reason,
        "too-few-qualified-first-stage-target-widths");
}

TEST(WidthSourceTest, KeepsAsymmetricAndJointInvaliditySeparateFromSeedAvailability)
{
    auto input=SourceFixture();
    for(std::size_t i:{0u,2u})
    {
        input.atoms[i].mdpde[1]=0.; input.atoms[i].original_tau=0.;
        input.atoms[i].reason="invalid-width"; input.tau[i]=0.;
    }
    input.atoms[1].ols[1]=0.; input.atoms[2].ols[1]=0.;
    d::ApplyJointSeedFallback(input);
    const auto m=width_research::SelectWidths(input,width_research::WidthSource::MDPDE);
    const auto o=width_research::SelectWidths(input,width_research::WidthSource::OLS);
    EXPECT_TRUE(o.atoms[0].width_valid); EXPECT_FALSE(m.atoms[0].width_valid);
    EXPECT_FALSE(o.atoms[1].width_valid); EXPECT_TRUE(m.atoms[1].width_valid);
    EXPECT_FALSE(o.atoms[2].width_valid); EXPECT_FALSE(m.atoms[2].width_valid);
    EXPECT_TRUE(o.atoms[3].width_valid); EXPECT_TRUE(m.atoms[3].width_valid);
    EXPECT_TRUE(o.available); EXPECT_TRUE(m.available);
    EXPECT_TRUE(o.atoms[0].calibration_donor); EXPECT_FALSE(m.atoms[0].calibration_donor);
    EXPECT_FALSE(o.atoms[1].calibration_donor); EXPECT_TRUE(m.atoms[1].calibration_donor);
    EXPECT_EQ(m.fallback_count,2u); EXPECT_EQ(o.fallback_count,2u);
    EXPECT_EQ(m.tau,input.tau);
}

TEST(WidthSourceTest, LocksProductionFirstStageAndMmBaseline)
{
    const auto path=std::filesystem::path(__FILE__).parent_path().parent_path().parent_path().parent_path()/
        "docs/research/first-stage-sampling/baseline.json";
    std::ifstream stream(path);
    ASSERT_TRUE(stream.good());
    j::parse_options parse_options; parse_options.numbers=j::number_precision::precise;
    const auto baseline=j::parse(std::string(std::istreambuf_iterator<char>(stream),{}),{},parse_options).as_object();
    auto model=second_stage_test::FirstStageStudyModel(5.5,0);
    auto map=second_stage_test::FirstStageStudyMap(*model,.4,1e-4,0);
    const auto problem=BuildPureHaloSupportExclusionProblem(*map,*model);
    const auto workset=d::MakeJointFittingWorkset(*model,problem);
    model->EditAnalysis().InitializeFromSelection();
    FitOptions options; options.quiet_mode=true;
    std::size_t fits=0;
    d::FirstStageObserverForTesting()=[&](int,std::string_view) { ++fits; };
    const auto init=d::PrepareJointWidthInitialization(*map,*model,workset,options);
    d::FirstStageObserverForTesting()={};
    EXPECT_EQ(fits,24u);
    d::FirstStageTimingReport timing;
    const auto timed=d::RunJointFirstStageInitializationFromPreparedSamples(*model,workset,options,&timing);
    EXPECT_GT(timing.alpha_training_seconds,0.); EXPECT_GT(timing.fitting_seconds,0.);
    EXPECT_EQ(timed.tau,init.tau);
    for(std::size_t i=0;i<init.atoms.size();++i)
    {
        EXPECT_EQ(timed.atoms[i].ols,init.atoms[i].ols);
        EXPECT_EQ(timed.atoms[i].mdpde,init.atoms[i].mdpde);
        EXPECT_EQ(timed.atoms[i].alpha,init.atoms[i].alpha);
        EXPECT_EQ(timed.atoms[i].first_stage_qualification,init.atoms[i].first_stage_qualification);
    }
    ASSERT_TRUE(init.valid);
    ASSERT_EQ(init.atoms.size(),24u);
    for(std::size_t i=0;i<init.atoms.size();++i)
    {
        const auto & atom=init.atoms[i];
        const auto & expected=baseline.at("atoms").as_array()[i].as_object();
        EXPECT_EQ(atom.id,j::value_to<std::string>(expected.at("id")));
        EXPECT_EQ(atom.alpha,j::value_to<double>(expected.at("alpha")));
        for(std::size_t k=0;k<3;++k)
        {
            EXPECT_EQ(atom.ols[k],j::value_to<double>(expected.at("ols").as_array()[k]));
            EXPECT_EQ(atom.mdpde[k],j::value_to<double>(expected.at("mdpde").as_array()[k]));
        }
        EXPECT_TRUE(atom.target);
        EXPECT_TRUE(atom.first_stage_fit_qualified);
        EXPECT_EQ(atom.first_stage_qualification,RHBMSolveQualification::NativeSuccess);
        EXPECT_EQ(atom.seed_source,"fitted");
        EXPECT_EQ(atom.original_tau,atom.mdpde[1]);
        EXPECT_EQ(atom.used_tau,atom.mdpde[1]);
        EXPECT_EQ(init.tau[i],atom.mdpde[1]);
    }
    const auto calibration=n::CalibrateGlobalWidth(init,"baseline");
    EXPECT_EQ(calibration.qualified_donor_count,24u);
    EXPECT_EQ(calibration.qualified_donor_coverage,1.);
    ASSERT_TRUE(calibration.tau_max4);
    EXPECT_EQ(*calibration.tau_max4,j::value_to<double>(baseline.at("calibration").as_object().at("bound")));
    const auto bounded=n::MakeBoundedInitialization(init,calibration);
    ASSERT_TRUE(bounded.available);
    EXPECT_EQ(bounded.tau_values,init.tau);
    EXPECT_EQ(bounded.projection_count,0u);
    const auto selected=width_research::SelectWidths(init,width_research::WidthSource::MDPDE);
    EXPECT_EQ(selected.tau,init.tau);
    EXPECT_EQ(width_research::CalibrateWidths(init,selected).tau_max4,calibration.tau_max4);
    EXPECT_EQ(n::MakeBoundedInitialization(width_research::ProjectionInput(init,selected),calibration).tau_values,
        bounded.tau_values);
    const auto fit=FitJointComponentsInitialized(problem,init);
    ASSERT_TRUE(fit.assembled_state);
    const auto & expected=baseline.at("joint").as_object();
    EXPECT_EQ(fit.objective,j::value_to<double>(expected.at("objective")));
    for(std::size_t i=0;i<init.atoms.size();++i)
        EXPECT_EQ(fit.assembled_state->tau[i],j::value_to<double>(expected.at("tau").as_array()[i]));

    const auto ols=width_research::SelectWidths(init,width_research::WidthSource::OLS);
    const auto ols_calibration=width_research::CalibrateWidths(init,ols);
    ASSERT_EQ(ols_calibration.status,"Available");
    // Reuse the exact immutable problem and estimates for the full factorial.
    for(const auto * seeds:{&selected,&ols}) for(const auto * bound:{&calibration,&ols_calibration})
    {
        const auto start=n::MakeBoundedInitialization(width_research::ProjectionInput(init,*seeds),*bound);
        n::BcdExecutionPolicy policy; policy.eta_upper_bound=bound->eta_max4;
        const auto candidate=n::FitBcdComponents(problem,start.tau_values,policy);
        ASSERT_EQ(candidate.RuntimeConvergence(),JointCheckStatus::Passed);
        ASSERT_TRUE(candidate.assembled_state); ASSERT_TRUE(candidate.objective);
        // Solver qualification alone must not forge Production calibration provenance.
        EXPECT_THROW(joint_result_io::Encode(CaptureJointAnalysisResult(candidate)),std::invalid_argument);
        EXPECT_NEAR(*candidate.objective,*fit.objective,1e-17);
        for(std::size_t i=0;i<init.atoms.size();++i)
            EXPECT_NEAR(candidate.assembled_state->tau[i],fit.assembled_state->tau[i],3e-9);
    }
    d::JointPostprocessingMetrics metrics;
    d::RunJointPostprocessing(*map,*model,problem,CaptureJointAnalysisResult(fit),options,metrics);
    EXPECT_EQ(metrics.available_peeling_sample_count,metrics.grid_peeling_sample_count);
    const auto snapshot=model->GetAnalysisView().GetJointResult();
    ASSERT_TRUE(snapshot);
    const auto restored=joint_result_io::Decode(joint_result_io::Encode(*snapshot));
    EXPECT_EQ(restored.objective,fit.objective);
    EXPECT_EQ(restored.initialization.tau,init.tau);
}
}
