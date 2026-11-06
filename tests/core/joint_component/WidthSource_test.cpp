#include <gtest/gtest.h>
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/FittingWorkset.hpp"
#include "core/detail/joint_component/WidthCalibration.hpp"
#include "support/FirstStageSamplingFixture.hpp"
#include "support/ResearchWidthInitialization.hpp"
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <boost/json.hpp>
#include <filesystem>
#include <fstream>
#include <iterator>

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
    const auto init=d::PrepareJointWidthInitialization(*map,*model,workset,options);
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
    const auto fit=FitJointComponentsInitialized(problem,init);
    ASSERT_TRUE(fit.assembled_state);
    const auto & expected=baseline.at("joint").as_object();
    EXPECT_EQ(fit.objective,j::value_to<double>(expected.at("objective")));
    for(std::size_t i=0;i<init.atoms.size();++i)
        EXPECT_EQ(fit.assembled_state->tau[i],j::value_to<double>(expected.at("tau").as_array()[i]));
}
}
