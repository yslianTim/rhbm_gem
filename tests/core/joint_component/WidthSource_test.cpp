#include <gtest/gtest.h>
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/FittingWorkset.hpp"
#include "core/detail/joint_component/WidthCalibration.hpp"
#include "support/FirstStageSamplingFixture.hpp"
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
