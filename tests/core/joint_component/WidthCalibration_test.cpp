#include <gtest/gtest.h>

#include "core/detail/joint_component/WidthCalibration.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {
using namespace rhbm_gem;
namespace joint=rhbm_gem::core::joint_component;

JointInitialization ExampleInitialization()
{
    JointInitialization initialization;
    initialization.valid=true;
    initialization.b={1.0,2.0,3.0,4.0,5.0,100.0,1.5};
    for(std::size_t i=0;i<initialization.b.size();++i)
    {
        JointInitializationAtom atom;
        atom.id=std::to_string(i+1);
        atom.target=i<6;
        atom.mdpde={1.0,initialization.b[i],0.0};
        atom.original_b=initialization.b[i];
        atom.seed_source="fitted";
        atom.native_status=static_cast<int>(RHBMEstimationStatus::SUCCESS);
        atom.first_stage_fit_qualified=true;
        atom.first_stage_qualification=RHBMSolveQualification::NativeSuccess;
        atom.reason="valid-width";
        initialization.atoms.push_back(std::move(atom));
    }
    return initialization;
}

TEST(WidthCalibrationTest, SelectsOnlyQualifiedTargetFirstStageWidths)
{
    auto input=ExampleInitialization();
    input.atoms[1].target=false;
    input.atoms[2].first_stage_fit_qualified=false;
    input.atoms[3].seed_source="median-fallback";
    input.atoms[4].native_status=static_cast<int>(RHBMEstimationStatus::MAX_ITERATIONS_REACHED);
    input.atoms[4].first_stage_qualification=RHBMSolveQualification::RefinedSuccess;
    input.atoms[4].first_stage_recovery=RHBMFailedOnlyRecoveryDiagnostics{};
    input.atoms[4].first_stage_recovery->accepted=true;
    input.atoms[5].mdpde[1]=std::numeric_limits<double>::infinity();

    const auto calibration=joint::CalibrateGlobalWidth(input,"unit-case");
    ASSERT_EQ(calibration.status,"Available");
    EXPECT_EQ(calibration.target_count,5u);
    EXPECT_EQ(calibration.qualified_donor_count,2u);
    EXPECT_EQ(calibration.native_success_count,3u);
    EXPECT_EQ(calibration.refined_success_count,1u);
    EXPECT_EQ(calibration.unqualified_count,1u);
    EXPECT_EQ(calibration.invalid_width_count,1u);
    EXPECT_EQ(calibration.fallback_count,1u);
    EXPECT_DOUBLE_EQ(*calibration.p95,4.8);
    EXPECT_EQ(calibration.excluded_donor_reasons.at("unqualified-first-stage"),1u);
    EXPECT_EQ(calibration.excluded_donor_reasons.at("median-fallback"),1u);
    EXPECT_EQ(calibration.excluded_donor_reasons.at("non-finite-width"),1u);
}

TEST(WidthCalibrationTest, UsesDeterministicType7QuantilesForSmallAndDuplicateSamples)
{
    auto input=ExampleInitialization();
    input.atoms.resize(4);
    input.b.resize(4);
    input.atoms[0].mdpde[1]=1.0;
    input.atoms[1].mdpde[1]=1.0;
    input.atoms[2].mdpde[1]=2.0;
    input.atoms[3].mdpde[1]=3.0;
    const auto first=joint::CalibrateGlobalWidth(input,"duplicates");
    std::reverse(input.atoms.begin(),input.atoms.end());
    std::reverse(input.b.begin(),input.b.end());
    const auto second=joint::CalibrateGlobalWidth(input,"duplicates");
    ASSERT_EQ(first.status,"Available");
    ASSERT_EQ(second.status,"Available");
    EXPECT_DOUBLE_EQ(*first.p95,2.85);
    EXPECT_DOUBLE_EQ(*second.p95,2.85);

    input.atoms.resize(1);
    input.b.resize(1);
    input.atoms[0].mdpde[1]=2.5;
    const auto singleton=joint::CalibrateGlobalWidth(input,"one-donor");
    ASSERT_EQ(singleton.status,"Available");
    EXPECT_DOUBLE_EQ(*singleton.p95,2.5);
}

TEST(WidthCalibrationTest, ReportsEmptyAndOverflowAsCalibrationUnavailable)
{
    auto input=ExampleInitialization();
    input.atoms.clear();
    input.b.clear();
    const auto empty=joint::CalibrateGlobalWidth(input,"empty");
    EXPECT_EQ(empty.status,"CalibrationUnavailable");
    EXPECT_EQ(empty.unavailable_reason,"no-qualified-first-stage-target-widths");
    EXPECT_FALSE(empty.p95);

    input=ExampleInitialization();
    for(auto & atom:input.atoms) if(atom.target) atom.mdpde[1]=std::numeric_limits<double>::max();
    const auto overflow=joint::CalibrateGlobalWidth(input,"overflow");
    EXPECT_EQ(overflow.status,"CalibrationUnavailable");
    EXPECT_EQ(overflow.unavailable_reason,"global-width-bound-overflow-or-invalid-log-bound");
    EXPECT_FALSE(overflow.bmax4);
}

TEST(WidthCalibrationTest, ProducesOneGlobalBoundAndProjectsWithoutChangingStoredInitialization)
{
    auto input=ExampleInitialization();
    input.b[5]=1000.0;
    const auto original=input;
    const auto calibration=joint::CalibrateGlobalWidth(input,"global");
    ASSERT_EQ(calibration.status,"Available");
    EXPECT_DOUBLE_EQ(*calibration.bmax2,2.0*(*calibration.p95));
    EXPECT_DOUBLE_EQ(*calibration.bmax4,4.0*(*calibration.p95));
    EXPECT_EQ(calibration.initial_violations2,1u);
    EXPECT_EQ(calibration.initial_violations4,1u);

    const auto projected=joint::MakeBoundedInitializationForResearch(
        input,calibration,joint::GlobalWidthBoundMultiplier::TwoP95);
    ASSERT_TRUE(projected.available);
    EXPECT_EQ(projected.projection_count,1u);
    EXPECT_DOUBLE_EQ(projected.maximum_projection_magnitude,1000.0-*calibration.bmax2);
    EXPECT_DOUBLE_EQ(projected.widths[5],*calibration.bmax2);
    EXPECT_DOUBLE_EQ(input.b[5],original.b[5]);
    EXPECT_DOUBLE_EQ(input.atoms[5].mdpde[1],original.atoms[5].mdpde[1]);
    EXPECT_DOUBLE_EQ(projected.widths[6],input.b[6]);
    ASSERT_EQ(projected.projections.size(),6u);
    EXPECT_DOUBLE_EQ(projected.projections[5].original_first_stage_width.value(),100.0);
}

TEST(WidthCalibrationTest, RejectsMalformedInitializationAndUnavailableCalibration)
{
    auto input=ExampleInitialization();
    input.b.pop_back();
    EXPECT_EQ(joint::CalibrateGlobalWidth(input,"mismatch").unavailable_reason,
        "initialization-record-size-mismatch");
    auto valid=ExampleInitialization();
    auto unavailable=joint::CalibrateGlobalWidth(JointInitialization{},"none");
    const auto projected=joint::MakeBoundedInitializationForResearch(
        valid,unavailable,joint::GlobalWidthBoundMultiplier::FourP95);
    EXPECT_FALSE(projected.available);
    EXPECT_EQ(projected.reason,"no-qualified-first-stage-target-widths");
}

} // namespace
