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
    initialization.b.resize(100);
    for(std::size_t i=0;i<initialization.b.size();++i)
    {
        initialization.b[i]=static_cast<double>(i+1);
        JointInitializationAtom atom;
        atom.id=std::to_string(i+1);
        atom.target=true;
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

JointInitialization QualifiedInitialization(std::size_t count)
{
    JointInitialization initialization;
    initialization.valid=true;
    for(std::size_t i=0;i<count;++i)
    {
        initialization.b.push_back(static_cast<double>(i+1));
        JointInitializationAtom atom;
        atom.id=std::to_string(i+1); atom.target=true;
        atom.mdpde={1.0,static_cast<double>(i+1),0.0};
        atom.original_b=initialization.b.back(); atom.seed_source="fitted";
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
    for(auto & atom:input.atoms) atom.mdpde[1]=1.0;
    input.atoms[4].mdpde[1]=5.0;
    input.atoms[5].mdpde[1]=std::numeric_limits<double>::infinity();

    const auto calibration=joint::CalibrateGlobalWidth(input,"unit-case");
    ASSERT_EQ(calibration.status,"Available");
    EXPECT_EQ(calibration.target_count,99u);
    EXPECT_EQ(calibration.qualified_donor_count,96u);
    EXPECT_EQ(calibration.native_success_count,97u);
    EXPECT_EQ(calibration.refined_success_count,1u);
    EXPECT_EQ(calibration.unqualified_count,1u);
    EXPECT_EQ(calibration.invalid_width_count,1u);
    EXPECT_EQ(calibration.fallback_count,1u);
    EXPECT_NEAR(calibration.qualified_donor_coverage,96.0/99.0,1e-15);
    EXPECT_DOUBLE_EQ(*calibration.p95,1.0);
    EXPECT_EQ(calibration.excluded_donor_reasons.at("unqualified-first-stage"),1u);
    EXPECT_EQ(calibration.excluded_donor_reasons.at("median-fallback"),1u);
    EXPECT_EQ(calibration.excluded_donor_reasons.at("non-finite-width"),1u);
}

TEST(WidthCalibrationTest, UsesDeterministicType7QuantilesForSmallAndDuplicateSamples)
{
    auto input=QualifiedInitialization(20);
    const auto first=joint::CalibrateGlobalWidth(input,"duplicates");
    std::reverse(input.atoms.begin(),input.atoms.end());
    std::reverse(input.b.begin(),input.b.end());
    const auto second=joint::CalibrateGlobalWidth(input,"duplicates");
    ASSERT_EQ(first.status,"Available");
    ASSERT_EQ(second.status,"Available");
    EXPECT_DOUBLE_EQ(*first.p95,19.05);
    EXPECT_DOUBLE_EQ(*second.p95,19.05);

    input=QualifiedInitialization(20);
    for(auto & atom:input.atoms) atom.mdpde[1]=1.0;
    input.atoms[18].mdpde[1]=2.0; input.atoms[19].mdpde[1]=3.0;
    const auto duplicates=joint::CalibrateGlobalWidth(input,"duplicates");
    ASSERT_EQ(duplicates.status,"Available");
    EXPECT_DOUBLE_EQ(*duplicates.p95,2.05);
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

    const auto projected=joint::MakeBoundedInitialization(
        input,calibration,joint::GlobalWidthBoundMultiplier::TwoP95);
    ASSERT_TRUE(projected.available);
    EXPECT_EQ(projected.projection_count,1u);
    EXPECT_DOUBLE_EQ(projected.maximum_projection_magnitude,1000.0-*calibration.bmax2);
    EXPECT_DOUBLE_EQ(projected.widths[5],*calibration.bmax2);
    EXPECT_DOUBLE_EQ(input.b[5],original.b[5]);
    EXPECT_DOUBLE_EQ(input.atoms[5].mdpde[1],original.atoms[5].mdpde[1]);
    EXPECT_DOUBLE_EQ(projected.widths[6],input.b[6]);
    ASSERT_EQ(projected.projections.size(),100u);
    EXPECT_DOUBLE_EQ(projected.projections[5].original_first_stage_width.value(),6.0);
}

TEST(WidthCalibrationTest, RejectsMalformedInitializationAndUnavailableCalibration)
{
    auto input=ExampleInitialization();
    input.b.pop_back();
    EXPECT_EQ(joint::CalibrateGlobalWidth(input,"mismatch").unavailable_reason,
        "initialization-record-size-mismatch");
    auto valid=ExampleInitialization();
    auto unavailable=joint::CalibrateGlobalWidth(JointInitialization{},"none");
    const auto projected=joint::MakeBoundedInitialization(
        valid,unavailable,joint::GlobalWidthBoundMultiplier::FourP95);
    EXPECT_FALSE(projected.available);
    EXPECT_EQ(projected.reason,"no-qualified-first-stage-target-widths");
}

TEST(WidthCalibrationTest, RejectsSmallOrPoorlyCoveredDonorPopulationsAndInvalidTargetStarts)
{
    auto too_small=QualifiedInitialization(19);
    EXPECT_EQ(joint::CalibrateGlobalWidth(too_small,"small").unavailable_reason,
        "too-few-qualified-first-stage-target-widths");

    auto low_coverage=QualifiedInitialization(22);
    low_coverage.atoms[0].first_stage_fit_qualified=false;
    low_coverage.atoms[1].seed_source="median-fallback";
    EXPECT_EQ(joint::CalibrateGlobalWidth(low_coverage,"coverage").unavailable_reason,
        "qualified-first-stage-target-coverage-below-0.95");

    auto invalid_start=QualifiedInitialization(20);
    invalid_start.b[0]=std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(joint::CalibrateGlobalWidth(invalid_start,"invalid-start").unavailable_reason,
        "invalid-target-initial-width");
}

} // namespace
