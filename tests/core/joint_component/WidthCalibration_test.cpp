#include <gtest/gtest.h>

#include "core/estimator/joint_component/WidthCalibration.hpp"

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
    initialization.tau.resize(100);
    for(std::size_t i=0;i<initialization.tau.size();++i)
    {
        initialization.tau[i]=static_cast<double>(i+1);
        JointInitializationAtom atom;
        atom.id=std::to_string(i+1);
        atom.target=true;
        atom.mdpde={1.0,initialization.tau[i],0.0};
        atom.original_tau=initialization.tau[i];
        atom.seed_source="fitted";
        atom.native_status=static_cast<int>(RHBMEstimationStatus::SUCCESS);
        atom.initial_stage_fit_qualified=true;
        atom.initial_stage_qualification=RHBMSolveQualification::NativeSuccess;
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
        initialization.tau.push_back(static_cast<double>(i+1));
        JointInitializationAtom atom;
        atom.id=std::to_string(i+1); atom.target=true;
        atom.mdpde={1.0,static_cast<double>(i+1),0.0};
        atom.original_tau=initialization.tau.back(); atom.seed_source="fitted";
        atom.native_status=static_cast<int>(RHBMEstimationStatus::SUCCESS);
        atom.initial_stage_fit_qualified=true;
        atom.initial_stage_qualification=RHBMSolveQualification::NativeSuccess;
        atom.reason="valid-width";
        initialization.atoms.push_back(std::move(atom));
    }
    return initialization;
}

TEST(WidthCalibrationTest, SelectsOnlyQualifiedTargetInitialStageWidths)
{
    auto input=ExampleInitialization();
    input.atoms[1].target=false;
    input.atoms[2].initial_stage_fit_qualified=false;
    input.atoms[3].seed_source="median-fallback";
    input.atoms[4].native_status=static_cast<int>(RHBMEstimationStatus::MAX_ITERATIONS_REACHED);
    input.atoms[4].initial_stage_qualification=RHBMSolveQualification::RefinedSuccess;
    input.atoms[4].initial_stage_recovery=RHBMFailedOnlyRecoveryDiagnostics{};
    input.atoms[4].initial_stage_recovery->accepted=true;
    for(auto & atom:input.atoms) atom.mdpde[1]=1.0;
    input.atoms[4].mdpde[1]=5.0;
    input.atoms[5].mdpde[1]=std::numeric_limits<double>::infinity();
    input.atoms[6].reason="initialization-exception: test failure";

    const auto calibration=joint::CalibrateGlobalWidth(input,"unit-case");
    ASSERT_EQ(calibration.status,"Available");
    EXPECT_EQ(calibration.target_count,99u);
    EXPECT_EQ(calibration.qualified_donor_count,95u);
    EXPECT_EQ(calibration.native_success_count,97u);
    EXPECT_EQ(calibration.refined_success_count,1u);
    EXPECT_EQ(calibration.unqualified_count,1u);
    EXPECT_EQ(calibration.invalid_width_count,1u);
    EXPECT_EQ(calibration.fallback_count,1u);
    EXPECT_EQ(calibration.initialization_failure_count,1u);
    EXPECT_NEAR(calibration.qualified_donor_coverage,95.0/99.0,1e-15);
    EXPECT_EQ(calibration.median_width,1.0);
    EXPECT_EQ(calibration.tau_max4,4.0);
    EXPECT_DOUBLE_EQ(*calibration.p95,1.0);
    EXPECT_EQ(calibration.excluded_donor_reasons.at("unqualified-initial-stage"),1u);
    EXPECT_EQ(calibration.excluded_donor_reasons.at("median-fallback"),1u);
    EXPECT_EQ(calibration.excluded_donor_reasons.at("non-finite-width"),1u);
    EXPECT_EQ(calibration.excluded_donor_reasons.at("initialization-failure"),1u);
}

TEST(WidthCalibrationTest, UsesDeterministicType7QuantilesForSmallAndDuplicateSamples)
{
    auto input=QualifiedInitialization(20);
    const auto first=joint::CalibrateGlobalWidth(input,"duplicates");
    std::reverse(input.atoms.begin(),input.atoms.end());
    std::reverse(input.tau.begin(),input.tau.end());
    const auto second=joint::CalibrateGlobalWidth(input,"duplicates");
    ASSERT_EQ(first.status,"Available");
    ASSERT_EQ(second.status,"Available");
    EXPECT_EQ(first.median_width,10.5);
    EXPECT_EQ(second.median_width,10.5);
    EXPECT_EQ(first.tau_max4,42.0);
    EXPECT_EQ(second.tau_max4,42.0);
    EXPECT_DOUBLE_EQ(*first.p95,19.05);
    EXPECT_DOUBLE_EQ(*second.p95,19.05);

    input=QualifiedInitialization(20);
    for(auto & atom:input.atoms) atom.mdpde[1]=1.0;
    input.atoms[18].mdpde[1]=2.0; input.atoms[19].mdpde[1]=3.0;
    const auto duplicates=joint::CalibrateGlobalWidth(input,"duplicates");
    ASSERT_EQ(duplicates.status,"Available");
    EXPECT_DOUBLE_EQ(*duplicates.p95,2.05);
}

TEST(WidthCalibrationTest, SingleExtremeDonorChangesP95WithoutChangingTheMedianBound)
{
    auto input=QualifiedInitialization(20);
    for(auto & atom:input.atoms) atom.mdpde[1]=1.0;
    input.atoms.back().mdpde[1]=1000.0;

    const auto calibration=joint::CalibrateGlobalWidth(input,"single-tail-donor");
    ASSERT_EQ(calibration.status,"Available");
    EXPECT_EQ(calibration.qualified_donor_count,20u);
    EXPECT_DOUBLE_EQ(calibration.qualified_donor_coverage,1.0);
    EXPECT_EQ(calibration.median_width,1.0);
    ASSERT_TRUE(calibration.p95);
    // The Type 7 position is 18.05, between 1 and 1000.
    EXPECT_NEAR(*calibration.p95,50.95,1e-12);
    EXPECT_EQ(calibration.tau_max4,4.0);
    EXPECT_EQ(calibration.eta_max4,std::log(4.0));
}

TEST(WidthCalibrationTest, ProjectsHighWidthDonorsToTheMedianBound)
{
    auto input=QualifiedInitialization(20);
    for(std::size_t i=0;i<input.tau.size();++i)
    {
        const double width=i<11 ? 1.0 : 10.0;
        input.tau[i]=width;
        input.atoms[i].mdpde[1]=width;
        input.atoms[i].original_tau=width;
    }
    const auto calibration=joint::CalibrateGlobalWidth(input,"multiple-tail-donors");
    ASSERT_EQ(calibration.status,"Available");
    EXPECT_EQ(calibration.qualified_donor_count,20u);
    EXPECT_EQ(calibration.median_width,1.0);
    EXPECT_EQ(calibration.p95,10.0);
    EXPECT_EQ(calibration.tau_max4,4.0);
    EXPECT_EQ(calibration.initial_violations4,9u);

    const auto projected=joint::MakeBoundedInitialization(input,calibration);
    ASSERT_TRUE(projected.available);
    EXPECT_EQ(projected.projection_count,9u);
    EXPECT_DOUBLE_EQ(projected.maximum_projection_magnitude,6.0);
    ASSERT_EQ(projected.tau_values.size(),20u);
    for(std::size_t i=0;i<input.tau.size();++i)
        EXPECT_DOUBLE_EQ(projected.tau_values[i],i<11 ? 1.0 : 4.0);
}

TEST(WidthCalibrationTest, EqualDonorsGiveFourTimesTheirWidthIncludingSubnormalWidths)
{
    auto input=QualifiedInitialization(20);
    for(double width:{0.5,std::numeric_limits<double>::denorm_min()})
    {
        for(auto & atom:input.atoms) atom.mdpde[1]=width;
        const auto calibration=joint::CalibrateGlobalWidth(input,"equal-donors");
        ASSERT_EQ(calibration.status,"Available");
        EXPECT_EQ(calibration.median_width,width);
        EXPECT_EQ(calibration.p95,width);
        EXPECT_EQ(calibration.tau_max4,4.0*width);
        EXPECT_EQ(calibration.eta_max4,std::log(4.0*width));
    }
}

TEST(WidthCalibrationTest, ExcludesInvalidMdpdeWidthsBeforeComputingTheMedian)
{
    auto input=QualifiedInitialization(20);
    for(double invalid:{0.0,-1.0,std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity()})
    {
        for(auto & atom:input.atoms) atom.mdpde[1]=invalid;
        const auto calibration=joint::CalibrateGlobalWidth(input,"invalid-donors");
        EXPECT_EQ(calibration.status,"CalibrationUnavailable");
        EXPECT_EQ(calibration.unavailable_reason,"no-qualified-initial-stage-target-widths");
        EXPECT_EQ(calibration.invalid_width_count,20u);
        EXPECT_FALSE(calibration.median_width);
        EXPECT_FALSE(calibration.eta_max4);
    }
}

TEST(WidthCalibrationTest, ReportsEmptyAndOverflowAsCalibrationUnavailable)
{
    auto input=ExampleInitialization();
    input.atoms.clear();
    input.tau.clear();
    const auto empty=joint::CalibrateGlobalWidth(input,"empty");
    EXPECT_EQ(empty.status,"CalibrationUnavailable");
    EXPECT_EQ(empty.unavailable_reason,"no-qualified-initial-stage-target-widths");
    EXPECT_FALSE(empty.p95);

    input=ExampleInitialization();
    for(auto & atom:input.atoms) if(atom.target) atom.mdpde[1]=std::numeric_limits<double>::max();
    const auto overflow=joint::CalibrateGlobalWidth(input,"overflow");
    EXPECT_EQ(overflow.status,"CalibrationUnavailable");
    EXPECT_EQ(overflow.unavailable_reason,"global-width-bound-overflow-or-invalid-log-bound");
    EXPECT_FALSE(overflow.tau_max4);
    EXPECT_FALSE(overflow.eta_max4);
}

TEST(WidthCalibrationTest, ProducesOneGlobalBoundAndProjectsWithoutChangingStoredInitialization)
{
    auto input=ExampleInitialization();
    input.tau[5]=1000.0;
    const auto original=input;
    const auto calibration=joint::CalibrateGlobalWidth(input,"global");
    ASSERT_EQ(calibration.status,"Available");
    EXPECT_EQ(calibration.median_width,50.5);
    EXPECT_EQ(calibration.p95,95.05);
    EXPECT_EQ(calibration.tau_max4,202.0);
    EXPECT_EQ(calibration.initial_violations4,1u);

    const auto projected=joint::MakeBoundedInitialization(input,calibration);
    ASSERT_TRUE(projected.available);
    EXPECT_DOUBLE_EQ(projected.tau_max,202.0);
    EXPECT_EQ(projected.projection_count,1u);
    EXPECT_DOUBLE_EQ(projected.maximum_projection_magnitude,1000.0-*calibration.tau_max4);
    EXPECT_DOUBLE_EQ(projected.tau_values[5],*calibration.tau_max4);
    EXPECT_DOUBLE_EQ(input.tau[5],original.tau[5]);
    EXPECT_DOUBLE_EQ(input.atoms[5].mdpde[1],original.atoms[5].mdpde[1]);
    EXPECT_DOUBLE_EQ(projected.tau_values[6],input.tau[6]);
    ASSERT_EQ(projected.projections.size(),100u);
    EXPECT_DOUBLE_EQ(projected.projections[5].original_initial_stage_tau.value(),6.0);
}

TEST(WidthCalibrationTest, RejectsMalformedInitializationAndUnavailableCalibration)
{
    auto input=ExampleInitialization();
    input.tau.pop_back();
    EXPECT_EQ(joint::CalibrateGlobalWidth(input,"mismatch").unavailable_reason,
        "initialization-record-size-mismatch");
    auto valid=ExampleInitialization();
    auto unavailable=joint::CalibrateGlobalWidth(JointInitialization{},"none");
    const auto projected=joint::MakeBoundedInitialization(valid,unavailable);
    EXPECT_FALSE(projected.available);
    EXPECT_EQ(projected.reason,"no-qualified-initial-stage-target-widths");
}

TEST(WidthCalibrationTest, RejectsSmallOrPoorlyCoveredDonorPopulationsAndInvalidTargetStarts)
{
    auto too_small=QualifiedInitialization(19);
    EXPECT_EQ(joint::CalibrateGlobalWidth(too_small,"small").unavailable_reason,
        "too-few-qualified-initial-stage-target-widths");

    auto low_coverage=QualifiedInitialization(22);
    low_coverage.atoms[0].initial_stage_fit_qualified=false;
    low_coverage.atoms[1].seed_source="median-fallback";
    EXPECT_EQ(joint::CalibrateGlobalWidth(low_coverage,"coverage").unavailable_reason,
        "qualified-initial-stage-target-coverage-below-0.95");

    auto invalid_start=QualifiedInitialization(20);
    invalid_start.tau[0]=std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(joint::CalibrateGlobalWidth(invalid_start,"invalid-start").unavailable_reason,
        "invalid-target-initial-width");
}

} // namespace
