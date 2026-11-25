#include <gtest/gtest.h>

#include "core/estimator/joint_component/WidthCalibration.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace {
using namespace rhbm_gem;
namespace joint=rhbm_gem::core::joint_component;

JointInitialization QualifiedInitialization(std::size_t count)
{
    JointInitialization initialization;
    initialization.valid=true;
    for(std::size_t i=0;i<count;++i)
    {
        const double width=static_cast<double>(i+1);
        initialization.tau.push_back(width);
        JointInitializationAtom atom;
        atom.id=std::to_string(i+1); atom.target=true;
        atom.mdpde={1.0,width,0.0};
        atom.original_tau=width; atom.seed_source="fitted";
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
    auto input=QualifiedInitialization(121);
    input.atoms[0].target=false;
    input.atoms[1].initial_stage_qualification=RHBMSolveQualification::Unqualified;
    input.atoms[2].seed_source="median-fallback";
    for(const auto i:{3u,7u})
    {
        auto & atom=input.atoms[i];
        atom.native_status=static_cast<int>(RHBMEstimationStatus::MAX_ITERATIONS_REACHED);
        atom.initial_stage_qualification=RHBMSolveQualification::RefinedSuccess;
        atom.initial_stage_recovery=RHBMFailedOnlyRecoveryDiagnostics{};
        atom.initial_stage_recovery->accepted=i==3;
    }
    input.atoms[3].mdpde[1]=1000.0;
    input.atoms[4].mdpde[1]=std::numeric_limits<double>::infinity();
    input.atoms[5].reason="initialization-exception: test failure";
    input.atoms[6].mdpde[1]=0.0;

    const auto calibration=joint::CalibrateGlobalWidth(input);
    ASSERT_EQ(calibration.status,"Available");
    EXPECT_EQ(calibration.target_count,120u);
    EXPECT_EQ(calibration.qualified_donor_count,114u);
    EXPECT_EQ(calibration.qualified_donor_coverage,0.95);
    EXPECT_EQ(calibration.median_width,65.5);
    EXPECT_EQ(calibration.tau_max4,262.0);
}

TEST(WidthCalibrationTest, ExplicitOutputDonorsExcludeBufferWidths)
{
    const auto input=QualifiedInitialization(22);
    std::vector<std::size_t> output(20);
    std::iota(output.begin(),output.end(),0);
    const auto calibration=joint::CalibrateGlobalWidth(input,output);
    ASSERT_EQ(calibration.status,"Available");
    EXPECT_EQ(calibration.target_count,20u);
    EXPECT_EQ(calibration.qualified_donor_count,20u);
    EXPECT_EQ(calibration.tau_max4,42.0);
    EXPECT_EQ(joint::CalibrateGlobalWidth(input).tau_max4,46.0);
    output.back()=output.front();
    EXPECT_EQ(joint::CalibrateGlobalWidth(input,output).unavailable_reason,
        "invalid-width-calibration-donor-mapping");
}

TEST(WidthCalibrationTest, FindsOddAndEvenMediansRegardlessOfDonorOrder)
{
    for(const auto count:{20u,21u})
    {
        auto input=QualifiedInitialization(count);
        const double median=count==20 ? 10.5 : 11.0;
        for(int order=0;order<2;++order)
        {
            const auto calibration=joint::CalibrateGlobalWidth(input);
            EXPECT_EQ(calibration.status,"Available");
            EXPECT_EQ(calibration.median_width,median);
            EXPECT_EQ(calibration.tau_max4,4.0*median);
            std::rotate(input.atoms.begin(),input.atoms.begin()+7,input.atoms.end());
            std::rotate(input.tau.begin(),input.tau.begin()+7,input.tau.end());
        }
    }
}

TEST(WidthCalibrationTest, SingleExtremeDonorDoesNotChangeTheMedianBound)
{
    auto input=QualifiedInitialization(20);
    for(auto & atom:input.atoms) atom.mdpde[1]=1.0;
    input.atoms.back().mdpde[1]=1000.0;

    const auto calibration=joint::CalibrateGlobalWidth(input);
    ASSERT_EQ(calibration.status,"Available");
    EXPECT_EQ(calibration.qualified_donor_count,20u);
    EXPECT_EQ(calibration.median_width,1.0);
    EXPECT_EQ(calibration.tau_max4,4.0);
    EXPECT_EQ(calibration.eta_max4,std::log(4.0));
}

TEST(WidthCalibrationTest, ProjectsHighWidthsWithoutChangingStoredInitialization)
{
    auto input=QualifiedInitialization(20);
    for(std::size_t i=0;i<input.tau.size();++i)
    {
        const double width=i<11 ? 1.0 : 10.0;
        input.tau[i]=width;
        input.atoms[i].mdpde[1]=width;
        input.atoms[i].original_tau=width;
    }
    const auto original=input;
    const auto calibration=joint::CalibrateGlobalWidth(input);
    ASSERT_EQ(calibration.status,"Available");
    EXPECT_EQ(calibration.median_width,1.0);
    EXPECT_EQ(calibration.tau_max4,4.0);

    const auto projected=joint::MakeBoundedInitialization(input,calibration);
    ASSERT_TRUE(projected.available);
    EXPECT_EQ(projected.tau_max,4.0);
    EXPECT_EQ(projected.projection_count,9u);
    EXPECT_EQ(projected.maximum_projection_magnitude,6.0);
    EXPECT_EQ(input.tau,original.tau);
    ASSERT_EQ(projected.tau_values.size(),input.tau.size());
    for(std::size_t i=0;i<input.tau.size();++i)
    {
        EXPECT_EQ(projected.tau_values[i],i<11 ? 1.0 : 4.0);
        EXPECT_EQ(input.atoms[i].original_tau,original.atoms[i].original_tau);
        EXPECT_EQ(input.atoms[i].mdpde,original.atoms[i].mdpde);
    }
}

TEST(WidthCalibrationTest, ChecksFiniteBoundsAtTheLimitsOfPositiveWidths)
{
    auto input=QualifiedInitialization(20);
    const double max=std::numeric_limits<double>::max();
    for(double width:{0.5,std::numeric_limits<double>::denorm_min(),max/4.0})
    {
        for(auto & atom:input.atoms) atom.mdpde[1]=width;
        const auto calibration=joint::CalibrateGlobalWidth(input);
        ASSERT_EQ(calibration.status,"Available");
        EXPECT_EQ(calibration.median_width,width);
        EXPECT_EQ(calibration.tau_max4,4.0*width);
        EXPECT_EQ(calibration.eta_max4,std::log(4.0*width));
    }
    for(double width:{std::nextafter(max/4.0,max),max})
    {
        for(auto & atom:input.atoms) atom.mdpde[1]=width;
        const auto calibration=joint::CalibrateGlobalWidth(input);
        EXPECT_EQ(calibration.status,"CalibrationUnavailable");
        EXPECT_EQ(calibration.unavailable_reason,"global-width-bound-overflow-or-invalid-log-bound");
        EXPECT_FALSE(calibration.median_width);
        EXPECT_FALSE(calibration.tau_max4);
        EXPECT_FALSE(calibration.eta_max4);
    }
}

TEST(WidthCalibrationTest, RejectsMissingInsufficientOrPoorlyCoveredDonors)
{
    const auto empty=joint::CalibrateGlobalWidth(JointInitialization{});
    EXPECT_EQ(empty.status,"CalibrationUnavailable");
    EXPECT_EQ(empty.unavailable_reason,"no-qualified-initial-stage-target-widths");
    EXPECT_EQ(joint::CalibrateGlobalWidth(QualifiedInitialization(19)).unavailable_reason,
        "too-few-qualified-initial-stage-target-widths");

    auto low_coverage=QualifiedInitialization(22);
    low_coverage.atoms[0].initial_stage_fit_qualified=false;
    low_coverage.atoms[1].seed_source="median-fallback";
    const auto calibration=joint::CalibrateGlobalWidth(low_coverage);
    EXPECT_EQ(calibration.unavailable_reason,"qualified-initial-stage-target-coverage-below-0.95");
    EXPECT_EQ(calibration.qualified_donor_count,20u);
    const auto projected=joint::MakeBoundedInitialization(low_coverage,calibration);
    EXPECT_FALSE(projected.available);
    EXPECT_EQ(projected.reason,calibration.unavailable_reason);
    EXPECT_TRUE(projected.tau_values.empty());
}

TEST(WidthCalibrationTest, RejectsMalformedInitializationAndInvalidWidths)
{
    auto input=QualifiedInitialization(20);
    for(double invalid:{0.0,-1.0,std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity()})
    {
        for(auto & atom:input.atoms) atom.mdpde[1]=invalid;
        const auto calibration=joint::CalibrateGlobalWidth(input);
        EXPECT_EQ(calibration.status,"CalibrationUnavailable");
        EXPECT_EQ(calibration.unavailable_reason,"no-qualified-initial-stage-target-widths");
        EXPECT_EQ(calibration.qualified_donor_count,0u);
        EXPECT_FALSE(calibration.median_width);

        input.tau[0]=invalid;
        EXPECT_EQ(joint::CalibrateGlobalWidth(input).unavailable_reason,"invalid-target-initial-width");
        input.tau[0]=1.0;
    }

    input.tau.pop_back();
    EXPECT_EQ(joint::CalibrateGlobalWidth(input).unavailable_reason,"initialization-record-size-mismatch");
}

} // namespace
