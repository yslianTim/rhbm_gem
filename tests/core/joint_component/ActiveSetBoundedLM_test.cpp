#include <gtest/gtest.h>

#include "core/detail/joint_component/ActiveSetBoundedLM.hpp"

#include <limits>

namespace {
namespace joint=rhbm_gem::core::joint_component;

TEST(ActiveSetBoundedLMTest, InteriorSolutionSatisfiesTrustRegionDampedKkt)
{
    Eigen::Matrix2d j=Eigen::Matrix2d::Identity();
    Eigen::Vector2d r(-1.0,2.0),d(1.0,2.0),upper(10.0,10.0);
    const double radius=std::hypot(2.0/3.0,4.0/3.0);
    const auto result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,radius,0.5);
    ASSERT_TRUE(result.valid)<<result.status;
    EXPECT_NEAR(result.damping,0.5,1e-12);
    EXPECT_FALSE(result.active[0]);
    EXPECT_FALSE(result.active[1]);
    EXPECT_NEAR(result.step(0),1.0/(1.0+0.5),1e-13);
    EXPECT_NEAR(result.step(1),-2.0/(1.0+0.5*4.0),1e-13);
    EXPECT_LE(result.kkt_inf_norm,1e-10);
}

TEST(ActiveSetBoundedLMTest, UpperBoundarySolutionHasCorrectMultiplierSign)
{
    Eigen::MatrixXd j(1,1); j<<1.0;
    Eigen::VectorXd r(1); r<<-2.0;
    Eigen::VectorXd d(1); d<<1.0;
    Eigen::VectorXd upper(1); upper<<0.5;
    const auto result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,10.0,0.0);
    ASSERT_TRUE(result.valid)<<result.status;
    EXPECT_TRUE(result.active[0]);
    EXPECT_DOUBLE_EQ(result.step(0),0.5);
    EXPECT_NEAR(result.gradient(0),-1.5,1e-14);
    EXPECT_NEAR(result.multipliers(0),1.5,1e-14);
}

TEST(ActiveSetBoundedLMTest, ReleasesAnUpperActiveSetCoordinateWhenKktRequiresAnInwardStep)
{
    Eigen::MatrixXd j(1,1); j<<1.0;
    Eigen::VectorXd r(1); r<<1.0;
    Eigen::VectorXd d(1); d<<1.0;
    Eigen::VectorXd upper(1); upper<<0.0;
    const auto result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,10.0,0.0);
    ASSERT_TRUE(result.valid)<<result.status;
    EXPECT_FALSE(result.active[0]);
    EXPECT_NEAR(result.step(0),-1.0,1e-14);
    EXPECT_EQ(result.releases,1u);
}

TEST(ActiveSetBoundedLMTest, EnforcesSeveralUpperBoundsAndChecksTheirMultipliers)
{
    Eigen::Matrix3d j=Eigen::Matrix3d::Identity();
    Eigen::Vector3d r(-2.0,-3.0,1.0),d=Eigen::Vector3d::Ones(),upper(0.25,0.5,0.0);
    const auto result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,10.0,0.0);
    ASSERT_TRUE(result.valid)<<result.status;
    EXPECT_TRUE(result.active[0]);
    EXPECT_TRUE(result.active[1]);
    EXPECT_FALSE(result.active[2]);
    EXPECT_NEAR(result.step(0),0.25,1e-14);
    EXPECT_NEAR(result.step(1),0.5,1e-14);
    EXPECT_NEAR(result.step(2),-1.0,1e-14);
    EXPECT_NEAR(result.multipliers(0),1.75,1e-14);
    EXPECT_NEAR(result.multipliers(1),2.5,1e-14);
}

TEST(ActiveSetBoundedLMTest, PreservesANearBoundaryStepAndVerifiesItsMultiplier)
{
    Eigen::MatrixXd j=Eigen::MatrixXd::Identity(1,1);
    Eigen::VectorXd r(1),d(1),upper(1);
    r<<-1.0; d<<1.0; upper<<1e-15;
    const auto result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,1e-14,0.25);
    ASSERT_TRUE(result.valid)<<result.status;
    EXPECT_TRUE(result.active[0]);
    EXPECT_DOUBLE_EQ(result.step(0),upper(0));
    EXPECT_LT(result.gradient(0),0.0);
    EXPECT_GT(result.multipliers(0),0.0);
    EXPECT_LE(result.kkt_inf_norm,1e-10);
}

TEST(ActiveSetBoundedLMTest, SolvesConstrainedTrustRegionWithTheReturnedDamping)
{
    Eigen::Matrix2d j=Eigen::Matrix2d::Identity();
    Eigen::Vector2d r(-4.0,1.0),d(1.0,2.0),upper(0.5,0.25);
    const auto result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,0.4,0.0);
    ASSERT_TRUE(result.valid)<<result.status;
    EXPECT_LE(result.scaled_step_norm,0.4+1e-13);
    EXPECT_LE(result.step(0),upper(0));
    EXPECT_LE(result.step(1),upper(1));
    Eigen::Vector2d gradient=j.transpose()*(j*result.step+r)+
        result.damping*d.cwiseProduct(d).cwiseProduct(result.step);
    EXPECT_NEAR((gradient-result.gradient).norm(),0.0,1e-12);
    EXPECT_LE(result.kkt_inf_norm,1e-9);
}

TEST(ActiveSetBoundedLMTest, TrustRegionDampingRecoversRankDeficiencyAndRejectsInfeasibleBounds)
{
    Eigen::MatrixXd j(1,2); j<<1.0,1.0;
    Eigen::VectorXd r(1); r<<-1.0;
    Eigen::VectorXd d=Eigen::VectorXd::Ones(2),upper=Eigen::VectorXd::Constant(2,10.0);
    const auto rank=joint::SolveBoundedTrustRegionStep(j,r,d,upper,0.5,0.0);
    ASSERT_TRUE(rank.valid)<<rank.status;
    EXPECT_GT(rank.damping,0.0);
    EXPECT_LE(rank.scaled_step_norm,0.5+1e-13);

    upper(0)=-1.0;
    const auto infeasible=joint::SolveBoundedTrustRegionStep(Eigen::Matrix2d::Identity(),
        Eigen::Vector2d::Zero(),Eigen::Vector2d::Ones(),upper,1.0,1.0);
    EXPECT_FALSE(infeasible.valid);
    EXPECT_EQ(infeasible.status,"invalid-subproblem-input");
}

TEST(ActiveSetBoundedLMTest, RejectsNonfiniteSubproblemCoefficients)
{
    Eigen::MatrixXd j=Eigen::MatrixXd::Identity(1,1);
    j(0,0)=std::numeric_limits<double>::quiet_NaN();
    Eigen::VectorXd r=Eigen::VectorXd::Zero(1),d=Eigen::VectorXd::Ones(1),
        upper=Eigen::VectorXd::Ones(1);
    const auto result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,1.0,0.5);
    EXPECT_FALSE(result.valid);
    EXPECT_EQ(result.status,"invalid-subproblem-input");
}

} // namespace
