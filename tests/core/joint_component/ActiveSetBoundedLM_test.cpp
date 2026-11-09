#include <gtest/gtest.h>

#include "core/detail/joint_component/ActiveSetBoundedLM.hpp"

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>

namespace {
namespace joint=rhbm_gem::core::joint_component;
std::uint64_t Bits(double value) {return std::bit_cast<std::uint64_t>(value);}

TEST(ActiveSetBoundedLMTest, InteriorSolutionSatisfiesTrustRegionDampedKkt)
{
    Eigen::Matrix2d j=Eigen::Matrix2d::Identity();
    Eigen::Vector2d r(-1.0,2.0),d(1.0,2.0),upper(10.0,10.0);
    const double radius=std::hypot(2.0/3.0,4.0/3.0);
    const auto result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,radius,0.5);
    ASSERT_TRUE(result.valid)<<result.status;
    EXPECT_NEAR(result.damping,0.5,1e-12);
    EXPECT_FALSE(result.working_bound_mask[0]);
    EXPECT_FALSE(result.working_bound_mask[1]);
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
    EXPECT_TRUE(result.working_bound_mask[0]);
    EXPECT_DOUBLE_EQ(result.step(0),0.5);
    EXPECT_NEAR(result.gradient(0),-1.5,1e-14);
    EXPECT_NEAR(result.multipliers(0),1.5,1e-14);
}

TEST(ActiveSetBoundedLMTest, ReleasesAnUpperWorkingSetCoordinateWhenKktRequiresAnInwardStep)
{
    Eigen::MatrixXd j(1,1); j<<1.0;
    Eigen::VectorXd r(1); r<<1.0;
    Eigen::VectorXd d(1); d<<1.0;
    Eigen::VectorXd upper(1); upper<<0.0;
    const auto result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,10.0,0.0);
    ASSERT_TRUE(result.valid)<<result.status;
    EXPECT_FALSE(result.working_bound_mask[0]);
    EXPECT_NEAR(result.step(0),-1.0,1e-14);
    EXPECT_EQ(result.working_set_releases,1u);
    EXPECT_EQ(result.iterations,2u);
    EXPECT_EQ(result.damping_solves,1u);
    EXPECT_EQ(result.damping_qr_factorizations,1u);
    EXPECT_DOUBLE_EQ(result.gradient(0),0.0);
    EXPECT_DOUBLE_EQ(result.multipliers(0),0.0);
    EXPECT_LE(result.kkt_inf_norm,1e-10);
}

TEST(ActiveSetBoundedLMTest, EnforcesSeveralUpperBoundsAndChecksTheirMultipliers)
{
    Eigen::Matrix3d j=Eigen::Matrix3d::Identity();
    Eigen::Vector3d r(-2.0,-3.0,1.0),d=Eigen::Vector3d::Ones(),upper(0.25,0.5,0.0);
    const auto result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,10.0,0.0);
    ASSERT_TRUE(result.valid)<<result.status;
    EXPECT_TRUE(result.working_bound_mask[0]);
    EXPECT_TRUE(result.working_bound_mask[1]);
    EXPECT_FALSE(result.working_bound_mask[2]);
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
    EXPECT_TRUE(result.working_bound_mask[0]);
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
    EXPECT_NEAR(rank.scaled_step_norm,0.5,1e-13);
    EXPECT_NEAR(rank.step(0),rank.step(1),1e-13);
    EXPECT_LE(rank.kkt_inf_norm,1e-10);
    EXPECT_GT(rank.bracket_solves,0u);
    EXPECT_GT(rank.bisection_solves,0u);
    EXPECT_EQ(rank.damping_solves,1+rank.bracket_solves+rank.bisection_solves);

    upper(0)=-1.0;
    const auto infeasible=joint::SolveBoundedTrustRegionStep(Eigen::Matrix2d::Identity(),
        Eigen::Vector2d::Zero(),Eigen::Vector2d::Ones(),upper,1.0,1.0);
    EXPECT_FALSE(infeasible.valid);
    EXPECT_EQ(infeasible.status,"invalid-subproblem-input");
}

TEST(ActiveSetBoundedLMTest, BoundaryEarlyStopPreservesTheDampedStepAndKkt)
{
    Eigen::Matrix2d j=Eigen::Matrix2d::Identity();
    Eigen::Vector2d r(-1.0,2.0),d(1.0,2.0),upper(10.0,10.0);
    const double radius=std::hypot(2.0/3.0,4.0/3.0);
    const auto result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,radius,0.5);

    ASSERT_TRUE(result.valid)<<result.status;
    EXPECT_NEAR(result.damping,.5,1e-12);
    EXPECT_NEAR(result.step(0),2.0/3.0,1e-13);
    EXPECT_NEAR(result.step(1),-2.0/3.0,1e-13);
    EXPECT_LE(result.kkt_inf_norm,1e-10);
    EXPECT_NEAR(result.scaled_step_norm,radius,1e-13*radius);
    EXPECT_EQ(result.bracket_solves,2u);
    EXPECT_LT(result.bisection_solves,54u);
    EXPECT_EQ(result.damping_solves,1+result.bracket_solves+result.bisection_solves);
    EXPECT_EQ(result.damping_active_set_iterations,result.damping_solves);
    EXPECT_EQ(result.damping_qr_factorizations,result.damping_solves);
    EXPECT_EQ(result.damping_working_set_releases,0u);
    EXPECT_EQ(Bits(result.multipliers(0)),Bits(0.0));
    EXPECT_EQ(Bits(result.multipliers(1)),Bits(0.0));
}

TEST(ActiveSetBoundedLMTest, CallWorkspaceResetsAcrossFullBoundaryFailureAndRecovery)
{
    const Eigen::Matrix2d j=Eigen::Matrix2d::Identity();
    const Eigen::Vector2d d=Eigen::Vector2d::Ones();
    const Eigen::Vector2d full_upper=Eigen::Vector2d::Constant(10.0);
    const Eigen::Vector2d full_residual(-1.0,2.0);
    const auto full=joint::SolveBoundedTrustRegionStep(j,full_residual,d,full_upper,10.0,0.0);
    ASSERT_TRUE(full.valid)<<full.status;

    const Eigen::Vector2d boundary_residual(-2.0,-3.0),boundary_upper(.25,.5);
    const auto boundary=joint::SolveBoundedTrustRegionStep(j,boundary_residual,d,
        boundary_upper,10.0,0.0);
    ASSERT_TRUE(boundary.valid)<<boundary.status;
    EXPECT_TRUE(boundary.working_bound_mask[0]);
    EXPECT_TRUE(boundary.working_bound_mask[1]);
    EXPECT_EQ(Bits(boundary.step(0)),Bits(.25));
    EXPECT_EQ(Bits(boundary.step(1)),Bits(.5));

    auto invalid_residual=Eigen::VectorXd::Zero(1);
    const auto invalid=joint::SolveBoundedTrustRegionStep(j,invalid_residual,d,
        full_upper,10.0,0.0);
    EXPECT_FALSE(invalid.valid);
    EXPECT_EQ(invalid.status,"invalid-subproblem-input");

    const auto recovered=joint::SolveBoundedTrustRegionStep(j,full_residual,d,full_upper,10.0,0.0);
    ASSERT_TRUE(recovered.valid)<<recovered.status;
    EXPECT_EQ(Bits(recovered.damping),Bits(full.damping));
    EXPECT_EQ(Bits(recovered.step(0)),Bits(full.step(0)));
    EXPECT_EQ(Bits(recovered.step(1)),Bits(full.step(1)));
    EXPECT_EQ(Bits(recovered.gradient(0)),Bits(full.gradient(0)));
    EXPECT_EQ(Bits(recovered.gradient(1)),Bits(full.gradient(1)));
    EXPECT_EQ(Bits(recovered.multipliers(0)),Bits(full.multipliers(0)));
    EXPECT_EQ(Bits(recovered.multipliers(1)),Bits(full.multipliers(1)));
    EXPECT_EQ(recovered.working_bound_mask,full.working_bound_mask);
}

TEST(ActiveSetBoundedLMTest, InvalidDimensionsKeepTheSubproblemFailureReason)
{
    Eigen::MatrixXd j=Eigen::MatrixXd::Identity(2,2);
    Eigen::VectorXd r=Eigen::VectorXd::Zero(1),d=Eigen::VectorXd::Ones(2),
        upper=Eigen::VectorXd::Ones(2);
    auto result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,0.0,0.0);
    EXPECT_FALSE(result.valid);
    EXPECT_EQ(result.status,"invalid-trust-region-input");
    result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,1.0,0.0);
    EXPECT_FALSE(result.valid);
    EXPECT_EQ(result.status,"invalid-subproblem-input");

    j.resize(2,0); r=Eigen::VectorXd::Zero(2); d.resize(0); upper.resize(0);
    result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,1.0,0.0);
    EXPECT_FALSE(result.valid);
    EXPECT_EQ(result.status,"invalid-subproblem-input");
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

TEST(ActiveSetBoundedLMTest, StationaryNullAtUpperBoundRemainsOutsideTheWorkingSet)
{
    const Eigen::Matrix2d j=Eigen::Vector2d(1.0,1e-161).asDiagonal();
    const Eigen::Vector2d r(-2.0,0.0),d(1.0,1e-161),upper(.5,0.0);
    const auto result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,1.0,0.0);
    ASSERT_TRUE(result.valid)<<result.status;
    ASSERT_EQ(result.working_bound_mask.size(),2u);
    EXPECT_DOUBLE_EQ(result.step(0),upper(0)); EXPECT_TRUE(result.working_bound_mask[0]);
    EXPECT_DOUBLE_EQ(result.step(1),upper(1)); EXPECT_FALSE(result.working_bound_mask[1]);
    EXPECT_DOUBLE_EQ(result.multipliers(1),0.0);
    EXPECT_EQ(result.iterations,2u);
    EXPECT_EQ(result.damping_qr_factorizations,1u);
    EXPECT_LE(result.kkt_inf_norm,1e-10);
}

// Four informative Widths must remain movable beside an unidentifiable Width.
TEST(ActiveSetBoundedLMTest, DegenerateFifthWidthDoesNotBlockInformativeUpdates)
{
    Eigen::MatrixXd j=Eigen::MatrixXd::Zero(5,5);
    j.diagonal()<<.02,.007,.006,.004,1.56536679e-161;
    Eigen::VectorXd r(5); r<<-.004,.002,-.001,.003,0.0;
    const Eigen::VectorXd d=j.colwise().blueNorm();
    const Eigen::VectorXd upper=Eigen::VectorXd::Constant(5,.5);
    ASSERT_GT((j.transpose()*r).head(4).cwiseAbs().maxCoeff(),1e-9);
    const auto result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,.1,0.0);
    ASSERT_TRUE(result.valid)<<result.status;
    EXPECT_TRUE(result.step.allFinite());
    EXPECT_TRUE((result.step.array()<=upper.array()).all());
    EXPECT_EQ(result.step(4),0.0);
    EXPECT_GT(result.step.head(4).norm(),.1);
    const Eigen::VectorXd gradient=j.transpose()*(j*result.step+r);
    EXPECT_LE(gradient.cwiseAbs().maxCoeff(),1e-9);
    EXPECT_LT((j*result.step+r).squaredNorm(),r.squaredNorm());
    EXPECT_LT(result.bracket_solves,80u);
    EXPECT_LE((gradient-result.gradient).norm(),1e-15);
    // A later profile can make the formerly null direction informative again.
    j(4,4)=.01; r(4)=-.002;
    const Eigen::VectorXd updated_d=j.colwise().blueNorm();
    const auto updated=joint::SolveBoundedTrustRegionStep(j,r,updated_d,upper,.1,0.0);
    ASSERT_TRUE(updated.valid)<<updated.status;
    EXPECT_GT(updated.step(4),.1);
}

} // namespace

TEST(ActiveSetBoundedLMTest, SmallJacobianWithNonstationaryGradientIsNotFrozen)
{
    Eigen::Matrix2d j=Eigen::Vector2d(1.0,1e-16).asDiagonal();
    const Eigen::Vector2d r(0.0,-1e10),d(1.0,1e-16),upper(1.0,.2);
    const auto result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,1.0,0.0);
    ASSERT_TRUE(result.valid)<<result.status;
    EXPECT_DOUBLE_EQ(result.step(1),upper(1));
    EXPECT_TRUE(result.working_bound_mask[1]);
    EXPECT_LT(result.gradient(1),-1e-9);
    EXPECT_GT(result.multipliers(1),0.0);
}

TEST(ActiveSetBoundedLMTest, SmallDiagonalDoesNotPermitAnAbsoluteToleranceRadiusOvershoot)
{
    Eigen::MatrixXd j=Eigen::MatrixXd::Constant(1,1,1e-16);
    const Eigen::VectorXd r=Eigen::VectorXd::Constant(1,-1e-16);
    const Eigen::VectorXd d=Eigen::VectorXd::Constant(1,1e-16);
    const double radius=1e-17;
    const auto step=joint::SolveBoundedTrustRegionStep(j,r,d,Eigen::VectorXd::Constant(1,10.0),radius,0.0);
    ASSERT_TRUE(step.valid)<<step.status;
    EXPECT_TRUE(step.step.allFinite());
    EXPECT_GT(step.step(0),.09);
    EXPECT_LE(step.step(0),.1+1e-14);
    EXPECT_LE(step.scaled_step_norm,radius*(1+64*std::numeric_limits<double>::epsilon()));
    EXPECT_LT((j*step.step+r).norm(),r.norm());
}
