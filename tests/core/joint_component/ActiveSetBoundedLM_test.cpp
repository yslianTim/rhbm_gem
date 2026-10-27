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
    joint::BoundedLmWorkForTesting()={};
    joint::BoundedLmWorkForTesting().capture_trace=true;
    Eigen::MatrixXd j(1,1); j<<1.0;
    Eigen::VectorXd r(1); r<<1.0;
    Eigen::VectorXd d(1); d<<1.0;
    Eigen::VectorXd upper(1); upper<<0.0;
    const auto result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,10.0,0.0);
    ASSERT_TRUE(result.valid)<<result.status;
    EXPECT_FALSE(result.active[0]);
    EXPECT_NEAR(result.step(0),-1.0,1e-14);
    EXPECT_EQ(result.releases,1u);
    const auto & trace=joint::BoundedLmWorkForTesting().damping_trace;
    ASSERT_EQ(trace.size(),1u);
    ASSERT_EQ(trace[0].active_set_sequence.size(),2u);
    ASSERT_EQ(trace[0].step_sequence.size(),2u);
    ASSERT_EQ(trace[0].gradient_sequence.size(),2u);
    ASSERT_EQ(trace[0].kkt_tolerance_sequence.size(),2u);
    ASSERT_EQ(trace[0].free_ranks.size(),1u);
    EXPECT_TRUE(trace[0].active_set_sequence[0][0]);
    EXPECT_FALSE(trace[0].active_set_sequence[1][0]);
    EXPECT_EQ(Bits(trace[0].step_sequence[0](0)),Bits(0.0));
    EXPECT_EQ(Bits(trace[0].step_sequence[1](0)),Bits(0.0));
    EXPECT_EQ(Bits(trace[0].gradient_sequence[0](0)),Bits(1.0));
    EXPECT_EQ(Bits(trace[0].gradient_sequence[1](0)),Bits(0.0));
    EXPECT_EQ(trace[0].free_ranks[0],1);
    EXPECT_TRUE(trace[0].within_trust_region);
    EXPECT_TRUE(trace[0].search_accepted);
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
    joint::BoundedLmWorkForTesting()={};
    joint::BoundedLmWorkForTesting().capture_trace=true;
    const auto rank=joint::SolveBoundedTrustRegionStep(j,r,d,upper,0.5,0.0);
    ASSERT_TRUE(rank.valid)<<rank.status;
    EXPECT_GT(rank.damping,0.0);
    EXPECT_LE(rank.scaled_step_norm,0.5+1e-13);
    const auto & work=joint::BoundedLmWorkForTesting();
    EXPECT_GT(work.bracket_solves,0u);
    EXPECT_GT(work.bisection_solves,0u);
    EXPECT_EQ(work.damping_solves,work.damping_trace.size());
    EXPECT_EQ(work.dense_qr_factorizations,work.damping_solves);
    EXPECT_GT(work.workspace_allocation_sites,0u);
    EXPECT_GT(work.workspace_peak_bytes,0u);
    ASSERT_FALSE(work.damping_trace.empty());
    EXPECT_EQ(work.damping_trace[0].status,"rank-deficient-free-subproblem");
    EXPECT_FALSE(work.damping_trace[0].valid);
    EXPECT_FALSE(work.damping_trace[0].search_accepted);

    upper(0)=-1.0;
    const auto infeasible=joint::SolveBoundedTrustRegionStep(Eigen::Matrix2d::Identity(),
        Eigen::Vector2d::Zero(),Eigen::Vector2d::Ones(),upper,1.0,1.0);
    EXPECT_FALSE(infeasible.valid);
    EXPECT_EQ(infeasible.status,"invalid-subproblem-input");
}

TEST(ActiveSetBoundedLMTest, RecordedInteriorDampedStepIsBitwiseStable)
{
    Eigen::Matrix2d j=Eigen::Matrix2d::Identity();
    Eigen::Vector2d r(-1.0,2.0),d(1.0,2.0),upper(10.0,10.0);
    const double radius=std::hypot(2.0/3.0,4.0/3.0);
    joint::BoundedLmWorkForTesting()={};
    joint::BoundedLmWorkForTesting().capture_trace=true;
    const auto result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,radius,0.5);

    ASSERT_TRUE(result.valid)<<result.status;
    EXPECT_EQ(Bits(result.damping),4602678819172646914ULL);
    EXPECT_EQ(Bits(result.step(0)),4604180019048437080ULL);
    EXPECT_EQ(Bits(result.step(1)),13827552055903212883ULL);
    EXPECT_EQ(Bits(result.gradient(0)),4379187687664386048ULL);
    EXPECT_EQ(Bits(result.gradient(1)),Bits(0.0));
    ASSERT_EQ(joint::BoundedLmWorkForTesting().damping_trace.size(),57u);
    const auto & trace=joint::BoundedLmWorkForTesting().damping_trace;
    EXPECT_EQ(Bits(trace[0].damping),Bits(0.0));
    EXPECT_EQ(Bits(trace[0].step(0)),Bits(1.0));
    EXPECT_EQ(Bits(trace[0].step(1)),Bits(-2.0));
    EXPECT_EQ(Bits(trace[1].damping),Bits(0.5));
    EXPECT_EQ(trace[1].status,"solved");
    EXPECT_TRUE(trace[1].valid);
    EXPECT_FALSE(trace[0].within_trust_region);
    EXPECT_FALSE(trace[0].search_accepted);
    ASSERT_EQ(trace[0].active_set_sequence.size(),1u);
    EXPECT_FALSE(trace[0].active_set_sequence[0][0]);
    EXPECT_FALSE(trace[0].active_set_sequence[0][1]);
    ASSERT_EQ(trace[0].step_sequence.size(),1u);
    EXPECT_EQ(Bits(trace[0].step_sequence[0](0)),Bits(0.0));
    EXPECT_EQ(Bits(trace[0].step_sequence[0](1)),Bits(0.0));
    ASSERT_EQ(trace[0].gradient_sequence.size(),1u);
    EXPECT_EQ(Bits(trace[0].gradient_sequence[0](0)),Bits(0.0));
    EXPECT_EQ(Bits(trace[0].gradient_sequence[0](1)),Bits(0.0));
    ASSERT_EQ(trace[0].free_ranks.size(),1u);
    EXPECT_EQ(trace[0].free_ranks[0],2);
    ASSERT_EQ(trace[0].kkt_tolerance_sequence.size(),1u);
    EXPECT_EQ(Bits(trace[0].kkt_tolerance_sequence[0]),
        Bits(2e-11*std::sqrt(2.0)*(std::sqrt(2.0)*std::sqrt(5.0))));
    EXPECT_EQ(Bits(result.multipliers(0)),Bits(0.0));
    EXPECT_EQ(Bits(result.multipliers(1)),Bits(0.0));
    EXPECT_EQ(Bits(trace[0].kkt_inf_norm),Bits(0.0));
    EXPECT_EQ(joint::BoundedLmWorkForTesting().bracket_solves,2u);
    EXPECT_EQ(joint::BoundedLmWorkForTesting().bisection_solves,54u);
}

TEST(ActiveSetBoundedLMTest, InvalidDimensionsKeepTheSubproblemFailureReason)
{
    Eigen::MatrixXd j=Eigen::MatrixXd::Identity(2,2);
    Eigen::VectorXd r=Eigen::VectorXd::Zero(1),d=Eigen::VectorXd::Ones(2),
        upper=Eigen::VectorXd::Ones(2);
    auto result=joint::SolveBoundedTrustRegionStep(j,r,d,upper,1.0,0.0);
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

} // namespace
