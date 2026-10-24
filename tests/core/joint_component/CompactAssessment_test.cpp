#include <gtest/gtest.h>
#include "core/detail/joint_component/Numerics.hpp"
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include "support/JointDenseReference.hpp"
#include <rhbm_gem/utils/math/EigenHelper.hpp>
#include <algorithm>

namespace {
namespace c=rhbm_gem::core;
namespace n=c::joint_component;

n::IndependentReplayEvidence PassingReplay()
{
    n::IndependentReplayEvidence replay;
    replay.available=true; replay.valid=true; replay.passed=true; replay.feasible=true;
    replay.primary_kkt_passed=true; replay.independent_kkt_passed=true;
    replay.prediction_passed=true; replay.residual_passed=true; replay.objective_passed=true;
    replay.kkt_agrees=true; replay.gradient_agrees=true;
    replay.primary_width_gradient=0.; replay.independent_width_gradient=0.;
    replay.primary_ac_kkt=0.; replay.independent_ac_kkt=0.;
    return replay;
}

TEST(JointCompactAssessmentTest, RankDeficientDesignDoesNotGateNumericalEndpoint)
{
    const n::Vector y=n::Vector::Zero(2);
    const n::Domain domain(2,{{{0,0.},{1,1.}},{{0,0.},{1,1.}}});
    auto context=n::CreateContext(y,2);
    n::Evaluation endpoint; endpoint.valid=true;
    endpoint.eta=n::Vector::Zero(2); endpoint.beta=n::Vector::Ones(4);
    endpoint.gradient=n::Vector::Zero(2);
    endpoint.x.resize(2,4);
    endpoint.x.insert(0,0)=1.; endpoint.x.insert(0,2)=1.;
    endpoint.x.insert(1,1)=1.; endpoint.x.insert(1,3)=1.; endpoint.x.makeCompressed();
    endpoint.residual=n::Vector::Zero(2);

    const auto assessment=n::VerifyJointEndpoint(domain,y,endpoint,PassingReplay(),context);
    EXPECT_TRUE(assessment.endpoint_verified);
    EXPECT_TRUE(assessment.numerical_qualified)<<assessment.failure;
}

TEST(JointCompactAssessmentTest, KktReplayAndParameterFailuresStillRejectEndpoint)
{
    const n::Vector y=n::Vector::Zero(1);
    const n::Domain domain(1,{{{0,0.}}});
    const auto context=n::CreateContext(y,1);
    n::Evaluation endpoint; endpoint.valid=true;
    endpoint.eta=n::Vector::Zero(1); endpoint.beta=n::Vector::Ones(2);
    endpoint.gradient=n::Vector::Zero(1); endpoint.residual=n::Vector::Zero(1);

    auto replay=PassingReplay(); replay.primary_kkt_passed=false;
    const auto kkt=n::VerifyJointEndpoint(domain,y,endpoint,replay,context);
    EXPECT_FALSE(kkt.numerical_qualified);
    EXPECT_EQ(kkt.failure,"global-ac-kkt-failed");

    replay=PassingReplay(); replay.passed=false; replay.reason="replay-disagreement";
    const auto disagreement=n::VerifyJointEndpoint(domain,y,endpoint,replay,context);
    EXPECT_FALSE(disagreement.numerical_qualified);
    EXPECT_EQ(disagreement.failure,"replay-disagreement");

    endpoint.eta(0)=std::numeric_limits<double>::quiet_NaN();
    const auto invalid=n::VerifyJointEndpoint(domain,y,endpoint,PassingReplay(),context);
    EXPECT_FALSE(invalid.numerical_qualified);
    EXPECT_EQ(invalid.failure,"invalid-accepted-parameter-state");
}

TEST(JointCompactAssessmentTest, PracticalCriteriaRelaxOnlyGlobalEndpointStationarity)
{
    const n::Vector y=n::Vector::Zero(1);
    const n::Domain domain(1,{{{0,0.}}});
    const auto context=n::CreateContext(y,1);
    n::Evaluation endpoint; endpoint.valid=true;
    endpoint.eta=n::Vector::Zero(1); endpoint.beta=n::Vector::Ones(2);
    endpoint.gradient=n::Vector::Constant(1,5e-10); endpoint.residual=n::Vector::Zero(1);
    endpoint.certificate.available=true; endpoint.certificate.feasible=true;
    endpoint.certificate.projected_kkt=5e-8; endpoint.certificate.kkt_passed=false;
    auto replay=PassingReplay();
    replay.primary_ac_kkt=5e-8; replay.independent_ac_kkt=5e-8;
    replay.primary_width_gradient=5e-10; replay.independent_width_gradient=5e-10;

    const auto strict=n::VerifyJointEndpoint(domain,y,endpoint,replay,context,
        n::StrictJointConvergenceCriteria);
    const auto practical=n::VerifyJointEndpoint(domain,y,endpoint,replay,context,
        n::PracticalJointConvergenceCriteria);
    EXPECT_FALSE(strict.numerical_qualified);
    EXPECT_TRUE(practical.numerical_qualified)<<practical.failure;
    EXPECT_FALSE(endpoint.certificate.kkt_passed);

    auto invalid_replay=replay;
    invalid_replay.kkt_agrees=false;
    const auto kkt_disagreement=n::VerifyJointEndpoint(domain,y,endpoint,invalid_replay,context,
        n::PracticalJointConvergenceCriteria);
    EXPECT_FALSE(kkt_disagreement.numerical_qualified);
    invalid_replay=replay; invalid_replay.gradient_agrees=false;
    const auto gradient_disagreement=n::VerifyJointEndpoint(domain,y,endpoint,invalid_replay,context,
        n::PracticalJointConvergenceCriteria);
    EXPECT_FALSE(gradient_disagreement.numerical_qualified);
    invalid_replay=replay; invalid_replay.prediction_passed=false;
    const auto prediction_disagreement=n::VerifyJointEndpoint(domain,y,endpoint,invalid_replay,context,
        n::PracticalJointConvergenceCriteria);
    EXPECT_FALSE(prediction_disagreement.numerical_qualified);
    invalid_replay=replay; invalid_replay.residual_passed=false;
    const auto residual_disagreement=n::VerifyJointEndpoint(domain,y,endpoint,invalid_replay,context,
        n::PracticalJointConvergenceCriteria);
    EXPECT_FALSE(residual_disagreement.numerical_qualified);
    invalid_replay=replay; invalid_replay.objective_passed=false;
    const auto objective_disagreement=n::VerifyJointEndpoint(domain,y,endpoint,invalid_replay,context,
        n::PracticalJointConvergenceCriteria);
    EXPECT_FALSE(objective_disagreement.numerical_qualified);

    replay.primary_ac_kkt=std::nextafter(n::PracticalJointConvergenceCriteria.global_ac_kkt,
        std::numeric_limits<double>::infinity());
    const auto above_kkt=n::VerifyJointEndpoint(domain,y,endpoint,replay,context,
        n::PracticalJointConvergenceCriteria);
    EXPECT_FALSE(above_kkt.numerical_qualified);
}

TEST(JointCompactAssessmentTest, ProjectedWidthRouteMatchesFullProfileDerivative)
{
    rhbm_gem::eigen_helper::ScopedEigenThreadCount threads{1};
    for(const std::string topology:{"chain","cube"}) for(const int atoms:{8,32})
    {
        SCOPED_TRACE(topology+"-"+std::to_string(atoms));
        const c::JointProblem problem(second_stage_test::SyntheticJointWorkload(topology,atoms));
        const auto & data=c::JointProblemAccess::Get(problem);
        const auto eta=n::Vector::Constant(atoms,std::log(.55));
        const auto endpoint=n::EvaluateProfile(data.domain,data.y,eta,&data.context);
        ASSERT_TRUE(endpoint.valid);
        const auto full=n::PrepareDerivative(endpoint,data.context.scale,&data.context);
        const auto projected=n::PrepareProjectedWidthDerivative(endpoint,data.context.scale,&data.context);
        ASSERT_TRUE(full.valid); ASSERT_TRUE(projected.valid);
        const auto dense=second_stage_test::matched::joint_abc::DenseDifferentiate(
            endpoint,data.context.scale,&data.context);
        ASSERT_TRUE(dense.valid);
        const auto a=n::ReduceDerivative(full,endpoint.residual);
        const auto b=n::ReduceProjectedWidthDerivative(projected);
        ASSERT_TRUE(a.valid); ASSERT_TRUE(b.valid);
        const auto dense_profile=n::ComputeSpectrum(dense.jacobian,data.context.rank,atoms,false);
        const auto reduced_profile=n::ComputeSpectrum(a.jacobian,data.context.rank,atoms,false);
        const auto dense_width=n::ComputeSpectrum(dense.projected,data.context.rank,atoms,false);
        const auto reduced_width=n::ComputeSpectrum(b.projected,data.context.rank,atoms,false);
        ASSERT_TRUE(dense_profile.available && reduced_profile.available);
        ASSERT_TRUE(dense_width.available && reduced_width.available);
        EXPECT_EQ(dense_profile.rank,reduced_profile.rank);
        EXPECT_LE((dense_profile.singular_values-reduced_profile.singular_values).lpNorm<Eigen::Infinity>()/
            dense_profile.singular_values(0),1e-10);
        EXPECT_EQ(dense_width.rank,reduced_width.rank);
        EXPECT_LE((dense_width.singular_values-reduced_width.singular_values).lpNorm<Eigen::Infinity>()/
            dense_width.singular_values(0),1e-10);
        EXPECT_LE((dense.projected.colwise().norm().transpose()-b.projected_norms).norm(),
            1e-10*(1.+b.projected_norms.norm()));
    }
}
}
