#include <gtest/gtest.h>
#include "core/detail/joint_component/Numerics.hpp"
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include <Eigen/SVD>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {
namespace c=rhbm_gem::core;
namespace n=c::joint_component;
using Matrix=Eigen::MatrixXd;

struct EigenThreads
{
    int saved{Eigen::nbThreads()};
    EigenThreads() {Eigen::setNbThreads(1);}
    ~EigenThreads() {Eigen::setNbThreads(saved);}
};

void CompareSpectrum(const Matrix & reference,const Matrix & compact)
{
    ASSERT_EQ(reference.rows(),compact.rows()); ASSERT_EQ(reference.cols(),compact.cols());
    const Eigen::JacobiSVD<Matrix> a(reference),b(compact);
    ASSERT_EQ(a.singularValues().size(),b.singularValues().size());
    EXPECT_LE((a.singularValues()-b.singularValues()).norm(),
        1e-10*(1.+a.singularValues().norm()));
}

void CompareDerivativeRoutes(const n::Evaluation & endpoint,const n::EvaluationContext & context)
{
    const auto reference_prepared=n::PrepareDerivative(endpoint,context.scale,&context);
    const auto compact_prepared=n::PrepareDerivativeCompact(endpoint,context.scale,&context);
    ASSERT_TRUE(reference_prepared.valid); ASSERT_TRUE(compact_prepared.valid);
    const auto reference=n::ReduceDerivative(reference_prepared,endpoint.residual);
    const auto compact=n::ReduceDerivativeCompact(compact_prepared,endpoint.residual);
    ASSERT_TRUE(reference.valid); ASSERT_TRUE(compact.valid);
    CompareSpectrum(reference.projected,compact.projected);
    CompareSpectrum(reference.jacobian,compact.jacobian);
    EXPECT_LE((reference.projected_norms-compact.projected_norms).norm(),
        1e-10*(1.+reference.projected_norms.norm()));
    EXPECT_LE((reference.jacobian_norms-compact.jacobian_norms).norm(),
        1e-10*(1.+reference.jacobian_norms.norm()));
    const Eigen::JacobiSVD<Matrix> a(reference.jacobian,Eigen::ComputeThinU|Eigen::ComputeThinV);
    const Eigen::JacobiSVD<Matrix> b(compact.jacobian,Eigen::ComputeThinU|Eigen::ComputeThinV);
    const auto reference_correction=a.solve(-reference.response);
    const auto compact_correction=b.solve(-compact.response);
    EXPECT_LE((reference_correction-compact_correction).norm(),
        1e-9*(1.+reference_correction.norm()));
}

n::Evaluation BoundaryEndpoint(double minimum,bool active_face=false)
{
    n::Evaluation endpoint;
    endpoint.valid=true; endpoint.eta=n::Vector::Zero(2); endpoint.beta=n::Vector::Ones(4);
    endpoint.gradient=n::Vector::Zero(2); endpoint.residual=n::Vector::Zero(6);
    endpoint.x.resize(6,4);
    std::vector<Eigen::Triplet<double>> design;
    for(Eigen::Index k=0;k<4;++k) design.emplace_back(k,k,1.);
    endpoint.x.setFromTriplets(design.begin(),design.end());
    if(active_face) {endpoint.beta(0)=0.; endpoint.beta(2)=0.;}
    const Eigen::Index first_column=active_face ? 1 : 0;
    const Eigen::Index second_column=active_face ? 3 : 2;
    endpoint.derivative.resize(6,4);
    endpoint.derivative.insert(4,first_column)=1.;
    endpoint.derivative.insert(5,second_column)=minimum;
    return endpoint;
}

TEST(JointCompactAssessmentTest, CompactStackMatchesObservationReferenceAcrossLattices)
{
    EigenThreads threads;
    for(const std::string topology:{"chain","cube"}) for(const int atoms:{8,32})
    {
        SCOPED_TRACE(topology+"-"+std::to_string(atoms));
        const c::JointProblem problem(second_stage_test::SyntheticJointWorkload(topology,atoms));
        const auto & data=c::JointProblemAccess::Get(problem);
        const auto eta=n::Vector::Constant(atoms,std::log(.55));
        const auto endpoint=n::EvaluateProfile(data.domain,data.y,eta,false,&data.context);
        ASSERT_TRUE(endpoint.valid);
        CompareDerivativeRoutes(endpoint,data.context);
    }
}

TEST(JointCompactAssessmentTest, RankBoundaryUsesObservationFallbackWithinOneUlpBand)
{
    EigenThreads threads;
    const n::Domain domain(6,{});
    const n::Vector y=n::Vector::Zero(6);
    auto context=n::CreateContext(y,2);
    context.rank={6,4,2};
    const double cutoff=context.rank.Relative(2);
    const std::array<std::pair<const char *,double>,5> cases{{
        {"below",.5*cutoff},
        {"just-below",std::nextafter(cutoff,0.)},
        {"boundary",cutoff},
        {"just-above",std::nextafter(cutoff,std::numeric_limits<double>::infinity())},
        {"above",2.*cutoff}}};
    for(const auto & [name,minimum]:cases)
    {
        SCOPED_TRACE(name);
        const auto endpoint=BoundaryEndpoint(minimum);
        n::AssessmentWorkForTesting()={};
        const auto assessment=n::AssessEvaluated(domain,y,endpoint,endpoint,context);
        ASSERT_TRUE(assessment.jacobian);
        const bool near=std::string(name)=="just-below" || std::string(name)=="boundary" ||
            std::string(name)=="just-above";
        EXPECT_EQ(n::AssessmentWorkForTesting().compact_boundary_fallbacks,near ? 1 : 0);
        EXPECT_EQ(n::AssessmentWorkForTesting().compact_accepted,near ? 0 : 1);
        EXPECT_EQ(n::AssessmentWorkForTesting().compact_attempts,1);
    }
}

TEST(JointCompactAssessmentTest, CancellationFallbackAndActiveFaceRemainValid)
{
    EigenThreads threads;
    {
        const c::JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
        const auto & data=c::JointProblemAccess::Get(problem);
        auto endpoint=n::EvaluateProfile(data.domain,data.y,n::Vector::Constant(8,std::log(.55)),false,&data.context);
        ASSERT_TRUE(endpoint.valid);
        endpoint.derivative=endpoint.x;
        const auto reference=n::EvaluateProfile(data.domain,data.y,endpoint.eta,true,&data.context);
        ASSERT_TRUE(reference.valid);
        n::AssessmentWorkForTesting()={};
        const auto assessment=n::AssessEvaluated(data.domain,data.y,endpoint,reference,data.context);
        EXPECT_TRUE(assessment.jacobian);
        EXPECT_EQ(n::AssessmentWorkForTesting().compact_other_fallbacks,1);
    }
    {
        const n::Domain domain(6,{});
        const n::Vector y=n::Vector::Zero(6);
        auto context=n::CreateContext(y,2); context.rank={6,4,2};
        const auto endpoint=BoundaryEndpoint(1e-3,true);
        const auto assessment=n::AssessEvaluated(domain,y,endpoint,endpoint,context);
        EXPECT_EQ(assessment.primary.beta(0),0.);
        EXPECT_EQ(assessment.primary.beta(2),0.);
        EXPECT_TRUE(assessment.jacobian);
    }
}

} // namespace
