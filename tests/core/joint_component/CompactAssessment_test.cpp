#include <gtest/gtest.h>
#include "core/detail/joint_component/Numerics.hpp"
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "core/detail/joint_component/CompactSvd.hpp"
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

void CompareProjectedWidthRoutes(const n::Evaluation & endpoint,const n::EvaluationContext & context)
{
    const auto reference_prepared=n::PrepareDerivative(endpoint,context.scale,&context);
    const auto projected_prepared=n::PrepareProjectedWidthDerivative(endpoint,context.scale,&context);
    ASSERT_TRUE(reference_prepared.valid); ASSERT_TRUE(projected_prepared.valid);
    const auto reference=n::ReduceDerivative(reference_prepared,endpoint.residual);
    const auto projected=n::ReduceProjectedWidthDerivative(projected_prepared);
    ASSERT_TRUE(reference.valid); ASSERT_TRUE(projected.valid);
    EXPECT_TRUE(projected.jacobian.size()==0);
    CompareSpectrum(reference.projected,projected.projected);
    EXPECT_LE((reference.projected_norms-projected.projected_norms).norm(),
        1e-10*(1.+reference.projected_norms.norm()));
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

n::Evaluation RankEndpoint(double separation,Eigen::Index atoms=2)
{
    n::Evaluation endpoint;
    const Eigen::Index columns=2*atoms,rows=columns+4;
    endpoint.valid=true; endpoint.eta=n::Vector::Zero(atoms); endpoint.beta=n::Vector::Ones(columns);
    endpoint.gradient=n::Vector::Zero(atoms); endpoint.residual=n::Vector::Zero(rows);
    endpoint.x.resize(rows,columns);
    endpoint.x.insert(0,0)=1.; endpoint.x.insert(0,1)=1.; endpoint.x.insert(1,1)=separation;
    for(Eigen::Index column=2;column<columns;++column) endpoint.x.insert(column,column)=1.;
    endpoint.x.makeCompressed();
    endpoint.derivative.resize(rows,columns);
    for(Eigen::Index column=0;column<columns;++column) endpoint.derivative.insert(columns+column%4,column)=1.;
    endpoint.derivative.makeCompressed();
    return endpoint;
}

struct CapturedFreeDesignRank
{
    n::TiledDifferential prepared;
    Matrix factor;
    n::CompactSvdResult svd;
    std::size_t calls{};
};

CapturedFreeDesignRank CaptureFreeDesignRank(const n::Evaluation & endpoint,const n::EvaluationContext & context)
{
    struct RestoreCapture
    {
        n::CompactSvdCapture saved{n::CompactSvdCaptureForTesting()};
        ~RestoreCapture() {n::CompactSvdCaptureForTesting()=std::move(saved);}
    } restore;
    CapturedFreeDesignRank out;
    n::CompactSvdCaptureForTesting()=[&](const Matrix & factor,double,double,const n::Vector *,
        const n::CompactSvdResult & svd) {
        ++out.calls; out.factor=factor; out.svd=svd;
    };
    out.prepared=n::PrepareProjectedWidthDerivative(endpoint,context.scale,&context);
    return out;
}

TEST(JointCompactAssessmentTest, ProjectedWidthOnlyRouteMatchesProfileDerivativeAcrossLattices)
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
        CompareProjectedWidthRoutes(endpoint,data.context);
    }
}

TEST(JointCompactAssessmentTest, ProjectedWidthRankKeepsThresholdSensitiveResultsUnavailable)
{
    EigenThreads threads;
    const n::Vector y=n::Vector::Zero(6);
    auto context=n::CreateContext(y,2);
    context.rank={6,4,2};
    const n::Domain domain(6,{{{0,0.0},{1,1.0}},{{2,2.0},{3,3.0}}});
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
        const auto prepared=n::PrepareProjectedWidthDerivative(endpoint,context.scale,&context);
        ASSERT_TRUE(prepared.valid);
        const auto differential=n::ReduceProjectedWidthDerivative(prepared);
        ASSERT_TRUE(differential.valid);
        const auto rank=n::EvaluateRank(differential.projected,
            {context.rank,endpoint.eta.size(),-1,n::RankBoundary::StrictGreater},
            nullptr,n::CompactSvdVectors::Right);
        ASSERT_TRUE(rank.valid);
        const bool near=std::string(name)=="just-below" || std::string(name)=="boundary" ||
            std::string(name)=="just-above";
        EXPECT_EQ(n::CompactRankDecisionSafe(rank),!near);
        EXPECT_EQ(rank.rank,minimum>cutoff ? 2 : 1);
        const auto assessment=n::AssessEvaluated(domain,y,endpoint,endpoint,context);
        const auto evidence=n::AssessmentEvidence(assessment,rhbm_gem::JointEvidenceScope::ComponentLocal);
        const auto identified=std::find_if(evidence.begin(),evidence.end(),[](const auto & check) {
            return check.name=="numerical-identifiability";
        });
        ASSERT_NE(identified,evidence.end());
        EXPECT_EQ(identified->status,near ? rhbm_gem::JointCheckStatus::Unavailable :
            rank.rank==2 ? rhbm_gem::JointCheckStatus::Passed : rhbm_gem::JointCheckStatus::Failed);
    }
}

TEST(JointCompactAssessmentTest, DenseDesignOracleMatchesInteriorFreeDesignFactor)
{
    EigenThreads threads;
    const n::Vector y=n::Vector::Zero(8);
    auto context=n::CreateContext(y,2); context.rank={8,4,2};
    for(const auto & [name,separation]:std::array<std::pair<const char *,double>,3>{{
        {"well-conditioned",.5},{"ill-conditioned-full-rank",1e-8},{"exactly-rank-deficient",0.}}})
    {
        SCOPED_TRACE(name);
        const auto endpoint=RankEndpoint(separation);
        const auto dense=n::DesignSpectrum(endpoint.x,n::Vector::Ones(y.size()),&context.rank);
        const auto free=CaptureFreeDesignRank(endpoint,context);
        ASSERT_TRUE(dense.available);
        ASSERT_TRUE(free.svd.valid);
        ASSERT_EQ(free.calls,1u);
        ASSERT_EQ(free.factor.rows(),4);
        ASSERT_EQ(free.factor.cols(),4);
        ASSERT_EQ(free.svd.singular_values.size(),dense.singular_values.size());
        ASSERT_TRUE(free.prepared.free_design_rank.has_value());
        EXPECT_EQ(free.prepared.free_design_rank->rows,context.rank.rows);
        EXPECT_EQ(free.prepared.free_design_rank->columns,4);
        EXPECT_EQ(free.prepared.free_design_rank->full_design_columns,4);
        EXPECT_EQ(free.prepared.free_design_rank->column_norms.size(),4);
        EXPECT_FALSE(free.prepared.free_design_rank->threshold_sensitive);
        EXPECT_EQ(free.svd.rank,dense.rank);
        EXPECT_NEAR(free.svd.threshold,dense.threshold,2e-14*dense.threshold);
        EXPECT_LE((free.svd.singular_values-dense.singular_values).norm(),
            2e-12*(1.+dense.singular_values.norm()));
        EXPECT_EQ(free.prepared.valid,dense.rank==4);
        EXPECT_EQ(endpoint.beta.minCoeff()>0.,true);
        for(Eigen::Index column=0;column<endpoint.x.cols();++column)
            EXPECT_GT(endpoint.x.col(column).norm(),0.);
    }
}

TEST(JointCompactAssessmentTest, ActiveAmplitudeFreeRankDoesNotRepresentFullDesign)
{
    EigenThreads threads;
    const n::Vector y=n::Vector::Zero(10);
    auto context=n::CreateContext(y,3); context.rank={10,6,3};
    for(const Eigen::Index active_amplitudes:{1,2,3})
    {
        auto endpoint=RankEndpoint(.5,3);
        for(Eigen::Index atom=0;atom<active_amplitudes;++atom) endpoint.beta(2*atom)=0.;
        const auto dense=n::DesignSpectrum(endpoint.x,n::Vector::Ones(y.size()),&context.rank);
        const auto free=CaptureFreeDesignRank(endpoint,context);
        ASSERT_TRUE(dense.available);
        ASSERT_TRUE(free.svd.valid);
        ASSERT_EQ(free.calls,1u);
        EXPECT_EQ(dense.columns,6);
        EXPECT_LT(free.svd.singular_values.size(),dense.columns);
        ASSERT_TRUE(free.prepared.free_design_rank.has_value());
        EXPECT_EQ(free.prepared.free_design_rank->rows,context.rank.rows);
        EXPECT_EQ(free.prepared.free_design_rank->columns,free.svd.singular_values.size());
        EXPECT_EQ(free.prepared.free_design_rank->full_design_columns,dense.columns);
        EXPECT_EQ(free.prepared.free_design_rank->column_norms.size(),free.svd.singular_values.size());
        EXPECT_EQ(free.factor.cols(),free.svd.singular_values.size());
        EXPECT_EQ(free.svd.rank,free.svd.singular_values.size());
        EXPECT_TRUE(free.prepared.valid);
    }
}

TEST(JointCompactAssessmentTest, StrictGreaterBoundaryControlsNearThresholdRank)
{
    EigenThreads threads;
    const n::RankPolicy policy{8,2,1};
    const std::array<std::pair<double,Eigen::Index>,3> cases{{
        {.5,1},{1.,1},{1.5,2}}};
    for(const auto & [minimum,expected_rank]:cases)
    {
        Matrix factor=Matrix::Zero(2,2); factor.diagonal()<<4.,minimum;
        const auto strict=n::EvaluateRank(factor,{policy,2,1.,n::RankBoundary::StrictGreater});
        ASSERT_TRUE(strict.valid);
        EXPECT_EQ(strict.rank,expected_rank);
        if(minimum==1.) EXPECT_FALSE(n::CompactRankDecisionSafe(strict));
        else EXPECT_TRUE(n::CompactRankDecisionSafe(strict));
    }
}

TEST(JointCompactAssessmentTest, SmallProjectedWidthEstimateCanBeComparedWithObservedRss)
{
    EigenThreads threads;
    const c::JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    const auto & data=c::JointProblemAccess::Get(problem);
    const auto endpoint=n::EvaluateProfile(data.domain,data.y,n::Vector::Constant(8,std::log(.55)),false,&data.context);
    ASSERT_TRUE(endpoint.valid);
    const auto before=n::CurrentRssBytes();
    const auto assessment=n::AssessEvaluated(data.domain,data.y,endpoint,endpoint,data.context);
    const auto after=n::CurrentRssBytes();
    ASSERT_EQ(assessment.rank_memory_estimate.stage,n::MemoryStage::RankAssessment);
    RecordProperty("small_projected_width_estimate_bytes",
        std::to_string(assessment.rank_memory_estimate.estimated_peak_bytes));
    if(before && after)
    {
        const auto growth=*after>*before ? *after-*before : 0;
        RecordProperty("small_projected_width_observed_rss_growth_bytes",std::to_string(growth));
        RecordProperty("small_projected_width_rss_status","available");
        EXPECT_LE(growth,assessment.rank_memory_estimate.estimated_peak_bytes);
    }
    else RecordProperty("small_projected_width_rss_status","unavailable");
}

TEST(JointCompactAssessmentTest, ProjectedWidthCancellationFallbackAndActiveFaceRemainValid)
{
    EigenThreads threads;
    {
        const c::JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
        const auto & data=c::JointProblemAccess::Get(problem);
        auto endpoint=n::EvaluateProfile(data.domain,data.y,n::Vector::Constant(8,std::log(.55)),false,&data.context);
        ASSERT_TRUE(endpoint.valid);
        endpoint.derivative=endpoint.x;
        const auto projected=n::PrepareProjectedWidthDerivative(endpoint,data.context.scale,&data.context);
        ASSERT_TRUE(projected.valid);
        EXPECT_TRUE(projected.reference_order);
        const auto differential=n::ReduceProjectedWidthDerivative(projected);
        ASSERT_TRUE(differential.valid);
        EXPECT_EQ(differential.projected.cols(),8);
    }
    {
        const n::Vector y=n::Vector::Zero(6);
        auto context=n::CreateContext(y,2); context.rank={6,4,2};
        const auto endpoint=BoundaryEndpoint(1e-3,true);
        EXPECT_EQ(endpoint.beta(0),0.); EXPECT_EQ(endpoint.beta(2),0.);
        const auto projected=n::PrepareProjectedWidthDerivative(endpoint,context.scale,&context);
        ASSERT_TRUE(projected.valid);
        const auto differential=n::ReduceProjectedWidthDerivative(projected);
        ASSERT_TRUE(differential.valid);
        EXPECT_EQ(differential.projected.cols(),2);
    }
}

} // namespace
