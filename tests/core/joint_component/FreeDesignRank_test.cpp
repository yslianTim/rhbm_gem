#include <gtest/gtest.h>
#include "core/detail/joint_component/FreeDesignRank.hpp"
#include <limits>
#include <numeric>

namespace {
namespace n=rhbm_gem::core::joint_component;
n::FreeDesignRankResult Rank(const n::Matrix & a,const n::RankRequest & request,const n::RankBudget & budget={})
{
    const n::Sparse z=a.sparseView(); n::Indices columns(static_cast<std::size_t>(a.cols())); std::iota(columns.begin(),columns.end(),0);
    if(a.rows()<a.cols()) return n::EvaluateFreeDesignRank(z,nullptr,request,budget);
    const auto factor=n::FreeDesignFactor::Fixed(z,columns);
    n::SparseWorkForTesting()={};
    const auto result=n::EvaluateFreeDesignRank(z,factor.get(),request,budget);
    EXPECT_EQ(n::SparseWorkForTesting().compact_extractions,0);
    EXPECT_EQ(n::SparseWorkForTesting().free_design_svds,0);
    return result;
}
}
TEST(JointFreeDesignRank,FullRankWithPermutedRowsColumnsAndNoDenseWork)
{
    n::Matrix a=n::Matrix::Zero(19,5);
    for(int k=0;k<5;++k) {a((k*3+2)%19,k)=1.; a((k*3+3)%19,k)=.2; a(18,k)=.01; a.col(k).normalize();}
    const auto out=Rank(a,{{1000,10,5},10});
    if(!n::SparseBackendEnabled()) {EXPECT_EQ(out.reason,"rank-backend-unavailable"); return;}
    EXPECT_EQ(out.status,n::FreeDesignRankStatus::FullRank) << out.reason;
    EXPECT_EQ(out.rank_lower,5); EXPECT_EQ(out.rank_upper,5);
    EXPECT_GT(out.minimum_lower,out.threshold_upper);
    const auto oracle=n::EvaluateRank(a,{{1000,10,5},10});
    EXPECT_LE(out.minimum_lower,oracle.singular_values.tail(1)(0));
    EXPECT_LE(out.maximum_lower,oracle.singular_values(0)); EXPECT_GE(out.maximum_upper,oracle.singular_values(0));
    EXPECT_LT(out.reconstruction_error,1e-10); EXPECT_GT(out.orthogonal_minimum,.99);
}
TEST(JointFreeDesignRank,StructuralDeficiencyDoesNotNeedAFactor)
{
    for(const auto boundary:{n::RankBoundary::SvdNative,n::RankBoundary::StrictGreater})
    {
        n::Matrix a=n::Matrix::Identity(5,5); a.col(3)=a.col(1);
        for(int kind=0;kind<3;++kind)
        {
            if(kind==1) a.col(3).setZero(); if(kind==2) a.setZero();
            const n::Sparse z=a.sparseView(); const auto out=n::EvaluateFreeDesignRank(z,nullptr,{{100,5,0},5,-1,boundary});
            EXPECT_EQ(out.status,n::FreeDesignRankStatus::Deficient); EXPECT_LE(out.rank_upper,4);
            if(kind==2) EXPECT_EQ(out.rank_upper,0);
        }
    }
    const n::Sparse wide=n::Matrix::Ones(2,3).sparseView();
    EXPECT_LE(n::EvaluateFreeDesignRank(wide,nullptr,{{2,3,0},3}).rank_upper,2);
}
TEST(JointFreeDesignRank,ThresholdRulesNeverMakeAnUnsupportedDecision)
{
    for(const auto boundary:{n::RankBoundary::SvdNative,n::RankBoundary::StrictGreater})
    for(const double absolute:{-1.,1e-8}) for(const double multiple:{.5,1.,2.})
    {
        n::RankRequest request{{1000,8,4},8,absolute,boundary};
        const double threshold=absolute<0 ? request.policy.Relative(request.columns) : absolute;
        n::Matrix a=n::Matrix::Identity(4,4); a(3,3)=multiple*threshold;
        const auto out=Rank(a,request); const auto oracle=n::EvaluateRank(a,request);
        if(out.status==n::FreeDesignRankStatus::FullRank) EXPECT_EQ(oracle.rank,4);
        if(out.status==n::FreeDesignRankStatus::Deficient) EXPECT_LT(oracle.rank,4);
        if(!n::SparseBackendEnabled()) continue;
        if(multiple==.5) EXPECT_EQ(out.status,n::FreeDesignRankStatus::Deficient) << out.reason;
        if(multiple==1.) EXPECT_EQ(out.status,n::FreeDesignRankStatus::Unavailable);
        EXPECT_LE(out.threshold_lower,oracle.threshold); EXPECT_GE(out.threshold_upper,oracle.threshold);
    }
}
TEST(JointFreeDesignRank,BudgetsInvalidValuesAndFactorIdentity)
{
    const n::Sparse z=n::Matrix::Identity(4,4).sparseView(); const n::RankRequest request{{100,4,0},4};
    EXPECT_EQ(n::EvaluateFreeDesignRank(z,nullptr,request,{0,1000,10000}).reason,"rank-time-budget");
    EXPECT_EQ(n::EvaluateFreeDesignRank(z,nullptr,request,{1,0,10000}).reason,"rank-work-budget");
    EXPECT_EQ(n::EvaluateFreeDesignRank(z,nullptr,request,{1,1000,0}).reason,"rank-memory-budget");
    n::Sparse invalid=z; invalid.coeffRef(0,0)=n::unavailable;
    EXPECT_EQ(n::EvaluateFreeDesignRank(invalid,nullptr,request).reason,"rank-nonfinite-design");
    n::Indices columns{0,1,2,3}; const auto factor=n::FreeDesignFactor::Fixed(z,columns);
    invalid=z; invalid.coeffRef(0,0)=2;
    const auto out=n::EvaluateFreeDesignRank(invalid,factor.get(),request);
    EXPECT_EQ(out.reason,n::SparseBackendEnabled() ? "rank-factor-mismatch" : "rank-backend-unavailable");
    invalid=z; invalid.coeffRef(0,0)=1e308;
    EXPECT_EQ(n::EvaluateFreeDesignRank(invalid,nullptr,request).reason,"rank-bound-overflow");
}
TEST(JointFreeDesignRank,ReportsExactReconstructionForecastAndBudgetStage)
{
    n::Matrix a=n::Matrix::Identity(12,4);
    for(Eigen::Index col=0;col<a.cols();++col) {a(6+col,col)=.15; a.col(col).normalize();}
    const n::RankRequest request{{1000,4,0},4};
    const n::RankBudget generous{120,std::numeric_limits<std::size_t>::max(),256*1024*1024};
    const auto full=Rank(a,request,generous);
    if(!n::SparseBackendEnabled()) {EXPECT_EQ(full.reason,"rank-backend-unavailable"); return;}
    ASSERT_EQ(full.status,n::FreeDesignRankStatus::FullRank)<<full.reason;
    ASSERT_TRUE(full.estimated_total_entries); ASSERT_TRUE(full.estimated_remaining_entries);
    ASSERT_TRUE(full.estimated_reconstruction_entries);
    EXPECT_EQ(*full.estimated_reconstruction_entries,
        static_cast<std::size_t>(a.cols())*(2*full.reflector_nonzeros+static_cast<std::size_t>(a.rows())));
    EXPECT_EQ(*full.estimated_total_entries,full.entries);
    EXPECT_EQ(*full.estimated_remaining_entries,0);
    const n::Sparse sparse=a.sparseView();
    EXPECT_EQ(full.design_nonzeros,static_cast<std::size_t>(sparse.nonZeros()));
    EXPECT_GT(full.factor_r_nonzeros,0);
    EXPECT_GT(full.reflector_count,0);
    EXPECT_EQ(full.work_stage,n::FreeDesignRankWorkStage::Reconstruction);

    const n::RankBudget just_short{120,*full.estimated_total_entries-1,256*1024*1024};
    const auto stopped=Rank(a,request,just_short);
    EXPECT_EQ(stopped.status,n::FreeDesignRankStatus::Unavailable);
    EXPECT_EQ(stopped.reason,"rank-work-budget");
    EXPECT_EQ(stopped.work_stage,n::FreeDesignRankWorkStage::Reconstruction);
    ASSERT_TRUE(stopped.estimated_total_entries); ASSERT_TRUE(stopped.estimated_remaining_entries);
    EXPECT_EQ(*stopped.estimated_total_entries,*full.estimated_total_entries);
    EXPECT_EQ(*stopped.estimated_total_entries,stopped.entries+*stopped.estimated_remaining_entries);
    EXPECT_GE(*stopped.estimated_total_entries,stopped.entries);
}
TEST(JointFreeDesignRank,DiagnosticsPreserveRankDecisionsAndProductionDefaults)
{
    const n::RankBudget defaults{};
    EXPECT_DOUBLE_EQ(defaults.seconds,120);
    EXPECT_EQ(defaults.entries,100000000);
    EXPECT_EQ(defaults.workspace_bytes,256*1024*1024);
    if(!n::SparseBackendEnabled()) GTEST_SKIP()<<"SPQR rank certificate unavailable";

    n::Matrix full=n::Matrix::Identity(9,4); full(6,0)=.1; full(7,1)=.2;
    const n::RankRequest request{{1000,4,0},4};
    const auto bounded=Rank(full,request,{120,1000000,256*1024*1024});
    const auto roomy=Rank(full,request,{120,10000000,256*1024*1024});
    EXPECT_EQ(roomy.status,bounded.status);
    EXPECT_EQ(roomy.reason,bounded.reason);
    EXPECT_EQ(bounded.status,n::EvaluateRank(full,request).rank==4 ?
        n::FreeDesignRankStatus::FullRank : n::FreeDesignRankStatus::Deficient);

    n::Matrix deficient=n::Matrix::Identity(9,4); deficient.col(3)=deficient.col(1);
    const auto deficient_small=Rank(deficient,request,{120,1000000,256*1024*1024});
    const auto deficient_roomy=Rank(deficient,request,{120,10000000,256*1024*1024});
    EXPECT_EQ(deficient_small.status,n::FreeDesignRankStatus::Deficient);
    EXPECT_EQ(deficient_roomy.status,deficient_small.status);
    EXPECT_EQ(deficient_roomy.reason,deficient_small.reason);
}
TEST(JointFreeDesignRank,TinyWorkBudgetNamesTheStageWithoutChangingReason)
{
    const n::Matrix a=n::Matrix::Identity(4,4); const n::RankRequest request{{100,4,0},4};
    const auto result=Rank(a,request,{120,0,256*1024*1024});
    EXPECT_EQ(result.status,n::FreeDesignRankStatus::Unavailable);
    EXPECT_EQ(result.reason,"rank-work-budget");
    EXPECT_EQ(result.work_stage,n::FreeDesignRankWorkStage::StructuralScan);
    EXPECT_FALSE(result.estimated_total_entries);
}
TEST(JointFreeDesignRank,TimeBudgetIncludesStructuralColumnScans)
{
    constexpr Eigen::Index rows=60000;
    n::Sparse z(rows,2); z.reserve(rows);
    for(Eigen::Index row=0;row<rows;++row) z.insert(row,0)=1;
    z.makeCompressed();
    const auto result=n::EvaluateFreeDesignRank(z,nullptr,{{rows,2,0},2},{1e-4,1000000,32*1024*1024});
    EXPECT_EQ(result.status,n::FreeDesignRankStatus::Unavailable);
    EXPECT_EQ(result.reason,"rank-time-budget");
}
