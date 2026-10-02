#include <gtest/gtest.h>
#include "core/detail/joint_component/FreeDesignRank.hpp"
#include "support/JointRankWitnessCensus.hpp"
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
    EXPECT_EQ(out.certificate,n::FreeDesignRankCertificate::LocalSupport);
    EXPECT_TRUE(out.local_witness.would_certify);
    EXPECT_EQ(out.work_stage,n::FreeDesignRankWorkStage::LocalWitness);
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
    n::Matrix dense(5,3);
    for(Eigen::Index row=0;row<dense.rows();++row)
    {
        const double x=static_cast<double>(row);
        dense(row,0)=1+.1*x; dense(row,1)=1+.01*x*x; dense(row,2)=1+.002*x*x*x;
    }
    const n::Sparse design=dense.sparseView(); n::Indices columns{0,1,2};
    const auto factor=n::FreeDesignFactor::Fixed(design,columns);
    invalid=design; invalid.coeffRef(0,0)+=.01;
    const auto out=n::EvaluateFreeDesignRank(invalid,factor.get(),{{100,5,0},3});
    EXPECT_EQ(out.reason,n::SparseBackendEnabled() ? "rank-factor-mismatch" : "rank-backend-unavailable");
    invalid=z; invalid.coeffRef(0,0)=1e308;
    EXPECT_EQ(n::EvaluateFreeDesignRank(invalid,nullptr,request).reason,"rank-bound-overflow");
}
TEST(JointFreeDesignRank,ReportsExactReconstructionForecastAndBudgetStage)
{
    n::Matrix a(12,3);
    for(Eigen::Index row=0;row<a.rows();++row)
    {
        const double x=static_cast<double>(row);
        a(row,0)=1+.1*x; a(row,1)=1+.01*x*x; a(row,2)=1+.002*x*x*x;
    }
    const n::RankRequest request{{1000,3,0},3};
    const n::RankBudget generous{120,std::numeric_limits<std::size_t>::max(),256*1024*1024};
    const auto full=Rank(a,request,generous);
    if(!n::SparseBackendEnabled()) {EXPECT_EQ(full.reason,"rank-backend-unavailable"); return;}
    ASSERT_EQ(full.status,n::FreeDesignRankStatus::FullRank)<<full.reason;
    EXPECT_EQ(full.certificate,n::FreeDesignRankCertificate::SpqrReconstruction);
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
    EXPECT_EQ(n::EvaluateRank(a,request).rank,a.cols());

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
TEST(JointFreeDesignRank,LocalWitnessBudgetExhaustionKeepsTheExistingReason)
{
    const n::Matrix a=n::Matrix::Identity(4,4); const n::RankRequest request{{100,4,0},4};
    const auto result=Rank(a,request,{120,4,256*1024*1024});
    EXPECT_EQ(result.status,n::FreeDesignRankStatus::Unavailable);
    EXPECT_EQ(result.reason,"rank-work-budget");
    EXPECT_EQ(result.work_stage,n::FreeDesignRankWorkStage::LocalWitness);
    EXPECT_EQ(result.certificate,n::FreeDesignRankCertificate::None);
}
TEST(JointLocalRankWitness,CertifiesDisjointTwoColumnSupportGroups)
{
    n::Matrix a=n::Matrix::Zero(4,4);
    a(0,0)=1; a(0,1)=.2; a(1,0)=.1; a(1,1)=1;
    a(2,2)=.9; a(2,3)=.15; a(3,2)=.1; a(3,3)=.8;
    const n::Sparse sparse=a.sparseView();
    const auto census=second_stage_test::DiagnoseLocalRankWitnesses(sparse,.01);
    EXPECT_EQ(census.groups,2); EXPECT_EQ(census.covered_columns,4); EXPECT_EQ(census.total_columns,4);
    EXPECT_DOUBLE_EQ(census.coverage_fraction,1); EXPECT_EQ(census.exclusive_rows,4);
    EXPECT_EQ(census.max_group_size,2); ASSERT_TRUE(census.minimum_lower);
    EXPECT_GT(*census.minimum_lower,census.threshold_upper);
    EXPECT_TRUE(census.exclusive_rows_disjoint); EXPECT_TRUE(census.would_certify);
    EXPECT_EQ(census.reason,"local-support-full-rank-witness");
    EXPECT_EQ(n::EvaluateRank(a,{{100,4,0},4}).rank,4);
}
TEST(JointLocalRankWitness,CertifiesExclusiveSingleColumnGroups)
{
    n::Matrix a=n::Matrix::Zero(3,2); a(0,0)=2; a(1,1)=3;
    const n::Sparse sparse=a.sparseView();
    const auto census=second_stage_test::DiagnoseLocalRankWitnesses(sparse,.5);
    EXPECT_EQ(census.groups,2); EXPECT_EQ(census.covered_columns,2);
    EXPECT_EQ(census.max_group_size,1); ASSERT_TRUE(census.minimum_lower);
    EXPECT_GT(*census.minimum_lower,.5); EXPECT_TRUE(census.would_certify);
    EXPECT_EQ(n::EvaluateRank(a,{{20,2,0},2}).rank,2);
}
TEST(JointLocalRankWitness,StrictThresholdAndDependentColumnsDoNotCertify)
{
    n::Matrix a=n::Matrix::Zero(4,4);
    a(0,0)=1; a(0,1)=.2; a(1,0)=.1; a(1,1)=1;
    a(2,2)=.9; a(2,3)=.15; a(3,2)=.1; a(3,3)=.8;
    const n::Sparse sparse=a.sparseView();
    const auto baseline=second_stage_test::DiagnoseLocalRankWitnesses(sparse,.01);
    ASSERT_TRUE(baseline.would_certify); ASSERT_TRUE(baseline.minimum_lower);
    const auto equality=second_stage_test::DiagnoseLocalRankWitnesses(sparse,*baseline.minimum_lower);
    EXPECT_FALSE(equality.would_certify);
    EXPECT_EQ(equality.reason,"local-lower-not-above-threshold");

    a.col(1)=a.col(0);
    const n::Sparse deficient_sparse=a.sparseView();
    const auto deficient=second_stage_test::DiagnoseLocalRankWitnesses(deficient_sparse,1e-12);
    EXPECT_FALSE(deficient.would_certify);
    EXPECT_EQ(n::EvaluateRank(a,{{100,4,0},4}).rank,3);
}
TEST(JointLocalRankWitness,RequiresTwoExclusiveRowsAndFallsBackForLargerGroups)
{
    n::Matrix shared=n::Matrix::Zero(3,3);
    shared(0,0)=1; shared(0,1)=.2; shared(1,0)=.1; shared(1,1)=1;
    shared(1,2)=.3; shared(2,2)=1;
    const n::Sparse shared_sparse=shared.sparseView();
    const auto missing=second_stage_test::DiagnoseLocalRankWitnesses(shared_sparse,.01);
    EXPECT_FALSE(missing.would_certify);
    EXPECT_EQ(missing.reason,"missing-exclusive-rows");

    n::Matrix larger=n::Matrix::Zero(6,5);
    larger(0,0)=1; larger(0,1)=.2; larger(1,0)=.1; larger(1,1)=1;
    for(Eigen::Index col=2;col<5;++col)
    {larger(2,col)=1; larger(3,col)=.2*static_cast<double>(col); larger(4,col)=.1*static_cast<double>(col*col);}
    const n::Sparse larger_sparse=larger.sparseView();
    const auto partial=second_stage_test::DiagnoseLocalRankWitnesses(larger_sparse,.01);
    EXPECT_EQ(partial.max_group_size,3); EXPECT_EQ(partial.covered_columns,2);
    EXPECT_DOUBLE_EQ(partial.coverage_fraction,.4); EXPECT_FALSE(partial.would_certify);
    EXPECT_EQ(partial.reason,"unsupported-for-local-witness");
}
TEST(JointLocalRankWitness,TwoByTwoLowerBoundUsesOutwardArithmetic)
{
    n::Matrix block(2,2); block<<1,.3,.2,.9;
    const auto lower=n::CertifiedSmallestSingularLowerBound2x2(block(0,0),block(0,1),block(1,0),block(1,1));
    ASSERT_TRUE(lower);
    const auto oracle=n::EvaluateRank(block,{{2,2,0},2});
    ASSERT_EQ(oracle.rank,2);
    EXPECT_LE(*lower,oracle.singular_values(1));
    EXPECT_GT(*lower,0);
}
TEST(JointFreeDesignRank,LocalSupportCertificateIsSufficientAndAvoidsReconstruction)
{
    n::Matrix a=n::Matrix::Zero(4,4);
    a(0,0)=1; a(0,1)=.2; a(1,0)=.1; a(1,1)=1;
    a(2,2)=.9; a(2,3)=.15; a(3,2)=.1; a(3,3)=.8;
    const n::RankRequest request{{1000,4,0},4};
    const auto result=Rank(a,request);
    ASSERT_EQ(result.status,n::FreeDesignRankStatus::FullRank)<<result.reason;
    EXPECT_EQ(result.certificate,n::FreeDesignRankCertificate::LocalSupport);
    EXPECT_EQ(result.reason,"rank-verified-full");
    EXPECT_TRUE(result.local_witness.would_certify);
    EXPECT_EQ(result.work_stage,n::FreeDesignRankWorkStage::LocalWitness);
    EXPECT_FALSE(result.estimated_reconstruction_entries);
    ASSERT_TRUE(result.estimated_total_entries); ASSERT_TRUE(result.estimated_remaining_entries);
    EXPECT_EQ(*result.estimated_total_entries,result.entries);
    EXPECT_EQ(*result.estimated_remaining_entries,0);
    EXPECT_GT(result.minimum_lower,result.threshold_upper);
    EXPECT_EQ(n::EvaluateRank(a,request).rank,4);
}
TEST(JointFreeDesignRank,LocalWitnessUsesStrictThresholdAndFallsBack)
{
    n::Matrix a=n::Matrix::Zero(4,4);
    a(0,0)=1; a(0,1)=.2; a(1,0)=.1; a(1,1)=1;
    a(2,2)=.9; a(2,3)=.15; a(3,2)=.1; a(3,3)=.8;
    const auto witness=n::DiagnoseFreeDesignLocalWitnesses(a.sparseView(),0);
    ASSERT_TRUE(witness.would_certify); ASSERT_TRUE(witness.minimum_lower);
    const n::RankRequest request{{1000,4,0},4,*witness.minimum_lower};
    const auto result=Rank(a,request);
    EXPECT_FALSE(result.local_witness.would_certify);
    EXPECT_EQ(result.local_witness.reason,"local-lower-not-above-threshold");
    EXPECT_NE(result.certificate,n::FreeDesignRankCertificate::LocalSupport);
    const auto oracle=n::EvaluateRank(a,request);
    if(result.status==n::FreeDesignRankStatus::FullRank) EXPECT_EQ(oracle.rank,4);
}
TEST(JointFreeDesignRank,StructuralDeficiencyStillPrecedesLocalWitness)
{
    n::Matrix a=n::Matrix::Identity(4,4); a.col(3)=a.col(1);
    const auto result=Rank(a,{{100,4,0},4});
    EXPECT_EQ(result.status,n::FreeDesignRankStatus::Deficient);
    EXPECT_EQ(result.certificate,n::FreeDesignRankCertificate::Structural);
    EXPECT_FALSE(result.local_witness.would_certify);
}
TEST(JointFreeDesignRank,UnsupportedAndNonexclusiveGroupsFallBackToSpqr)
{
    n::Matrix shared=n::Matrix::Zero(3,3);
    shared(0,0)=1; shared(0,1)=.2; shared(1,0)=.1; shared(1,1)=1;
    shared(1,2)=.3; shared(2,2)=1;
    const auto nonexclusive=Rank(shared,{{100,3,0},3});
    EXPECT_EQ(nonexclusive.local_witness.reason,"missing-exclusive-rows");
    EXPECT_NE(nonexclusive.certificate,n::FreeDesignRankCertificate::LocalSupport);

    n::Matrix larger=n::Matrix::Zero(6,5);
    larger(0,0)=1; larger(0,1)=.2; larger(1,0)=.1; larger(1,1)=1;
    for(Eigen::Index col=2;col<5;++col)
    {larger(2,col)=1; larger(3,col)=.2*static_cast<double>(col); larger(4,col)=.1*static_cast<double>(col*col);}
    const auto unsupported=Rank(larger,{{100,5,0},5});
    EXPECT_EQ(unsupported.local_witness.reason,"unsupported-for-local-witness");
    EXPECT_DOUBLE_EQ(unsupported.local_witness.coverage_fraction,.4);
    EXPECT_NE(unsupported.certificate,n::FreeDesignRankCertificate::LocalSupport);
}
TEST(JointFreeDesignRank,RowColumnPermutationsAndScaleKeepLocalConclusion)
{
    n::Matrix a=n::Matrix::Zero(4,4);
    a(0,0)=1; a(0,1)=.2; a(1,0)=.1; a(1,1)=1;
    a(2,2)=.9; a(2,3)=.15; a(3,2)=.1; a(3,3)=.8;
    const n::RankRequest request{{1000,4,0},4};
    const auto original=Rank(a,request);
    ASSERT_EQ(original.status,n::FreeDesignRankStatus::FullRank);
    ASSERT_EQ(original.certificate,n::FreeDesignRankCertificate::LocalSupport);
    n::Matrix permuted(a.rows(),a.cols());
    for(Eigen::Index row=0;row<a.rows();++row)
        for(Eigen::Index col=0;col<a.cols();++col)
            permuted(row,col)=a(a.rows()-1-row,a.cols()-1-col);
    const auto reordered=Rank(permuted,request);
    EXPECT_EQ(reordered.status,n::FreeDesignRankStatus::FullRank);
    EXPECT_EQ(reordered.certificate,n::FreeDesignRankCertificate::LocalSupport);
    for(const double scale:{1e-4,1e4})
    {
        const auto scaled=Rank(scale*a,request);
        EXPECT_EQ(scaled.status,n::FreeDesignRankStatus::FullRank)<<scaled.reason;
        EXPECT_EQ(scaled.certificate,n::FreeDesignRankCertificate::LocalSupport);
        EXPECT_EQ(n::EvaluateRank(scale*a,request).rank,4);
    }
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
