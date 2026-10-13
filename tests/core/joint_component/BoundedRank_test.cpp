#include <gtest/gtest.h>
#include "core/detail/joint_component/CompactSvd.hpp"
#include "core/detail/joint_component/Numerics.hpp"
#include "support/BoundedRankResearch.hpp"
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {
namespace n=rhbm_gem::core::joint_component;
using Matrix=Eigen::MatrixXd;
using Vector=Eigen::VectorXd;

struct RankFixture
{
    std::string name;
    Matrix design;
    std::vector<std::string> column_roles;
};

std::vector<RankFixture> ReferenceFixtures()
{
    Matrix ill(2,2); ill<<1,1,0,1e-7;
    Matrix duplicate(3,2); duplicate<<1,1,0,0,1,1;
    Matrix zero_column=Matrix::Zero(4,3); zero_column(0,0)=1; zero_column(1,1)=1;
    Matrix structural(3,3); structural<<1,1,0,0,0,1,0,0,0;
    Matrix numeric_deficient(3,2); numeric_deficient<<1,2,0,0,1,2;
    Matrix coupled(3,2); coupled<<1,1,0,1e-10,0,0;
    Matrix connected_full=Matrix::Zero(5,5);
    for(Eigen::Index k=0;k<5;++k)
    {
        connected_full(k,k)=1;
        if(k) connected_full(k,k-1)=.25;
    }
    Matrix connected_deficient(3,3); connected_deficient<<1,0,1,1,1,2,0,1,1;
    Matrix weak(3,2); weak<<1,1,0,1e-12,0,0;
    Matrix active_boundary(4,4); active_boundary<<0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1;
    Matrix mixed(4,3); mixed<<1,0,1,0,1,0,0,0,1,0,0,0;
    return {
        {"full-rank-well-conditioned",Matrix::Identity(4,4),{}},
        {"full-rank-ill-conditioned",ill,{}},
        {"exact-duplicate-columns",duplicate,{}},
        {"zero-support-column",zero_column,{}},
        {"structural-rank-deficiency",structural,{}},
        {"structurally-full-numerically-deficient",numeric_deficient,{}},
        {"near-threshold",weak,{}},
        {"active-amplitude-boundary",active_boundary,{"active-A","free-C","free-A","free-C"}},
        {"mixed-target-halo",mixed,{"target","halo","target"}},
        {"strongly-coupled-A-C",coupled,{"A","C"}},
        {"connected-full-rank",connected_full,{}},
        {"connected-rank-deficient",connected_deficient,{}},
        {"weak-signal-cancellation-sensitive",weak,{"A","C"}},
    };
}

Matrix NormalizeColumns(Matrix x,Vector & norms)
{
    norms.resize(x.cols());
    for(Eigen::Index column=0;column<x.cols();++column)
    {
        norms(column)=x.col(column).norm();
        if(norms(column)>0) x.col(column)/=norms(column);
    }
    return x;
}

n::Sparse SparseFrom(const Matrix & matrix)
{return matrix.sparseView();}

TEST(JointComponentBoundedRankTest, DenseSvdOracleCoversReferenceControls)
{
    for(const auto & fixture:ReferenceFixtures())
    {
        SCOPED_TRACE(fixture.name);
        const n::RankPolicy policy{fixture.design.rows(),fixture.design.cols(),0};
        Vector norms;
        const Matrix normalized=NormalizeColumns(fixture.design,norms);
        const auto reference=n::ComputeSpectrum(fixture.design,policy,fixture.design.cols(),true);
        ASSERT_TRUE(reference.available)<<reference.reason;
        Eigen::JacobiSVD<Matrix> svd(normalized,Eigen::ComputeThinU|Eigen::ComputeThinV);
        ASSERT_EQ(reference.singular_values.size(),svd.singularValues().size());
        EXPECT_LE((reference.singular_values-svd.singularValues()).lpNorm<Eigen::Infinity>(),
            2e-12*std::max(1.0,svd.singularValues()(0)));
        EXPECT_LE((reference.column_norms-norms).lpNorm<Eigen::Infinity>(),1e-14);
        const double expected_threshold=policy.Relative(fixture.design.cols())*svd.singularValues()(0);
        EXPECT_DOUBLE_EQ(reference.threshold,expected_threshold);
        EXPECT_EQ(reference.rank,(reference.singular_values.array()>reference.threshold).count());
        EXPECT_EQ(reference.rows,policy.rows);
        EXPECT_EQ(reference.columns,fixture.design.cols());
    }
}

TEST(JointComponentBoundedRankTest, FullDesignAndFreeFaceUseDifferentColumnDomains)
{
    const auto fixtures=ReferenceFixtures();
    const auto found=std::find_if(fixtures.begin(),fixtures.end(),[](const auto & f) {
        return f.name=="active-amplitude-boundary";
    });
    ASSERT_NE(found,fixtures.end());
    const auto full=n::ComputeSpectrum(found->design,{100,4,2},4,true);
    Matrix free_face=found->design(Eigen::all,Eigen::seq(1,3));
    const auto free=n::ComputeSpectrum(free_face,{100,4,2},3,true);
    ASSERT_TRUE(full.available && free.available);
    EXPECT_EQ(full.columns,4);
    EXPECT_EQ(free.columns,3);
    EXPECT_EQ(full.rank,3);
    EXPECT_EQ(free.rank,3);
}

TEST(JointComponentBoundedRankTest, StrictRankBoundaryUsesOriginalRowsAndStrictGreater)
{
    constexpr Eigen::Index original_rows=1000;
    const n::RankPolicy policy{original_rows,2,0};
    const double threshold=policy.Absolute(2,1.0);
    for(const auto & item:std::vector<std::pair<const char *,double>>{
        {"below",std::nextafter(threshold,0.0)},
        {"equal",threshold},
        {"above",std::nextafter(threshold,std::numeric_limits<double>::infinity())}})
    {
        SCOPED_TRACE(item.first);
        Matrix x=Matrix::Zero(2,2); x(0,0)=1; x(1,1)=item.second;
        const auto spectrum=n::ComputeSpectrum(x,policy,2,false);
        ASSERT_TRUE(spectrum.available)<<spectrum.reason;
        EXPECT_DOUBLE_EQ(spectrum.threshold,threshold);
        EXPECT_EQ(spectrum.rows,original_rows);
        EXPECT_EQ(spectrum.rank,item.second>threshold ? 2 : 1);
    }
}

TEST(JointComponentBoundedRankTest, StructuralObstructionsAndMatchingDoNotClaimNumericalRank)
{
    n::Sparse zero_column(3,3);
    for(Eigen::Index row=0;row<3;++row) {zero_column.insert(row,0)=1; zero_column.insert(row,1)=1;}
    auto zero=n::AssessStructuralRank(zero_column);
    ASSERT_TRUE(zero.available);
    EXPECT_TRUE(zero.proven_deficient);
    EXPECT_EQ(zero.supported_columns,2);
    EXPECT_EQ(zero.rank_upper_bound,2);
    EXPECT_EQ(zero.reason,"zero-design-column");

    n::Sparse row_obstruction(2,3); row_obstruction.insert(0,0)=1;
    row_obstruction.insert(1,1)=1; row_obstruction.insert(1,2)=1;
    const auto rows=n::AssessStructuralRank(row_obstruction);
    ASSERT_TRUE(rows.available);
    EXPECT_TRUE(rows.proven_deficient);
    EXPECT_EQ(rows.supported_rows,2);
    EXPECT_EQ(rows.rank_upper_bound,2);

    n::Sparse hall(3,3); hall.insert(0,0)=1; hall.insert(0,1)=1;
    hall.insert(1,2)=1; hall.insert(2,2)=1;
    const auto matching=n::AssessStructuralRank(hall);
    ASSERT_TRUE(matching.available);
    ASSERT_TRUE(matching.matching_available);
    EXPECT_TRUE(matching.proven_deficient);
    EXPECT_FALSE(matching.structurally_full_rank);
    EXPECT_EQ(matching.supported_rows,3);
    EXPECT_EQ(matching.matching_size,2);
    EXPECT_EQ(matching.rank_upper_bound,2);

    Matrix proportional(2,2); proportional<<1,2,2,4;
    const auto proportional_structure=n::AssessStructuralRank(SparseFrom(proportional));
    ASSERT_TRUE(proportional_structure.matching_available);
    EXPECT_TRUE(proportional_structure.structurally_full_rank);
    EXPECT_FALSE(proportional_structure.proven_deficient);
    EXPECT_EQ(proportional_structure.identical_support_pattern_groups,1);
    const auto proportional_numeric=n::ComputeSpectrum(proportional,{2,2,0},2,false);
    ASSERT_TRUE(proportional_numeric.available);
    EXPECT_EQ(proportional_numeric.rank,1);

    const auto full=n::AssessStructuralRank(SparseFrom(Matrix::Identity(3,3)));
    ASSERT_TRUE(full.matching_available);
    EXPECT_TRUE(full.structurally_full_rank);
    EXPECT_FALSE(full.proven_deficient);
    EXPECT_EQ(full.matching_size,3);
}

TEST(JointComponentBoundedRankTest, MixedRolesAndActiveFaceKeepSeparateColumnDomains)
{
    Matrix design=Matrix::Zero(4,4);
    design(0,1)=1; design(1,2)=1; design(2,3)=1;
    const auto full=n::AssessStructuralRank(SparseFrom(design));
    ASSERT_TRUE(full.available);
    EXPECT_TRUE(full.proven_deficient);
    EXPECT_EQ(full.columns,4);
    EXPECT_EQ(full.rank_upper_bound,3);
    Matrix free_face=design(Eigen::all,Eigen::seq(1,3));
    const auto free=n::AssessStructuralRank(SparseFrom(free_face));
    ASSERT_TRUE(free.matching_available);
    EXPECT_TRUE(free.structurally_full_rank);
    EXPECT_EQ(free.columns,3);

    Matrix mixed(4,3); mixed<<1,0,0,0,1,0,0,0,1,0,0,0;
    const auto roles=n::AssessStructuralRank(SparseFrom(mixed));
    ASSERT_TRUE(roles.matching_available);
    EXPECT_TRUE(roles.structurally_full_rank);
}

TEST(JointComponentBoundedRankTest, StructuralMemoryPreflightRejectsBeforeGraphAllocation)
{
    n::Sparse x(100,4); x.insert(0,0)=1; x.insert(1,1)=1; x.insert(2,2)=1; x.insert(3,3)=1;
    const auto rejected=n::AssessStructuralRank(x,0,0);
    EXPECT_FALSE(rejected.available);
    EXPECT_FALSE(rejected.matching_available);
    EXPECT_EQ(rejected.reason,"structural-rank-memory-limit");

    n::MemoryShape overflow;
    overflow.component_observations=std::numeric_limits<std::size_t>::max();
    overflow.component_design_columns=std::numeric_limits<std::size_t>::max();
    overflow.component_structural_nonzeros=std::numeric_limits<std::size_t>::max();
    overflow.component_structural_matching=true;
    EXPECT_EQ(n::PlanMemory(overflow,n::MemoryStage::StructuralRank).status,n::MemoryStatus::Overflow);
}

TEST(JointComponentBoundedRankTest, StreamingGershgorinControlsNeverOverstateTheirEvidence)
{
    const auto limited=second_stage_test::bounded_rank::ScalarGershgorin(
        SparseFrom(Matrix::Identity(4,4)),{4,4,0},0,0);
    EXPECT_FALSE(limited.available);
    EXPECT_EQ(limited.reason,"streaming-memory-limit");

    std::size_t exact_full_rank_cases{},certified{},unavailable{},false_certifications{};
    for(const auto & fixture:ReferenceFixtures())
    {
        SCOPED_TRACE(fixture.name);
        const n::RankPolicy policy{fixture.design.rows(),fixture.design.cols(),0};
        const auto oracle=n::ComputeSpectrum(fixture.design,policy,fixture.design.cols(),true);
        ASSERT_TRUE(oracle.available)<<oracle.reason;
        const auto sparse=SparseFrom(fixture.design);
        const auto estimate=second_stage_test::bounded_rank::ScalarGershgorin(sparse,policy);
        if(oracle.rank==oracle.columns) ++exact_full_rank_cases;
        if(!estimate.certified_full_rank) ++unavailable;
        if(estimate.certified_full_rank)
        {
            ++certified;
            if(oracle.rank!=oracle.columns) ++false_certifications;
        }
        EXPECT_FALSE(estimate.certified_full_rank && oracle.rank!=oracle.columns);
        if(fixture.name=="full-rank-well-conditioned")
        {
            EXPECT_TRUE(estimate.available);
            EXPECT_TRUE(estimate.gershgorin_suggests_full_rank);
            EXPECT_TRUE(estimate.certified_full_rank);
        }
        if(fixture.name=="structurally-full-numerically-deficient" ||
            fixture.name=="connected-rank-deficient" || fixture.name=="near-threshold")
            EXPECT_FALSE(estimate.certified_full_rank);
        std::cout<<"streaming-control fixture="<<fixture.name
            <<" rows="<<fixture.design.rows()<<" columns="<<fixture.design.cols()
            <<" nonzeros="<<sparse.nonZeros()<<" oracle_rank="<<oracle.rank
            <<" sigma_max="<<oracle.singular_values(0)<<" sigma_min="<<oracle.minimum
            <<" rank_threshold="<<oracle.threshold
            <<" lambda_min_lower_estimate="<<estimate.lambda_min_lower
            <<" lambda_max_upper_estimate="<<estimate.lambda_max_upper
            <<" bound_gap="<<estimate.bound_gap
            <<" estimated_full_rank="<<estimate.gershgorin_suggests_full_rank
            <<" certified_full_rank="<<estimate.certified_full_rank
            <<" status="<<estimate.reason
            <<" estimated_peak_bytes="<<estimate.estimated_peak_bytes
            <<" rss_before_bytes="<<(estimate.current_rss_before_bytes ? std::to_string(*estimate.current_rss_before_bytes) : "unavailable")
            <<" rss_after_bytes="<<(estimate.current_rss_after_bytes ? std::to_string(*estimate.current_rss_after_bytes) : "unavailable")
            <<" process_peak_rss_bytes="<<(estimate.process_peak_rss_bytes ? std::to_string(*estimate.process_peak_rss_bytes) : "unavailable")
            <<" seconds="<<estimate.seconds<<"\n";
    }
    constexpr Eigen::Index original_rows=1000;
    const n::RankPolicy near_policy{original_rows,2,0};
    const double tau=near_policy.Absolute(2,1.0);
    for(const auto & item:std::vector<std::pair<const char *,double>>{{"below",.5},{"near",1.0},{"above",2.0}})
    {
        const double angle=std::sqrt(2.0)*item.second*tau;
        Matrix x(2,2); x<<1,std::cos(angle),0,std::sin(angle);
        const auto oracle=n::ComputeSpectrum(x,near_policy,2,true);
        const auto estimate=second_stage_test::bounded_rank::ScalarGershgorin(SparseFrom(x),near_policy);
        ASSERT_TRUE(oracle.available);
        EXPECT_FALSE(estimate.certified_full_rank);
        std::cout<<"streaming-near-threshold case="<<item.first
            <<" rows="<<original_rows<<" columns=2 sigma_min="<<oracle.minimum
            <<" rank_threshold="<<oracle.threshold
            <<" lambda_min_lower_estimate="<<estimate.lambda_min_lower
            <<" lambda_max_upper_estimate="<<estimate.lambda_max_upper
            <<" status="<<estimate.reason<<" certified_full_rank="<<estimate.certified_full_rank<<"\n";
    }
    std::cout<<"streaming-summary full_rank_cases="<<exact_full_rank_cases
        <<" certified_full_rank_cases="<<certified<<" unavailable_cases="<<unavailable
        <<" near_threshold_unavailable_cases=3 false_certifications="<<false_certifications<<"\n";
    EXPECT_GT(exact_full_rank_cases,0u);
    EXPECT_EQ(certified,1u);
    EXPECT_EQ(false_certifications,0u);
    RecordProperty("full_rank_cases",std::to_string(exact_full_rank_cases));
    RecordProperty("certified_full_rank_cases",std::to_string(certified));
    RecordProperty("unavailable_cases",std::to_string(unavailable));
    RecordProperty("false_certifications",std::to_string(false_certifications));
}

} // namespace
