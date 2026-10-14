#include <gtest/gtest.h>
#include "support/ProjectedWidthResearch.hpp"
#include "core/detail/joint_component/CompactSvd.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>

namespace {
namespace pw=second_stage_test::projected_width;
namespace n=rhbm_gem::core::joint_component;

std::string Values(const Eigen::VectorXd & x)
{
    std::string out="[";
    for(Eigen::Index k=0;k<x.size();++k)
    {if(k) out+=","; out+=std::to_string(x(k));}
    return out+"]";
}

TEST(JointComponentProjectedWidthResearchTest, DenseOracleBaselinesAreReproducible)
{
    std::size_t full_rank{},deficient{},threshold_sensitive{},false_full_rank{};
    for(const auto & fixture:pw::Fixtures())
    {
        SCOPED_TRACE(fixture.name);
        const auto started=std::chrono::steady_clock::now();
        const auto oracle=pw::EvaluateDenseOracle(fixture);
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        ASSERT_TRUE(oracle.valid)<<oracle.reason;
        EXPECT_EQ(oracle.free_columns,fixture.free_design.cols());
        EXPECT_EQ(oracle.width_columns,fixture.width_derivative.cols());
        EXPECT_EQ(oracle.design.full_rank,oracle.design.rank==oracle.free_columns);
        EXPECT_EQ(oracle.width.full_rank,oracle.width.rank==oracle.width_columns);
        EXPECT_GE(oracle.projection_residual,0);
        if(oracle.design.full_rank && oracle.width.full_rank) ++full_rank;
        if(!oracle.design.full_rank || !oracle.width.full_rank) ++deficient;
        if(oracle.width.threshold_sensitive) ++threshold_sensitive;
        if(oracle.width.full_rank && fixture.name=="exactly-rank-deficient") ++false_full_rank;
        std::cout<<"projected-width-baseline fixture="<<fixture.name
            <<" rows="<<oracle.rows<<" free_columns="<<oracle.free_columns
            <<" full_design_columns="<<fixture.full_design_columns
            <<" width_columns="<<oracle.width_columns
            <<" policy_rows="<<fixture.policy.rows
            <<" structural_nonzeros="<<oracle.structural_nonzeros
            <<" supported_rows="<<oracle.supported_rows
            <<" column_norms="<<Values(oracle.design_column_norms)
            <<" design_singular_values="<<Values(oracle.design.singular_values)
            <<" design_rank="<<oracle.design.rank<<" design_threshold="<<oracle.design.threshold
            <<" design_threshold_sensitive="<<oracle.design.threshold_sensitive
            <<" width_singular_values="<<Values(oracle.width.singular_values)
            <<" width_rank="<<oracle.width.rank<<" width_threshold="<<oracle.width.threshold
            <<" width_threshold_sensitive="<<oracle.width.threshold_sensitive
            <<" projection_residual="<<oracle.projection_residual
            <<" full_rank="<<(oracle.width.full_rank && oracle.design.full_rank)
            <<" target_columns="<<fixture.target_columns
            <<" active_amplitude_boundary="<<fixture.active_amplitude_boundary
            <<" seconds="<<seconds<<"\n";
    }
    EXPECT_GE(full_rank,3u);
    EXPECT_GE(deficient,1u);
    EXPECT_GE(threshold_sensitive,1u);
    EXPECT_EQ(false_full_rank,0u);
    const auto peak=n::ProcessPeakRssBytes();
    std::cout<<"projected-width-baseline-summary full_rank_controls="<<full_rank
        <<" deficient_controls="<<deficient<<" threshold_sensitive_controls="<<threshold_sensitive
        <<" false_full_rank="<<false_full_rank
        <<" process_peak_rss_bytes="<<(peak ? std::to_string(*peak) : "unavailable")<<"\n";
    RecordProperty("full_rank_controls",std::to_string(full_rank));
    RecordProperty("deficient_controls",std::to_string(deficient));
    RecordProperty("threshold_sensitive_controls",std::to_string(threshold_sensitive));
    RecordProperty("false_full_rank",std::to_string(false_full_rank));
}

TEST(JointComponentProjectedWidthResearchTest, ActiveFaceAndHaloDimensionsRemainExplicit)
{
    const auto fixtures=pw::Fixtures();
    const auto active=std::find_if(fixtures.begin(),fixtures.end(),[](const auto & f) {
        return f.name=="active-amplitude-boundary";
    });
    const auto mixed=std::find_if(fixtures.begin(),fixtures.end(),[](const auto & f) {
        return f.name=="mixed-target-halo";
    });
    ASSERT_NE(active,fixtures.end()); ASSERT_NE(mixed,fixtures.end());
    EXPECT_TRUE(active->active_amplitude_boundary);
    EXPECT_LT(active->free_design.cols(),active->full_design_columns);
    EXPECT_EQ(active->width_derivative.cols(),active->target_columns);
    EXPECT_EQ(mixed->target_columns,3);
    EXPECT_EQ(mixed->free_design.cols(),4);
}

TEST(JointComponentProjectedWidthResearchTest, ValuesOnlyAndRightVectorModesPreserveRankEvidence)
{
    for(const auto & fixture:pw::Fixtures())
    {
        SCOPED_TRACE(fixture.name);
        const auto oracle=pw::EvaluateDenseOracle(fixture);
        ASSERT_TRUE(oracle.valid)<<oracle.reason;
        const n::RankRequest request{fixture.policy,fixture.width_derivative.cols(),-1,n::RankBoundary::SvdNative};
        const auto started=std::chrono::steady_clock::now();
        const auto values=n::EvaluateRank(oracle.projected_width,request,nullptr,n::CompactSvdVectors::None);
        const double values_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        const auto vectors_started=std::chrono::steady_clock::now();
        const auto with_vectors=n::EvaluateRank(oracle.projected_width,request,nullptr,n::CompactSvdVectors::Right);
        const double vectors_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-vectors_started).count();
        ASSERT_TRUE(values.valid && with_vectors.valid);
        EXPECT_EQ(values.rank,with_vectors.rank);
        EXPECT_EQ(values.rank,oracle.width.rank);
        EXPECT_DOUBLE_EQ(values.threshold,with_vectors.threshold);
        EXPECT_EQ(n::CompactRankDecisionSafe(values),n::CompactRankDecisionSafe(with_vectors));
        EXPECT_LE((values.singular_values-with_vectors.singular_values).lpNorm<Eigen::Infinity>(),
            2e-12*std::max(1.0,with_vectors.singular_values(0)));
        EXPECT_EQ(values.right_vectors.size(),0);
        EXPECT_EQ(with_vectors.right_vectors.cols(),std::min(oracle.rows,oracle.width_columns));
        const auto vector_bytes=static_cast<std::size_t>(with_vectors.right_vectors.size())*sizeof(double);
        std::cout<<"projected-width-svd-mode fixture="<<fixture.name
            <<" rows="<<oracle.rows<<" width_columns="<<oracle.width_columns
            <<" rank="<<values.rank<<" threshold="<<values.threshold
            <<" threshold_sensitive="<<!n::CompactRankDecisionSafe(values)
            <<" right_vector_bytes="<<vector_bytes
            <<" values_only_seconds="<<values_seconds<<" right_vector_seconds="<<vectors_seconds
            <<" process_peak_rss_bytes="<<(n::ProcessPeakRssBytes() ?
                std::to_string(*n::ProcessPeakRssBytes()) : "unavailable")<<"\n";
    }
    constexpr Eigen::Index original_rows=1000;
    const n::RankPolicy policy{original_rows,3,3};
    const double threshold=policy.Relative(3);
    for(const auto & item:std::vector<std::pair<const char *,double>>{
        {"below",std::nextafter(threshold,0.0)},
        {"equal",threshold},
        {"above",std::nextafter(threshold,std::numeric_limits<double>::infinity())}})
    {
        Matrix boundary=Matrix::Zero(3,3); boundary(0,0)=1; boundary(1,1)=.5; boundary(2,2)=item.second;
        const n::RankRequest request{policy,3,-1,n::RankBoundary::StrictGreater};
        const auto values=n::EvaluateRank(boundary,request,nullptr,n::CompactSvdVectors::None);
        const auto with_vectors=n::EvaluateRank(boundary,request,nullptr,n::CompactSvdVectors::Right);
        ASSERT_TRUE(values.valid && with_vectors.valid)<<item.first;
        EXPECT_DOUBLE_EQ(values.threshold,threshold);
        EXPECT_EQ(values.rank,with_vectors.rank);
        EXPECT_EQ(values.rank,item.second>threshold ? 3 : 2);
        EXPECT_EQ(n::CompactRankDecisionSafe(values),n::CompactRankDecisionSafe(with_vectors));
        std::cout<<"projected-width-svd-boundary case="<<item.first
            <<" rows="<<original_rows<<" columns=3 threshold="<<values.threshold
            <<" minimum="<<values.singular_values(2)<<" rank="<<values.rank
            <<" threshold_sensitive="<<!n::CompactRankDecisionSafe(values)<<"\n";
    }
}
}
