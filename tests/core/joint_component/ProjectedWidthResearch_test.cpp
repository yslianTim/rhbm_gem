#include <gtest/gtest.h>
#include "support/ProjectedWidthResearch.hpp"
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
}
