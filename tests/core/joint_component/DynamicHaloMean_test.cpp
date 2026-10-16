#include <gtest/gtest.h>

#include "support/DynamicHaloMean.hpp"

using namespace joint_test::dynamic_halo;

namespace {
rhbm_gem::core::JointProblemInput Input()
{
    rhbm_gem::core::JointProblemInput input;
    input.atom_ids={"z-target","a-target","halo"};
    input.row_ids={"r0","r1","r2","r3"};
    input.observations.assign(4,0.0);
    input.support={{{0,0.1},{1,0.2}},{{1,0.3},{2,0.35}},{{0,0.4},{1,0.5},{2,0.6},{3,0.7}}};
    input.selection_domain.emplace();
    input.selection_domain->target_indices={0,1};
    return input;
}
}

TEST(DynamicHaloMeanTest, AssignsOneOwnerAndExpandsRowsToOwnedHaloSupport)
{
    auto input=Input();
    const auto result=BuildOwnership(input,{{"block-a",0,{1}},{"block-z",0,{0}}});
    ASSERT_EQ(result.halos.size(),1u);
    EXPECT_EQ(result.blocks[result.halos[0].owner_block].id,"block-a");
    EXPECT_EQ(result.halos[0].owner_overlap,2u);
    EXPECT_EQ(result.halos[0].intersecting_blocks.size(),2u);
    EXPECT_EQ(result.halos[0].overlapping_targets.size(),2u);
    EXPECT_EQ(result.block_dynamic_rows[result.halos[0].owner_block],
        (std::vector<std::size_t>{0,1,2,3}));
}

TEST(DynamicHaloMeanTest, TieBreakDoesNotDependOnBlockTraversalOrder)
{
    const auto input=Input();
    const auto forward=BuildOwnership(input,{{"z-block",0,{0}},{"a-block",0,{1}}});
    const auto reverse=BuildOwnership(input,{{"a-block",0,{1}},{"z-block",0,{0}}});
    EXPECT_EQ(forward.blocks[forward.halos[0].owner_block].id,"a-block");
    EXPECT_EQ(reverse.blocks[reverse.halos[0].owner_block].id,"a-block");
}

TEST(DynamicHaloMeanTest, FlagsHaloThatCouplesSeparateComponents)
{
    auto input=Input();
    const auto result=BuildOwnership(input,{{"left",0,{0}},{"right",1,{1}}});
    EXPECT_EQ(result.cross_component_halos,1u);
}

TEST(DynamicHaloMeanTest, ComputesPhysicalMeansAndArithmeticWidthChainRule)
{
    const std::array values{Parameters{2.0,0.4,-0.3},Parameters{4.0,0.8,0.5}};
    const auto ab=Mean(values,false),abc=Mean(values,true);
    EXPECT_DOUBLE_EQ(ab.a,3.0); EXPECT_DOUBLE_EQ(ab.b,0.6); EXPECT_DOUBLE_EQ(ab.c,0.0);
    EXPECT_DOUBLE_EQ(abc.c,0.1);
    EXPECT_DOUBLE_EQ(MeanWidthLogDerivative(0.8,2),0.4);
}
