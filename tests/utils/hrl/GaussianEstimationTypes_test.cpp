#include <gtest/gtest.h>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>

namespace rg = rhbm_gem;

TEST(GaussianEstimationTypesTest, ResultTypesDefaultToZeroValues)
{
    const rg::LocalGaussianResult local_result;
    EXPECT_DOUBLE_EQ(0.0, local_result.alpha_r);
    EXPECT_DOUBLE_EQ(0.0, local_result.ols.GetModel().GetAmplitude());
    EXPECT_DOUBLE_EQ(0.0, local_result.ols.GetModel().GetWidth());
    EXPECT_DOUBLE_EQ(0.0, local_result.mdpde.GetModel().GetAmplitude());
    EXPECT_DOUBLE_EQ(0.0, local_result.mdpde.GetModel().GetWidth());
    EXPECT_FALSE(local_result.is_outlier);
    EXPECT_DOUBLE_EQ(0.0, local_result.statistical_distance);

    const rg::GroupGaussianResult group_result;
    EXPECT_DOUBLE_EQ(0.0, group_result.alpha_g);
    EXPECT_DOUBLE_EQ(0.0, group_result.mean.GetAmplitude());
    EXPECT_DOUBLE_EQ(0.0, group_result.mean.GetWidth());
    EXPECT_DOUBLE_EQ(0.0, group_result.mdpde.GetAmplitude());
    EXPECT_DOUBLE_EQ(0.0, group_result.mdpde.GetWidth());
    EXPECT_DOUBLE_EQ(0.0, group_result.prior.GetModel().GetAmplitude());
    EXPECT_DOUBLE_EQ(0.0, group_result.prior.GetModel().GetWidth());
    EXPECT_TRUE(group_result.member_results.empty());
}
