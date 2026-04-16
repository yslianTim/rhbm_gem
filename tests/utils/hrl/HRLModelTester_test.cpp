#include <gtest/gtest.h>

#include <initializer_list>
#include <vector>

#include <rhbm_gem/utils/hrl/HRLModelTester.hpp>

namespace {

Eigen::VectorXd MakeVector(std::initializer_list<double> values)
{
    Eigen::VectorXd result(static_cast<Eigen::Index>(values.size()));
    Eigen::Index index{ 0 };
    for (double value : values)
    {
        result(index++) = value;
    }
    return result;
}

} // namespace

TEST(HRLModelTesterTest, RunDataEntryWithNeighborhoodTestReturnsSamplingEntries)
{
    HRLModelTester tester(3, 2, 2);
    tester.SetFittingRange(0.0, 1.0);

    const auto sampling_entries{
        tester.RunDataEntryWithNeighborhoodTest(
            MakeVector({ 1.0, 0.5, 0.0 }),
            4,
            0.0,
            1.0,
            2.0,
            1,
            0.0)
    };

    EXPECT_EQ(sampling_entries.size(), 44u);
}

TEST(HRLModelTesterTest, RunBetaMDPDETestPopulatesResidualOutputs)
{
    HRLModelTester tester(3, 2, 2);
    tester.SetFittingRange(0.0, 1.0);

    std::vector<Eigen::VectorXd> residual_mean_ols_list;
    std::vector<Eigen::VectorXd> residual_mean_mdpde_list;
    std::vector<Eigen::VectorXd> residual_sigma_ols_list;
    std::vector<Eigen::VectorXd> residual_sigma_mdpde_list;

    const std::vector<double> alpha_r_list{ 0.0, 0.5 };
    const bool result{
        tester.RunBetaMDPDETest(
            alpha_r_list,
            residual_mean_ols_list,
            residual_mean_mdpde_list,
            residual_sigma_ols_list,
            residual_sigma_mdpde_list,
            MakeVector({ 1.0, 0.5, 0.0 }),
            10,
            0.01,
            0.0,
            1)
    };

    ASSERT_TRUE(result);
    ASSERT_EQ(residual_mean_ols_list.size(), alpha_r_list.size() + 1);
    ASSERT_EQ(residual_mean_mdpde_list.size(), alpha_r_list.size() + 1);
    ASSERT_EQ(residual_sigma_ols_list.size(), alpha_r_list.size() + 1);
    ASSERT_EQ(residual_sigma_mdpde_list.size(), alpha_r_list.size() + 1);

    for (const auto & residual : residual_mean_ols_list)
    {
        EXPECT_EQ(residual.size(), 3);
    }
    for (const auto & residual : residual_mean_mdpde_list)
    {
        EXPECT_EQ(residual.size(), 3);
    }
    for (const auto & residual : residual_sigma_ols_list)
    {
        EXPECT_EQ(residual.size(), 3);
    }
    for (const auto & residual : residual_sigma_mdpde_list)
    {
        EXPECT_EQ(residual.size(), 3);
    }
}
