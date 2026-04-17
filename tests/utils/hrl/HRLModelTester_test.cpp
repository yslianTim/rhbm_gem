#include <gtest/gtest.h>

#include <initializer_list>
#include <vector>

#include <rhbm_gem/utils/hrl/HRLModelTestDataFactory.hpp>
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

TEST(HRLModelTesterTest, RunBetaMDPDETestPopulatesResidualOutputs)
{
    HRLModelTestDataFactory factory(
        3,
        rhbm_gem::GaussianLinearizationSpec::DefaultDataset());
    factory.SetFittingRange(0.0, 1.0);
    HRLModelTester tester(3);

    const auto test_input{
        factory.BuildBetaTestInput(HRLModelTestDataFactory::BetaScenario{
            MakeVector({ 1.0, 0.5, 0.0 }),
            10,
            0.01,
            0.0,
            2,
            42
        })
    };

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
            test_input,
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

TEST(HRLModelTesterTest, RunMuMDPDETestPopulatesResidualOutputs)
{
    HRLModelTestDataFactory factory(
        3,
        rhbm_gem::GaussianLinearizationSpec::DefaultDataset());
    HRLModelTester tester(3);

    const auto test_input{
        factory.BuildMuTestInput(HRLModelTestDataFactory::MuScenario{
            12,
            MakeVector({ 1.0, 0.5, 0.1 }),
            MakeVector({ 0.05, 0.025, 0.01 }),
            MakeVector({ 1.5, 0.5, 0.1 }),
            MakeVector({ 0.05, 0.025, 0.01 }),
            0.1,
            2,
            33
        })
    };

    std::vector<Eigen::VectorXd> residual_mean_median_list;
    std::vector<Eigen::VectorXd> residual_mean_mdpde_list;
    std::vector<Eigen::VectorXd> residual_sigma_median_list;
    std::vector<Eigen::VectorXd> residual_sigma_mdpde_list;

    const std::vector<double> alpha_g_list{ 0.2 };
    const bool result{
        tester.RunMuMDPDETest(
            alpha_g_list,
            residual_mean_median_list,
            residual_mean_mdpde_list,
            residual_sigma_median_list,
            residual_sigma_mdpde_list,
            test_input,
            1)
    };

    ASSERT_TRUE(result);
    ASSERT_EQ(residual_mean_median_list.size(), alpha_g_list.size() + 1);
    ASSERT_EQ(residual_mean_mdpde_list.size(), alpha_g_list.size() + 1);
    ASSERT_EQ(residual_sigma_median_list.size(), alpha_g_list.size() + 1);
    ASSERT_EQ(residual_sigma_mdpde_list.size(), alpha_g_list.size() + 1);
}

TEST(HRLModelTesterTest, RunBetaMDPDEWithNeighborhoodTestConsumesPreparedInputs)
{
    HRLModelTestDataFactory factory(
        3,
        rhbm_gem::GaussianLinearizationSpec::DefaultDataset());
    factory.SetFittingRange(0.0, 1.0);
    HRLModelTester tester(3);

    const auto test_input{
        factory.BuildNeighborhoodTestInput(HRLModelTestDataFactory::NeighborhoodScenario{
            MakeVector({ 1.0, 0.5, 0.0 }),
            8,
            0.01,
            0.0,
            1.0,
            2.0,
            1,
            45.0,
            true,
            0.0,
            4.0,
            2,
            101
        })
    };

    std::vector<Eigen::VectorXd> residual_mean_list;
    std::vector<Eigen::VectorXd> residual_sigma_list;
    double training_alpha_r_average{ 0.0 };

    const bool result{
        tester.RunBetaMDPDEWithNeighborhoodTest(
            residual_mean_list,
            residual_sigma_list,
            test_input,
            training_alpha_r_average,
            1,
            45.0)
    };

    ASSERT_TRUE(result);
    ASSERT_EQ(residual_mean_list.size(), 3u);
    ASSERT_EQ(residual_sigma_list.size(), 3u);
    EXPECT_GT(training_alpha_r_average, 0.0);
}
