#include <gtest/gtest.h>

#include <limits>
#include <stdexcept>
#include <vector>

#include <rhbm_gem/utils/hrl/TestDataFactory.hpp>
#include <rhbm_gem/utils/hrl/RHBMTester.hpp>

namespace {
namespace rt = rhbm_gem::rhbm_tester;
namespace tdf = rhbm_gem::test_data_factory;
namespace rg = rhbm_gem;

tdf::GaussianParameterDistribution MakeDistribution(
    const rg::GaussianModel3D & mean,
    const rg::GaussianModel3DUncertainty & sigma = rg::GaussianModel3DUncertainty{ 0.05, 0.025, 0.01 })
{
    return tdf::GaussianParameterDistribution{ mean, sigma };
}

rt::BetaMDPDETestOptions MakeBetaOptions(
    const std::vector<double> & alpha_r_list,
    bool alpha_training)
{
    rt::BetaMDPDETestOptions options;
    options.requested_alpha_r_list = alpha_r_list;
    options.alpha_training = alpha_training;
    options.thread_size = 1;
    return options;
}

rt::MuMDPDETestOptions MakeMuOptions(
    const std::vector<double> & alpha_g_list,
    bool alpha_training)
{
    rt::MuMDPDETestOptions options;
    options.requested_alpha_g_list = alpha_g_list;
    options.alpha_training = alpha_training;
    options.thread_size = 1;
    return options;
}

void ExpectBiasStatisticSize(const rt::BiasStatistics & bias)
{
    EXPECT_EQ(bias.mean.size(), rg::GaussianModel3D::ParameterSize());
    EXPECT_EQ(bias.sigma.size(), rg::GaussianModel3D::ParameterSize());
}

} // namespace

TEST(RHBMTesterTest, RunBetaMDPDETestPopulatesBiasOutputs)
{
    const std::vector<double> alpha_r_list{ 0.0, 0.5 };
    const auto test_input{
        tdf::BuildLocalTestData(tdf::LocalScenario{
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42
        })
    };

    const auto bias{
        rt::RunBetaMDPDETest(test_input, MakeBetaOptions(alpha_r_list, true))
    };

    ASSERT_EQ(bias.mdpde.requested_alpha.size(), alpha_r_list.size());
    ExpectBiasStatisticSize(bias.ols);
    for (const auto & requested_alpha_bias : bias.mdpde.requested_alpha)
    {
        ExpectBiasStatisticSize(requested_alpha_bias);
    }
    ASSERT_TRUE(bias.mdpde.trained_alpha.has_value());
    ExpectBiasStatisticSize(bias.mdpde.trained_alpha.value());
    ASSERT_TRUE(bias.mdpde.trained_alpha_average.has_value());
    EXPECT_GE(bias.mdpde.trained_alpha_average.value(), 0.0);
}

TEST(RHBMTesterTest, RunMuMDPDETestPopulatesBiasOutputs)
{
    const std::vector<double> alpha_g_list{ 0.2 };
    const auto test_input{
        tdf::BuildGroupTestData(tdf::GroupScenario{
            12,
            10,
            MakeDistribution(rg::GaussianModel3D{ 1.0, 0.5, 0.1 }),
            MakeDistribution(rg::GaussianModel3D{ 1.5, 0.5, 0.1 }),
            0.1,
            2,
            33
        })
    };

    const auto bias{
        rt::RunMuMDPDETest(test_input, MakeMuOptions(alpha_g_list, true))
    };

    ASSERT_EQ(bias.mdpde.requested_alpha.size(), alpha_g_list.size());
    ExpectBiasStatisticSize(bias.median);
    for (const auto & requested_alpha_bias : bias.mdpde.requested_alpha)
    {
        ExpectBiasStatisticSize(requested_alpha_bias);
    }
    ASSERT_TRUE(bias.mdpde.trained_alpha.has_value());
    ExpectBiasStatisticSize(bias.mdpde.trained_alpha.value());
    ASSERT_TRUE(bias.mdpde.trained_alpha_average.has_value());
    EXPECT_GE(bias.mdpde.trained_alpha_average.value(), 0.0);
}

TEST(RHBMTesterTest, RunBetaMDPDETestSkipsTrainedAlphaWhenDisabled)
{
    const std::vector<double> alpha_r_list{ 0.5 };
    const auto test_input{
        tdf::BuildLocalTestData(tdf::LocalScenario{
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42
        })
    };

    const auto bias{
        rt::RunBetaMDPDETest(test_input, MakeBetaOptions(alpha_r_list, false))
    };

    ASSERT_EQ(bias.mdpde.requested_alpha.size(), alpha_r_list.size());
    ExpectBiasStatisticSize(bias.ols);
    ExpectBiasStatisticSize(bias.mdpde.requested_alpha.front());
    EXPECT_FALSE(bias.mdpde.trained_alpha.has_value());
    EXPECT_FALSE(bias.mdpde.trained_alpha_average.has_value());
}

TEST(RHBMTesterTest, RunMuMDPDETestSkipsTrainedAlphaWhenDisabled)
{
    const std::vector<double> alpha_g_list{ 0.2 };
    const auto test_input{
        tdf::BuildGroupTestData(tdf::GroupScenario{
            12,
            10,
            MakeDistribution(rg::GaussianModel3D{ 1.0, 0.5, 0.1 }),
            MakeDistribution(rg::GaussianModel3D{ 1.5, 0.5, 0.1 }),
            0.1,
            2,
            33
        })
    };

    const auto bias{
        rt::RunMuMDPDETest(test_input, MakeMuOptions(alpha_g_list, false))
    };

    ASSERT_EQ(bias.mdpde.requested_alpha.size(), alpha_g_list.size());
    ExpectBiasStatisticSize(bias.median);
    ExpectBiasStatisticSize(bias.mdpde.requested_alpha.front());
    EXPECT_FALSE(bias.mdpde.trained_alpha.has_value());
    EXPECT_FALSE(bias.mdpde.trained_alpha_average.has_value());
}

TEST(RHBMTesterTest, RunBetaMDPDETestAllowsEmptyAlphaListWithoutTraining)
{
    const auto test_input{
        tdf::BuildLocalTestData(tdf::LocalScenario{
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42
        })
    };

    const auto bias{
        rt::RunBetaMDPDETest(test_input, MakeBetaOptions({}, false))
    };

    EXPECT_TRUE(bias.mdpde.requested_alpha.empty());
    ExpectBiasStatisticSize(bias.ols);
    EXPECT_FALSE(bias.mdpde.trained_alpha.has_value());
    EXPECT_FALSE(bias.mdpde.trained_alpha_average.has_value());
}

TEST(RHBMTesterTest, RunBetaMDPDETestAllowsEmptyAlphaListWithTraining)
{
    const auto test_input{
        tdf::BuildLocalTestData(tdf::LocalScenario{
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42
        })
    };

    const auto bias{
        rt::RunBetaMDPDETest(test_input, MakeBetaOptions({}, true))
    };

    EXPECT_TRUE(bias.mdpde.requested_alpha.empty());
    ExpectBiasStatisticSize(bias.ols);
    ASSERT_TRUE(bias.mdpde.trained_alpha.has_value());
    ExpectBiasStatisticSize(bias.mdpde.trained_alpha.value());
    ASSERT_TRUE(bias.mdpde.trained_alpha_average.has_value());
    EXPECT_GE(bias.mdpde.trained_alpha_average.value(), 0.0);
}

TEST(RHBMTesterTest, RunBetaMDPDETestRejectsNonFiniteTruth)
{
    auto test_input{
        tdf::BuildLocalTestData(tdf::LocalScenario{
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42
        })
    };
    test_input.gaus_true = rg::GaussianModel3D{
        std::numeric_limits<double>::quiet_NaN(),
        0.5,
        0.0
    };

    EXPECT_THROW(
        rt::RunBetaMDPDETest(test_input, MakeBetaOptions({}, true)),
        std::invalid_argument
    );
}
