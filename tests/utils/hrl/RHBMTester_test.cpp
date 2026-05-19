#include <gtest/gtest.h>

#include <cmath>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <vector>

#include <rhbm_gem/utils/hrl/TestDataFactory.hpp>
#include <rhbm_gem/utils/hrl/RHBMTester.hpp>

namespace {
namespace rt = rhbm_gem::rhbm_tester;
namespace tdf = rhbm_gem::test_data_factory;


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

void ExpectBiasStatisticSize(const rt::BiasStatistics & bias)
{
    EXPECT_EQ(bias.mean.size(), rhbm_gem::GaussianModel3D::ParameterSize());
    EXPECT_EQ(bias.sigma.size(), rhbm_gem::GaussianModel3D::ParameterSize());
}

void ExpectBiasStatisticFinite(const rt::BiasStatistics & bias)
{
    for (Eigen::Index i = 0; i < bias.mean.size(); i++)
    {
        EXPECT_TRUE(std::isfinite(bias.mean(i)));
        EXPECT_TRUE(std::isfinite(bias.sigma(i)));
    }
}

} // namespace

TEST(RHBMTesterTest, RunBetaMDPDETestPopulatesBiasOutputs)
{
    const std::vector<double> alpha_r_list{ 0.0, 0.5 };
    const auto test_input{
        tdf::BuildBetaTestInput(tdf::BetaScenario{
            rhbm_gem::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42,
            alpha_r_list,
            true
        })
    };

    const auto bias{
        rt::RunBetaMDPDETest(test_input, 1)
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
        tdf::BuildMuTestInput(tdf::MuScenario{
            12,
            MakeVector({ 1.0, 0.5, 0.1 }),
            MakeVector({ 0.05, 0.025, 0.01 }),
            MakeVector({ 1.5, 0.5, 0.1 }),
            MakeVector({ 0.05, 0.025, 0.01 }),
            0.1,
            2,
            33,
            alpha_g_list,
            true
        })
    };

    const auto bias{
        rt::RunMuMDPDETest(test_input, 1)
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
        tdf::BuildBetaTestInput(tdf::BetaScenario{
            rhbm_gem::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42,
            alpha_r_list,
            false
        })
    };

    const auto bias{
        rt::RunBetaMDPDETest(test_input, 1)
    };

    ASSERT_EQ(bias.mdpde.requested_alpha.size(), alpha_r_list.size());
    ExpectBiasStatisticSize(bias.ols);
    ExpectBiasStatisticSize(bias.mdpde.requested_alpha.front());
    EXPECT_FALSE(bias.mdpde.trained_alpha.has_value());
    EXPECT_FALSE(bias.mdpde.trained_alpha_average.has_value());
}

TEST(RHBMTesterTest, RunBetaMDPDETestUsesOffsetTauGridFitModel)
{
    auto scenario{ tdf::BetaScenario{
        rhbm_gem::GaussianModel3D{ 1.0, 0.5, -0.1 },
        8,
        0.0,
        0.0,
        1,
        42,
        { 0.0 },
        false
    } };
    scenario.local_fit_model = rhbm_gem::LocalGaussianFitModel::OffsetTauGrid;
    const auto test_input{ tdf::BuildBetaTestInput(scenario) };

    const auto bias{ rt::RunBetaMDPDETest(test_input, 1) };

    ASSERT_EQ(bias.mdpde.requested_alpha.size(), 1u);
    ExpectBiasStatisticSize(bias.ols);
    ExpectBiasStatisticSize(bias.mdpde.requested_alpha.front());
    ExpectBiasStatisticFinite(bias.ols);
    ExpectBiasStatisticFinite(bias.mdpde.requested_alpha.front());
    EXPECT_FALSE(bias.mdpde.trained_alpha.has_value());
    EXPECT_FALSE(bias.mdpde.trained_alpha_average.has_value());
}

TEST(RHBMTesterTest, RunMuMDPDETestSkipsTrainedAlphaWhenDisabled)
{
    const std::vector<double> alpha_g_list{ 0.2 };
    const auto test_input{
        tdf::BuildMuTestInput(tdf::MuScenario{
            12,
            MakeVector({ 1.0, 0.5, 0.1 }),
            MakeVector({ 0.05, 0.025, 0.01 }),
            MakeVector({ 1.5, 0.5, 0.1 }),
            MakeVector({ 0.05, 0.025, 0.01 }),
            0.1,
            2,
            33,
            alpha_g_list,
            false
        })
    };

    const auto bias{
        rt::RunMuMDPDETest(test_input, 1)
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
        tdf::BuildBetaTestInput(tdf::BetaScenario{
            rhbm_gem::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42,
            {},
            false
        })
    };

    const auto bias{
        rt::RunBetaMDPDETest(test_input, 1)
    };

    EXPECT_TRUE(bias.mdpde.requested_alpha.empty());
    ExpectBiasStatisticSize(bias.ols);
    EXPECT_FALSE(bias.mdpde.trained_alpha.has_value());
    EXPECT_FALSE(bias.mdpde.trained_alpha_average.has_value());
}

TEST(RHBMTesterTest, RunBetaMDPDETestAllowsEmptyAlphaListWithTraining)
{
    const auto test_input{
        tdf::BuildBetaTestInput(tdf::BetaScenario{
            rhbm_gem::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42,
            {},
            true
        })
    };

    const auto bias{
        rt::RunBetaMDPDETest(test_input, 1)
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
        tdf::BuildBetaTestInput(tdf::BetaScenario{
            rhbm_gem::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42
        })
    };
    test_input.gaus_true = rhbm_gem::GaussianModel3D{
        std::numeric_limits<double>::quiet_NaN(),
        0.5,
        0.0
    };

    EXPECT_THROW(
        rt::RunBetaMDPDETest(test_input, 1),
        std::invalid_argument
    );
}
