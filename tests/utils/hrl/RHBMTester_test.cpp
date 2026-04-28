#include <gtest/gtest.h>

#include <initializer_list>
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

} // namespace

TEST(RHBMTesterTest, RunBetaMDPDETestPopulatesBiasOutputs)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::AtomDecode());
    factory.SetFittingRange(0.0, 1.0);
    const std::vector<double> alpha_r_list{ 0.0, 0.5 };
    const auto test_input{
        factory.BuildBetaTestInput(tdf::TestDataFactory::BetaScenario{
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

    rt::BetaMDPDETestBias bias;
    const bool result{
        rt::RunBetaMDPDETest(
            bias,
            test_input,
            1)
    };

    ASSERT_TRUE(result);
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
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::AtomDecode());
    const std::vector<double> alpha_g_list{ 0.2 };
    const auto test_input{
        factory.BuildMuTestInput(tdf::TestDataFactory::MuScenario{
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

    rt::MuMDPDETestBias bias;
    const bool result{
        rt::RunMuMDPDETest(
            bias,
            test_input,
            1)
    };

    ASSERT_TRUE(result);
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
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::AtomDecode());
    factory.SetFittingRange(0.0, 1.0);
    const std::vector<double> alpha_r_list{ 0.5 };
    const auto test_input{
        factory.BuildBetaTestInput(tdf::TestDataFactory::BetaScenario{
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

    rt::BetaMDPDETestBias bias;
    const bool result{
        rt::RunBetaMDPDETest(
            bias,
            test_input,
            1)
    };

    ASSERT_TRUE(result);
    ASSERT_EQ(bias.mdpde.requested_alpha.size(), alpha_r_list.size());
    ExpectBiasStatisticSize(bias.ols);
    ExpectBiasStatisticSize(bias.mdpde.requested_alpha.front());
    EXPECT_FALSE(bias.mdpde.trained_alpha.has_value());
    EXPECT_FALSE(bias.mdpde.trained_alpha_average.has_value());
}

TEST(RHBMTesterTest, RunMuMDPDETestSkipsTrainedAlphaWhenDisabled)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::AtomDecode());
    const std::vector<double> alpha_g_list{ 0.2 };
    const auto test_input{
        factory.BuildMuTestInput(tdf::TestDataFactory::MuScenario{
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

    rt::MuMDPDETestBias bias;
    const bool result{
        rt::RunMuMDPDETest(
            bias,
            test_input,
            1)
    };

    ASSERT_TRUE(result);
    ASSERT_EQ(bias.mdpde.requested_alpha.size(), alpha_g_list.size());
    ExpectBiasStatisticSize(bias.median);
    ExpectBiasStatisticSize(bias.mdpde.requested_alpha.front());
    EXPECT_FALSE(bias.mdpde.trained_alpha.has_value());
    EXPECT_FALSE(bias.mdpde.trained_alpha_average.has_value());
}

TEST(RHBMTesterTest, RunBetaMDPDETestAllowsEmptyAlphaListWithoutTraining)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::AtomDecode());
    factory.SetFittingRange(0.0, 1.0);
    const auto test_input{
        factory.BuildBetaTestInput(tdf::TestDataFactory::BetaScenario{
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

    rt::BetaMDPDETestBias bias;
    const bool result{
        rt::RunBetaMDPDETest(
            bias,
            test_input,
            1)
    };

    ASSERT_TRUE(result);
    EXPECT_TRUE(bias.mdpde.requested_alpha.empty());
    ExpectBiasStatisticSize(bias.ols);
    EXPECT_FALSE(bias.mdpde.trained_alpha.has_value());
    EXPECT_FALSE(bias.mdpde.trained_alpha_average.has_value());
}

TEST(RHBMTesterTest, RunBetaMDPDETestAllowsEmptyAlphaListWithTraining)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::AtomDecode());
    factory.SetFittingRange(0.0, 1.0);
    const auto test_input{
        factory.BuildBetaTestInput(tdf::TestDataFactory::BetaScenario{
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

    rt::BetaMDPDETestBias bias;
    const bool result{
        rt::RunBetaMDPDETest(
            bias,
            test_input,
            1)
    };

    ASSERT_TRUE(result);
    EXPECT_TRUE(bias.mdpde.requested_alpha.empty());
    ExpectBiasStatisticSize(bias.ols);
    ASSERT_TRUE(bias.mdpde.trained_alpha.has_value());
    ExpectBiasStatisticSize(bias.mdpde.trained_alpha.value());
    ASSERT_TRUE(bias.mdpde.trained_alpha_average.has_value());
    EXPECT_GE(bias.mdpde.trained_alpha_average.value(), 0.0);
}

TEST(RHBMTesterTest, RunBetaMDPDETestRejectsWrongSizedTruth)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::AtomDecode());
    factory.SetFittingRange(0.0, 1.0);
    auto test_input{
        factory.BuildBetaTestInput(tdf::TestDataFactory::BetaScenario{
            rhbm_gem::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42
        })
    };
    test_input.gaus_true = MakeVector({ 1.0, 0.5 });

    rt::BetaMDPDETestBias bias;

    EXPECT_THROW(
        rt::RunBetaMDPDETest(
            bias,
            test_input,
            1),
        std::invalid_argument
    );
}
