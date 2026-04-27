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

void ExpectResidualStatisticSize(const rt::ResidualStatistics & residual)
{
    EXPECT_EQ(residual.mean.size(), rhbm_gem::GaussianModel3D::kParameterSize);
    EXPECT_EQ(residual.sigma.size(), rhbm_gem::GaussianModel3D::kParameterSize);
}

} // namespace

TEST(RHBMTesterTest, RunBetaMDPDETestPopulatesResidualOutputs)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());
    factory.SetFittingRange(0.0, 1.0);
    const std::vector<double> alpha_r_list{ 0.0, 0.5 };
    const auto test_input{
        factory.BuildBetaTestInput(tdf::TestDataFactory::BetaScenario{
            MakeVector({ 1.0, 0.5, 0.0 }),
            10,
            0.01,
            0.0,
            2,
            42,
            alpha_r_list,
            true
        })
    };

    rt::BetaMDPDETestResidual residual;
    const bool result{
        rt::RunBetaMDPDETest(
            residual,
            test_input,
            1)
    };

    ASSERT_TRUE(result);
    ASSERT_EQ(residual.mdpde.requested_alpha.size(), alpha_r_list.size());
    ExpectResidualStatisticSize(residual.ols);
    for (const auto & requested_alpha_residual : residual.mdpde.requested_alpha)
    {
        ExpectResidualStatisticSize(requested_alpha_residual);
    }
    ASSERT_TRUE(residual.mdpde.trained_alpha.has_value());
    ExpectResidualStatisticSize(residual.mdpde.trained_alpha.value());
}

TEST(RHBMTesterTest, RunMuMDPDETestPopulatesResidualOutputs)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());
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

    rt::MuMDPDETestResidual residual;
    const bool result{
        rt::RunMuMDPDETest(
            residual,
            test_input,
            1)
    };

    ASSERT_TRUE(result);
    ASSERT_EQ(residual.mdpde.requested_alpha.size(), alpha_g_list.size());
    ExpectResidualStatisticSize(residual.median);
    for (const auto & requested_alpha_residual : residual.mdpde.requested_alpha)
    {
        ExpectResidualStatisticSize(requested_alpha_residual);
    }
    ASSERT_TRUE(residual.mdpde.trained_alpha.has_value());
    ExpectResidualStatisticSize(residual.mdpde.trained_alpha.value());
}

TEST(RHBMTesterTest, RunBetaMDPDETestSkipsTrainedAlphaWhenDisabled)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());
    factory.SetFittingRange(0.0, 1.0);
    const std::vector<double> alpha_r_list{ 0.5 };
    const auto test_input{
        factory.BuildBetaTestInput(tdf::TestDataFactory::BetaScenario{
            MakeVector({ 1.0, 0.5, 0.0 }),
            10,
            0.01,
            0.0,
            2,
            42,
            alpha_r_list,
            false
        })
    };

    rt::BetaMDPDETestResidual residual;
    const bool result{
        rt::RunBetaMDPDETest(
            residual,
            test_input,
            1)
    };

    ASSERT_TRUE(result);
    ASSERT_EQ(residual.mdpde.requested_alpha.size(), alpha_r_list.size());
    ExpectResidualStatisticSize(residual.ols);
    ExpectResidualStatisticSize(residual.mdpde.requested_alpha.front());
    EXPECT_FALSE(residual.mdpde.trained_alpha.has_value());
}

TEST(RHBMTesterTest, RunMuMDPDETestSkipsTrainedAlphaWhenDisabled)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());
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

    rt::MuMDPDETestResidual residual;
    const bool result{
        rt::RunMuMDPDETest(
            residual,
            test_input,
            1)
    };

    ASSERT_TRUE(result);
    ASSERT_EQ(residual.mdpde.requested_alpha.size(), alpha_g_list.size());
    ExpectResidualStatisticSize(residual.median);
    ExpectResidualStatisticSize(residual.mdpde.requested_alpha.front());
    EXPECT_FALSE(residual.mdpde.trained_alpha.has_value());
}

TEST(RHBMTesterTest, RunBetaMDPDETestAllowsEmptyAlphaListWithoutTraining)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());
    factory.SetFittingRange(0.0, 1.0);
    const auto test_input{
        factory.BuildBetaTestInput(tdf::TestDataFactory::BetaScenario{
            MakeVector({ 1.0, 0.5, 0.0 }),
            10,
            0.01,
            0.0,
            2,
            42,
            {},
            false
        })
    };

    rt::BetaMDPDETestResidual residual;
    const bool result{
        rt::RunBetaMDPDETest(
            residual,
            test_input,
            1)
    };

    ASSERT_TRUE(result);
    EXPECT_TRUE(residual.mdpde.requested_alpha.empty());
    ExpectResidualStatisticSize(residual.ols);
    EXPECT_FALSE(residual.mdpde.trained_alpha.has_value());
}

TEST(RHBMTesterTest, RunBetaMDPDEWithNeighborhoodTestConsumesPreparedInputs)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());
    factory.SetFittingRange(0.0, 1.0);
    const auto test_input{
        factory.BuildNeighborhoodTestInput(tdf::TestDataFactory::NeighborhoodScenario{
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

    rt::NeighborhoodMDPDETestResidual residual;

    const bool result{
        rt::RunBetaMDPDEWithNeighborhoodTest(
            residual,
            test_input,
            1)
    };

    ASSERT_TRUE(result);
    ExpectResidualStatisticSize(residual.no_cut_ols);
    ExpectResidualStatisticSize(residual.no_cut_mdpde);
    ExpectResidualStatisticSize(residual.cut_mdpde);
    EXPECT_GT(residual.trained_alpha_r_average, 0.0);
}

TEST(RHBMTesterTest, RunBetaMDPDETestRejectsWrongSizedTruth)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());
    factory.SetFittingRange(0.0, 1.0);
    auto test_input{
        factory.BuildBetaTestInput(tdf::TestDataFactory::BetaScenario{
            MakeVector({ 1.0, 0.5, 0.0 }),
            10,
            0.01,
            0.0,
            2,
            42
        })
    };
    test_input.gaus_true = MakeVector({ 1.0, 0.5 });

    rt::BetaMDPDETestResidual residual;

    EXPECT_THROW(
        rt::RunBetaMDPDETest(
            residual,
            test_input,
            1),
        std::invalid_argument
    );
}
