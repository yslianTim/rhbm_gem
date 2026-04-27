#include <gtest/gtest.h>

#include <cmath>
#include <initializer_list>
#include <limits>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/hrl/TestDataFactory.hpp>

namespace {
namespace tdf = rhbm_gem::test_data_factory;
namespace rg = rhbm_gem;

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

double ComputeExpectedGaussianResponseAtDistance3D(double distance, double width)
{
    const auto width_square{ width * width };
    return 1.0 / std::pow(Constants::two_pi * width_square, 1.5) *
        std::exp(-0.5 * distance * distance / width_square);
}

void ExpectDatasetEquals(
    const rg::RHBMMemberDataset & lhs,
    const rg::RHBMMemberDataset & rhs)
{
    ASSERT_EQ(lhs.X.rows(), rhs.X.rows());
    ASSERT_EQ(lhs.X.cols(), rhs.X.cols());
    ASSERT_EQ(lhs.y.rows(), rhs.y.rows());
    ASSERT_EQ(lhs.score.rows(), rhs.score.rows());
    for (Eigen::Index row = 0; row < lhs.X.rows(); row++)
    {
        for (Eigen::Index col = 0; col < lhs.X.cols(); col++)
        {
            EXPECT_DOUBLE_EQ(lhs.X(row, col), rhs.X(row, col));
        }
        EXPECT_DOUBLE_EQ(lhs.y(row), rhs.y(row));
        EXPECT_DOUBLE_EQ(lhs.score(row), rhs.score(row));
    }
}

void ExpectMatrixEquals(const Eigen::MatrixXd & lhs, const Eigen::MatrixXd & rhs)
{
    ASSERT_EQ(lhs.rows(), rhs.rows());
    ASSERT_EQ(lhs.cols(), rhs.cols());
    for (Eigen::Index row = 0; row < lhs.rows(); row++)
    {
        for (Eigen::Index col = 0; col < lhs.cols(); col++)
        {
            EXPECT_DOUBLE_EQ(lhs(row, col), rhs(row, col));
        }
    }
}

void ExpectSamplingEntriesEquals(
    const LocalPotentialSampleList & lhs,
    const LocalPotentialSampleList & rhs)
{
    ASSERT_EQ(lhs.size(), rhs.size());
    for (size_t i = 0; i < lhs.size(); i++)
    {
        EXPECT_FLOAT_EQ(lhs.at(i).distance, rhs.at(i).distance);
        EXPECT_FLOAT_EQ(lhs.at(i).response, rhs.at(i).response);
        EXPECT_FLOAT_EQ(lhs.at(i).score, rhs.at(i).score);
        EXPECT_EQ(lhs.at(i).position.has_value(), rhs.at(i).position.has_value());
        if (lhs.at(i).position.has_value() && rhs.at(i).position.has_value())
        {
            EXPECT_EQ(lhs.at(i).position.value(), rhs.at(i).position.value());
        }
    }
}

} // namespace

TEST(TestDataFactoryTest, BuildBetaTestInputIsReproducibleWithFixedSeed)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());
    factory.SetFittingRange(0.0, 1.0);

    const auto scenario{ tdf::TestDataFactory::BetaScenario{
        MakeVector({ 1.0, 0.5, 0.0 }),
        8,
        0.05,
        0.1,
        3,
        42
    } };

    const auto first_input{ factory.BuildBetaTestInput(scenario) };
    const auto second_input{ factory.BuildBetaTestInput(scenario) };

    ASSERT_EQ(first_input.replica_datasets.size(), second_input.replica_datasets.size());
    for (size_t i = 0; i < first_input.replica_datasets.size(); i++)
    {
        ExpectDatasetEquals(
            first_input.replica_datasets.at(i),
            second_input.replica_datasets.at(i)
        );
    }
}

TEST(TestDataFactoryTest, BuildBetaTestInputChangesWhenOutlierPolicyChanges)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());
    factory.SetFittingRange(0.0, 1.0);

    const auto base_scenario{ tdf::TestDataFactory::BetaScenario{
        MakeVector({ 1.0, 0.5, 0.0 }),
        8,
        0.0,
        0.0,
        1,
        42
    } };

    auto outlier_scenario{ base_scenario };
    outlier_scenario.outlier_ratio = 1.0;

    const auto baseline_input{ factory.BuildBetaTestInput(base_scenario) };
    const auto outlier_input{ factory.BuildBetaTestInput(outlier_scenario) };

    ASSERT_EQ(baseline_input.replica_datasets.size(), 1u);
    ASSERT_EQ(outlier_input.replica_datasets.size(), 1u);
    EXPECT_NE(
        baseline_input.replica_datasets.front().y(0),
        outlier_input.replica_datasets.front().y(0)
    );
}

TEST(TestDataFactoryTest, BuildBetaTestInputChangesWhenNoisePolicyChanges)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());
    factory.SetFittingRange(0.0, 1.0);

    const auto noiseless_scenario{ tdf::TestDataFactory::BetaScenario{
        MakeVector({ 1.0, 0.5, 0.0 }),
        8,
        0.0,
        0.0,
        1,
        42
    } };

    auto noisy_scenario{ noiseless_scenario };
    noisy_scenario.data_error_sigma = 0.1;

    const auto noiseless_input{ factory.BuildBetaTestInput(noiseless_scenario) };
    const auto noisy_input{ factory.BuildBetaTestInput(noisy_scenario) };

    ASSERT_EQ(noiseless_input.replica_datasets.size(), 1u);
    ASSERT_EQ(noisy_input.replica_datasets.size(), 1u);
    EXPECT_NE(
        noiseless_input.replica_datasets.front().y(0),
        noisy_input.replica_datasets.front().y(0)
    );
}

TEST(TestDataFactoryTest, BuildBetaTestInputUsesExpectedZeroDistanceGaussianResponse)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());
    factory.SetFittingRange(0.0, 0.0);

    constexpr double amplitude{ 2.0 };
    constexpr double width{ 0.5 };
    constexpr double intercept{ 0.0 };
    const auto input{
        factory.BuildBetaTestInput(tdf::TestDataFactory::BetaScenario{
            MakeVector({ amplitude, width, intercept }),
            1,
            0.0,
            0.0,
            1,
            42
        })
    };

    ASSERT_EQ(input.replica_datasets.size(), 1u);
    const auto & dataset{ input.replica_datasets.front() };
    ASSERT_EQ(dataset.X.rows(), 1);
    ASSERT_EQ(dataset.X.cols(), 2);
    ASSERT_EQ(dataset.y.rows(), 1);
    EXPECT_DOUBLE_EQ(dataset.X(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(dataset.X(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(dataset.score(0), 1.0);

    const auto expected_response{
        amplitude * ComputeExpectedGaussianResponseAtDistance3D(0.0, width) + intercept
    };
    EXPECT_NEAR(std::log(static_cast<float>(expected_response)), dataset.y(0), 1.0e-7);
}

TEST(TestDataFactoryTest, BuildMuTestInputIsReproducibleWithFixedSeed)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());

    const auto scenario{ tdf::TestDataFactory::MuScenario{
        4,
        MakeVector({ 1.0, 0.5, 0.1 }),
        MakeVector({ 0.05, 0.025, 0.01 }),
        MakeVector({ 1.5, 0.5, 0.1 }),
        MakeVector({ 0.05, 0.025, 0.01 }),
        0.2,
        3,
        77
    } };

    const auto first_input{ factory.BuildMuTestInput(scenario) };
    const auto second_input{ factory.BuildMuTestInput(scenario) };

    ASSERT_EQ(first_input.replica_beta_matrices.size(), second_input.replica_beta_matrices.size());
    for (size_t i = 0; i < first_input.replica_beta_matrices.size(); i++)
    {
        ExpectMatrixEquals(
            first_input.replica_beta_matrices.at(i),
            second_input.replica_beta_matrices.at(i)
        );
    }
}

TEST(TestDataFactoryTest, BuildNeighborhoodTestInputProvidesPairedDatasetsAndSamplingSummary)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());
    factory.SetFittingRange(0.0, 1.0);

    const auto input{
        factory.BuildNeighborhoodTestInput(tdf::TestDataFactory::NeighborhoodScenario{
            MakeVector({ 1.0, 0.5, 0.0 }),
            8,
            0.05,
            0.0,
            1.0,
            2.0,
            1,
            120.0,
            true,
            0.0,
            4.0,
            2,
            11
        })
    };

    ASSERT_EQ(input.no_cut_input.replica_datasets.size(), 2u);
    ASSERT_EQ(input.cut_input.replica_datasets.size(), 2u);
    EXPECT_TRUE(input.no_cut_input.gaus_true.isApprox(MakeVector({ 1.0, 0.5, 0.0 })));
    EXPECT_TRUE(input.cut_input.gaus_true.isApprox(MakeVector({ 1.0, 0.5, 0.0 })));
    EXPECT_TRUE(input.no_cut_input.alpha_training);
    EXPECT_TRUE(input.cut_input.alpha_training);
    ASSERT_EQ(input.sampling_summaries.size(), 1u);
    ASSERT_FALSE(input.sampling_summaries.front().empty());
    for (const auto & sample : input.sampling_summaries.front())
    {
        EXPECT_TRUE(sample.position.has_value());
    }

    for (size_t i = 0; i < input.no_cut_input.replica_datasets.size(); i++)
    {
        EXPECT_LE(
            input.cut_input.replica_datasets.at(i).y.rows(),
            input.no_cut_input.replica_datasets.at(i).y.rows()
        );
    }
}

TEST(TestDataFactoryTest, BuildNeighborhoodTestInputSamplingSummaryIncludesNeighborContribution)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());
    factory.SetFittingRange(0.0, 1.0);

    constexpr double amplitude{ 2.0 };
    constexpr double width{ 0.5 };
    constexpr double intercept{ 0.25 };
    constexpr double neighbor_distance{ 1.0 };
    const auto input{
        factory.BuildNeighborhoodTestInput(tdf::TestDataFactory::NeighborhoodScenario{
            MakeVector({ amplitude, width, intercept }),
            1,
            0.0,
            0.0,
            1.0,
            neighbor_distance,
            1,
            0.0,
            true,
            0.0,
            0.0,
            1,
            11
        })
    };

    ASSERT_EQ(input.sampling_summaries.size(), 1u);
    const auto & summary{ input.sampling_summaries.front() };
    ASSERT_EQ(summary.size(), 1u);
    EXPECT_FLOAT_EQ(summary.front().distance, 0.0f);
    ASSERT_TRUE(summary.front().position.has_value());

    const auto expected_response{
        amplitude * (
            ComputeExpectedGaussianResponseAtDistance3D(0.0, width) +
            ComputeExpectedGaussianResponseAtDistance3D(neighbor_distance, width)) +
        intercept
    };
    EXPECT_NEAR(expected_response, summary.front().response, 1.0e-5);
}

TEST(TestDataFactoryTest, BuildNeighborhoodTestInputIsReproducibleWithFixedSeed)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());
    factory.SetFittingRange(0.0, 1.0);

    const auto scenario{ tdf::TestDataFactory::NeighborhoodScenario{
        MakeVector({ 1.0, 0.5, 0.0 }),
        8,
        0.05,
        0.0,
        1.0,
        2.0,
        1,
        120.0,
        true,
        0.0,
        4.0,
        2,
        11
    } };

    const auto first_input{ factory.BuildNeighborhoodTestInput(scenario) };
    const auto second_input{ factory.BuildNeighborhoodTestInput(scenario) };

    ASSERT_EQ(
        first_input.no_cut_input.replica_datasets.size(),
        second_input.no_cut_input.replica_datasets.size());
    ASSERT_EQ(
        first_input.cut_input.replica_datasets.size(),
        second_input.cut_input.replica_datasets.size());
    ASSERT_EQ(first_input.sampling_summaries.size(), second_input.sampling_summaries.size());
    for (size_t i = 0; i < first_input.no_cut_input.replica_datasets.size(); i++)
    {
        ExpectDatasetEquals(
            first_input.no_cut_input.replica_datasets.at(i),
            second_input.no_cut_input.replica_datasets.at(i));
        ExpectDatasetEquals(
            first_input.cut_input.replica_datasets.at(i),
            second_input.cut_input.replica_datasets.at(i));
    }
    for (size_t i = 0; i < first_input.sampling_summaries.size(); i++)
    {
        ExpectSamplingEntriesEquals(
            first_input.sampling_summaries.at(i),
            second_input.sampling_summaries.at(i));
    }
}

TEST(TestDataFactoryTest, BuildTestInputsRejectNonPositiveGaussianWidth)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());

    EXPECT_THROW(
        factory.BuildBetaTestInput(tdf::TestDataFactory::BetaScenario{
            MakeVector({ 1.0, 0.0, 0.0 }),
            8,
            0.05,
            0.1,
            1,
            42
        }),
        std::invalid_argument);
    EXPECT_THROW(
        factory.BuildNeighborhoodTestInput(tdf::TestDataFactory::NeighborhoodScenario{
            MakeVector({ 1.0, -0.5, 0.0 }),
            8,
            0.05,
            0.0,
            1.0,
            2.0,
            1,
            120.0,
            false,
            0.0,
            4.0,
            1,
            11
        }),
        std::invalid_argument);
}

TEST(TestDataFactoryTest, ConstructorRejectsInvalidNumericInputs)
{
    auto spec{ rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset() };
    spec.basis_size = 0;
    EXPECT_THROW(
        {
            tdf::TestDataFactory factory{ spec };
        },
        std::invalid_argument);
}

TEST(TestDataFactoryTest, BuildTestInputsRejectInvalidGaussianVectorSize)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());

    EXPECT_THROW(
        factory.BuildBetaTestInput(tdf::TestDataFactory::BetaScenario{
            MakeVector({ 1.0, 0.5 }),
            8,
            0.05,
            0.1,
            1,
            42
        }),
        std::invalid_argument);
    EXPECT_THROW(
        factory.BuildMuTestInput(tdf::TestDataFactory::MuScenario{
            4,
            MakeVector({ 1.0, 0.5 }),
            MakeVector({ 0.05, 0.025, 0.01 }),
            MakeVector({ 1.5, 0.5, 0.1 }),
            MakeVector({ 0.05, 0.025, 0.01 }),
            0.2,
            1,
            77
        }),
        std::invalid_argument);
    EXPECT_THROW(
        factory.BuildNeighborhoodTestInput(tdf::TestDataFactory::NeighborhoodScenario{
            MakeVector({ 1.0, 0.5 }),
            8,
            0.05,
            0.0,
            1.0,
            2.0,
            1,
            120.0,
            false,
            0.0,
            4.0,
            1,
            11
        }),
        std::invalid_argument);
}

TEST(TestDataFactoryTest, SetFittingRangeRejectsInvalidRange)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());

    EXPECT_THROW(factory.SetFittingRange(-1.0, 1.0), std::invalid_argument);
    EXPECT_THROW(factory.SetFittingRange(2.0, 1.0), std::invalid_argument);
    EXPECT_THROW(
        factory.SetFittingRange(0.0, std::numeric_limits<double>::infinity()),
        std::invalid_argument);
}

TEST(TestDataFactoryTest, BuildNeighborhoodTestInputRejectsInvalidSamplingInputs)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());

    EXPECT_THROW(
        factory.BuildNeighborhoodTestInput(tdf::TestDataFactory::NeighborhoodScenario{
            MakeVector({ 1.0, 0.5, 0.0 }),
            8,
            0.05,
            0.0,
            1.0,
            std::numeric_limits<double>::infinity(),
            1,
            120.0,
            false,
            0.0,
            4.0,
            1,
            11
        }),
        std::invalid_argument);
    EXPECT_THROW(
        factory.BuildNeighborhoodTestInput(tdf::TestDataFactory::NeighborhoodScenario{
            MakeVector({ 1.0, 0.5, 0.0 }),
            8,
            0.05,
            0.0,
            1.0,
            2.0,
            5,
            120.0,
            false,
            0.0,
            4.0,
            1,
            11
        }),
        std::invalid_argument);
    EXPECT_THROW(
        factory.BuildNeighborhoodTestInput(tdf::TestDataFactory::NeighborhoodScenario{
            MakeVector({ 1.0, 0.5, 0.0 }),
            8,
            0.05,
            2.0,
            1.0,
            2.0,
            1,
            120.0,
            false,
            0.0,
            4.0,
            1,
            11
        }),
        std::invalid_argument);
}

TEST(TestDataFactoryTest, BuildBetaTestInputRejectsNonPositiveScenarioSizes)
{
    tdf::TestDataFactory factory(
        rhbm_gem::linearization_service::LinearizationSpec::DefaultDataset());

    EXPECT_THROW(
        factory.BuildBetaTestInput(tdf::TestDataFactory::BetaScenario{
            MakeVector({ 1.0, 0.5, 0.0 }),
            0,
            0.05,
            0.1,
            3,
            42
        }),
        std::invalid_argument);
    EXPECT_THROW(
        factory.BuildBetaTestInput(tdf::TestDataFactory::BetaScenario{
            MakeVector({ 1.0, 0.5, 0.0 }),
            8,
            0.05,
            0.1,
            0,
            42
        }),
        std::invalid_argument);
}
