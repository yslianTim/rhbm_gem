#include <gtest/gtest.h>

#include <initializer_list>
#include <limits>

#include <rhbm_gem/utils/hrl/HRLModelTestDataFactory.hpp>

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

void ExpectDatasetEquals(
    const HRLMemberDataset & lhs,
    const HRLMemberDataset & rhs)
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

} // namespace

TEST(HRLModelTestDataFactoryTest, BuildBetaTestInputIsReproducibleWithFixedSeed)
{
    HRLModelTestDataFactory factory(
        3,
        rhbm_gem::GaussianLinearizationSpec::DefaultDataset());
    factory.SetFittingRange(0.0, 1.0);

    const auto scenario{ HRLModelTestDataFactory::BetaScenario{
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

TEST(HRLModelTestDataFactoryTest, BuildBetaTestInputChangesWhenOutlierPolicyChanges)
{
    HRLModelTestDataFactory factory(
        3,
        rhbm_gem::GaussianLinearizationSpec::DefaultDataset());
    factory.SetFittingRange(0.0, 1.0);

    const auto base_scenario{ HRLModelTestDataFactory::BetaScenario{
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

TEST(HRLModelTestDataFactoryTest, BuildBetaTestInputChangesWhenNoisePolicyChanges)
{
    HRLModelTestDataFactory factory(
        3,
        rhbm_gem::GaussianLinearizationSpec::DefaultDataset());
    factory.SetFittingRange(0.0, 1.0);

    const auto noiseless_scenario{ HRLModelTestDataFactory::BetaScenario{
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

TEST(HRLModelTestDataFactoryTest, BuildMuTestInputIsReproducibleWithFixedSeed)
{
    HRLModelTestDataFactory factory(
        3,
        rhbm_gem::GaussianLinearizationSpec::DefaultDataset());

    const auto scenario{ HRLModelTestDataFactory::MuScenario{
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

TEST(HRLModelTestDataFactoryTest, BuildNeighborhoodTestInputProvidesPairedDatasetsAndSamplingSummary)
{
    HRLModelTestDataFactory factory(
        3,
        rhbm_gem::GaussianLinearizationSpec::DefaultDataset());
    factory.SetFittingRange(0.0, 1.0);

    const auto input{
        factory.BuildNeighborhoodTestInput(HRLModelTestDataFactory::NeighborhoodScenario{
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

    ASSERT_EQ(input.no_cut_datasets.size(), 2u);
    ASSERT_EQ(input.cut_datasets.size(), 2u);
    ASSERT_EQ(input.sampling_summaries.size(), 1u);
    ASSERT_FALSE(input.sampling_summaries.front().empty());

    for (size_t i = 0; i < input.no_cut_datasets.size(); i++)
    {
        EXPECT_LE(
            input.cut_datasets.at(i).y.rows(),
            input.no_cut_datasets.at(i).y.rows()
        );
    }
}

TEST(HRLModelTestDataFactoryTest, ConstructorRejectsInvalidNumericInputs)
{
    EXPECT_THROW(
        HRLModelTestDataFactory(
            0,
            rhbm_gem::GaussianLinearizationSpec::DefaultDataset()),
        std::invalid_argument);

    auto spec{ rhbm_gem::GaussianLinearizationSpec::DefaultDataset() };
    spec.basis_size = 0;
    EXPECT_THROW(HRLModelTestDataFactory(3, spec), std::invalid_argument);
}

TEST(HRLModelTestDataFactoryTest, SetFittingRangeRejectsInvalidRange)
{
    HRLModelTestDataFactory factory(
        3,
        rhbm_gem::GaussianLinearizationSpec::DefaultDataset());

    EXPECT_THROW(factory.SetFittingRange(-1.0, 1.0), std::invalid_argument);
    EXPECT_THROW(factory.SetFittingRange(2.0, 1.0), std::invalid_argument);
    EXPECT_THROW(
        factory.SetFittingRange(0.0, std::numeric_limits<double>::infinity()),
        std::invalid_argument);
}

TEST(HRLModelTestDataFactoryTest, BuildBetaTestInputRejectsNonPositiveScenarioSizes)
{
    HRLModelTestDataFactory factory(
        3,
        rhbm_gem::GaussianLinearizationSpec::DefaultDataset());

    EXPECT_THROW(
        factory.BuildBetaTestInput(HRLModelTestDataFactory::BetaScenario{
            MakeVector({ 1.0, 0.5, 0.0 }),
            0,
            0.05,
            0.1,
            3,
            42
        }),
        std::invalid_argument);
    EXPECT_THROW(
        factory.BuildBetaTestInput(HRLModelTestDataFactory::BetaScenario{
            MakeVector({ 1.0, 0.5, 0.0 }),
            8,
            0.05,
            0.1,
            0,
            42
        }),
        std::invalid_argument);
}
