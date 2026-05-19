#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <map>
#include <vector>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/TestDataFactory.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>

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

double ComputeDistanceFromOrigin(const std::array<float, 3> & position)
{
    return std::sqrt(
        static_cast<double>(position[0]) * static_cast<double>(position[0]) +
        static_cast<double>(position[1]) * static_cast<double>(position[1]) +
        static_cast<double>(position[2]) * static_cast<double>(position[2]));
}

void ExpectDatasetEquals(
    const rg::RHBMMemberDataset & lhs,
    const rg::RHBMMemberDataset & rhs)
{
    ASSERT_EQ(lhs.X.rows(), rhs.X.rows());
    ASSERT_EQ(lhs.X.cols(), rhs.X.cols());
    ASSERT_EQ(lhs.y.rows(), rhs.y.rows());
    for (Eigen::Index row = 0; row < lhs.X.rows(); row++)
    {
        for (Eigen::Index col = 0; col < lhs.X.cols(); col++)
        {
            EXPECT_DOUBLE_EQ(lhs.X(row, col), rhs.X(row, col));
        }
        EXPECT_DOUBLE_EQ(lhs.y(row), rhs.y(row));
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
        EXPECT_FLOAT_EQ(lhs.at(i).point.distance, rhs.at(i).point.distance);
        EXPECT_FLOAT_EQ(lhs.at(i).response, rhs.at(i).response);
        EXPECT_EQ(lhs.at(i).point.position, rhs.at(i).point.position);
    }
}

void ExpectDatasetMatchesSamplingEntries(
    const rg::RHBMMemberDataset & dataset,
    const LocalPotentialSampleList & sampling_entries,
    double fit_range_min = 0.0,
    double fit_range_max = 1.0)
{
    const auto expected_dataset{
        rg::rhbm_helper::BuildMemberDataset(sampling_entries, fit_range_min, fit_range_max)
    };
    ExpectDatasetEquals(dataset, expected_dataset);
}

} // namespace

TEST(TestDataFactoryTest, BuildBetaTestInputIsReproducibleWithFixedSeed)
{
    const auto scenario{ tdf::BetaScenario{
        rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
        8,
        0.05,
        0.1,
        3,
        42
    } };

    const auto first_input{ tdf::BuildBetaTestInput(scenario) };
    const auto second_input{ tdf::BuildBetaTestInput(scenario) };

    ASSERT_EQ(
        first_input.replica_sampling_entries.size(),
        second_input.replica_sampling_entries.size());
    ASSERT_EQ(first_input.replica_datasets.size(), second_input.replica_datasets.size());
    for (size_t i = 0; i < first_input.replica_datasets.size(); i++)
    {
        ExpectSamplingEntriesEquals(
            first_input.replica_sampling_entries.at(i),
            second_input.replica_sampling_entries.at(i)
        );
        ExpectDatasetEquals(
            first_input.replica_datasets.at(i),
            second_input.replica_datasets.at(i)
        );
        ExpectDatasetMatchesSamplingEntries(
            first_input.replica_datasets.at(i),
            first_input.replica_sampling_entries.at(i)
        );
    }
}

TEST(TestDataFactoryTest, BuildBetaTestInputChangesWhenOutlierPolicyChanges)
{
    const auto base_scenario{ tdf::BetaScenario{
        rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
        8,
        0.0,
        0.0,
        1,
        42
    } };

    auto outlier_scenario{ base_scenario };
    outlier_scenario.outlier_ratio = 1.0;

    const auto baseline_input{ tdf::BuildBetaTestInput(base_scenario) };
    const auto outlier_input{ tdf::BuildBetaTestInput(outlier_scenario) };

    ASSERT_EQ(baseline_input.replica_datasets.size(), 1u);
    ASSERT_EQ(outlier_input.replica_datasets.size(), 1u);
    ASSERT_EQ(baseline_input.replica_sampling_entries.size(), 1u);
    ASSERT_EQ(outlier_input.replica_sampling_entries.size(), 1u);
    ASSERT_FALSE(baseline_input.replica_sampling_entries.front().empty());
    ASSERT_FALSE(outlier_input.replica_sampling_entries.front().empty());
    EXPECT_NE(
        baseline_input.replica_sampling_entries.front().front().response,
        outlier_input.replica_sampling_entries.front().front().response
    );
    EXPECT_NE(
        baseline_input.replica_datasets.front().y(0),
        outlier_input.replica_datasets.front().y(0)
    );
}

TEST(TestDataFactoryTest, BuildBetaTestInputChangesWhenNoisePolicyChanges)
{
    const auto noiseless_scenario{ tdf::BetaScenario{
        rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
        8,
        0.0,
        0.0,
        1,
        42
    } };

    auto noisy_scenario{ noiseless_scenario };
    noisy_scenario.data_error_sigma = 0.1;

    const auto noiseless_input{ tdf::BuildBetaTestInput(noiseless_scenario) };
    const auto noisy_input{ tdf::BuildBetaTestInput(noisy_scenario) };

    ASSERT_EQ(noiseless_input.replica_datasets.size(), 1u);
    ASSERT_EQ(noisy_input.replica_datasets.size(), 1u);
    ASSERT_EQ(noiseless_input.replica_sampling_entries.size(), 1u);
    ASSERT_EQ(noisy_input.replica_sampling_entries.size(), 1u);
    ASSERT_FALSE(noiseless_input.replica_sampling_entries.front().empty());
    ASSERT_FALSE(noisy_input.replica_sampling_entries.front().empty());
    EXPECT_NE(
        noiseless_input.replica_sampling_entries.front().front().response,
        noisy_input.replica_sampling_entries.front().front().response
    );
    EXPECT_NE(
        noiseless_input.replica_datasets.front().y(0),
        noisy_input.replica_datasets.front().y(0)
    );
}

TEST(TestDataFactoryTest, BuildBetaTestInputUsesExpectedZeroDistanceGaussianResponse)
{
    const auto options{ tdf::TestDataBuildOptions{
        0.0,
        0.0
    } };

    constexpr double amplitude{ 2.0 };
    constexpr double width{ 0.5 };
    constexpr double intercept{ 0.0 };
    const auto input{
        tdf::BuildBetaTestInput(tdf::BetaScenario{
            rg::GaussianModel3D{ amplitude, width, intercept },
            1,
            0.0,
            0.0,
            1,
            42
        }, options)
    };

    ASSERT_EQ(input.replica_datasets.size(), 1u);
    ASSERT_EQ(input.replica_sampling_entries.size(), 1u);
    ASSERT_EQ(input.replica_sampling_entries.front().size(), 1u);
    const auto & dataset{ input.replica_datasets.front() };
    ASSERT_EQ(dataset.X.rows(), 1);
    ASSERT_EQ(dataset.X.cols(), 2);
    ASSERT_EQ(dataset.y.rows(), 1);
    EXPECT_DOUBLE_EQ(dataset.X(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(dataset.X(0, 1), 0.0);

    const auto expected_response{
        amplitude * ComputeExpectedGaussianResponseAtDistance3D(0.0, width) + intercept
    };
    EXPECT_NEAR(expected_response, input.replica_sampling_entries.front().front().response, 1.0e-7);
    EXPECT_NEAR(std::log(static_cast<float>(expected_response)), dataset.y(0), 1.0e-7);
    ExpectDatasetMatchesSamplingEntries(
        dataset,
        input.replica_sampling_entries.front(),
        options.fit_range_min,
        options.fit_range_max);
}

TEST(TestDataFactoryTest, BuildMuTestInputIsReproducibleWithFixedSeed)
{
    const auto scenario{ tdf::MuScenario{
        4,
        MakeVector({ 1.0, 0.5, 0.1 }),
        MakeVector({ 0.05, 0.025, 0.01 }),
        MakeVector({ 1.5, 0.5, 0.1 }),
        MakeVector({ 0.05, 0.025, 0.01 }),
        0.2,
        3,
        77
    } };

    const auto first_input{ tdf::BuildMuTestInput(scenario) };
    const auto second_input{ tdf::BuildMuTestInput(scenario) };

    ASSERT_EQ(
        first_input.replica_member_sampling_entries.size(),
        second_input.replica_member_sampling_entries.size());
    ASSERT_EQ(first_input.replica_beta_matrices.size(), second_input.replica_beta_matrices.size());
    for (size_t i = 0; i < first_input.replica_beta_matrices.size(); i++)
    {
        ASSERT_EQ(
            first_input.replica_member_sampling_entries.at(i).size(),
            static_cast<size_t>(first_input.replica_beta_matrices.at(i).cols()));
        ASSERT_EQ(
            second_input.replica_member_sampling_entries.at(i).size(),
            static_cast<size_t>(second_input.replica_beta_matrices.at(i).cols()));
        for (size_t j = 0; j < first_input.replica_member_sampling_entries.at(i).size(); j++)
        {
            ExpectSamplingEntriesEquals(
                first_input.replica_member_sampling_entries.at(i).at(j),
                second_input.replica_member_sampling_entries.at(i).at(j));
        }
        ExpectMatrixEquals(
            first_input.replica_beta_matrices.at(i),
            second_input.replica_beta_matrices.at(i)
        );
    }
}

TEST(TestDataFactoryTest, BuildAtomNeighborhoodTestInputProvidesPairedDatasetsAndIsReproducible)
{
    const auto scenario{ tdf::AtomNeighborhoodScenario{
        tdf::AtomNeighborType::O,
        rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
        8,
        0.05,
        0.0,
        1.0,
        15.0,
        true,
        0.0,
        4.0,
        2,
        11,
        { 0.25 },
        false
    } };

    const auto input{ tdf::BuildAtomNeighborhoodTestInput(scenario) };
    const auto repeated_input{ tdf::BuildAtomNeighborhoodTestInput(scenario) };

    ASSERT_EQ(input.no_cut_input.replica_datasets.size(), 2u);
    ASSERT_EQ(input.cut_input.replica_datasets.size(), 2u);
    ASSERT_EQ(input.no_cut_input.replica_sampling_entries.size(), 2u);
    ASSERT_EQ(input.cut_input.replica_sampling_entries.size(), 2u);
    EXPECT_TRUE(input.no_cut_input.gaus_true.ToVector().isApprox(MakeVector({ 1.0, 0.5, 0.0 })));
    EXPECT_TRUE(input.cut_input.gaus_true.ToVector().isApprox(MakeVector({ 1.0, 0.5, 0.0 })));
    EXPECT_EQ(input.no_cut_input.requested_alpha_r_list, std::vector<double>{ 0.25 });
    EXPECT_EQ(input.cut_input.requested_alpha_r_list, std::vector<double>{ 0.25 });
    EXPECT_FALSE(input.no_cut_input.alpha_training);
    EXPECT_FALSE(input.cut_input.alpha_training);
    ASSERT_EQ(input.sampling_summaries.size(), 1u);
    ASSERT_FALSE(input.sampling_summaries.front().empty());
    for (const auto & sample : input.sampling_summaries.front())
    {
        EXPECT_NEAR(
            sample.point.distance,
            ComputeDistanceFromOrigin(sample.point.position),
            1.0e-5);
    }

    for (size_t i = 0; i < input.no_cut_input.replica_datasets.size(); i++)
    {
        EXPECT_LE(
            input.cut_input.replica_datasets.at(i).y.rows(),
            input.no_cut_input.replica_datasets.at(i).y.rows()
        );
        ExpectDatasetMatchesSamplingEntries(
            input.no_cut_input.replica_datasets.at(i),
            input.no_cut_input.replica_sampling_entries.at(i));
        ExpectDatasetMatchesSamplingEntries(
            input.cut_input.replica_datasets.at(i),
            input.cut_input.replica_sampling_entries.at(i));
        ExpectSamplingEntriesEquals(
            input.no_cut_input.replica_sampling_entries.at(i),
            repeated_input.no_cut_input.replica_sampling_entries.at(i));
        ExpectSamplingEntriesEquals(
            input.cut_input.replica_sampling_entries.at(i),
            repeated_input.cut_input.replica_sampling_entries.at(i));
        ExpectDatasetEquals(
            input.no_cut_input.replica_datasets.at(i),
            repeated_input.no_cut_input.replica_datasets.at(i));
        ExpectDatasetEquals(
            input.cut_input.replica_datasets.at(i),
            repeated_input.cut_input.replica_datasets.at(i));
    }
    for (size_t i = 0; i < input.sampling_summaries.size(); i++)
    {
        ExpectSamplingEntriesEquals(
            input.sampling_summaries.at(i),
            repeated_input.sampling_summaries.at(i));
    }
}

TEST(TestDataFactoryTest, BuildTestInputsRejectNonPositiveGaussianWidth)
{
    EXPECT_THROW(
        tdf::BuildBetaTestInput(tdf::BetaScenario{
            rg::GaussianModel3D{ 1.0, 0.0, 0.0 },
            8,
            0.05,
            0.1,
            1,
            42
        }),
        std::invalid_argument);
    EXPECT_THROW(
        tdf::BuildAtomNeighborhoodTestInput(tdf::AtomNeighborhoodScenario{
            tdf::AtomNeighborType::O,
            rg::GaussianModel3D{ 1.0, -0.5, 0.0 },
            8,
            0.05,
            0.0,
            1.0,
            15.0,
            false,
            0.0,
            4.0,
            1,
            11
        }),
        std::invalid_argument);
}

TEST(TestDataFactoryTest, BuildMuTestInputRejectsInvalidGaussianVectorSize)
{
    EXPECT_THROW(
        tdf::BuildMuTestInput(tdf::MuScenario{
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
}

TEST(TestDataFactoryTest, BuildBetaTestInputRejectsInvalidFittingRange)
{
    const auto scenario{ tdf::BetaScenario{
        rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
        8,
        0.05,
        0.1,
        1,
        42
    } };

    EXPECT_THROW(
        tdf::BuildBetaTestInput(
            scenario,
            tdf::TestDataBuildOptions{
                -1.0,
                1.0
            }),
        std::invalid_argument);
    EXPECT_THROW(
        tdf::BuildBetaTestInput(
            scenario,
            tdf::TestDataBuildOptions{
                2.0,
                1.0
            }),
        std::invalid_argument);
    EXPECT_THROW(
        tdf::BuildBetaTestInput(
            scenario,
            tdf::TestDataBuildOptions{
                0.0,
                std::numeric_limits<double>::infinity()
            }),
        std::invalid_argument);
}

TEST(TestDataFactoryTest, BuildBetaTestInputRejectsNonPositiveScenarioSizes)
{
    EXPECT_THROW(
        tdf::BuildBetaTestInput(tdf::BetaScenario{
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
            0,
            0.05,
            0.1,
            3,
            42
        }),
        std::invalid_argument);
    EXPECT_THROW(
        tdf::BuildBetaTestInput(tdf::BetaScenario{
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
            8,
            0.05,
            0.1,
            0,
            42
        }),
        std::invalid_argument);
}
