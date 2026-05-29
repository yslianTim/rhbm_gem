#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/core/TestDataFactory.hpp>

namespace {
namespace tdf = rhbm_gem::core;
namespace rg = rhbm_gem;

tdf::GaussianParameterDistribution MakeDistribution(
    const rg::GaussianModel3D & mean,
    const rg::GaussianModel3DUncertainty & sigma = rg::GaussianModel3DUncertainty{ 0.05, 0.025, 0.01 })
{
    return tdf::GaussianParameterDistribution{ mean, sigma };
}

double ComputeExpectedGaussianResponseAtDistance3D(double distance, double width)
{
    const auto width_square{ width * width };
    return 1.0 / std::pow(Constants::two_pi * width_square, 1.5) *
        std::exp(-0.5 * distance * distance / width_square);
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

} // namespace

TEST(TestDataFactoryTest, BuildLocalTestDataIsReproducibleWithFixedSeed)
{
    const auto scenario{ tdf::LocalScenario{
        rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
        8,
        0.05,
        0.1,
        3,
        42
    } };

    const auto first_input{ tdf::BuildLocalTestData(scenario) };
    const auto second_input{ tdf::BuildLocalTestData(scenario) };

    ASSERT_EQ(
        first_input.replica_sampling_entries.size(),
        second_input.replica_sampling_entries.size());
    for (size_t i = 0; i < first_input.replica_sampling_entries.size(); i++)
    {
        ExpectSamplingEntriesEquals(
            first_input.replica_sampling_entries.at(i),
            second_input.replica_sampling_entries.at(i)
        );
    }
}

TEST(TestDataFactoryTest, BuildLocalTestDataChangesWhenOutlierPolicyChanges)
{
    const auto base_scenario{ tdf::LocalScenario{
        rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
        8,
        0.0,
        0.0,
        1,
        42
    } };

    auto outlier_scenario{ base_scenario };
    outlier_scenario.outlier_ratio = 1.0;

    const auto baseline_input{ tdf::BuildLocalTestData(base_scenario) };
    const auto outlier_input{ tdf::BuildLocalTestData(outlier_scenario) };

    ASSERT_EQ(baseline_input.replica_sampling_entries.size(), 1u);
    ASSERT_EQ(outlier_input.replica_sampling_entries.size(), 1u);
    ASSERT_FALSE(baseline_input.replica_sampling_entries.front().empty());
    ASSERT_FALSE(outlier_input.replica_sampling_entries.front().empty());
    EXPECT_NE(
        baseline_input.replica_sampling_entries.front().front().response,
        outlier_input.replica_sampling_entries.front().front().response
    );
}

TEST(TestDataFactoryTest, BuildLocalTestDataChangesWhenNoisePolicyChanges)
{
    const auto noiseless_scenario{ tdf::LocalScenario{
        rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
        8,
        0.0,
        0.0,
        1,
        42
    } };

    auto noisy_scenario{ noiseless_scenario };
    noisy_scenario.data_error_sigma = 0.1;

    const auto noiseless_input{ tdf::BuildLocalTestData(noiseless_scenario) };
    const auto noisy_input{ tdf::BuildLocalTestData(noisy_scenario) };

    ASSERT_EQ(noiseless_input.replica_sampling_entries.size(), 1u);
    ASSERT_EQ(noisy_input.replica_sampling_entries.size(), 1u);
    ASSERT_FALSE(noiseless_input.replica_sampling_entries.front().empty());
    ASSERT_FALSE(noisy_input.replica_sampling_entries.front().empty());
    EXPECT_NE(
        noiseless_input.replica_sampling_entries.front().front().response,
        noisy_input.replica_sampling_entries.front().front().response
    );
}

TEST(TestDataFactoryTest, BuildLocalTestDataUsesDefaultSamplingDistanceRange)
{
    constexpr double amplitude{ 2.0 };
    constexpr double width{ 0.5 };
    constexpr double intercept{ 0.0 };
    const auto input{
        tdf::BuildLocalTestData(tdf::LocalScenario{
            rg::GaussianModel3D{ amplitude, width, intercept },
            1,
            0.0,
            0.0,
            1,
            42
        })
    };

    ASSERT_EQ(input.replica_sampling_entries.size(), 1u);
    ASSERT_EQ(input.replica_sampling_entries.front().size(), 1u);
    const auto & sample{ input.replica_sampling_entries.front().front() };
    EXPECT_GE(sample.point.distance, 0.0f);
    EXPECT_LE(sample.point.distance, 1.0f);

    const auto expected_response{
        amplitude * ComputeExpectedGaussianResponseAtDistance3D(sample.point.distance, width) +
            intercept
    };
    EXPECT_NEAR(expected_response, sample.response, 1.0e-7);
}

TEST(TestDataFactoryTest, BuildGroupTestDataIsReproducibleWithFixedSeed)
{
    const auto scenario{ tdf::GroupScenario{
        4,
        10,
        MakeDistribution(rg::GaussianModel3D{ 1.0, 0.5, 0.1 }),
        MakeDistribution(rg::GaussianModel3D{ 1.5, 0.5, 0.1 }),
        0.2,
        3,
        77
    } };

    const auto first_input{ tdf::BuildGroupTestData(scenario) };
    const auto second_input{ tdf::BuildGroupTestData(scenario) };

    ASSERT_EQ(
        first_input.replica_member_sampling_entries.size(),
        second_input.replica_member_sampling_entries.size());
    for (size_t i = 0; i < first_input.replica_member_sampling_entries.size(); i++)
    {
        ASSERT_EQ(
            first_input.replica_member_sampling_entries.at(i).size(),
            static_cast<size_t>(scenario.member_size));
        ASSERT_EQ(
            second_input.replica_member_sampling_entries.at(i).size(),
            static_cast<size_t>(scenario.member_size));
        for (size_t j = 0; j < first_input.replica_member_sampling_entries.at(i).size(); j++)
        {
            ExpectSamplingEntriesEquals(
                first_input.replica_member_sampling_entries.at(i).at(j),
                second_input.replica_member_sampling_entries.at(i).at(j));
        }
    }
}

TEST(TestDataFactoryTest, BuildGroupTestDataUsesDefaultSamplingEntrySize)
{
    tdf::GroupScenario scenario;
    scenario.member_size = 2;
    scenario.inlier_distribution = MakeDistribution(rg::GaussianModel3D{ 1.0, 0.5, 0.1 });
    scenario.outlier_distribution = MakeDistribution(rg::GaussianModel3D{ 1.5, 0.5, 0.1 });
    scenario.replica_size = 1;
    scenario.random_seed = 77;

    const auto input{ tdf::BuildGroupTestData(scenario) };

    ASSERT_EQ(input.replica_member_sampling_entries.size(), 1u);
    ASSERT_EQ(input.replica_member_sampling_entries.front().size(), 2u);
    for (const auto & member_samples : input.replica_member_sampling_entries.front())
    {
        EXPECT_EQ(member_samples.size(), 10u);
    }
}

TEST(TestDataFactoryTest, BuildTestDataRejectNonPositiveGaussianWidth)
{
    EXPECT_THROW(
        tdf::BuildLocalTestData(tdf::LocalScenario{
            rg::GaussianModel3D{ 1.0, 0.0, 0.0 },
            8,
            0.05,
            0.1,
            1,
            42
        }),
        std::invalid_argument);
    EXPECT_THROW(
        tdf::BuildLocalTestData(tdf::AtomModelScenario{
            Spot::O,
            rg::GaussianModel3D{ 1.0, -0.5, 0.0 },
            0.05,
            1,
            11
        }),
        std::invalid_argument);
}

TEST(TestDataFactoryTest, BuildGroupTestDataRejectsInvalidGaussianDistribution)
{
    EXPECT_THROW(
        tdf::BuildGroupTestData(tdf::GroupScenario{
            4,
            10,
            MakeDistribution(rg::GaussianModel3D{ 1.0, 0.0, 0.1 }),
            MakeDistribution(rg::GaussianModel3D{ 1.5, 0.5, 0.1 }),
            0.2,
            1,
            77
        }),
        std::invalid_argument);
    EXPECT_THROW(
        tdf::BuildGroupTestData(tdf::GroupScenario{
            4,
            10,
            MakeDistribution(
                rg::GaussianModel3D{ 1.0, 0.5, 0.1 },
                rg::GaussianModel3DUncertainty{ 0.05, -0.025, 0.01 }),
            MakeDistribution(rg::GaussianModel3D{ 1.5, 0.5, 0.1 }),
            0.2,
            1,
            77
        }),
        std::invalid_argument);
}

TEST(TestDataFactoryTest, BuildLocalTestDataRejectsNonPositiveScenarioSizes)
{
    EXPECT_THROW(
        tdf::BuildLocalTestData(tdf::LocalScenario{
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
            0,
            0.05,
            0.1,
            3,
            42
        }),
        std::invalid_argument);
    EXPECT_THROW(
        tdf::BuildLocalTestData(tdf::LocalScenario{
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
            8,
            0.05,
            0.1,
            0,
            42
        }),
        std::invalid_argument);
}
