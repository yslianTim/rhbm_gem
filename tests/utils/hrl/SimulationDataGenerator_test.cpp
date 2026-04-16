#include <gtest/gtest.h>

#include <initializer_list>
#include <random>

#include <rhbm_gem/utils/hrl/SimulationDataGenerator.hpp>

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

TEST(SimulationDataGeneratorTest, NeighborhoodSamplingKeepsPointCountWhenAngleIsZero)
{
    SimulationDataGenerator generator(3);
    generator.SetFittingRange(0.0, 1.0);

    const auto sampling_entries{
        generator.GenerateGaussianSamplingWithNeighborhood(
            4,
            MakeVector({ 1.0, 0.5, 0.0 }),
            SimulationDataGenerator::NeighborhoodOptions{
                0.0,
                1.0,
                2.0,
                1,
                0.0
            })
    };

    EXPECT_EQ(sampling_entries.size(), 44u);
}

TEST(SimulationDataGeneratorTest, NeighborhoodSamplingRemovesPointsWhenAngleFilteringEnabled)
{
    SimulationDataGenerator generator(3);
    generator.SetFittingRange(0.0, 1.0);

    const auto unfiltered_sampling_entries{
        generator.GenerateGaussianSamplingWithNeighborhood(
            4,
            MakeVector({ 1.0, 0.5, 0.0 }),
            SimulationDataGenerator::NeighborhoodOptions{
                0.0,
                1.0,
                2.0,
                1,
                0.0
            })
    };
    const auto filtered_sampling_entries{
        generator.GenerateGaussianSamplingWithNeighborhood(
            4,
            MakeVector({ 1.0, 0.5, 0.0 }),
            SimulationDataGenerator::NeighborhoodOptions{
                0.0,
                1.0,
                2.0,
                1,
                120.0
            })
    };

    ASSERT_FALSE(unfiltered_sampling_entries.empty());
    EXPECT_LT(filtered_sampling_entries.size(), unfiltered_sampling_entries.size());
}

TEST(SimulationDataGeneratorTest, NeighborhoodSamplingRejectsTooManyNeighbors)
{
    SimulationDataGenerator generator(3);

    EXPECT_THROW(
        generator.GenerateGaussianSamplingWithNeighborhood(
            4,
            MakeVector({ 1.0, 0.5, 0.0 }),
            SimulationDataGenerator::NeighborhoodOptions{
                0.0,
                1.0,
                2.0,
                5,
                0.0
            }),
        std::invalid_argument);
}

TEST(SimulationDataGeneratorTest, LinearDatasetGenerationIsReproducibleWithFixedSeed)
{
    SimulationDataGenerator generator(3);
    generator.SetFittingRange(0.0, 1.0);

    std::mt19937 first_rng(42);
    std::mt19937 second_rng(42);
    const auto gaus_par{ MakeVector({ 1.0, 0.5, 0.0 }) };

    const auto first_dataset{
        generator.GenerateLinearDataset(8, gaus_par, 0.05, 0.1, first_rng)
    };
    const auto second_dataset{
        generator.GenerateLinearDataset(8, gaus_par, 0.05, 0.1, second_rng)
    };

    ASSERT_EQ(first_dataset.size(), second_dataset.size());
    for (size_t i = 0; i < first_dataset.size(); i++)
    {
        EXPECT_DOUBLE_EQ(first_dataset.at(i).response, second_dataset.at(i).response);
        EXPECT_DOUBLE_EQ(first_dataset.at(i).weight, second_dataset.at(i).weight);
        EXPECT_EQ(
            first_dataset.at(i).basis_values.size(),
            second_dataset.at(i).basis_values.size()
        );
        for (size_t j = 0; j < first_dataset.at(i).basis_values.size(); j++)
        {
            EXPECT_DOUBLE_EQ(
                first_dataset.at(i).basis_values.at(j),
                second_dataset.at(i).basis_values.at(j)
            );
        }
    }
}
