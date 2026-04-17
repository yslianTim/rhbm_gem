#include <gtest/gtest.h>

#include <random>

#include <rhbm_gem/utils/math/GaussianPotentialSampler.hpp>

TEST(GaussianPotentialSamplerTest, NeighborhoodSamplingKeepsPointCountWhenRejectAngleIsZero)
{
    GaussianPotentialSampler sampler;

    const auto sampling_entries{
        sampler.GenerateNeighborhoodSamples(
            4,
            GaussianModel3D{ 1.0, 0.5, 0.0 },
            NeighborhoodSamplingOptions{
                0.0,
                1.0,
                2.0,
                1,
                0.0
            })
    };

    EXPECT_EQ(sampling_entries.size(), 44u);
}

TEST(GaussianPotentialSamplerTest, NeighborhoodSamplingRemovesPointsWhenAngleFilteringEnabled)
{
    GaussianPotentialSampler sampler;

    const auto unfiltered_sampling_entries{
        sampler.GenerateNeighborhoodSamples(
            4,
            GaussianModel3D{ 1.0, 0.5, 0.0 },
            NeighborhoodSamplingOptions{
                0.0,
                1.0,
                2.0,
                1,
                0.0
            })
    };
    const auto filtered_sampling_entries{
        sampler.GenerateNeighborhoodSamples(
            4,
            GaussianModel3D{ 1.0, 0.5, 0.0 },
            NeighborhoodSamplingOptions{
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

TEST(GaussianPotentialSamplerTest, NeighborhoodSamplingRejectsTooManyNeighbors)
{
    GaussianPotentialSampler sampler;

    EXPECT_THROW(
        sampler.GenerateNeighborhoodSamples(
            4,
            GaussianModel3D{ 1.0, 0.5, 0.0 },
            NeighborhoodSamplingOptions{
                0.0,
                1.0,
                2.0,
                5,
                0.0
            }),
        std::invalid_argument);
}

TEST(GaussianPotentialSamplerTest, RadialSamplingIsReproducibleWithFixedSeed)
{
    GaussianPotentialSampler sampler;
    std::mt19937 first_rng(42);
    std::mt19937 second_rng(42);

    const auto first_samples{
        sampler.GenerateRadialSamples(8, GaussianModel3D{ 1.0, 0.5, 0.0 }, 0.0, 1.0, first_rng)
    };
    const auto second_samples{
        sampler.GenerateRadialSamples(8, GaussianModel3D{ 1.0, 0.5, 0.0 }, 0.0, 1.0, second_rng)
    };

    ASSERT_EQ(first_samples.size(), second_samples.size());
    for (size_t i = 0; i < first_samples.size(); i++)
    {
        EXPECT_FLOAT_EQ(first_samples.at(i).distance, second_samples.at(i).distance);
        EXPECT_FLOAT_EQ(first_samples.at(i).response, second_samples.at(i).response);
    }
}

TEST(GaussianPotentialSamplerTest, InterceptContributesToRadialResponse)
{
    GaussianPotentialSampler sampler;
    std::mt19937 base_rng(7);
    std::mt19937 shifted_rng(7);

    const auto base_samples{
        sampler.GenerateRadialSamples(8, GaussianModel3D{ 1.0, 0.5, 0.0 }, 0.0, 1.0, base_rng)
    };
    const auto shifted_samples{
        sampler.GenerateRadialSamples(8, GaussianModel3D{ 1.0, 0.5, 0.25 }, 0.0, 1.0, shifted_rng)
    };

    ASSERT_EQ(base_samples.size(), shifted_samples.size());
    for (size_t i = 0; i < base_samples.size(); i++)
    {
        EXPECT_NEAR(
            shifted_samples.at(i).response - base_samples.at(i).response,
            0.25,
            1.0e-6
        );
    }
}
