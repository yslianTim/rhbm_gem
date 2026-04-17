#include <gtest/gtest.h>

#include <algorithm>
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

    EXPECT_EQ(sampling_entries.size(), 40u);
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

TEST(GaussianPotentialSamplerTest, NeighborhoodSamplingCanKeepRejectedPointsWithZeroScore)
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
    const auto removed_sampling_entries{
        sampler.GenerateNeighborhoodSamples(
            4,
            GaussianModel3D{ 1.0, 0.5, 0.0 },
            NeighborhoodSamplingOptions{
                0.0,
                1.0,
                2.0,
                1,
                120.0,
                NeighborhoodRejectedPointPolicy::RemoveRejectedPoints
            })
    };
    const auto zero_score_sampling_entries{
        sampler.GenerateNeighborhoodSamples(
            4,
            GaussianModel3D{ 1.0, 0.5, 0.0 },
            NeighborhoodSamplingOptions{
                0.0,
                1.0,
                2.0,
                1,
                120.0,
                NeighborhoodRejectedPointPolicy::KeepRejectedPointsWithZeroScore
            })
    };

    ASSERT_EQ(zero_score_sampling_entries.size(), unfiltered_sampling_entries.size());
    EXPECT_FALSE(std::none_of(
        zero_score_sampling_entries.begin(),
        zero_score_sampling_entries.end(),
        [](const LocalPotentialSample & sample)
        {
            return sample.score == 0.0f;
        }));
    EXPECT_EQ(
        static_cast<size_t>(std::count_if(
            zero_score_sampling_entries.begin(),
            zero_score_sampling_entries.end(),
            [](const LocalPotentialSample & sample)
            {
                return sample.score > 0.0f;
            })),
        removed_sampling_entries.size());

    for (size_t i = 0; i < zero_score_sampling_entries.size(); i++)
    {
        const auto & scored_sample{ zero_score_sampling_entries.at(i) };
        const auto & unfiltered_sample{ unfiltered_sampling_entries.at(i) };
        EXPECT_FLOAT_EQ(scored_sample.distance, unfiltered_sample.distance);
        ASSERT_TRUE(scored_sample.position.has_value());
        ASSERT_TRUE(unfiltered_sample.position.has_value());
        EXPECT_EQ(scored_sample.position.value(), unfiltered_sample.position.value());
    }
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
