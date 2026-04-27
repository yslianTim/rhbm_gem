#include <gtest/gtest.h>

#include <cmath>

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>

#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/LocalPotentialEntry.hpp"

namespace rg = rhbm_gem;

TEST(RHBMTypesTest, LocalAndGroupIntensityStayEquivalent)
{
    const rg::GaussianEstimate estimate{ 9.0, 1.5 };
    const rg::GaussianParameterUncertainty standard_deviation{ 0.5, 0.2 };

    rg::LocalPotentialEntry local_entry;
    local_entry.SetEstimateMDPDE(estimate);
    local_entry.SetAnnotation(
        "component",
        rg::LocalPotentialAnnotation{
            rg::GaussianEstimateWithUncertainty{ estimate, standard_deviation }, false, 0.0 });

    rg::AtomGroupPotentialEntry group_entry;
    group_entry.SetGroupStatistics(42, {}, {}, estimate, standard_deviation, 0.0);
    const auto group_gaussian{ group_entry.GetPriorWithUncertainty(42) };

    EXPECT_DOUBLE_EQ(local_entry.GetEstimateMDPDE().Intensity(), estimate.Intensity());
    EXPECT_DOUBLE_EQ(local_entry.GetEstimateMDPDE().GetParameter(2), estimate.Intensity());
    EXPECT_DOUBLE_EQ(group_entry.GetPrior(42).Intensity(), estimate.Intensity());
    EXPECT_DOUBLE_EQ(group_entry.GetPrior(42).GetParameter(2), estimate.Intensity());
    EXPECT_DOUBLE_EQ(
        local_entry.FindAnnotation("component")->gaussian.GetStandardDeviation(2),
        group_gaussian.GetStandardDeviation(2));
}

TEST(RHBMTypesTest, PublicRHBMTypesHeaderExposesStableGaussianValueMath)
{
    const rg::GaussianEstimate estimate{ 9.0, 1.5 };
    const rg::GaussianParameterUncertainty standard_deviation{ 0.5, 0.2 };
    const rg::GaussianEstimateWithUncertainty gaussian{ estimate, standard_deviation };

    const auto expected_intensity{
        estimate.amplitude
        * std::pow(Constants::two_pi * estimate.width * estimate.width, -1.5) };

    EXPECT_DOUBLE_EQ(estimate.GetParameter(2), expected_intensity);
    EXPECT_DOUBLE_EQ(estimate.Intensity(), expected_intensity);
    EXPECT_DOUBLE_EQ(gaussian.GetEstimate(2), expected_intensity);
    EXPECT_DOUBLE_EQ(
        gaussian.GetStandardDeviation(2),
        std::sqrt(
            std::pow(expected_intensity * standard_deviation.amplitude / estimate.amplitude, 2)
            + std::pow(
                -3.0 * estimate.amplitude * std::pow(Constants::two_pi, -1.5)
                    * std::pow(estimate.width, -4) * standard_deviation.width,
                2)));
}
