#include <gtest/gtest.h>

#include <cmath>

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>

#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/LocalPotentialEntry.hpp"

namespace rg = rhbm_gem;

TEST(RHBMTypesTest, LocalAndGroupIntensityStayEquivalent)
{
    const rg::GaussianModel3D estimate{ 9.0, 1.5 };
    const rg::GaussianModel3DUncertainty standard_deviation{ 0.5, 0.2 };

    rg::LocalPotentialEntry local_entry;
    local_entry.SetEstimateMDPDE(estimate);
    local_entry.SetAnnotation(
        "component",
        rg::LocalPotentialAnnotation{
            rg::GaussianModel3DWithUncertainty{ estimate, standard_deviation }, false, 0.0 });

    rg::AtomGroupPotentialEntry group_entry;
    group_entry.SetGroupStatistics(42, {}, {}, estimate, standard_deviation, 0.0);
    const auto group_gaussian{ group_entry.GetPriorWithUncertainty(42) };

    EXPECT_DOUBLE_EQ(local_entry.GetEstimateMDPDE().Intensity(), estimate.Intensity());
    EXPECT_DOUBLE_EQ(local_entry.GetEstimateMDPDE().GetDisplayParameter(2), estimate.Intensity());
    EXPECT_DOUBLE_EQ(group_entry.GetPrior(42).Intensity(), estimate.Intensity());
    EXPECT_DOUBLE_EQ(group_entry.GetPrior(42).GetDisplayParameter(2), estimate.Intensity());
    EXPECT_DOUBLE_EQ(
        local_entry.FindAnnotation("component")->gaussian.GetDisplayStandardDeviation(2),
        group_gaussian.GetDisplayStandardDeviation(2));
}

TEST(RHBMTypesTest, PublicRHBMTypesHeaderExposesStableGaussianValueMath)
{
    const rg::GaussianModel3D estimate{ 9.0, 1.5 };
    const rg::GaussianModel3DUncertainty standard_deviation{ 0.5, 0.2 };
    const rg::GaussianModel3DWithUncertainty gaussian{ estimate, standard_deviation };

    const auto expected_intensity{
        estimate.GetAmplitude()
        * std::pow(Constants::two_pi * estimate.GetWidth() * estimate.GetWidth(), -1.5) };

    EXPECT_DOUBLE_EQ(estimate.GetDisplayParameter(2), expected_intensity);
    EXPECT_DOUBLE_EQ(estimate.GetModelParameter(2), 0.0);
    EXPECT_DOUBLE_EQ(estimate.Intensity(), expected_intensity);
    EXPECT_DOUBLE_EQ(gaussian.GetDisplayParameter(2), expected_intensity);
    EXPECT_DOUBLE_EQ(
        gaussian.GetStandardDeviationModel().GetModelParameter(
            rg::GaussianModel3D::InterceptIndex()),
        0.0);
    EXPECT_DOUBLE_EQ(
        gaussian.GetDisplayStandardDeviation(2),
        std::sqrt(
            std::pow(expected_intensity * standard_deviation.GetAmplitude() / estimate.GetAmplitude(), 2)
            + std::pow(
                -3.0 * estimate.GetAmplitude() * std::pow(Constants::two_pi, -1.5)
                    * std::pow(estimate.GetWidth(), -4) * standard_deviation.GetWidth(),
                2)));
}
