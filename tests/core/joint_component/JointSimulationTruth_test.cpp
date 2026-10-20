#include <gtest/gtest.h>

#include "core/detail/joint_component/Numerics.hpp"
#include "support/JointSimulationTruth.hpp"

#include <cmath>

namespace {
namespace joint=rhbm_gem::core::joint_component;
namespace simulation=joint_test::simulation;

TEST(JointSimulationTruthTest, MapScalingDividesACAndPreservesPhysicalWidth)
{
    const simulation::Truth original{12.0,1.25,-0.75};
    const auto normalized=simulation::ToNormalizedDomain(original,4.0);
    EXPECT_DOUBLE_EQ(normalized.amplitude,3.0);
    EXPECT_DOUBLE_EQ(normalized.width,original.width);
    EXPECT_DOUBLE_EQ(normalized.charge,-0.1875);
    const auto reconstructed=simulation::ToOriginalDomain(normalized,4.0);
    EXPECT_DOUBLE_EQ(reconstructed.amplitude,original.amplitude);
    EXPECT_DOUBLE_EQ(reconstructed.width,original.width);
    EXPECT_DOUBLE_EQ(reconstructed.charge,original.charge);
}

TEST(JointSimulationTruthTest, PredictionsReconstructInOriginalAndNormalizedDomains)
{
    const simulation::Truth original{6.0,0.8,-0.3};
    const double divisor=3.0;
    const auto normalized=simulation::ToNormalizedDomain(original,divisor);
    double original_prediction{},normalized_prediction{};
    for(const double square:{0.0,0.25,1.0,6.25})
    {
        const auto original_basis=joint::EvaluateKernel(square,original.width,2.5);
        const auto normalized_basis=joint::EvaluateKernel(square,normalized.width,2.5);
        original_prediction+=original.amplitude*original_basis.gaussian+original.charge*original_basis.charge;
        normalized_prediction+=normalized.amplitude*normalized_basis.gaussian+
            normalized.charge*normalized_basis.charge;
    }
    EXPECT_NEAR(normalized_prediction,original_prediction/divisor,1e-14);
    const auto restored=simulation::ToOriginalDomain(normalized,divisor);
    double restored_prediction{};
    for(const double square:{0.0,0.25,1.0,6.25})
    {
        const auto basis=joint::EvaluateKernel(square,restored.width,2.5);
        restored_prediction+=restored.amplitude*basis.gaussian+restored.charge*basis.charge;
    }
    EXPECT_NEAR(restored_prediction,original_prediction,1e-14);
}

TEST(JointSimulationTruthTest, ParameterRmseUsesTheDeclaredComparisonScale)
{
    const simulation::Truth sidecar_truth{9.0,1.1,-1.5};
    const auto comparison_truth=simulation::ToNormalizedDomain(sidecar_truth,3.0);
    const simulation::Truth estimate{3.0,1.1,-0.5};
    const double rmse_a=std::sqrt(std::pow(estimate.amplitude-comparison_truth.amplitude,2));
    const double rmse_b=std::sqrt(std::pow(estimate.width-comparison_truth.width,2));
    const double rmse_c=std::sqrt(std::pow(estimate.charge-comparison_truth.charge,2));
    EXPECT_DOUBLE_EQ(rmse_a,0.0);
    EXPECT_DOUBLE_EQ(rmse_b,0.0);
    EXPECT_DOUBLE_EQ(rmse_c,0.0);
    EXPECT_NE(estimate.amplitude,sidecar_truth.amplitude);
}

} // namespace
