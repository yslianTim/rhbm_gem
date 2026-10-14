#include <gtest/gtest.h>

#include "core/detail/LocalConditionalCovariance.hpp"
#include "core/detail/joint_component/Numerics.hpp"

#include <rhbm_gem/utils/hrl/LinearizationService.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>

namespace {
using namespace rhbm_gem;
using namespace rhbm_gem::core::detail;

struct Fixture
{
    LocalStageEstimate stage;
    LocalPotentialSampleList raw;
    PostFitPeelingResult peeling;
    GaussianModel3D point{ 1.8, 0.65, -0.3 };

    explicit Fixture(const std::vector<double> & distances)
    {
        stage.point = point;
        stage.source = {EstimateMethod::JointComponents, "17", "component-2", "run-3", FittingRole::Target};
        stage.convergence = JointCheckStatus::Passed;
        peeling.source = stage.source;
        for (std::size_t i = 0; i < distances.size(); ++i)
        {
            LocalPotentialSample sample;
            sample.point.distance = distances[i];
            const auto basis = joint_component::EvaluateKernel(distances[i] * distances[i], point.GetWidth(), 2.5);
            sample.response = point.GetAmplitude() * basis.gaussian + point.GetOffset() * basis.charge +
                0.01 * std::sin(static_cast<double>(i + 1));
            raw.push_back(sample);
            peeling.samples.push_back({sample.response, {}});
        }
    }
};

Eigen::Matrix2d FiniteDifferenceGroupCovariance(const Fixture & fixture)
{
    const auto & point = fixture.point;
    Eigen::MatrixXd jacobian(static_cast<Eigen::Index>(fixture.raw.size()), 2);
    double sum_squares = 0.0;
    const auto evaluate = [&](double amplitude, double width, double distance) {
        const auto basis = joint_component::EvaluateKernel(distance * distance, width, 2.5);
        return amplitude * basis.gaussian + point.GetOffset() * basis.charge;
    };
    for (std::size_t i = 0; i < fixture.raw.size(); ++i)
    {
        const auto distance = fixture.raw[i].point.distance;
        const double h_a = point.GetAmplitude() * 1.0e-6;
        const double h_b = point.GetWidth() * 1.0e-6;
        jacobian(static_cast<Eigen::Index>(i), 0) =
            (evaluate(point.GetAmplitude() + h_a, point.GetWidth(), distance) -
             evaluate(point.GetAmplitude() - h_a, point.GetWidth(), distance)) / (2.0 * h_a);
        jacobian(static_cast<Eigen::Index>(i), 1) =
            (evaluate(point.GetAmplitude(), point.GetWidth() + h_b, distance) -
             evaluate(point.GetAmplitude(), point.GetWidth() - h_b, distance)) / (2.0 * h_b);
        const double residual = fixture.peeling.samples[i].response.value() -
            evaluate(point.GetAmplitude(), point.GetWidth(), distance);
        sum_squares += residual * residual;
    }
    const double variance = sum_squares / static_cast<double>(fixture.raw.size() - 2);
    const Eigen::Matrix2d covariance_ab = variance *
        (jacobian.transpose() * jacobian).ldlt().solve(Eigen::Matrix2d::Identity());
    Eigen::Matrix2d transform;
    transform << 1.0 / point.GetAmplitude(), -3.0 / point.GetWidth(),
        0.0, -2.0 / std::pow(point.GetWidth(), 3);
    return transform * covariance_ab * transform.transpose();
}

TEST(JointLocalConditionalCovarianceTest, MatchesDenseFiniteDifferenceReference)
{
    Fixture fixture({0.08, 0.21, 0.39, 0.61, 0.88, 1.17, 1.49, 1.82, 2.16, 2.42});
    const auto covariance = ComputeLocalConditionalCovariance(
        fixture.stage, fixture.raw, fixture.peeling, true);
    ASSERT_EQ(covariance.status, EvidenceStatus::Available) << covariance.reason;
    ASSERT_TRUE(covariance.covariance);
    EXPECT_EQ(covariance.observation_count, fixture.raw.size());
    EXPECT_EQ(covariance.degrees_of_freedom, fixture.raw.size() - 2);
    EXPECT_TRUE(covariance.covariance->allFinite());
    EXPECT_TRUE(covariance.covariance->isApprox(covariance.covariance->transpose(), 1.0e-12));
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> spectrum(*covariance.covariance);
    ASSERT_EQ(spectrum.info(), Eigen::Success);
    EXPECT_GT(spectrum.eigenvalues()(0), 0.0);
    EXPECT_TRUE(covariance.covariance->isApprox(FiniteDifferenceGroupCovariance(fixture), 2.0e-5));

    const auto evidence = BuildLocalGroupEvidence(fixture.stage, &fixture.peeling, covariance);
    EXPECT_EQ(evidence.status, EvidenceStatus::Available);
    EXPECT_EQ(evidence.uncertainty_method, "local-conditional-iid-ls");
    ASSERT_TRUE(evidence.estimate);
    ASSERT_TRUE(evidence.covariance);
    EXPECT_TRUE(evidence.estimate->isApprox(
        linearization_service::EncodeGaussianToParameterVector(fixture.point), 1.0e-12));
    EXPECT_TRUE(evidence.covariance->isApprox(*covariance.covariance, 1.0e-12));
}

TEST(JointLocalConditionalCovarianceTest, RejectsRankDeficientAndInsufficientSamples)
{
    Fixture duplicate({0.4, 0.4, 0.4, 0.4, 0.4});
    EXPECT_EQ(ComputeLocalConditionalCovariance(duplicate.stage, duplicate.raw, duplicate.peeling, true).reason,
        "rank-deficient-local-design");

    Fixture short_samples({0.2, 0.7});
    EXPECT_EQ(ComputeLocalConditionalCovariance(short_samples.stage, short_samples.raw, short_samples.peeling, true).reason,
        "insufficient-observations");
}

TEST(JointLocalConditionalCovarianceTest, ExcludesUnavailableSamplesAndRejectsStaleSources)
{
    Fixture fixture({0.1, 0.3, 0.6, 0.9, 1.3, 1.8});
    fixture.peeling.samples[0].response.reset();
    fixture.peeling.samples[1].response.reset();
    fixture.peeling.samples[2].response.reset();
    const auto covariance = ComputeLocalConditionalCovariance(
        fixture.stage, fixture.raw, fixture.peeling, true);
    ASSERT_EQ(covariance.status, EvidenceStatus::Available) << covariance.reason;
    EXPECT_EQ(covariance.observation_count, 3u);

    fixture.peeling.source.run_id = "stale-run";
    const auto stale = ComputeLocalConditionalCovariance(fixture.stage, fixture.raw, fixture.peeling, true);
    EXPECT_EQ(stale.reason, "peeling-source-mismatch");

    fixture.peeling.source = fixture.stage.source;
    const auto missing_geometry = ComputeLocalConditionalCovariance(
        fixture.stage, fixture.raw, fixture.peeling, false);
    EXPECT_EQ(missing_geometry.reason, "sample-geometry-unavailable");
}

} // namespace
