#include <rhbm_gem/utils/hrl/TestDataFactory.hpp>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/LocalPotentialSampleScoring.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{
using namespace rhbm_gem;

struct NeighborhoodSamplingOptions
{
    double radius_min{ 0.0 };
    double radius_max{ 1.0 };
    double neighbor_distance{ 2.0 };
    size_t neighbor_count{ 1 };
    double reject_angle_deg{ 0.0 };
};

std::mt19937 BuildReplicaGenerator(
    int replica_index,
    const std::optional<std::uint32_t> & random_seed)
{
    if (random_seed.has_value())
    {
        std::seed_seq seed_sequence{
            *random_seed,
            static_cast<std::uint32_t>(replica_index)
        };
        return std::mt19937(seed_sequence);
    }

    std::random_device random_device;
    std::seed_seq seed_sequence{
        random_device(),
        random_device(),
        random_device(),
        static_cast<std::uint32_t>(replica_index)
    };
    return std::mt19937(seed_sequence);
}

void ValidateGaussianModel(const GaussianModel3D & model)
{
    numeric_validation::RequireFinite(model.amplitude, "GaussianModel3D amplitude");
    numeric_validation::RequireFinitePositive(model.width, "GaussianModel3D width");
    numeric_validation::RequireFinite(model.intercept, "GaussianModel3D intercept");
}

double ComputeGaussianResponseAtPoint3D(
    const Eigen::VectorXd & point,
    const Eigen::VectorXd & center,
    double width)
{
    const auto width_square{ width * width };
    return 1.0 / std::pow(Constants::two_pi * width_square, 1.5) *
        std::exp(-0.5 * (point - center).squaredNorm() / width_square);
}

double ComputeGaussianResponseWithNeighborhood3D(
    const Eigen::VectorXd & point,
    const Eigen::VectorXd & center,
    const std::vector<Eigen::VectorXd> & neighbor_center_list,
    double width)
{
    auto response{ ComputeGaussianResponseAtPoint3D(point, center, width) };
    for (const auto & neighbor_center : neighbor_center_list)
    {
        response += ComputeGaussianResponseAtPoint3D(point, neighbor_center, width);
    }
    return response;
}

std::vector<Eigen::VectorXd> BuildNeighborCenterList(
    const NeighborhoodSamplingOptions & options)
{
    constexpr size_t max_neighbor_count{ 4 };
    numeric_validation::RequireAtMost(
        options.neighbor_count,
        max_neighbor_count,
        "neighbor_count");

    std::vector<Eigen::VectorXd> neighbor_center_list;
    neighbor_center_list.reserve(options.neighbor_count);
    if (options.neighbor_count == 0)
    {
        return neighbor_center_list;
    }

    Eigen::VectorXd neighbor_center[max_neighbor_count];
    for (size_t i = 0; i < max_neighbor_count; i++)
    {
        neighbor_center[i] = Eigen::VectorXd::Zero(3);
    }

    if (options.neighbor_count <= 2)
    {
        neighbor_center[0] << 1.0, 0.0, 0.0;
        if (options.neighbor_count == 2)
        {
            neighbor_center[1] << -1.0, 0.0, 0.0;
        }
    }

    if (options.neighbor_count == 3)
    {
        neighbor_center[0] << 1.0, 0.0, 0.0;
        neighbor_center[1] << -0.5, std::sqrt(3) / 2.0, 0.0;
        neighbor_center[2] << -0.5, -std::sqrt(3) / 2.0, 0.0;
    }

    if (options.neighbor_count == 4)
    {
        neighbor_center[0] << 0.0, 0.0, 1.0;
        neighbor_center[1] << 0.0, 2.0 * std::sqrt(2) / 3.0, -1.0 / 3.0;
        neighbor_center[2] << -std::sqrt(6) / 3.0, -std::sqrt(2) / 3.0, -1.0 / 3.0;
        neighbor_center[3] << std::sqrt(6) / 3.0, -std::sqrt(2) / 3.0, -1.0 / 3.0;
    }

    for (size_t i = 0; i < options.neighbor_count; i++)
    {
        neighbor_center[i] *= options.neighbor_distance;
        neighbor_center_list.emplace_back(neighbor_center[i]);
    }

    return neighbor_center_list;
}

LocalPotentialSampleList GenerateRadialSamples(
    size_t sample_count,
    const GaussianModel3D & model,
    double distance_min,
    double distance_max,
    std::mt19937 & generator)
{
    numeric_validation::RequirePositive(sample_count, "sample_count");
    ValidateGaussianModel(model);
    numeric_validation::RequireFiniteNonNegativeRange(
        distance_min,
        distance_max,
        "distance range");

    std::uniform_real_distribution<> dist_distance(distance_min, distance_max);
    const auto sampling_scores{
        BuildDefaultLocalPotentialSampleScoreList(sample_count)
    };
    LocalPotentialSampleList sample_list;
    sample_list.reserve(sample_count);
    for (size_t i = 0; i < sample_count; i++)
    {
        const auto distance{ dist_distance(generator) };
        const auto response{
            model.ResponseAtDistance(distance)
        };
        sample_list.emplace_back(LocalPotentialSample{
            static_cast<float>(distance),
            static_cast<float>(response),
            sampling_scores.at(i)
        });
    }

    return sample_list;
}

LocalPotentialSampleList GenerateNeighborhoodSamples(
    size_t samples_per_radius,
    const GaussianModel3D & model,
    const NeighborhoodSamplingOptions & options)
{
    numeric_validation::RequirePositive(samples_per_radius, "samples_per_radius");
    ValidateGaussianModel(model);
    numeric_validation::RequireFiniteNonNegativeRange(
        options.radius_min,
        options.radius_max,
        "radius range");
    numeric_validation::RequireFinitePositive(
        options.neighbor_distance,
        "neighbor_distance");

    auto neighbor_center_list{ BuildNeighborCenterList(options) };
    const Eigen::VectorXd atom_center{ Eigen::VectorXd::Zero(3) };

    SphereSampler sampler;
    sampler.SetSamplingProfile(
        SphereSamplingProfile::FibonacciDeterministic(
            SphereDistanceRange{ options.radius_min, options.radius_max },
            0.1,
            static_cast<unsigned int>(samples_per_radius),
            false
        )
    );
    const auto sampling_points{ sampler.GenerateSamplingPoints({ 0.0f, 0.0f, 0.0f }) };
    const auto sampling_scores{
        BuildLocalPotentialSampleScoreList(
            sampling_points, neighbor_center_list, options.reject_angle_deg)
    };

    LocalPotentialSampleList sample_list;
    sample_list.reserve(static_cast<size_t>(std::count_if(
        sampling_scores.begin(),
        sampling_scores.end(),
        [](float score)
        {
            return score > 0.0f;
        })));

    for (size_t i = 0; i < sampling_points.size(); i++)
    {
        const auto & sampling_point{ sampling_points.at(i) };
        const auto sampling_score{ sampling_scores.at(i) };
        if (sampling_score <= 0.0f)
        {
            continue;
        }

        Eigen::VectorXd point{ Eigen::VectorXd::Zero(3) };
        point(0) = sampling_point.position[0];
        point(1) = sampling_point.position[1];
        point(2) = sampling_point.position[2];

        const auto response{
            model.amplitude * ComputeGaussianResponseWithNeighborhood3D(
                point,
                atom_center,
                neighbor_center_list,
                model.width
            ) +
            model.intercept
        };
        sample_list.emplace_back(LocalPotentialSample{
            sampling_point.distance,
            static_cast<float>(response),
            sampling_score,
            sampling_point.position
        });
    }

    return sample_list;
}

LocalPotentialSampleList BuildGaussianSampling(
    size_t sampling_entry_size,
    const GaussianModel3D & model,
    double outlier_ratio,
    double fit_range_min,
    double fit_range_max,
    std::mt19937 & generator)
{
    auto sampling_entries{
        GenerateRadialSamples(
            sampling_entry_size,
            model,
            fit_range_min,
            fit_range_max,
            generator
        )
    };
    std::uniform_real_distribution<> dist_outlier(0.0, 1.0);
    const auto outlier_response{ 0.5 * model.Intensity() };
    for (auto & sampling_entry : sampling_entries)
    {
        if (dist_outlier(generator) < outlier_ratio)
        {
            sampling_entry.response = static_cast<float>(outlier_response);
        }
    }
    return sampling_entries;
}

SeriesPointList BuildLinearDataset(
    const LocalPotentialSampleList & sampling_entries,
    const GaussianModel3D & model,
    double error_sigma,
    const linearization_service::LinearizationSpec & linearization_spec,
    double fit_range_min,
    double fit_range_max,
    std::mt19937 & generator)
{
    auto linear_data_entry_list{
        linearization_service::BuildDatasetSeries(
            linearization_spec,
            sampling_entries,
            fit_range_min,
            fit_range_max,
            linearization_service::LinearizationContext::FromModel(model)
        )
    };
    const auto max_response{ model.Intensity() };
    std::normal_distribution<> dist_error(0.0, error_sigma * max_response);
    for (auto & data_entry : linear_data_entry_list)
    {
        data_entry.response += dist_error(generator);
    }
    return linear_data_entry_list;
}

SeriesPointList BuildLinearDataset(
    size_t sampling_entry_size,
    const GaussianModel3D & model,
    double error_sigma,
    double outlier_ratio,
    const linearization_service::LinearizationSpec & linearization_spec,
    double fit_range_min,
    double fit_range_max,
    std::mt19937 & generator)
{
    const auto sampling_entries{
        BuildGaussianSampling(
            sampling_entry_size,
            model,
            outlier_ratio,
            fit_range_min,
            fit_range_max,
            generator)
    };
    return BuildLinearDataset(
        sampling_entries,
        model,
        error_sigma,
        linearization_spec,
        fit_range_min,
        fit_range_max,
        generator);
}

SeriesPointList BuildLinearDatasetWithNeighborhood(
    size_t samples_per_radius,
    const GaussianModel3D & model,
    double error_sigma,
    const NeighborhoodSamplingOptions & options,
    const linearization_service::LinearizationSpec & linearization_spec,
    double fit_range_min,
    double fit_range_max,
    std::mt19937 & generator)
{
    const auto sampling_entries{ GenerateNeighborhoodSamples(samples_per_radius, model, options) };
    return BuildLinearDataset(
        sampling_entries,
        model,
        error_sigma,
        linearization_spec,
        fit_range_min,
        fit_range_max,
        generator);
}

Eigen::MatrixXd BuildRandomGausParameters(
    int member_size,
    const Eigen::VectorXd & gaus_prior,
    const Eigen::VectorXd & gaus_sigma,
    const Eigen::VectorXd & outlier_prior,
    const Eigen::VectorXd & outlier_sigma,
    double outlier_ratio,
    std::mt19937 & generator)
{
    std::uniform_real_distribution<> dist_outlier(0.0, 1.0);
    std::vector<std::normal_distribution<>> dist_gaus_list;
    std::vector<std::normal_distribution<>> dist_outlier_list;
    for (int p = 0; p < GaussianModel3D::kParameterSize; p++)
    {
        std::normal_distribution<> dist_gaus_par(gaus_prior(p), gaus_sigma(p));
        std::normal_distribution<> dist_outlier_par(outlier_prior(p), outlier_sigma(p));
        dist_gaus_list.emplace_back(dist_gaus_par);
        dist_outlier_list.emplace_back(dist_outlier_par);
    }

    Eigen::MatrixXd gaus_par_matrix{
        Eigen::MatrixXd::Zero(GaussianModel3D::kParameterSize, member_size)
    };
    for (int i = 0; i < member_size; i++)
    {
        Eigen::VectorXd gaus_par{ Eigen::VectorXd::Zero(GaussianModel3D::kParameterSize) };
        Eigen::VectorXd outlier_par{ Eigen::VectorXd::Zero(GaussianModel3D::kParameterSize) };
        const bool outlier_flag{ dist_outlier(generator) < outlier_ratio };
        for (int p = 0; p < GaussianModel3D::kParameterSize; p++)
        {
            gaus_par(p) = dist_gaus_list.at(static_cast<size_t>(p))(generator);
            outlier_par(p) = dist_outlier_list.at(static_cast<size_t>(p))(generator);
            if (outlier_flag)
            {
                gaus_par(p) = outlier_par(p);
            }
        }
        gaus_par_matrix.col(i) = gaus_par;
    }
    return gaus_par_matrix;
}

Eigen::MatrixXd BuildBetaMatrix(
    const Eigen::MatrixXd & gaus_array,
    const linearization_service::LinearizationSpec & linearization_spec)
{
    const auto member_size{ static_cast<int>(gaus_array.cols()) };
    Eigen::MatrixXd beta_matrix{ Eigen::MatrixXd::Zero(linearization_spec.basis_size, member_size) };
    for (int i = 0; i < member_size; i++)
    {
        beta_matrix.col(i) = linearization_service::EncodeGaussianToBeta(
            linearization_spec, GaussianModel3D::FromVector(gaus_array.col(i)));
    }
    return beta_matrix;
}
} // namespace

namespace rhbm_gem::test_data_factory
{

TestDataFactory::TestDataFactory(linearization_service::LinearizationSpec linearization_spec) :
    m_linearization_spec{ std::move(linearization_spec) },
    m_fit_range_min{ 0.0 },
    m_fit_range_max{ 1.0 }
{
    numeric_validation::RequirePositive(m_linearization_spec.basis_size, "linearization_spec basis_size");
}

void TestDataFactory::SetFittingRange(double x_min, double x_max)
{
    numeric_validation::RequireFiniteNonNegativeRange(x_min, x_max, "fitting range");
    m_fit_range_min = x_min;
    m_fit_range_max = x_max;
}

RHBMBetaTestInput TestDataFactory::BuildBetaTestInput(const BetaScenario & scenario) const
{
    numeric_validation::RequirePositive(scenario.sampling_entry_size, "sampling_entry_size");
    numeric_validation::RequirePositive(scenario.replica_size, "replica_size");
    eigen_validation::RequireVectorSize(
        scenario.gaus_true,
        GaussianModel3D::kParameterSize,
        "scenario.gaus_true");
    const auto gaussian_model{ GaussianModel3D::FromVector(scenario.gaus_true) };

    RHBMBetaTestInput input;
    input.gaus_true = scenario.gaus_true;
    input.requested_alpha_r_list = scenario.requested_alpha_r_list;
    input.alpha_training = scenario.alpha_training;
    input.replica_datasets.reserve(static_cast<size_t>(scenario.replica_size));

    for (int i = 0; i < scenario.replica_size; i++)
    {
        auto generator{ BuildReplicaGenerator(i, scenario.random_seed) };
        const auto data_entry_list{
            BuildLinearDataset(
                static_cast<size_t>(scenario.sampling_entry_size),
                gaussian_model,
                scenario.data_error_sigma,
                scenario.outlier_ratio,
                m_linearization_spec,
                m_fit_range_min,
                m_fit_range_max,
                generator
            )
        };
        input.replica_datasets.emplace_back(rhbm_helper::BuildMemberDataset(data_entry_list));
    }

    return input;
}

RHBMMuTestInput TestDataFactory::BuildMuTestInput(const MuScenario & scenario) const
{
    numeric_validation::RequirePositive(scenario.member_size, "member_size");
    numeric_validation::RequirePositive(scenario.replica_size, "replica_size");
    eigen_validation::RequireVectorSize(
        scenario.gaus_prior,
        GaussianModel3D::kParameterSize,
        "scenario.gaus_prior");
    eigen_validation::RequireVectorSize(
        scenario.gaus_sigma,
        GaussianModel3D::kParameterSize,
        "scenario.gaus_sigma");
    eigen_validation::RequireVectorSize(
        scenario.outlier_prior,
        GaussianModel3D::kParameterSize,
        "scenario.outlier_prior");
    eigen_validation::RequireVectorSize(
        scenario.outlier_sigma,
        GaussianModel3D::kParameterSize,
        "scenario.outlier_sigma");

    RHBMMuTestInput input;
    input.gaus_true = scenario.gaus_prior;
    input.requested_alpha_g_list = scenario.requested_alpha_g_list;
    input.alpha_training = scenario.alpha_training;
    input.replica_beta_matrices.reserve(static_cast<size_t>(scenario.replica_size));

    for (int i = 0; i < scenario.replica_size; i++)
    {
        auto generator{ BuildReplicaGenerator(i, scenario.random_seed) };
        const auto random_gaus_array{
            BuildRandomGausParameters(
                scenario.member_size,
                scenario.gaus_prior,
                scenario.gaus_sigma,
                scenario.outlier_prior,
                scenario.outlier_sigma,
                scenario.outlier_ratio,
                generator
            )
        };
        input.replica_beta_matrices.emplace_back(
            BuildBetaMatrix(random_gaus_array, m_linearization_spec));
    }

    return input;
}

RHBMNeighborhoodTestInput TestDataFactory::BuildNeighborhoodTestInput(
    const NeighborhoodScenario & scenario) const
{
    numeric_validation::RequirePositive(scenario.sampling_entry_size, "sampling_entry_size");
    numeric_validation::RequirePositive(scenario.replica_size, "replica_size");
    eigen_validation::RequireVectorSize(
        scenario.gaus_true,
        GaussianModel3D::kParameterSize,
        "scenario.gaus_true");
    const auto gaussian_model{ GaussianModel3D::FromVector(scenario.gaus_true) };

    const NeighborhoodSamplingOptions no_cut_options{
        scenario.radius_min,
        scenario.radius_max,
        scenario.neighbor_distance,
        scenario.neighbor_count,
        0.0
    };
    const NeighborhoodSamplingOptions cut_options{
        scenario.radius_min,
        scenario.radius_max,
        scenario.neighbor_distance,
        scenario.neighbor_count,
        scenario.rejected_angle
    };

    RHBMNeighborhoodTestInput input;
    input.gaus_true = scenario.gaus_true;
    input.no_cut_datasets.reserve(static_cast<size_t>(scenario.replica_size));
    input.cut_datasets.reserve(static_cast<size_t>(scenario.replica_size));
    if (scenario.include_sampling_summary)
    {
        input.sampling_summaries.reserve(1);
        input.sampling_summaries.emplace_back(
            GenerateNeighborhoodSamples(
                static_cast<size_t>(scenario.sampling_entry_size),
                gaussian_model,
                NeighborhoodSamplingOptions{
                    scenario.summary_radius_min,
                    scenario.summary_radius_max,
                    scenario.neighbor_distance,
                    scenario.neighbor_count,
                    scenario.rejected_angle
                }
            )
        );
    }

    for (int i = 0; i < scenario.replica_size; i++)
    {
        auto generator{ BuildReplicaGenerator(i, scenario.random_seed) };
        const auto no_cut_data_entry_list{
            BuildLinearDatasetWithNeighborhood(
                static_cast<size_t>(scenario.sampling_entry_size),
                gaussian_model,
                scenario.data_error_sigma,
                no_cut_options,
                m_linearization_spec,
                m_fit_range_min,
                m_fit_range_max,
                generator
            )
        };
        const auto cut_data_entry_list{
            BuildLinearDatasetWithNeighborhood(
                static_cast<size_t>(scenario.sampling_entry_size),
                gaussian_model,
                scenario.data_error_sigma,
                cut_options,
                m_linearization_spec,
                m_fit_range_min,
                m_fit_range_max,
                generator
            )
        };
        input.no_cut_datasets.emplace_back(
            rhbm_helper::BuildMemberDataset(no_cut_data_entry_list)
        );
        input.cut_datasets.emplace_back(
            rhbm_helper::BuildMemberDataset(cut_data_entry_list)
        );
    }

    return input;
}

} // namespace rhbm_gem::test_data_factory
