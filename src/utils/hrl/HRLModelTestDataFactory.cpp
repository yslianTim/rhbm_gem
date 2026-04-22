#include <rhbm_gem/utils/hrl/HRLModelTestDataFactory.hpp>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/hrl/GaussianLinearizationService.hpp>
#include <rhbm_gem/utils/hrl/HRLDataTransform.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/GaussianResponseMath.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{
int ValidateGaussianParameterSize(int value)
{
    if (value != 3)
    {
        throw std::invalid_argument(
            "gaus_par_size must be 3 for GaussianPotentialSampler-backed test data");
    }
    return value;
}

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
} // namespace

HRLModelTestDataFactory::HRLModelTestDataFactory(
    int gaus_par_size,
    rhbm_gem::GaussianLinearizationSpec linearization_spec) :
    m_gaus_par_size{
        ValidateGaussianParameterSize(
            rhbm_gem::numeric_validation::RequirePositive(gaus_par_size, "gaus_par_size"))
    },
    m_linearization_spec{ std::move(linearization_spec) },
    m_fit_range_min{ 0.0 },
    m_fit_range_max{ 1.0 },
    m_potential_sampler{}
{
    rhbm_gem::numeric_validation::RequirePositive(
        m_linearization_spec.basis_size,
        "linearization_spec basis_size");
}

void HRLModelTestDataFactory::SetFittingRange(double x_min, double x_max)
{
    rhbm_gem::numeric_validation::RequireFiniteNonNegativeRange(x_min, x_max, "fitting range");
    m_fit_range_min = x_min;
    m_fit_range_max = x_max;
}

HRLBetaTestInput HRLModelTestDataFactory::BuildBetaTestInput(const BetaScenario & scenario) const
{
    rhbm_gem::numeric_validation::RequirePositive(
        scenario.sampling_entry_size,
        "sampling_entry_size");
    rhbm_gem::numeric_validation::RequirePositive(scenario.replica_size, "replica_size");
    ValidateGausParametersDimension(scenario.gaus_true);
    const auto gaussian_model{ BuildGaussianModel(scenario.gaus_true) };

    HRLBetaTestInput input;
    input.gaus_true = scenario.gaus_true;
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
                generator
            )
        };
        input.replica_datasets.emplace_back(
            HRLDataTransform::BuildMemberDataset(data_entry_list)
        );
    }

    return input;
}

HRLMuTestInput HRLModelTestDataFactory::BuildMuTestInput(const MuScenario & scenario) const
{
    rhbm_gem::numeric_validation::RequirePositive(scenario.member_size, "member_size");
    rhbm_gem::numeric_validation::RequirePositive(scenario.replica_size, "replica_size");
    ValidateGausParametersDimension(scenario.gaus_prior);
    ValidateGausParametersDimension(scenario.gaus_sigma);
    ValidateGausParametersDimension(scenario.outlier_prior);
    ValidateGausParametersDimension(scenario.outlier_sigma);

    HRLMuTestInput input;
    input.gaus_true = scenario.gaus_prior;
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
        input.replica_beta_matrices.emplace_back(BuildBetaMatrix(random_gaus_array));
    }

    return input;
}

HRLNeighborhoodTestInput HRLModelTestDataFactory::BuildNeighborhoodTestInput(
    const NeighborhoodScenario & scenario) const
{
    rhbm_gem::numeric_validation::RequirePositive(
        scenario.sampling_entry_size,
        "sampling_entry_size");
    rhbm_gem::numeric_validation::RequirePositive(scenario.replica_size, "replica_size");
    ValidateGausParametersDimension(scenario.gaus_true);
    const auto gaussian_model{ BuildGaussianModel(scenario.gaus_true) };

    const NeighborhoodSamplingOptions no_cut_options{
        scenario.radius_min,
        scenario.radius_max,
        scenario.neighbor_distance,
        scenario.neighbor_count,
        0.0,
        NeighborhoodRejectedPointPolicy::RemoveRejectedPoints
    };
    const NeighborhoodSamplingOptions cut_options{
        scenario.radius_min,
        scenario.radius_max,
        scenario.neighbor_distance,
        scenario.neighbor_count,
        scenario.rejected_angle,
        NeighborhoodRejectedPointPolicy::RemoveRejectedPoints
    };

    HRLNeighborhoodTestInput input;
    input.gaus_true = scenario.gaus_true;
    input.no_cut_datasets.reserve(static_cast<size_t>(scenario.replica_size));
    input.cut_datasets.reserve(static_cast<size_t>(scenario.replica_size));
    if (scenario.include_sampling_summary)
    {
        input.sampling_summaries.reserve(1);
        input.sampling_summaries.emplace_back(
            m_potential_sampler.GenerateNeighborhoodSamples(
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
                generator
            )
        };
        const auto cut_data_entry_list{
            BuildLinearDatasetWithNeighborhood(
                static_cast<size_t>(scenario.sampling_entry_size),
                gaussian_model,
                scenario.data_error_sigma,
                cut_options,
                generator
            )
        };
        input.no_cut_datasets.emplace_back(
            HRLDataTransform::BuildMemberDataset(no_cut_data_entry_list)
        );
        input.cut_datasets.emplace_back(
            HRLDataTransform::BuildMemberDataset(cut_data_entry_list)
        );
    }

    return input;
}

void HRLModelTestDataFactory::ValidateGausParametersDimension(const Eigen::VectorXd & gaus_par) const
{
    try
    {
        rhbm_gem::eigen_validation::RequireVectorSize(
            gaus_par,
            m_gaus_par_size,
            "gaus_par");
    }
    catch (const std::invalid_argument &)
    {
        throw std::invalid_argument(
            "model parameters size invalid, must be : " + std::to_string(m_gaus_par_size)
        );
    }
}

GaussianModel3D HRLModelTestDataFactory::BuildGaussianModel(const Eigen::VectorXd & gaus_par) const
{
    ValidateGausParametersDimension(gaus_par);
    return GaussianModel3D{ gaus_par(0), gaus_par(1), gaus_par(2) };
}

LocalPotentialSampleList HRLModelTestDataFactory::BuildGaussianSampling(
    size_t sampling_entry_size,
    const GaussianModel3D & model,
    double outlier_ratio,
    std::mt19937 & generator
) const
{
    auto sampling_entries{
        m_potential_sampler.GenerateRadialSamples(
            sampling_entry_size,
            model,
            m_fit_range_min,
            m_fit_range_max,
            generator
        )
    };
    std::uniform_real_distribution<> dist_outlier(0.0, 1.0);
    const auto outlier_response{
        0.5 * model.amplitude * std::pow(Constants::two_pi * std::pow(model.width, 2), -1.5)
    };
    for (auto & sampling_entry : sampling_entries)
    {
        if (dist_outlier(generator) < outlier_ratio)
        {
            sampling_entry.response = static_cast<float>(outlier_response);
        }
    }
    return sampling_entries;
}

SeriesPointList HRLModelTestDataFactory::BuildLinearDataset(
    const LocalPotentialSampleList & sampling_entries,
    const GaussianModel3D & model,
    double error_sigma,
    std::mt19937 & generator
) const
{
    const rhbm_gem::GaussianLinearizationService linearization_service{ m_linearization_spec };
    Eigen::VectorXd model_parameters{ Eigen::VectorXd::Zero(3) };
    model_parameters(0) = model.amplitude;
    model_parameters(1) = model.width;
    model_parameters(2) = model.intercept;
    auto linear_data_entry_list{
        linearization_service.BuildDatasetSeries(
            sampling_entries,
            m_fit_range_min,
            m_fit_range_max,
            rhbm_gem::GaussianLinearizationContext::FromModelParameters(model_parameters)
        )
    };
    const auto max_response{
        model.amplitude *
        rhbm_gem::GaussianResponseMath::GetGaussianResponseAtDistance(0.0, model.width)
    };
    std::normal_distribution<> dist_error(0.0, error_sigma * max_response);
    for (auto & data_entry : linear_data_entry_list)
    {
        data_entry.response += dist_error(generator);
    }
    return linear_data_entry_list;
}

SeriesPointList HRLModelTestDataFactory::BuildLinearDataset(
    size_t sampling_entry_size,
    const GaussianModel3D & model,
    double error_sigma,
    double outlier_ratio,
    std::mt19937 & generator
) const
{
    const auto sampling_entries{
        BuildGaussianSampling(sampling_entry_size, model, outlier_ratio, generator)
    };
    return BuildLinearDataset(sampling_entries, model, error_sigma, generator);
}

SeriesPointList HRLModelTestDataFactory::BuildLinearDatasetWithNeighborhood(
    size_t samples_per_radius,
    const GaussianModel3D & model,
    double error_sigma,
    const NeighborhoodSamplingOptions & options,
    std::mt19937 & generator
) const
{
    const auto sampling_entries{
        m_potential_sampler.GenerateNeighborhoodSamples(samples_per_radius, model, options)
    };
    return BuildLinearDataset(sampling_entries, model, error_sigma, generator);
}

Eigen::MatrixXd HRLModelTestDataFactory::BuildRandomGausParameters(
    int member_size,
    const Eigen::VectorXd & gaus_prior,
    const Eigen::VectorXd & gaus_sigma,
    const Eigen::VectorXd & outlier_prior,
    const Eigen::VectorXd & outlier_sigma,
    double outlier_ratio,
    std::mt19937 & generator) const
{
    std::uniform_real_distribution<> dist_outlier(0.0, 1.0);
    std::vector<std::normal_distribution<>> dist_gaus_list;
    std::vector<std::normal_distribution<>> dist_outlier_list;
    for (int p = 0; p < m_gaus_par_size; p++)
    {
        std::normal_distribution<> dist_gaus_par(gaus_prior(p), gaus_sigma(p));
        std::normal_distribution<> dist_outlier_par(outlier_prior(p), outlier_sigma(p));
        dist_gaus_list.emplace_back(dist_gaus_par);
        dist_outlier_list.emplace_back(dist_outlier_par);
    }

    Eigen::MatrixXd gaus_par_matrix{ Eigen::MatrixXd::Zero(m_gaus_par_size, member_size) };
    for (int i = 0; i < member_size; i++)
    {
        Eigen::VectorXd gaus_par{ Eigen::VectorXd::Zero(m_gaus_par_size) };
        Eigen::VectorXd outlier_par{ Eigen::VectorXd::Zero(m_gaus_par_size) };
        const bool outlier_flag{ dist_outlier(generator) < outlier_ratio };
        for (int p = 0; p < m_gaus_par_size; p++)
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

Eigen::MatrixXd HRLModelTestDataFactory::BuildBetaMatrix(const Eigen::MatrixXd & gaus_array) const
{
    const auto member_size{ static_cast<int>(gaus_array.cols()) };
    const rhbm_gem::GaussianLinearizationService linearization_service{ m_linearization_spec };
    Eigen::MatrixXd beta_matrix{
        Eigen::MatrixXd::Zero(m_linearization_spec.basis_size, member_size)
    };
    for (int i = 0; i < member_size; i++)
    {
        beta_matrix.col(i) = linearization_service.EncodeGaussianToBeta(gaus_array.col(i));
    }
    return beta_matrix;
}
