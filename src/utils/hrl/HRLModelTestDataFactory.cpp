#include <rhbm_gem/utils/hrl/HRLModelTestDataFactory.hpp>

#include <rhbm_gem/utils/hrl/HRLDataTransform.hpp>
#include <rhbm_gem/utils/math/GausLinearTransformHelper.hpp>

#include <stdexcept>
#include <string>

namespace
{
int ValidatePositive(int value, const char * name)
{
    if (value <= 0)
    {
        throw std::invalid_argument(std::string(name) + " must be positive value");
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

HRLModelTestDataFactory::HRLModelTestDataFactory(int gaus_par_size, int linear_basis_size) :
    m_gaus_par_size{ ValidatePositive(gaus_par_size, "gaus_par_size") },
    m_linear_basis_size{ ValidatePositive(linear_basis_size, "linear_basis_size") },
    m_data_generator{ m_gaus_par_size }
{
}

void HRLModelTestDataFactory::SetFittingRange(double x_min, double x_max)
{
    m_data_generator.SetFittingRange(x_min, x_max);
}

HRLBetaTestInput HRLModelTestDataFactory::BuildBetaTestInput(const BetaScenario & scenario) const
{
    ValidatePositive(scenario.sampling_entry_size, "sampling_entry_size");
    ValidatePositive(scenario.replica_size, "replica_size");
    ValidateGausParametersDimension(scenario.gaus_true);

    HRLBetaTestInput input;
    input.gaus_true = scenario.gaus_true;
    input.replica_datasets.reserve(static_cast<size_t>(scenario.replica_size));

    for (int i = 0; i < scenario.replica_size; i++)
    {
        auto generator{ BuildReplicaGenerator(i, scenario.random_seed) };
        const auto data_entry_list{
            m_data_generator.GenerateLinearDataset(
                static_cast<size_t>(scenario.sampling_entry_size),
                scenario.gaus_true,
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
    ValidatePositive(scenario.member_size, "member_size");
    ValidatePositive(scenario.replica_size, "replica_size");
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
    ValidatePositive(scenario.sampling_entry_size, "sampling_entry_size");
    ValidatePositive(scenario.replica_size, "replica_size");
    ValidateGausParametersDimension(scenario.gaus_true);

    const SimulationDataGenerator::NeighborhoodOptions no_cut_options{
        scenario.radius_min,
        scenario.radius_max,
        scenario.neighbor_distance,
        scenario.neighbor_count,
        0.0
    };
    const SimulationDataGenerator::NeighborhoodOptions cut_options{
        scenario.radius_min,
        scenario.radius_max,
        scenario.neighbor_distance,
        scenario.neighbor_count,
        scenario.rejected_angle
    };

    HRLNeighborhoodTestInput input;
    input.gaus_true = scenario.gaus_true;
    input.no_cut_datasets.reserve(static_cast<size_t>(scenario.replica_size));
    input.cut_datasets.reserve(static_cast<size_t>(scenario.replica_size));
    if (scenario.include_sampling_summary)
    {
        input.sampling_summaries.reserve(1);
        input.sampling_summaries.emplace_back(
            m_data_generator.GenerateGaussianSamplingWithNeighborhood(
                static_cast<size_t>(scenario.sampling_entry_size),
                scenario.gaus_true,
                SimulationDataGenerator::NeighborhoodOptions{
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
            m_data_generator.GenerateLinearDatasetWithNeighborhood(
                static_cast<size_t>(scenario.sampling_entry_size),
                scenario.gaus_true,
                scenario.data_error_sigma,
                no_cut_options,
                generator
            )
        };
        const auto cut_data_entry_list{
            m_data_generator.GenerateLinearDatasetWithNeighborhood(
                static_cast<size_t>(scenario.sampling_entry_size),
                scenario.gaus_true,
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
    if (gaus_par.rows() != m_gaus_par_size)
    {
        throw std::invalid_argument(
            "model parameters size invalid, must be : " + std::to_string(m_gaus_par_size)
        );
    }
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
    Eigen::MatrixXd beta_matrix{ Eigen::MatrixXd::Zero(m_linear_basis_size, member_size) };
    for (int i = 0; i < member_size; i++)
    {
        const Eigen::VectorXd gaus_par{ gaus_array.col(i) };
        const Eigen::VectorXd beta{
            GausLinearTransformHelper::BuildLinearModelCoefficentVector(gaus_par(0), gaus_par(1))
        };
        for (int j = 0; j < m_linear_basis_size; j++)
        {
            beta_matrix.col(i)(j) = beta(j);
        }
    }
    return beta_matrix;
}
