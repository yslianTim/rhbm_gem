#include <rhbm_gem/utils/hrl/HRLModelTester.hpp>
#include <rhbm_gem/utils/hrl/HRLAlphaTrainer.hpp>
#include <rhbm_gem/utils/hrl/HRLDataTransform.hpp>
#include <rhbm_gem/utils/hrl/HRLModelAlgorithms.hpp>
#include <rhbm_gem/utils/math/GausLinearTransformHelper.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <cmath>
#include <limits>
#include <utility>
#include <memory>
#include <stdexcept>
#include <sstream>

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
} // namespace

HRLModelTester::HRLModelTester(int gaus_par_size, int linear_basis_size, int replica_size) :
    m_gaus_par_size{ ValidatePositive(gaus_par_size, "gaus_par_size") },
    m_linear_basis_size{ ValidatePositive(linear_basis_size, "linear_basis_size") },
    m_replica_size{ replica_size },
    m_x_min{ 0.0 }, m_x_max{ 1.0 }
{
}

Eigen::MatrixXd HRLModelTester::BuildRandomGausParameters(
    int member_size,
    const Eigen::VectorXd & gaus_prior,
    const Eigen::VectorXd & gaus_sigma,
    const Eigen::VectorXd & outlier_prior,
    const Eigen::VectorXd & outlier_sigma,
    double outlier_ratio)
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
        bool outlier_flag{ false };
        if (dist_outlier(m_generator) < outlier_ratio) outlier_flag = true;
        for (int p = 0; p < m_gaus_par_size; p++)
        {
            gaus_par(p) = dist_gaus_list.at(static_cast<size_t>(p))(m_generator);
            outlier_par(p) = dist_outlier_list.at(static_cast<size_t>(p))(m_generator);
            if (outlier_flag == true) gaus_par(p) = outlier_par(p);
        }
        gaus_par_matrix.col(i) = gaus_par;
    }
    return gaus_par_matrix;
}

Eigen::MatrixXd HRLModelTester::BuildBetaMatrix(const Eigen::MatrixXd & gaus_array)
{
    int member_size{ static_cast<int>(gaus_array.cols()) };
    Eigen::MatrixXd beta_matrix{ Eigen::MatrixXd::Zero(m_linear_basis_size, member_size) };
    for (int i = 0; i < member_size; i++)
    {
        Eigen::VectorXd gaus_par{ gaus_array.col(i) };
        Eigen::VectorXd beta{
            GausLinearTransformHelper::BuildLinearModelCoefficentVector(gaus_par(0), gaus_par(1))
        };
        for (int j = 0; j < m_linear_basis_size; j++) beta_matrix.col(i)(j) = beta(j);
    }
    return beta_matrix;
}

LocalPotentialSampleList HRLModelTester::BuildRandomGausSamplingEntry(
    size_t sampling_entry_size, const Eigen::VectorXd & gaus_par, double outlier_ratio)
{
    CheckGausParametersDimension(gaus_par);
    auto amplitude{ gaus_par(0) };
    auto width{ gaus_par(1) };
    //auto intersect{ gaus_par(2) };
    auto intersect{ 0.0 };
    std::uniform_real_distribution<> dist_distance(m_x_min, m_x_max);
    std::uniform_real_distribution<> dist_outlier(0.0, 1.0);
    LocalPotentialSampleList sampling_entry_list;
    sampling_entry_list.reserve(sampling_entry_size);
    for (size_t i = 0; i < sampling_entry_size; i++)
    {
        auto r{ dist_distance(m_generator) };
        auto y{
            amplitude * GausLinearTransformHelper::GetGaussianResponseAtDistance(r, width) + intersect
        };
        auto y_outlier{ 0.5 * amplitude * std::pow(Constants::two_pi * std::pow(width, 2), -1.5) };
        if (dist_outlier(m_generator) < outlier_ratio)
        {
            y = y_outlier;
        }
        sampling_entry_list.emplace_back(LocalPotentialSample{
            static_cast<float>(r),
            static_cast<float>(y)
        });
    }
    return sampling_entry_list;
}

LocalPotentialSampleList HRLModelTester::BuildRandomGausSamplingEntryWithNeighborhood(
    size_t sampling_entry_size,
    const Eigen::VectorXd & gaus_par,
    double neighbor_distance,
    size_t neighbor_count)
{
    CheckGausParametersDimension(gaus_par);
    const size_t max_neighbor_count{ 4 };
    if (neighbor_count > max_neighbor_count)
    {
        throw std::invalid_argument("neighbor_count should be less than or equal to 4");
    }
    auto amplitude{ gaus_par(0) };
    auto width{ gaus_par(1) };
    auto intersect{ 0.0 };
    Eigen::VectorXd atom_center{ Eigen::VectorXd::Zero(3) };
    Eigen::VectorXd neighbor_center[max_neighbor_count];
    for (size_t i = 0; i < max_neighbor_count; i++)
    {
        neighbor_center[i] = Eigen::VectorXd::Zero(3);
    }

    // 1 & 2 neighbor atoms
    if (neighbor_count <= 2)
    {
        neighbor_center[0] << 1.0, 0.0, 0.0;
        neighbor_center[1] << -1.0, 0.0, 0.0;
    }
    
    // 3 neighbor atoms (angle = 120 degree)
    if (neighbor_count == 3)
    {
        neighbor_center[0] << 1.0, 0.0, 0.0;
        neighbor_center[1] << -0.5, std::sqrt(3) / 2.0, 0.0;
        neighbor_center[2] << -0.5, -std::sqrt(3) / 2.0, 0.0;
    }

    // 4 neighbor atoms in tetrahedral arrangement
    if (neighbor_count == 4)
    {
        neighbor_center[0] << 0.0, 0.0, 1.0;
        neighbor_center[1] << 0.0, 2.0 * std::sqrt(2) / 3.0, -1.0 / 3.0;
        neighbor_center[2] << -std::sqrt(6) / 3.0, -std::sqrt(2) / 3.0, -1.0 / 3.0;
        neighbor_center[3] << std::sqrt(6) / 3.0, -std::sqrt(2) / 3.0, -1.0 / 3.0;
    }

    std::vector<Eigen::VectorXd> neighbor_center_list;
    neighbor_center_list.reserve(neighbor_count);
    for (size_t i = 0; i < neighbor_count; i++)
    {
        neighbor_center[i] *= neighbor_distance;
        neighbor_center_list.emplace_back(neighbor_center[i]);
    }

    SphereSampler sampler;
    sampler.SetSamplingProfile(
        //SphereSamplingProfile::RadiusUniformRandom(
        //    SphereDistanceRange{ m_x_min, m_x_max },
        //    static_cast<unsigned int>(sampling_entry_size))
        SphereSamplingProfile::FibonacciDeterministic(
            SphereDistanceRange{ m_x_min, m_x_max },
            0.1,
            static_cast<unsigned int>(sampling_entry_size))
    );
    const auto sampling_points{ sampler.GenerateSamplingPoints({ 0.0f, 0.0f, 0.0f }) };
    LocalPotentialSampleList sampling_entry_list;
    sampling_entry_list.reserve(sampling_entry_size);
    for (const auto & sampling_point : sampling_points)
    {
        Eigen::VectorXd point{ Eigen::VectorXd::Zero(3) };
        point(0) = sampling_point.position[0];
        point(1) = sampling_point.position[1];
        point(2) = sampling_point.position[2];
        auto y{
            amplitude * GausLinearTransformHelper::GetGaussianPesponseAtPointWithNeighborhood(
                point, atom_center, neighbor_center_list, width) + intersect
        };
        sampling_entry_list.emplace_back(LocalPotentialSample{
            sampling_point.distance,
            static_cast<float>(y),
            1.0f,
            sampling_point.position
        });
    }
    return sampling_entry_list;
}

std::vector<Eigen::VectorXd> HRLModelTester::BuildRandomLinearDataEntry(
    size_t sampling_entry_size,
    const Eigen::VectorXd & gaus_par,
    double error_sigma,
    double outlier_ratio)
{
    auto sampling_entries{
        BuildRandomGausSamplingEntry(static_cast<size_t>(sampling_entry_size), gaus_par, outlier_ratio)
    };
    auto linear_data_entry_list{
        GausLinearTransformHelper::MapValueTransform(sampling_entries, m_x_min, m_x_max)
    };
    auto max_response{
        gaus_par(0) * GausLinearTransformHelper::GetGaussianResponseAtDistance(0.0, gaus_par(1))
    };
    std::normal_distribution<> dist_error(0.0, error_sigma * max_response);
    for (auto & data_entry : linear_data_entry_list)
    {
        data_entry(m_linear_basis_size) += dist_error(m_generator);
    }
    return linear_data_entry_list;
}

std::vector<Eigen::VectorXd> HRLModelTester::BuildRandomLinearDataEntryWithNeighborhood(
    size_t sampling_entry_size,
    const Eigen::VectorXd & gaus_par,
    double error_sigma,
    double neighbor_distance)
{
    auto sampling_entries{
        BuildRandomGausSamplingEntryWithNeighborhood(
            static_cast<size_t>(sampling_entry_size), gaus_par, neighbor_distance, 1)
    };
    auto linear_data_entry_list{
        GausLinearTransformHelper::MapValueTransform(sampling_entries, m_x_min, m_x_max)
    };
    auto max_response{
        gaus_par(0) * GausLinearTransformHelper::GetGaussianResponseAtDistance(0.0, gaus_par(1))
    };
    std::normal_distribution<> dist_error(0.0, error_sigma * max_response);
    for (auto & data_entry : linear_data_entry_list)
    {
        data_entry(m_linear_basis_size) += dist_error(m_generator);
    }
    return linear_data_entry_list;
}

bool HRLModelTester::RunBetaMDPDETest(
    const std::vector<double> & alpha_r_list,
    std::vector<Eigen::VectorXd> & residual_mean_ols_list,
    std::vector<Eigen::VectorXd> & residual_mean_mdpde_list,
    std::vector<Eigen::VectorXd> & residual_sigma_ols_list,
    std::vector<Eigen::VectorXd> & residual_sigma_mdpde_list,
    const Eigen::VectorXd & gaus_true,
    int sampling_entry_size,
    double data_error_sigma,
    double outlier_ratio,
    int thread_size)
{
#ifndef USE_OPENMP
    (void)thread_size;
#endif

    auto local_alpha_r_list{ alpha_r_list };
    auto alpha_size{ local_alpha_r_list.size() + 1 }; // add one for training alpha_r
    std::vector<Eigen::MatrixXd> residual_matrix_ols_list(alpha_size);
    std::vector<Eigen::MatrixXd> residual_matrix_mdpde_list(alpha_size);
    residual_matrix_ols_list.assign(
        alpha_size, Eigen::MatrixXd::Zero(m_gaus_par_size, m_replica_size)
    );
    residual_matrix_mdpde_list.assign(
        alpha_size, Eigen::MatrixXd::Zero(m_gaus_par_size, m_replica_size)
    );
    residual_mean_ols_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));
    residual_mean_mdpde_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));
    residual_sigma_ols_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));
    residual_sigma_mdpde_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));

    const size_t subset_size_alpha_r{ 5 };
    std::vector<double> train_alpha_r_list{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5 };
    local_alpha_r_list.emplace_back(0.0);

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(thread_size)
#endif
    for (int i = 0; i < m_replica_size; i++)
    {
        auto data_entry_list{
            BuildRandomLinearDataEntry(
                static_cast<size_t>(sampling_entry_size), gaus_true, data_error_sigma, outlier_ratio
            )
        };

        auto error_array{
            HRLAlphaTrainer::EvaluateAlphaR(
                data_entry_list,
                subset_size_alpha_r,
                train_alpha_r_list
            )
        };
        int error_min_id;
        error_array.minCoeff(&error_min_id);
        auto alpha_r_train{ train_alpha_r_list.at(static_cast<size_t>(error_min_id)) };
        local_alpha_r_list.back() = alpha_r_train; // update training alpha_r for each replica

        const auto dataset{ HRLDataTransform::BuildMemberDataset(data_entry_list) };
        HRLExecutionOptions options;
        options.quiet_mode = true;
        options.thread_size = thread_size;

        for (size_t j = 0; j < alpha_size; j++)
        {
            const auto beta_result = HRLModelAlgorithms::EstimateBetaMDPDE(
                local_alpha_r_list.at(j),
                dataset.X,
                dataset.y,
                options
            );

            Eigen::VectorXd gaus_0{ Eigen::VectorXd::Zero(m_gaus_par_size) };
            auto gaus_ols{
                GausLinearTransformHelper::BuildGaus3DModel(beta_result.beta_ols, gaus_0)
            };
            auto gaus_mdpde{
                GausLinearTransformHelper::BuildGaus3DModel(beta_result.beta_mdpde, gaus_0)
            };
            residual_matrix_ols_list.at(j).col(i) = CalculateNormalizedResidual(gaus_ols, gaus_true);
            residual_matrix_mdpde_list.at(j).col(i) = CalculateNormalizedResidual(gaus_mdpde, gaus_true);
        }
    }

    for (size_t j = 0; j < alpha_size; j++)
    {
        residual_mean_ols_list.at(j) = residual_matrix_ols_list.at(j).rowwise().mean();
        residual_mean_mdpde_list.at(j) = residual_matrix_mdpde_list.at(j).rowwise().mean();
        residual_sigma_ols_list.at(j) =
            (residual_matrix_ols_list.at(j).colwise() - residual_mean_ols_list.at(j)).rowwise().norm()
            / std::sqrt(m_replica_size - 1);
        residual_sigma_mdpde_list.at(j) =
            (residual_matrix_mdpde_list.at(j).colwise() - residual_mean_mdpde_list.at(j)).rowwise().norm()
            / std::sqrt(m_replica_size - 1);
    }

    return true;
}

bool HRLModelTester::RunMuMDPDETest(
    const std::vector<double> & alpha_g_list,
    std::vector<Eigen::VectorXd> & residual_mean_median_list,
    std::vector<Eigen::VectorXd> & residual_mean_mdpde_list,
    std::vector<Eigen::VectorXd> & residual_sigma_median_list,
    std::vector<Eigen::VectorXd> & residual_sigma_mdpde_list,
    int member_size,
    const Eigen::VectorXd & gaus_prior,
    const Eigen::VectorXd & gaus_sigma,
    const Eigen::VectorXd & outlier_prior,
    const Eigen::VectorXd & outlier_sigma,
    double outlier_ratio,
    int thread_size)
{
#ifndef USE_OPENMP
    (void)thread_size;
#endif

    auto local_alpha_g_list{ alpha_g_list };
    auto alpha_size{ local_alpha_g_list.size() + 1 }; // add one for training alpha_g
    std::vector<Eigen::MatrixXd> residual_matrix_median_list(alpha_size);
    std::vector<Eigen::MatrixXd> residual_matrix_mdpde_list(alpha_size);
    residual_matrix_median_list.assign(
        alpha_size, Eigen::MatrixXd::Zero(m_gaus_par_size, m_replica_size)
    );
    residual_matrix_mdpde_list.assign(
        alpha_size, Eigen::MatrixXd::Zero(m_gaus_par_size, m_replica_size)
    );
    residual_mean_median_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));
    residual_mean_mdpde_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));
    residual_sigma_median_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));
    residual_sigma_mdpde_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));

    const size_t subset_size_alpha_g{ 10 };
    std::vector<double> train_alpha_g_list{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    local_alpha_g_list.emplace_back(0.0);

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(thread_size)
#endif
    for (int i = 0; i < m_replica_size; i++)
    {
        Eigen::VectorXd gaus_true{ gaus_prior };
        auto random_gaus_array{
            BuildRandomGausParameters(
                member_size, gaus_prior, gaus_sigma,
                outlier_prior, outlier_sigma, outlier_ratio
            )
        };

        auto beta_matrix{ BuildBetaMatrix(random_gaus_array) };
        std::vector<Eigen::VectorXd> train_data_entry_list;
        for (int m = 0; m < beta_matrix.cols(); m++) train_data_entry_list.emplace_back(beta_matrix.col(m));

        auto error_array{
            HRLAlphaTrainer::EvaluateAlphaG(
                train_data_entry_list,
                subset_size_alpha_g,
                train_alpha_g_list
            )
        };
        int error_min_id;
        error_array.minCoeff(&error_min_id);
        local_alpha_g_list.back() = train_alpha_g_list.at(static_cast<size_t>(error_min_id));
        
        HRLExecutionOptions options;
        options.quiet_mode = true;
        options.thread_size = thread_size;

        for (size_t j = 0; j < alpha_size; j++)
        {
            const auto mdpde_result = HRLModelAlgorithms::EstimateMuMDPDE(
                local_alpha_g_list.at(j),
                beta_matrix,
                options
            );
            const auto ols_result = HRLModelAlgorithms::EstimateMuMDPDE(
                0.0,
                beta_matrix,
                options
            );

            Eigen::VectorXd gaus_0{ Eigen::VectorXd::Zero(m_gaus_par_size) };
            auto gaus_ols{
                GausLinearTransformHelper::BuildGaus3DModel(ols_result.mu_mdpde, gaus_0)
            };
            auto gaus_mdpde{
                GausLinearTransformHelper::BuildGaus3DModel(mdpde_result.mu_mdpde, gaus_0)
            };
            residual_matrix_median_list.at(j).col(i) = CalculateNormalizedResidual(gaus_ols, gaus_true);
            residual_matrix_mdpde_list.at(j).col(i) = CalculateNormalizedResidual(gaus_mdpde, gaus_true);
        }
    }

    for (size_t j = 0; j < alpha_size; j++)
    {
        residual_mean_median_list.at(j) = residual_matrix_median_list.at(j).rowwise().mean();
        residual_mean_mdpde_list.at(j) = residual_matrix_mdpde_list.at(j).rowwise().mean();
        residual_sigma_median_list.at(j) =
            (residual_matrix_median_list.at(j).colwise() - residual_mean_median_list.at(j)).rowwise().norm()
            / std::sqrt(m_replica_size - 1);
        residual_sigma_mdpde_list.at(j) =
            (residual_matrix_mdpde_list.at(j).colwise() - residual_mean_mdpde_list.at(j)).rowwise().norm()
            / std::sqrt(m_replica_size - 1);
    }

    return true;
}

bool HRLModelTester::RunBetaMDPDEWithNeighborhoodTest(
    const std::vector<double> & alpha_r_list,
    std::vector<Eigen::VectorXd> & residual_mean_ols_list,
    std::vector<Eigen::VectorXd> & residual_mean_mdpde_list,
    std::vector<Eigen::VectorXd> & residual_sigma_ols_list,
    std::vector<Eigen::VectorXd> & residual_sigma_mdpde_list,
    const Eigen::VectorXd & gaus_true,
    int sampling_entry_size,
    double data_error_sigma,
    double neighbor_distance,
    int thread_size)
{
#ifndef USE_OPENMP
    (void)thread_size;
#endif

    auto local_alpha_r_list{ alpha_r_list };
    auto alpha_size{ local_alpha_r_list.size() + 1 }; // add one for training alpha_r
    std::vector<Eigen::MatrixXd> residual_matrix_ols_list(alpha_size);
    std::vector<Eigen::MatrixXd> residual_matrix_mdpde_list(alpha_size);
    residual_matrix_ols_list.assign(
        alpha_size, Eigen::MatrixXd::Zero(m_gaus_par_size, m_replica_size)
    );
    residual_matrix_mdpde_list.assign(
        alpha_size, Eigen::MatrixXd::Zero(m_gaus_par_size, m_replica_size)
    );
    residual_mean_ols_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));
    residual_mean_mdpde_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));
    residual_sigma_ols_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));
    residual_sigma_mdpde_list.assign(alpha_size, Eigen::VectorXd::Zero(m_gaus_par_size));

    const size_t subset_size_alpha_r{ 5 };
    std::vector<double> train_alpha_r_list{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5 };
    local_alpha_r_list.emplace_back(0.0);

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(thread_size)
#endif
    for (int i = 0; i < m_replica_size; i++)
    {
        auto data_entry_list{
            BuildRandomLinearDataEntryWithNeighborhood(
                static_cast<size_t>(sampling_entry_size), gaus_true, data_error_sigma, neighbor_distance
            )
        };

        auto error_array{
            HRLAlphaTrainer::EvaluateAlphaR(
                data_entry_list,
                subset_size_alpha_r,
                train_alpha_r_list
            )
        };
        int error_min_id;
        error_array.minCoeff(&error_min_id);
        auto alpha_r_train{ train_alpha_r_list.at(static_cast<size_t>(error_min_id)) };
        local_alpha_r_list.back() = alpha_r_train; // update training alpha_r for each replica

        const auto dataset{ HRLDataTransform::BuildMemberDataset(data_entry_list) };
        HRLExecutionOptions options;
        options.quiet_mode = true;
        options.thread_size = thread_size;

        for (size_t j = 0; j < alpha_size; j++)
        {
            const auto beta_result = HRLModelAlgorithms::EstimateBetaMDPDE(
                local_alpha_r_list.at(j),
                dataset.X,
                dataset.y,
                options
            );

            Eigen::VectorXd gaus_0{ Eigen::VectorXd::Zero(m_gaus_par_size) };
            auto gaus_ols{
                GausLinearTransformHelper::BuildGaus3DModel(beta_result.beta_ols, gaus_0)
            };
            auto gaus_mdpde{
                GausLinearTransformHelper::BuildGaus3DModel(beta_result.beta_mdpde, gaus_0)
            };
            residual_matrix_ols_list.at(j).col(i) = CalculateNormalizedResidual(gaus_ols, gaus_true);
            residual_matrix_mdpde_list.at(j).col(i) = CalculateNormalizedResidual(gaus_mdpde, gaus_true);
        }
    }

    for (size_t j = 0; j < alpha_size; j++)
    {
        residual_mean_ols_list.at(j) = residual_matrix_ols_list.at(j).rowwise().mean();
        residual_mean_mdpde_list.at(j) = residual_matrix_mdpde_list.at(j).rowwise().mean();
        residual_sigma_ols_list.at(j) =
            (residual_matrix_ols_list.at(j).colwise() - residual_mean_ols_list.at(j)).rowwise().norm()
            / std::sqrt(m_replica_size - 1);
        residual_sigma_mdpde_list.at(j) =
            (residual_matrix_mdpde_list.at(j).colwise() - residual_mean_mdpde_list.at(j)).rowwise().norm()
            / std::sqrt(m_replica_size - 1);
    }

    return true;
}

bool HRLModelTester::CheckGausParametersDimension(const Eigen::VectorXd & gaus_par)
{
    if (gaus_par.rows() != m_gaus_par_size)
    {
        throw std::invalid_argument("model parameters size invalid, must be : " +
            std::to_string(m_gaus_par_size));
    }
    return true;
}

Eigen::VectorXd HRLModelTester::CalculateNormalizedResidual(
    const Eigen::VectorXd & estimate, const Eigen::VectorXd & truth)
{
    if (estimate.rows() != truth.rows())
    {
        Logger::Log(LogLevel::Error,
            "estimate size " + std::to_string(estimate.rows()) +
            " != model size " + std::to_string(truth.rows()));
        throw std::invalid_argument("model parameters size inconsistant.");
    }
    return ((estimate - truth).array() / truth.array());
}
