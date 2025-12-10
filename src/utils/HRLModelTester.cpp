#include "HRLModelTester.hpp"
#include "HRLModelHelper.hpp"
#include "GausLinearTransformHelper.hpp"
#include "Constants.hpp"
#include "Logger.hpp"

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
    m_member_size{ 1 },
    m_replica_size{ replica_size },
    m_sampling_entry_size{ 1000 },
    m_x_min{ 0.0 }, m_x_max{ 1.0 },
    m_data_error_sigma{ 0.05 },
    m_outlier_ratio{ 0.0 },
    m_gaus_par_prior{ Eigen::VectorXd::Zero(m_gaus_par_size) },
    m_gaus_par_sigma{ Eigen::VectorXd::Ones(m_gaus_par_size) }
{
}

void HRLModelTester::SetGausParametersPrior(const Eigen::VectorXd & gaus_par_prior)
{
    CheckGausParametersDimension(gaus_par_prior);
    m_gaus_par_prior = gaus_par_prior;
}

void HRLModelTester::SetGausParametersSigma(const Eigen::VectorXd & gaus_par_sigma)
{
    CheckGausParametersDimension(gaus_par_sigma);
    m_gaus_par_sigma = gaus_par_sigma;
}

Eigen::MatrixXd HRLModelTester::BuildRandomGausParameters(int replica_size)
{
    std::vector<std::normal_distribution<>> dist_gaus_par_list;
    for (int p = 0; p < m_gaus_par_size; p++)
    {
        std::normal_distribution<> dist_gaus_par(m_gaus_par_prior(p), m_gaus_par_sigma(p));
        dist_gaus_par_list.emplace_back(dist_gaus_par);
    }

    Eigen::MatrixXd gaus_par_matrix{ Eigen::MatrixXd::Zero(m_gaus_par_size, replica_size) };
    for (int i = 0; i < replica_size; i++)
    {
        Eigen::VectorXd gaus_par{ Eigen::VectorXd::Zero(m_gaus_par_size) };
        for (int p = 0; p < m_gaus_par_size; p++)
        {
            gaus_par(p) = dist_gaus_par_list.at(static_cast<size_t>(p))(m_generator);
        }
        gaus_par_matrix.col(i) = gaus_par;
    }
    return gaus_par_matrix;
}

std::vector<std::tuple<float, float>> HRLModelTester::BuildRandomGausSamplingEntry(
    size_t sampling_entry_size, const Eigen::VectorXd & gaus_par, double outlier_ratio)
{
    CheckGausParametersDimension(gaus_par);
    auto amplitude{ gaus_par(0) };
    auto width{ gaus_par(1) };
    //auto intersect{ gaus_par(2) };
    auto intersect{ 0.0 };
    std::uniform_real_distribution<> dist_distance(m_x_min, m_x_max);
    std::uniform_real_distribution<> dist_outlier(0.0, 1.0);
    std::vector<std::tuple<float, float>> sampling_entry_list;
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
        sampling_entry_list.emplace_back(r, y);
    }
    return sampling_entry_list;
}

std::vector<Eigen::VectorXd> HRLModelTester::BuildRandomLinearDataEntry(
    size_t sampling_entry_size,
    const Eigen::VectorXd & gaus_par,
    double error_sigma,
    double outlier_ratio)
{
    auto sampling_entries{ BuildRandomGausSamplingEntry(sampling_entry_size, gaus_par, outlier_ratio) };
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
    double alpha_r,
    Eigen::VectorXd & residual_mean_ols,
    Eigen::VectorXd & residual_mean_mdpde,
    Eigen::VectorXd & residual_sigma_ols,
    Eigen::VectorXd & residual_sigma_mdpde,
    int thread_size)
{
    auto gaus_par_local_matrix{ BuildRandomGausParameters(m_replica_size) };
    Eigen::MatrixXd gaus_residual_matrix_ols(m_gaus_par_size, m_replica_size);
    Eigen::MatrixXd gaus_residual_matrix_mdpde(m_gaus_par_size, m_replica_size);

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(thread_size)
#endif
    for (int i = 0; i < m_replica_size; i++)
    {
        Eigen::VectorXd gaus_true{ gaus_par_local_matrix.col(i) };
        auto data_entry_list{
            BuildRandomLinearDataEntry(
                m_sampling_entry_size, gaus_true, m_data_error_sigma, m_outlier_ratio
            )
        };
        auto data_array{ HRLModelHelper::BuildBasisVectorAndResponseArray(data_entry_list) };
        const auto & X{ std::get<0>(data_array) };
        const auto & y{ std::get<1>(data_array) };

        Eigen::VectorXd beta_ols;
        Eigen::VectorXd beta_mdpde;
        double sigma_square;
        Eigen::DiagonalMatrix<double, Eigen::Dynamic> W;
        Eigen::DiagonalMatrix<double, Eigen::Dynamic> capital_sigma;
        HRLModelHelper::AlgorithmBetaMDPDE(
            alpha_r, X, y,
            beta_ols, beta_mdpde, sigma_square, W, capital_sigma, true
        );

        Eigen::VectorXd gaus_0{ Eigen::VectorXd::Zero(m_gaus_par_size) };
        auto gaus_ols{ GausLinearTransformHelper::BuildGaus3DModel(beta_ols, gaus_0) };
        auto gaus_mdpde{ GausLinearTransformHelper::BuildGaus3DModel(beta_mdpde, gaus_0) };
        gaus_residual_matrix_ols.col(i) = CalculateResidual(gaus_ols, gaus_true);
        gaus_residual_matrix_mdpde.col(i) = CalculateResidual(gaus_mdpde, gaus_true);
    }
    residual_mean_ols = gaus_residual_matrix_ols.rowwise().mean();
    residual_mean_mdpde = gaus_residual_matrix_mdpde.rowwise().mean();
    residual_sigma_ols = (gaus_residual_matrix_ols.colwise() - residual_mean_ols).rowwise().norm() / std::sqrt(m_replica_size - 1);
    residual_sigma_mdpde = (gaus_residual_matrix_mdpde.colwise() - residual_mean_mdpde).rowwise().norm() / std::sqrt(m_replica_size - 1);
    
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

Eigen::VectorXd HRLModelTester::CalculateResidual(
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
