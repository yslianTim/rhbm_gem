#include "HRLModelTester.hpp"
#include "HRLModelHelper.hpp"
#include "GausLinearTransformHelper.hpp"
#include "Constants.hpp"
#include "Logger.hpp"

#include <cmath>
#include <random>
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

HRLModelTester::HRLModelTester(int model_par_size, int linear_basis_size, int member_size) :
    m_model_par_size{ ValidatePositive(model_par_size, "model_par_size") },
    m_linear_basis_size{ ValidatePositive(linear_basis_size, "linear_basis_size") },
    m_member_size{ ValidatePositive(member_size, "member_size") },
    m_sampling_entry_size{ 1000 },
    m_x_min{ 0.0 }, m_x_max{ 1.0 },
    m_data_error_sigma{ 0.05 },
    m_model_par_prior{ Eigen::VectorXd::Zero(m_model_par_size) },
    m_model_par_sigma{ Eigen::VectorXd::Ones(m_model_par_size) },
    m_sampling_entries_list{},
    m_data_array{},
    m_local_estimate_mean_ols{ Eigen::VectorXd::Zero(m_model_par_size) },
    m_local_estimate_mean_mdpde{ Eigen::VectorXd::Zero(m_model_par_size) },
    m_local_estimate_mean_posterior{ Eigen::VectorXd::Zero(m_model_par_size) },
    m_group_estimate_residual_mean{ Eigen::VectorXd::Zero(m_model_par_size) },
    m_group_estimate_residual_mdpde{ Eigen::VectorXd::Zero(m_model_par_size) },
    m_group_estimate_residual_hrl{ Eigen::VectorXd::Zero(m_model_par_size) }
{
    auto size{ static_cast<size_t>(m_member_size) };
    m_sampling_entries_list.reserve(size);
    m_data_array.reserve(size);
    m_model_par_local_list.reserve(size);
    m_model_par_0_list.reserve(size);
    m_local_estimate_residual_ols_list.reserve(size);
    m_local_estimate_residual_mdpde_list.reserve(size);
    m_local_estimate_residual_posterior_list.reserve(size);
    m_local_estimate_deviation_ols_list.reserve(size);
    m_local_estimate_deviation_mdpde_list.reserve(size);
    m_local_estimate_deviation_posterior_list.reserve(size);
}

void HRLModelTester::SetModelParametersPrior(const Eigen::VectorXd & model_par_prior)
{
    CheckModelParametersDimension(model_par_prior);
    m_model_par_prior = model_par_prior;
}

void HRLModelTester::SetModelParametersSigma(const Eigen::VectorXd & model_par_sigma)
{
    CheckModelParametersDimension(model_par_sigma);
    m_model_par_sigma = model_par_sigma;
}

void HRLModelTester::BuildDataArray(void)
{
    std::random_device random_device;
    std::mt19937 generator{ random_device() };
    std::vector<std::normal_distribution<>> dist_model_par_list;
    for (int p = 0; p < m_model_par_size; p++)
    {
        std::normal_distribution<> dist_model_par(m_model_par_prior(p), m_model_par_sigma(p));
        dist_model_par_list.emplace_back(dist_model_par);
    }

    m_data_array.clear();
    for (int i = 0; i < m_member_size; i++)
    {
        Eigen::VectorXd model_par_local{ Eigen::VectorXd::Zero(m_model_par_size) };
        for (int p = 0; p < m_model_par_size; p++)
        {
            model_par_local(p) = dist_model_par_list.at(static_cast<size_t>(p))(generator);
        }

        auto sampling_entries{ BuildRandomSamplingEntry(m_sampling_entry_size, model_par_local) };
        Eigen::VectorXd model_par_0{ CalculateMoment(sampling_entries) };
        m_sampling_entries_list.emplace_back(sampling_entries);
        std::vector<Eigen::VectorXd> data_vector_list;
        data_vector_list.reserve(sampling_entries.size());
        for (auto & [r, y] : sampling_entries)
        {
            auto data_vector{
                GausLinearTransformHelper::BuildLinearModelDataVector(
                    r, y, model_par_0, m_linear_basis_size)
            };
            data_vector_list.emplace_back(data_vector);
        }
        m_data_array.emplace_back(data_vector_list, "member_"+std::to_string(i));
        m_model_par_local_list.emplace_back(model_par_local);
        m_model_par_0_list.emplace_back(model_par_0);
    }
}

std::vector<std::tuple<double, double>> HRLModelTester::BuildRandomSamplingEntry(
    size_t sampling_entry_size, const Eigen::VectorXd & model_par)
{
    CheckModelParametersDimension(model_par);
    auto gaus_amplitude{ model_par(0) };
    auto gaus_width{ model_par(1) };
    //auto intersect{ model_par(2) };
    auto intersect{ 0.0 };
    std::random_device random_device;
    std::mt19937 generator{ random_device() };
    std::uniform_real_distribution<> dist_distance(m_x_min, m_x_max);
    std::normal_distribution<> dist_error(0.0, m_data_error_sigma);
    std::vector<std::tuple<double, double>> member_sampling_entry;
    member_sampling_entry.reserve(sampling_entry_size);
    for (size_t i = 0; i < sampling_entry_size; i++)
    {
        auto r{ dist_distance(generator) };
        auto error{ dist_error(generator) };
        auto y{
            gaus_amplitude * GausLinearTransformHelper::GetGaussianResponseAtDistance(r, gaus_width)
            + intersect + error
        };
        if (y <= 0.0) continue;
        member_sampling_entry.emplace_back(r, y);
    }
    return member_sampling_entry;
}

bool HRLModelTester::RunTest(double alpha_r, double alpha_g)
{
    BuildDataArray();

    auto estimator{ std::make_unique<HRLModelHelper>(m_linear_basis_size, m_member_size) };
    estimator->SetQuietMode();
    estimator->SetThreadSize(1);
    estimator->SetDataArray(std::move(m_data_array));
    estimator->RunEstimation(alpha_g);

    auto group_mean{
        GausLinearTransformHelper::BuildGaus3DModel(estimator->GetMuVectorMean(), m_model_par_prior)
    };
    auto group_mdpde{
        GausLinearTransformHelper::BuildGaus3DModel(estimator->GetMuVectorMDPDE(), m_model_par_prior)
    };
    auto group_hrl{
        GausLinearTransformHelper::BuildGaus3DModel(estimator->GetMuVectorPrior(), m_model_par_prior)
    };
    m_group_estimate_residual_mean = CalculateResidual(group_mean, m_model_par_prior);
    m_group_estimate_residual_mdpde = CalculateResidual(group_mdpde, m_model_par_prior);
    m_group_estimate_residual_hrl = CalculateResidual(group_hrl, m_model_par_prior);

    for (int i = 0; i < m_member_size; i++)
    {
        auto model_par_local{ m_model_par_local_list.at(static_cast<size_t>(i)) };
        auto model_par_0{ m_model_par_0_list.at(static_cast<size_t>(i)) };
        auto local_ols{
            GausLinearTransformHelper::BuildGaus3DModel(
                estimator->GetBetaMatrixMDPDE(i), model_par_0) // to be update
        };
        auto local_mdpde{
            GausLinearTransformHelper::BuildGaus3DModel(
                estimator->GetBetaMatrixMDPDE(i), model_par_0)
        };
        auto local_posterior{
            GausLinearTransformHelper::BuildGaus3DModel(
                estimator->GetBetaMatrixPosterior(i), model_par_0)
        };
        m_local_estimate_residual_ols_list.emplace_back(
            CalculateResidual(local_ols, model_par_local));
        m_local_estimate_residual_mdpde_list.emplace_back(
            CalculateResidual(local_mdpde, model_par_local));
        m_local_estimate_residual_posterior_list.emplace_back(
            CalculateResidual(local_posterior, model_par_local));
        
        m_local_estimate_deviation_ols_list.emplace_back(local_ols);
        m_local_estimate_deviation_mdpde_list.emplace_back(local_mdpde);
        m_local_estimate_deviation_posterior_list.emplace_back(local_posterior);

        m_local_estimate_mean_ols += local_ols;
        m_local_estimate_mean_mdpde += local_mdpde;
        m_local_estimate_mean_posterior += local_posterior;
    }
    m_local_estimate_mean_ols /= m_member_size;
    m_local_estimate_mean_mdpde /= m_member_size;
    m_local_estimate_mean_posterior /= m_member_size;

    for (size_t i = 0; i < static_cast<size_t>(m_member_size); i++)
    {
        m_local_estimate_deviation_ols_list.at(i) -= m_local_estimate_mean_ols;
        m_local_estimate_deviation_mdpde_list.at(i) -= m_local_estimate_mean_mdpde;
        m_local_estimate_deviation_posterior_list.at(i) -= m_local_estimate_mean_posterior;
    }

    return true;
}

bool HRLModelTester::CheckModelParametersDimension(const Eigen::VectorXd & model_par)
{
    if (model_par.rows() != m_model_par_size)
    {
        throw std::invalid_argument("model parameters size invalid, must be : " +
            std::to_string(m_model_par_size));
    }
    return true;
}

Eigen::VectorXd HRLModelTester::CalculateResidual(
    const Eigen::VectorXd & estimate, const Eigen::VectorXd & true_value)
{
    if (estimate.rows() != true_value.rows())
    {
        Logger::Log(LogLevel::Error,
            "estimate size " + std::to_string(estimate.rows()) +
            " != model size " + std::to_string(true_value.rows()));
        throw std::invalid_argument("model parameters size inconsistant.");
    }
    return ((estimate - true_value).array() / true_value.array());
}

Eigen::VectorXd HRLModelTester::CalculateMoment(
    const std::vector<std::tuple<double, double>> & sampling_entries) const
{
    Eigen::VectorXd model_par{ Eigen::VectorXd::Zero(3) };
    double m_0{ 0.0 }, m_2{ 0.0 };
    double y_max{ 0.0 };
    for (auto & [r, y] : sampling_entries)
    {
        auto y_new{ y - model_par(2) };
        m_0 += y_new;
        m_2 += y_new * r * r;
        y_max = std::max(y_max, y_new);
    }
    model_par(1) = std::sqrt(m_2 / m_0 / 3.0);
    model_par(0) = y_max / GausLinearTransformHelper::GetGaussianResponseAtDistance(0.0, model_par(1));
    return model_par;
}
