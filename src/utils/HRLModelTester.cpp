#include "HRLModelTester.hpp"
#include "HRLModelHelper.hpp"
#include "GausLinearTransformHelper.hpp"
#include "Constants.hpp"
#include "Logger.hpp"

#include <cmath>
#include <random>
#include <limits>
#include <utility>
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

HRLModelTester::HRLModelTester(int basis_size, int member_size) :
    m_basis_size{ ValidatePositive(basis_size, "basis_size") },
    m_member_size{ ValidatePositive(member_size, "member_size") },
    m_x_min{ 0.0 }, m_x_max{ 1.0 },
    m_data_error_sigma{ 0.05 },
    m_model_par_prior{ Eigen::VectorXd::Zero(3) },
    m_model_par_sigma{ Eigen::VectorXd::Ones(3) },
    m_data_array{}
{
    m_data_array.reserve(static_cast<size_t>(m_member_size));
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
    std::normal_distribution<> dist_amplitude(m_model_par_prior(0), m_model_par_sigma(0));
    std::normal_distribution<> dist_width(m_model_par_prior(1), m_model_par_sigma(1));
    std::normal_distribution<> dist_intersect(m_model_par_prior(2), m_model_par_sigma(2));

    m_data_array.clear();
    const size_t sampling_entry_size{ 100 };
    for (int i = 0; i < m_member_size; i++)
    {
        Eigen::VectorXd model_par{ Eigen::VectorXd::Zero(3) };
        model_par(0) = dist_amplitude(generator);
        model_par(1) = dist_width(generator);
        model_par(2) = dist_intersect(generator);

        Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
        auto sampling_entries{ BuildRandomSamplingEntry(sampling_entry_size, model_par) };
        std::vector<Eigen::VectorXd> data_vector_list;
        data_vector_list.reserve(sampling_entry_size);
        for (auto & [r, y] : sampling_entries)
        {
            auto data_vector{
                GausLinearTransformHelper::BuildLinearModelDataVector(r, y, model_par_init, m_basis_size)
            };
            data_vector_list.emplace_back(data_vector);
        }
        m_data_array.emplace_back(data_vector_list, "member_"+std::to_string(i));
    }
}

std::vector<std::tuple<double, double>> HRLModelTester::BuildRandomSamplingEntry(
    size_t sampling_entry_size, const Eigen::VectorXd & model_par)
{
    CheckModelParametersDimension(model_par);
    auto gaus_amplitude{ model_par(0) };
    auto gaus_width{ model_par(1) };
    auto intersect{ model_par(2) };
    std::random_device random_device;
    std::mt19937 generator{ random_device() };
    std::uniform_real_distribution<> dist_distance(m_x_min, m_x_max);
    auto data_max_value{ gaus_amplitude * std::pow(Constants::two_pi * std::pow(gaus_width, 2), -1.5) };
    std::normal_distribution<> dist_error(0.0, m_data_error_sigma * data_max_value);
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
        member_sampling_entry.emplace_back(r, y);
    }
    return member_sampling_entry;
}

bool HRLModelTester::RunTest(double alpha_r, double alpha_g)
{
    BuildDataArray();

    auto model_estimator{ std::make_unique<HRLModelHelper>(m_basis_size, m_member_size) };
    model_estimator->SetThreadSize(1);
    model_estimator->SetDataArray(std::move(m_data_array));
    model_estimator->RunEstimation(alpha_r, alpha_g);
    return true;
}

bool HRLModelTester::CheckModelParametersDimension(const Eigen::VectorXd & model_par)
{
    if (model_par.rows() != 3)
    {
        throw std::invalid_argument("model parameters size invalid, must be 3");
    }
    return true;
}
