#include "LocalPotentialEntry.hpp"
#include "GausLinearTransformHelper.hpp"
#include "Constants.hpp"
#include "Logger.hpp"

#include <cmath>

namespace rhbm_gem {

LocalPotentialEntry::LocalPotentialEntry(void) :
    m_alpha_r{ 0.0 }, m_basis_and_response_entry_list_tmp{}
{

}

LocalPotentialEntry::~LocalPotentialEntry()
{

}

void LocalPotentialEntry::AddDistanceAndMapValueList(std::vector<std::tuple<float, float>> && list)
{
    m_distance_and_map_value_list = std::move(list);
}

void LocalPotentialEntry::AddBasisAndResponseEntryList(std::vector<Eigen::VectorXd> && list)
{
    m_basis_and_response_entry_list_tmp = std::move(list);
}

void LocalPotentialEntry::ClearDistanceAndMapValueList(void)
{
    m_distance_and_map_value_list.clear();
    m_distance_and_map_value_list.shrink_to_fit();
}

void LocalPotentialEntry::AddGausEstimateOLS(double v0, double v1)
{
    m_gaus_estimate_ols = std::make_tuple(v0, v1);
}

void LocalPotentialEntry::AddGausEstimateMDPDE(double v0, double v1)
{
    m_gaus_estimate_mdpde = std::make_tuple(v0, v1);
}

void LocalPotentialEntry::AddGausEstimatePosterior(const std::string & key, double v0, double v1)
{
    m_gaus_estimate_posterior_map[key] = std::make_tuple(v0, v1);
}

void LocalPotentialEntry::AddGausVariancePosterior(const std::string & key, double v0, double v1)
{
    m_gaus_variance_posterior_map[key] = std::make_tuple(v0, v1);
}

void LocalPotentialEntry::AddOutlierTag(const std::string & key, bool value)
{
    m_outlier_tag_map[key] = value;
}

void LocalPotentialEntry::AddStatisticalDistance(const std::string & key, double value)
{
    m_statistical_distance_map[key] = value;
}

int LocalPotentialEntry::GetDistanceAndMapValueListSize(void) const
{
    return static_cast<int>(m_distance_and_map_value_list.size());
}

int LocalPotentialEntry::GetBasisAndResponseEntryListSize(void) const
{
    return static_cast<int>(m_basis_and_response_entry_list_tmp.size());
}

const std::vector<std::tuple<float, float>> & LocalPotentialEntry::GetDistanceAndMapValueList(void) const
{
    return m_distance_and_map_value_list;
}

const std::vector<Eigen::VectorXd> & LocalPotentialEntry::GetBasisAndResponseEntryList(void) const
{
    if (m_basis_and_response_entry_list_tmp.empty())
    {
        throw std::runtime_error("m_basis_and_response_entry_list_tmp is empty");
    }
    return m_basis_and_response_entry_list_tmp;
}

double LocalPotentialEntry::GetMapValueNearCenter(void) const
{
    const auto & data_list{ GetDistanceAndMapValueList() };
    if (data_list.empty()) return 0.0;
    int count{ 0 };
    double sum{ 0.0 };
    for (const auto & [distance, map_value] : data_list)
    {
        if (distance > 0.05f) continue;
        sum += map_value;
        count++;
    }
    if (count == 0) return 0.0;
    return sum/static_cast<double>(count);
}

double LocalPotentialEntry::GetMomentZeroEstimate(void) const
{
    auto data_size{ GetDistanceAndMapValueListSize() };
    if (data_size == 0) return 0.0;
    auto y_sum{ 0.0 };
    auto count{ 0 };
    for (const auto & [distance, map_value] : m_distance_and_map_value_list)
    {
        auto y_value{ static_cast<double>(map_value) };
        if (y_value <= 0.0) continue;
        if (distance > 1.0f) continue;
        y_sum += y_value;
        count++;
    }
    return y_sum/static_cast<double>(count);
}

double LocalPotentialEntry::GetMomentTwoEstimate(void) const
{
    auto data_size{ GetDistanceAndMapValueListSize() };
    if (data_size == 0) return 0.0;
    auto m_0{ GetMomentZeroEstimate() };
    auto y_sum{ 0.0 };
    auto count{ 0 };
    for (const auto & [distance, map_value] : m_distance_and_map_value_list)
    {
        auto y_value{ static_cast<double>(map_value) };
        if (y_value <= 0.0) continue;
        if (distance > 1.0f) continue;
        y_sum += y_value * distance * distance;
        count++;
    }
    y_sum /= static_cast<double>(count);
    return std::sqrt(y_sum/m_0/3.0);
}

double LocalPotentialEntry::GetBetaEstimateOLS(int par_id) const
{
    auto amplitude{ GetAmplitudeEstimateOLS() };
    auto width{ GetWidthEstimateOLS() };
    Eigen::VectorXd coefficient_vector{
        GausLinearTransformHelper::BuildLinearModelCoefficentVector(amplitude, width)
    };
    return coefficient_vector(par_id);
}

double LocalPotentialEntry::GetBetaEstimateMDPDE(int par_id) const
{
    auto amplitude{ GetAmplitudeEstimateMDPDE() };
    auto width{ GetWidthEstimateMDPDE() };
    Eigen::VectorXd coefficient_vector{
        GausLinearTransformHelper::BuildLinearModelCoefficentVector(amplitude, width)
    };
    return coefficient_vector(par_id);
}

double LocalPotentialEntry::GetGausEstimateOLS(int par_id) const
{
    switch (par_id)
    {
        case 0: return GetAmplitudeEstimateOLS();
        case 1: return GetWidthEstimateOLS();
        case 2: return GetIntensityEstimateOLS();
        default:
            Logger::Log(LogLevel::Error, "Invalid parameter index: " + std::to_string(par_id));
            return 0.0;
    }
}

double LocalPotentialEntry::GetAmplitudeEstimateOLS(void) const
{
    return std::get<0>(m_gaus_estimate_ols);
}

double LocalPotentialEntry::GetWidthEstimateOLS(void) const
{
    return std::get<1>(m_gaus_estimate_ols);
}

double LocalPotentialEntry::GetIntensityEstimateOLS(void) const
{
    return CalculateIntensityEstimate(GetAmplitudeEstimateOLS(), GetWidthEstimateOLS());
}

double LocalPotentialEntry::GetGausEstimateMDPDE(int par_id) const
{
    switch (par_id)
    {
        case 0: return GetAmplitudeEstimateMDPDE();
        case 1: return GetWidthEstimateMDPDE();
        case 2: return GetIntensityEstimateMDPDE();
        default:
            Logger::Log(LogLevel::Error, "Invalid parameter index: " + std::to_string(par_id));
            return 0.0;
    }
}

double LocalPotentialEntry::GetAmplitudeEstimateMDPDE(void) const
{
    return std::get<0>(m_gaus_estimate_mdpde);
}

double LocalPotentialEntry::GetWidthEstimateMDPDE(void) const
{
    return std::get<1>(m_gaus_estimate_mdpde);
}

double LocalPotentialEntry::GetIntensityEstimateMDPDE(void) const
{
    return CalculateIntensityEstimate(GetAmplitudeEstimateMDPDE(), GetWidthEstimateMDPDE());
}

double LocalPotentialEntry::GetGausEstimatePosterior(const std::string & key, int par_id) const
{
    switch (par_id)
    {
        case 0: return GetAmplitudeEstimatePosterior(key);
        case 1: return GetWidthEstimatePosterior(key);
        case 2: return GetIntensityEstimatePosterior(key);
        default:
            Logger::Log(LogLevel::Error, "Invalid parameter index: " + std::to_string(par_id));
            return 0.0;
    }
}

double LocalPotentialEntry::GetAmplitudeEstimatePosterior(const std::string & key) const
{
    return std::get<0>(m_gaus_estimate_posterior_map.at(key));
}

double LocalPotentialEntry::GetWidthEstimatePosterior(const std::string & key) const
{
    return std::get<1>(m_gaus_estimate_posterior_map.at(key));
}

double LocalPotentialEntry::GetIntensityEstimatePosterior(const std::string & key) const
{
    return CalculateIntensityEstimate(
        GetAmplitudeEstimatePosterior(key), GetWidthEstimatePosterior(key));
}

double LocalPotentialEntry::GetGausVariancePosterior(const std::string & key, int par_id) const
{
    switch (par_id)
    {
        case 0: return GetAmplitudeVariancePosterior(key);
        case 1: return GetWidthVariancePosterior(key);
        case 2: return GetIntensityVariancePosterior(key);
        default:
            Logger::Log(LogLevel::Error, "Invalid parameter index: " + std::to_string(par_id));
            return 0.0;
    }
}

double LocalPotentialEntry::GetAmplitudeVariancePosterior(const std::string & key) const
{
    return std::get<0>(m_gaus_variance_posterior_map.at(key));
}

double LocalPotentialEntry::GetWidthVariancePosterior(const std::string & key) const
{
    return std::get<1>(m_gaus_variance_posterior_map.at(key));
}

double LocalPotentialEntry::GetIntensityVariancePosterior(const std::string & key) const
{
    return CalculateIntensityVariance(
        GetAmplitudeEstimatePosterior(key), GetAmplitudeVariancePosterior(key),
        GetWidthEstimatePosterior(key), GetWidthVariancePosterior(key));
}

bool LocalPotentialEntry::GetOutlierTag(const std::string & key) const
{
    return m_outlier_tag_map.at(key);
}

double LocalPotentialEntry::GetStatisticalDistance(const std::string & key) const
{
    return m_statistical_distance_map.at(key);
}

double LocalPotentialEntry::CalculateIntensityEstimate(double amplitude, double width) const
{
    if (width == 0.0) return 0.0;
    return amplitude * std::pow(Constants::two_pi * width * width, -1.5);
}

double LocalPotentialEntry::CalculateIntensityVariance(
    double amplitude, double sigma_amplitude, double width, double sigma_width) const
{
    return std::sqrt(
        std::pow(std::pow(Constants::two_pi * width * width, -1.5) * sigma_amplitude, 2) +
        std::pow(-3.0 * amplitude * std::pow(Constants::two_pi, -1.5) * std::pow(width, -4) * sigma_width, 2)
    );
}

} // namespace rhbm_gem
