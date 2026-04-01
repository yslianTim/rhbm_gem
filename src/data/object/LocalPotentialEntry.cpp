#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <cmath>

namespace rhbm_gem {

LocalPotentialEntry::LocalPotentialEntry() :
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

void LocalPotentialEntry::ClearDistanceAndMapValueList()
{
    m_distance_and_map_value_list.clear();
    m_distance_and_map_value_list.shrink_to_fit();
}

void LocalPotentialEntry::SetPosterior(
    const std::string & key,
    const GaussianPosterior & posterior)
{
    m_gaus_posterior_map[key] = posterior;
}

void LocalPotentialEntry::AddOutlierTag(const std::string & key, bool value)
{
    m_outlier_tag_map[key] = value;
}

void LocalPotentialEntry::AddStatisticalDistance(const std::string & key, double value)
{
    m_statistical_distance_map[key] = value;
}

int LocalPotentialEntry::GetDistanceAndMapValueListSize() const
{
    return static_cast<int>(m_distance_and_map_value_list.size());
}

int LocalPotentialEntry::GetBasisAndResponseEntryListSize() const
{
    return static_cast<int>(m_basis_and_response_entry_list_tmp.size());
}

const GaussianPosterior & LocalPotentialEntry::GetPosterior(const std::string & key) const
{
    return m_gaus_posterior_map.at(key);
}

const std::vector<std::tuple<float, float>> & LocalPotentialEntry::GetDistanceAndMapValueList() const
{
    return m_distance_and_map_value_list;
}

const std::vector<Eigen::VectorXd> & LocalPotentialEntry::GetBasisAndResponseEntryList() const
{
    if (m_basis_and_response_entry_list_tmp.empty())
    {
        throw std::runtime_error("m_basis_and_response_entry_list_tmp is empty");
    }
    return m_basis_and_response_entry_list_tmp;
}

double LocalPotentialEntry::GetMapValueNearCenter() const
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

double LocalPotentialEntry::GetMomentZeroEstimate() const
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

double LocalPotentialEntry::GetMomentTwoEstimate() const
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
    return GetEstimateOLS().ToBeta()(par_id);
}

double LocalPotentialEntry::GetBetaEstimateMDPDE(int par_id) const
{
    return GetEstimateMDPDE().ToBeta()(par_id);
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

double LocalPotentialEntry::GetAmplitudeEstimateOLS() const
{
    return m_gaus_estimate_ols.amplitude;
}

double LocalPotentialEntry::GetWidthEstimateOLS() const
{
    return m_gaus_estimate_ols.width;
}

double LocalPotentialEntry::GetIntensityEstimateOLS() const
{
    return m_gaus_estimate_ols.Intensity();
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

double LocalPotentialEntry::GetAmplitudeEstimateMDPDE() const
{
    return m_gaus_estimate_mdpde.amplitude;
}

double LocalPotentialEntry::GetWidthEstimateMDPDE() const
{
    return m_gaus_estimate_mdpde.width;
}

double LocalPotentialEntry::GetIntensityEstimateMDPDE() const
{
    return m_gaus_estimate_mdpde.Intensity();
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
    return GetPosterior(key).estimate.amplitude;
}

double LocalPotentialEntry::GetWidthEstimatePosterior(const std::string & key) const
{
    return GetPosterior(key).estimate.width;
}

double LocalPotentialEntry::GetIntensityEstimatePosterior(const std::string & key) const
{
    return GetPosterior(key).estimate.Intensity();
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
    return GetPosterior(key).variance.amplitude;
}

double LocalPotentialEntry::GetWidthVariancePosterior(const std::string & key) const
{
    return GetPosterior(key).variance.width;
}

double LocalPotentialEntry::GetIntensityVariancePosterior(const std::string & key) const
{
    return GetPosterior(key).IntensityVariance();
}

bool LocalPotentialEntry::GetOutlierTag(const std::string & key) const
{
    return m_outlier_tag_map.at(key);
}

double LocalPotentialEntry::GetStatisticalDistance(const std::string & key) const
{
    return m_statistical_distance_map.at(key);
}

double LocalPotentialEntry::CalculateQScore(int par_choice) const
{
    if (m_distance_and_map_value_list.empty())
    {
        return 0.0;
    }

    auto q_score{ 0.0 };
    auto amplitude{ 0.0 };
    auto width{ 0.0 };
    auto intersect{ 0.0 };
    if (par_choice == 0)
    {
        amplitude = 0.05; // TODO
        width = 0.6; // TODO
        intersect = -0.005; // TODO
    }
    else if (par_choice == 1)
    {
        amplitude = GetAmplitudeEstimateMDPDE() * std::pow(2.0*Constants::pi*GetWidthEstimateMDPDE()*GetWidthEstimateMDPDE(), -1.5);
        width = GetWidthEstimateMDPDE();
        intersect = 0.0;
    }

    if (std::isfinite(amplitude) == false || std::isfinite(width) == false ||
        std::isfinite(intersect) == false || width <= 0.0)
    {
        return 0.0;
    }

    std::vector<float> distance_list;
    std::vector<float> map_value_list;
    std::vector<float> reference_value_list;
    distance_list.reserve(m_distance_and_map_value_list.size());
    map_value_list.reserve(m_distance_and_map_value_list.size());
    reference_value_list.reserve(m_distance_and_map_value_list.size());
    for (auto & [distance, map_value] : m_distance_and_map_value_list)
    {
        distance_list.emplace_back(distance);
        map_value_list.emplace_back(map_value);
        reference_value_list.emplace_back(amplitude * std::exp(-0.5 * std::pow(distance/width, 2)) + intersect);
    }

    float map_value_mean{ ArrayStats<float>::ComputeMean(map_value_list.data(), map_value_list.size()) };
    float reference_value_mean{ ArrayStats<float>::ComputeMean(reference_value_list.data(), reference_value_list.size()) };
    
    auto numerator{ 0.0 };
    auto denominator{ 0.0 };
    auto map_value_norm_square{ 0.0 };
    auto reference_value_norm_square{ 0.0 };
    for (size_t i = 0; i < map_value_list.size(); i++)
    {
        auto map_value_diff{ map_value_list.at(i) - map_value_mean };
        auto reference_value_diff{ reference_value_list.at(i) - reference_value_mean };
        numerator += map_value_diff * reference_value_diff;
        map_value_norm_square += map_value_diff * map_value_diff;
        reference_value_norm_square += reference_value_diff * reference_value_diff;
    }
    denominator = std::sqrt(map_value_norm_square) * std::sqrt(reference_value_norm_square);
    if (denominator == 0.0 || std::isfinite(denominator) == false || std::isfinite(numerator) == false)
    {
        return 0.0;
    }
    q_score = numerator/denominator;

    if (std::isfinite(q_score) == false)
    {
        return 0.0;
    }

    return q_score;
}

} // namespace rhbm_gem
