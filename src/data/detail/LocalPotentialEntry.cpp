#include "data/detail/LocalPotentialEntry.hpp"
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <cmath>
#include <utility>

namespace rhbm_gem {

LocalPotentialEntry::LocalPotentialEntry() :
    m_alpha_r{ 0.0 }
{

}

LocalPotentialEntry::~LocalPotentialEntry()
{

}

void LocalPotentialEntry::SetDistanceAndMapValueList(std::vector<std::tuple<float, float>> value)
{
    m_distance_and_map_value_list = std::move(value);
}

void LocalPotentialEntry::SetAnnotation(
    const std::string & key,
    LocalPotentialAnnotation annotation)
{
    m_annotation_map[key] = std::move(annotation);
}

LocalPotentialAnnotation * LocalPotentialEntry::FindAnnotation(const std::string & key)
{
    const auto iter{ m_annotation_map.find(key) };
    return (iter == m_annotation_map.end()) ? nullptr : &iter->second;
}

const LocalPotentialAnnotation * LocalPotentialEntry::FindAnnotation(const std::string & key) const
{
    const auto iter{ m_annotation_map.find(key) };
    return (iter == m_annotation_map.end()) ? nullptr : &iter->second;
}

int LocalPotentialEntry::GetDistanceAndMapValueListSize() const
{
    return static_cast<int>(m_distance_and_map_value_list.size());
}

const std::unordered_map<std::string, LocalPotentialAnnotation> &
LocalPotentialEntry::Annotations() const
{
    return m_annotation_map;
}

const std::vector<std::tuple<float, float>> & LocalPotentialEntry::GetDistanceAndMapValueList() const
{
    return m_distance_and_map_value_list;
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
        const auto & estimate{ GetEstimateMDPDE() };
        amplitude = estimate.Intensity();
        width = estimate.width;
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
