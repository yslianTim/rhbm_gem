#include "data/detail/LocalPotentialEntry.hpp"
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/utils/hrl/GaussianLinearizationService.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <cmath>
#include <stdexcept>
#include <utility>

namespace rhbm_gem {
namespace {

inline bool IsEffectiveSample(const LocalPotentialSample & sample)
{
    return sample.weight > 0.0f;
}

const GaussianLinearizationService & LinearizationService()
{
    static const GaussianLinearizationService service{
        GaussianLinearizationSpec::DefaultDataset()
    };
    return service;
}

} // namespace

LocalPotentialEntry::LocalPotentialEntry() :
    m_alpha_r{ 0.0 }
{

}

LocalPotentialEntry::~LocalPotentialEntry()
{

}

void LocalPotentialEntry::SetSamplingEntries(LocalPotentialSampleList value)
{
    m_sampling_entries = std::move(value);
}

void LocalPotentialEntry::SetDataset(
    HRLMemberDataset dataset)
{
    m_dataset = std::move(dataset);
}

void LocalPotentialEntry::SetFitResult(HRLBetaEstimateResult value)
{
    m_fit_result = std::move(value);
}

void LocalPotentialEntry::SetAnnotation(
    const std::string & key,
    LocalPotentialAnnotation annotation)
{
    m_annotation_map[key] = std::move(annotation);
}

void LocalPotentialEntry::ClearTransientFitState()
{
    m_dataset.reset();
    m_fit_result.reset();
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

int LocalPotentialEntry::GetSamplingEntryCount() const
{
    return static_cast<int>(m_sampling_entries.size());
}

const std::unordered_map<std::string, LocalPotentialAnnotation> &
LocalPotentialEntry::Annotations() const
{
    return m_annotation_map;
}

const HRLMemberDataset & LocalPotentialEntry::GetDataset() const
{
    if (!m_dataset.has_value())
    {
        throw std::runtime_error("LocalPotentialEntry dataset is not available");
    }
    return *m_dataset;
}

const HRLBetaEstimateResult & LocalPotentialEntry::GetFitResult() const
{
    if (!m_fit_result.has_value())
    {
        throw std::runtime_error("LocalPotentialEntry fit result is not available");
    }
    return *m_fit_result;
}

const LocalPotentialSampleList & LocalPotentialEntry::GetSamplingEntries() const
{
    return m_sampling_entries;
}

std::tuple<float, float> LocalPotentialEntry::GetDistanceRange(double margin_rate) const
{
    std::vector<float> distance_array;
    distance_array.reserve(m_sampling_entries.size());
    for (const auto & sample : m_sampling_entries)
    {
        if (!IsEffectiveSample(sample)) continue;
        distance_array.emplace_back(sample.distance);
    }
    return ArrayStats<float>::ComputeScalingRangeTuple(
        distance_array, static_cast<float>(margin_rate));
}

std::tuple<float, float> LocalPotentialEntry::GetResponseRange(double margin_rate) const
{
    std::vector<float> map_value_array;
    map_value_array.reserve(m_sampling_entries.size());
    for (const auto & sample : m_sampling_entries)
    {
        if (!IsEffectiveSample(sample)) continue;
        map_value_array.emplace_back(sample.response);
    }
    return ArrayStats<float>::ComputeScalingRangeTuple(
        map_value_array, static_cast<float>(margin_rate));
}

SeriesPointList LocalPotentialEntry::GetBinnedDistanceResponseSeries(
    int bin_size,
    double x_min,
    double x_max) const
{
    const auto bin_spacing{ (x_max - x_min) / static_cast<double>(bin_size) };
    std::vector<std::vector<float>> bin_map(static_cast<size_t>(bin_size));
    std::vector<std::vector<float>> bin_weight_map(static_cast<size_t>(bin_size));
    for (const auto & sample : m_sampling_entries)
    {
        if (!IsEffectiveSample(sample)) continue;
        const auto distance{ sample.distance };
        if (distance < x_min || distance >= x_max)
        {
            continue;
        }
        const auto shifted_distance{ static_cast<double>(distance) - x_min };
        const auto bin_index{ static_cast<int>(std::floor(shifted_distance / bin_spacing)) };
        if (bin_index < 0 || bin_index >= bin_size)
        {
            continue;
        }
        bin_map.at(static_cast<size_t>(bin_index)).emplace_back(sample.response);
        bin_weight_map.at(static_cast<size_t>(bin_index)).emplace_back(sample.weight);
    }

    SeriesPointList binned_distance_response_series;
    binned_distance_response_series.reserve(static_cast<size_t>(bin_size));
    for (int i = 0; i < bin_size; i++)
    {
        const auto x_value{ static_cast<float>(x_min + (i + 0.5) * bin_spacing) };
        const auto & bin_values{ bin_map.at(static_cast<size_t>(i)) };
        const auto & bin_weights{ bin_weight_map.at(static_cast<size_t>(i)) };
        const auto y_value{ bin_values.empty() ? 0.0f : ArrayStats<float>::ComputeMedian(bin_values) };
        const auto weight_value{
            bin_weights.empty()
                ? 0.0f
                : ArrayStats<float>::ComputeMean(bin_weights.data(), bin_weights.size())
        };
        binned_distance_response_series.emplace_back(SeriesPoint{ { x_value }, y_value, weight_value });
    }
    return binned_distance_response_series;
}

SeriesPointList LocalPotentialEntry::GetLinearModelSeries() const
{
    return LinearizationService().BuildLinearModelSeries(m_sampling_entries);
}

double LocalPotentialEntry::GetMapValueNearCenter() const
{
    const auto & data_list{ GetSamplingEntries() };
    if (data_list.empty()) return 0.0;
    int count{ 0 };
    double sum{ 0.0 };
    for (const auto & sample : data_list)
    {
        if (!IsEffectiveSample(sample)) continue;
        if (sample.distance > 0.05f) continue;
        sum += sample.response;
        count++;
    }
    if (count == 0) return 0.0;
    return sum/static_cast<double>(count);
}

double LocalPotentialEntry::GetMomentZeroEstimate() const
{
    auto data_size{ GetSamplingEntryCount() };
    if (data_size == 0) return 0.0;
    auto y_sum{ 0.0 };
    auto count{ 0 };
    for (const auto & sample : m_sampling_entries)
    {
        if (!IsEffectiveSample(sample)) continue;
        auto y_value{ static_cast<double>(sample.response) };
        if (y_value <= 0.0) continue;
        if (sample.distance > 1.0f) continue;
        y_sum += y_value;
        count++;
    }
    if (count == 0) return 0.0;
    return y_sum/static_cast<double>(count);
}

double LocalPotentialEntry::GetMomentTwoEstimate() const
{
    auto data_size{ GetSamplingEntryCount() };
    if (data_size == 0) return 0.0;
    auto m_0{ GetMomentZeroEstimate() };
    auto y_sum{ 0.0 };
    auto count{ 0 };
    for (const auto & sample : m_sampling_entries)
    {
        if (!IsEffectiveSample(sample)) continue;
        auto y_value{ static_cast<double>(sample.response) };
        if (y_value <= 0.0) continue;
        if (sample.distance > 1.0f) continue;
        y_sum += y_value * sample.distance * sample.distance;
        count++;
    }
    if (count == 0 || m_0 == 0.0) return 0.0;
    y_sum /= static_cast<double>(count);
    return std::sqrt(y_sum/m_0/3.0);
}

double LocalPotentialEntry::CalculateQScore(int par_choice) const
{
    if (m_sampling_entries.empty())
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
    distance_list.reserve(m_sampling_entries.size());
    map_value_list.reserve(m_sampling_entries.size());
    reference_value_list.reserve(m_sampling_entries.size());
    for (const auto & sample : m_sampling_entries)
    {
        if (!IsEffectiveSample(sample)) continue;
        distance_list.emplace_back(sample.distance);
        map_value_list.emplace_back(sample.response);
        reference_value_list.emplace_back(
            amplitude * std::exp(-0.5 * std::pow(sample.distance / width, 2)) + intersect);
    }
    if (map_value_list.empty())
    {
        return 0.0;
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
