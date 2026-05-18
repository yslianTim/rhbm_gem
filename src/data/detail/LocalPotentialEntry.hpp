#pragma once

#include <string>
#include <unordered_map>
#include <utility>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem {

class LocalPotentialEntry
{
    LocalPotentialSampleList m_sampling_entries;
    LocalGaussianResult m_gaussian_result;
    std::unordered_map<std::string, LocalPotentialAnnotation> m_annotation_map;

public:
    LocalPotentialEntry() = default;
    ~LocalPotentialEntry() = default;

    void SetAlphaR(double value) { m_gaussian_result.alpha_r = value; }
    void SetSamplingEntries(LocalPotentialSampleList value)
    {
        m_sampling_entries = std::move(value);
    }
    void SetGaussianResult(LocalGaussianResult value)
    {
        m_gaussian_result = std::move(value);
    }
    void SetAnnotation(const std::string & key, LocalPotentialAnnotation annotation)
    {
        m_annotation_map[key] = std::move(annotation);
    }
    void ClearTransientFitState()
    {
        m_gaussian_result.fit_result.reset();
    }

    int SamplingEntryCount() const
    {
        return static_cast<int>(m_sampling_entries.size());
    }
    const LocalGaussianResult & GaussianResult() const { return m_gaussian_result; }
    LocalPotentialAnnotation * FindAnnotation(const std::string & key)
    {
        const auto iter{ m_annotation_map.find(key) };
        return (iter == m_annotation_map.end()) ? nullptr : &iter->second;
    }
    const LocalPotentialAnnotation * FindAnnotation(const std::string & key) const
    {
        const auto iter{ m_annotation_map.find(key) };
        return (iter == m_annotation_map.end()) ? nullptr : &iter->second;
    }
    const LocalPotentialSampleList & SamplingEntries() const { return m_sampling_entries; }
};

} // namespace rhbm_gem
