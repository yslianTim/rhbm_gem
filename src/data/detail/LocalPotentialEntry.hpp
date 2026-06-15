#pragma once

#include <optional>
#include <utility>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem {

class LocalPotentialEntry
{
    LocalPotentialSampleList m_sampling_entries;
    LocalGaussianResult m_gaussian_result;
    std::optional<LocalPotentialAnnotation> m_annotation;

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
    void SetAnnotation(LocalPotentialAnnotation annotation)
    {
        m_annotation = std::move(annotation);
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
    LocalPotentialAnnotation * FindAnnotation()
    {
        return m_annotation.has_value() ? &m_annotation.value() : nullptr;
    }
    const LocalPotentialAnnotation * FindAnnotation() const
    {
        return m_annotation.has_value() ? &m_annotation.value() : nullptr;
    }
    const LocalPotentialSampleList & SamplingEntries() const { return m_sampling_entries; }
};

} // namespace rhbm_gem
