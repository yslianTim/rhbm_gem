#pragma once

#include <utility>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem {

class LocalPotentialEntry
{
    LocalPotentialSampleList m_sampling_entries;
    LocalPotentialSampleList m_updated_sampling_entries;
    LocalGaussianResult m_gaussian_result;

public:
    LocalPotentialEntry() = default;
    ~LocalPotentialEntry() = default;

    void SetAlphaR(double value) { m_gaussian_result.alpha_r = value; }
    void SetSamplingEntries(LocalPotentialSampleList value)
    {
        m_sampling_entries = std::move(value);
    }
    void SetUpdatedSamplingEntries(LocalPotentialSampleList value)
    {
        m_updated_sampling_entries = std::move(value);
    }
    void SetGaussianResult(LocalGaussianResult value)
    {
        m_gaussian_result = std::move(value);
    }
    void SetPosteriorResult(
        GaussianModel3DWithUncertainty posterior,
        bool is_outlier,
        double statistical_distance)
    {
        m_gaussian_result.posterior = std::move(posterior);
        m_gaussian_result.is_outlier = is_outlier;
        m_gaussian_result.statistical_distance = statistical_distance;
    }
    void ClearTransientFitState()
    {
        m_gaussian_result.fit_result.reset();
        m_updated_sampling_entries.clear();
    }

    int SamplingEntryCount() const
    {
        return static_cast<int>(m_sampling_entries.size());
    }
    const LocalGaussianResult & GaussianResult() const { return m_gaussian_result; }
    const LocalPotentialSampleList & SamplingEntries() const { return m_sampling_entries; }
    const LocalPotentialSampleList & UpdatedSamplingEntries() const { return m_updated_sampling_entries; }
};

} // namespace rhbm_gem
