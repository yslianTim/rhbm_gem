#pragma once

#include <utility>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem {

class LocalPotentialEntry
{
    LocalPotentialSampleList m_raw_sampling_entries;
    LocalPotentialSampleList m_peeling_sampling_entries;
    LocalGaussianResult m_gaussian_result;
    int m_neighbor_count_for_peeling{ 0 };

public:
    LocalPotentialEntry() = default;
    ~LocalPotentialEntry() = default;

    void SetAlphaR(double value) { m_gaussian_result.alpha_r = value; }
    void SetRawSamplingEntries(LocalPotentialSampleList value)
    {
        m_raw_sampling_entries = std::move(value);
    }
    void SetPeelingSamplingEntries(LocalPotentialSampleList value)
    {
        m_peeling_sampling_entries = std::move(value);
    }
    void SetNeighborCountForPeeling(int value)
    {
        m_neighbor_count_for_peeling = value;
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
    }

    int RawSamplingEntryCount() const
    {
        return static_cast<int>(m_raw_sampling_entries.size());
    }
    int PeelingSamplingEntryCount() const
    {
        return static_cast<int>(m_peeling_sampling_entries.size());
    }
    int NeighborCountForPeeling() const { return m_neighbor_count_for_peeling; }
    const LocalGaussianResult & GaussianResult() const { return m_gaussian_result; }
    const LocalPotentialSampleList & RawSamplingEntries() const { return m_raw_sampling_entries; }
    const LocalPotentialSampleList & PeelingSamplingEntries() const { return m_peeling_sampling_entries; }
};

} // namespace rhbm_gem
