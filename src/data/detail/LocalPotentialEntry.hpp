#pragma once

#include <stdexcept>
#include <utility>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem {

class LocalPotentialEntry
{
    LocalPotentialSampleList m_raw_sampling_entries;
    LocalPotentialSampleList m_peeling_sampling_entries;
    LocalGaussianResult m_gaussian_result_1st;
    LocalGaussianResult m_gaussian_result_2nd;
    LocalGaussianResult m_gaussian_result_3rd;
    int m_neighbor_count_for_peeling{ 0 };

public:
    LocalPotentialEntry() = default;
    ~LocalPotentialEntry() = default;

    void SetAlphaR(LocalFittingStage stage, double value)
    {
        GaussianResult(stage).alpha_r = value;
    }
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
    void SetGaussianResult(LocalFittingStage stage, LocalGaussianResult value)
    {
        GaussianResult(stage) = std::move(value);
    }
    void SetPosteriorResult(
        LocalFittingStage stage,
        GaussianModel3DWithUncertainty posterior,
        bool is_outlier,
        double statistical_distance)
    {
        auto & result{ GaussianResult(stage) };
        result.posterior = std::move(posterior);
        result.is_outlier = is_outlier;
        result.statistical_distance = statistical_distance;
    }
    void ClearTransientFitState(LocalFittingStage stage)
    {
        GaussianResult(stage).fit_result.reset();
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
    LocalGaussianResult & GaussianResult(LocalFittingStage stage)
    {
        switch (stage)
        {
            case LocalFittingStage::First:  return m_gaussian_result_1st;
            case LocalFittingStage::Second: return m_gaussian_result_2nd;
            case LocalFittingStage::Third:  return m_gaussian_result_3rd;
        }
        throw std::invalid_argument("Unknown local fitting stage.");
    }
    const LocalGaussianResult & GaussianResult(LocalFittingStage stage) const
    {
        switch (stage)
        {
            case LocalFittingStage::First:  return m_gaussian_result_1st;
            case LocalFittingStage::Second: return m_gaussian_result_2nd;
            case LocalFittingStage::Third:  return m_gaussian_result_3rd;
        }
        throw std::invalid_argument("Unknown local fitting stage.");
    }
    const LocalPotentialSampleList & RawSamplingEntries() const { return m_raw_sampling_entries; }
    const LocalPotentialSampleList & PeelingSamplingEntries() const { return m_peeling_sampling_entries; }
};

} // namespace rhbm_gem
