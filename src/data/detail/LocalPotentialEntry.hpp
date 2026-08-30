#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>
#include <utility>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem {

class LocalPotentialEntry
{
    LocalPotentialSampleList m_raw_sampling_entries;
    LocalPotentialSampleList m_peeling_sampling_entries;
    std::array<LocalGaussianResult, 3> m_gaussian_results{};
    int m_neighbor_count_for_peeling{ 0 };

public:
    LocalPotentialEntry() = default;
    ~LocalPotentialEntry() = default;

    void SetAlphaR(FittingStage stage, double value)
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
    void SetGaussianResult(FittingStage stage, LocalGaussianResult value)
    {
        GaussianResult(stage) = std::move(value);
    }
    void SetPosteriorResult(
        FittingStage stage,
        GaussianModel3DWithUncertainty posterior,
        bool is_outlier,
        double statistical_distance)
    {
        auto & result{ GaussianResult(stage) };
        result.posterior = std::move(posterior);
        result.is_outlier = is_outlier;
        result.statistical_distance = statistical_distance;
    }
    void ClearTransientFitState(FittingStage stage)
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
    LocalGaussianResult & GaussianResult(FittingStage stage)
    {
        return m_gaussian_results.at(StageIndex(stage));
    }
    const LocalGaussianResult & GaussianResult(FittingStage stage) const
    {
        return m_gaussian_results.at(StageIndex(stage));
    }
    const LocalPotentialSampleList & RawSamplingEntries() const { return m_raw_sampling_entries; }
    const LocalPotentialSampleList & PeelingSamplingEntries() const { return m_peeling_sampling_entries; }

private:
    static std::size_t StageIndex(FittingStage stage)
    {
        const auto index{ static_cast<std::size_t>(stage) };
        if (index >= 3)
        {
            throw std::invalid_argument("Unknown local fitting stage.");
        }
        return index;
    }
};

} // namespace rhbm_gem
