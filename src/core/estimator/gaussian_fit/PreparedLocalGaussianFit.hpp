#pragma once

#include <cstddef>
#include <vector>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem::core::detail {

// Both operations consume the identical finite, positive, radius-selected design.
double TrainLocalGaussianAlpha(const RHBMMemberDataset &, int thread_size);
LocalGaussianResult EstimateLocalGaussianDataset(const RHBMMemberDataset &, double alpha_r,
    int thread_size, double offset);

class PreparedLocalGaussianDesign
{
    struct Row
    {
        std::size_t source_sample_index{ 0 };
        double distance{ 0.0 };
    };

    std::size_t m_source_sample_count{ 0 };
    std::vector<Row> m_row_list{};
    RHBMDesignMatrix m_design_matrix{};

public:
    PreparedLocalGaussianDesign() = default;

    PreparedLocalGaussianDesign(
        const LocalPotentialSampleList & sample_entries,
        double range_min,
        double range_max);

    // Exact retained rows, including an empty dataset when no response is positive.
    RHBMMemberDataset BuildSelectedDataset(
        const std::vector<double> & sample_response_list,
        const GaussianModel3D & offset_model) const;

    RHBMMemberDataset BuildDataset(
        const std::vector<double> & sample_response_list,
        const GaussianModel3D & offset_model) const;

    LocalGaussianResult Estimate(
        const std::vector<double> & sample_response_list,
        double alpha_r,
        int thread_size,
        const GaussianModel3D & offset_model) const;
};

} // namespace rhbm_gem::core::detail
