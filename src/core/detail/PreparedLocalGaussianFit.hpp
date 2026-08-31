#pragma once

#include <cstddef>
#include <vector>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem::core::detail {

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

    RHBMMemberDataset BuildDataset(
        const std::vector<double> & sample_response_list,
        const GaussianModel3D & offset_model) const;

    LocalGaussianResult Estimate(
        const std::vector<double> & sample_response_list,
        double alpha_r,
        int thread_size,
        const GaussianModel3D & offset_model) const;
};

LocalPotentialSampleList BuildSamplesForZeroOffsetGaussianFit(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & model);

} // namespace rhbm_gem::core::detail
