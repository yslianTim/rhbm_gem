#pragma once

#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/HRLModelTypes.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace HRLDataTransform
{
HRLMemberDataset BuildMemberDataset(const SeriesPointList & series_point_list);

HRLBetaMatrix BuildBetaMatrix(const std::vector<HRLBetaVector> & beta_list);

HRLGroupEstimationInput BuildGroupInput(
    const std::vector<HRLMemberDataset> & member_datasets,
    const std::vector<HRLBetaEstimateResult> & member_fit_results
);
} // namespace HRLDataTransform
