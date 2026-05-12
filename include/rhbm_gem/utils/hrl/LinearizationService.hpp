#pragma once

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem::linearization_service
{

struct LinearizationRange
{
    double min{ 0.0 };
    double max{ 1.0 };
};

SeriesPointList BuildDatasetSeries(
    const LocalPotentialSampleList & sampling_entries,
    const LinearizationRange & fit_range);

RHBMParameterVector EncodeGaussianToParameterVector(const GaussianModel3D & gaussian_model);
GaussianModel3D DecodeParameterVector(const RHBMParameterVector & parameter_vector);
GaussianModel3DWithUncertainty DecodeParameterVector(
    const RHBMParameterVector & parameter_vector,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix);

} // namespace rhbm_gem::linearization_service
