#pragma once

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem::linearization_service
{

SeriesPointList BuildDatasetSeries(
    const LocalPotentialSampleList & sampling_entries,
    double range_min,
    double range_max);

SeriesPointList BuildOffsetGaussianDatasetSeries(
    const LocalPotentialSampleList & sampling_entries,
    double range_min,
    double range_max,
    double tau);

RHBMParameterVector EncodeGaussianToParameterVector(const GaussianModel3D & gaussian_model);
GaussianModel3D DecodeParameterVector(const RHBMParameterVector & parameter_vector);
GaussianModel3DWithUncertainty DecodeParameterVector(
    const RHBMParameterVector & parameter_vector,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix);

} // namespace rhbm_gem::linearization_service
