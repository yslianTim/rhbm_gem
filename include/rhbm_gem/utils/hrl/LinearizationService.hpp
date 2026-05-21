#pragma once

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem::linearization_service
{

SeriesPointList BuildDatasetSeries(
    const LocalPotentialSampleList & sampling_entries,
    double range_min,
    double range_max,
    LocalGaussianFitModel fit_model);

RHBMParameterVector EncodeGaussianToParameterVector(const GaussianModel3D & gaussian_model);
RHBMParameterVector EncodeGaussianToParameterVector(
    const GaussianModel3D & gaussian_model,
    LocalGaussianFitModel fit_model);
GaussianModel3D DecodeParameterVector(const RHBMParameterVector & parameter_vector);
GaussianModel3D DecodeParameterVector(
    const RHBMParameterVector & parameter_vector,
    LocalGaussianFitModel fit_model);
GaussianModel3DWithUncertainty DecodeParameterVector(
    const RHBMParameterVector & parameter_vector,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix);
GaussianModel3DWithUncertainty DecodeParameterVector(
    const RHBMParameterVector & parameter_vector,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix,
    LocalGaussianFitModel fit_model);

} // namespace rhbm_gem::linearization_service
