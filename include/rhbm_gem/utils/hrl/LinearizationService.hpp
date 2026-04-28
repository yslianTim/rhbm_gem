#pragma once

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem::linearization_service
{

enum class GaussianModelKind
{
    MODEL_2D,
    MODEL_3D
};

struct LinearizationSpec
{
    GaussianModelKind model_kind{ GaussianModelKind::MODEL_3D };

    static LinearizationSpec AtomDecode() { return LinearizationSpec{}; }
    static LinearizationSpec BondDecode();
};

struct LinearizationRange
{
    double min{ 0.0 };
    double max{ 1.0 };
};

SeriesPointList BuildDatasetSeries(
    const LocalPotentialSampleList & sampling_entries,
    const LinearizationRange & fit_range);

RHBMParameterVector EncodeGaussianToParameterVector(
    const LinearizationSpec & spec,
    const GaussianModel3D & gaussian_model);

GaussianModel3D DecodeParameterVector(
    const LinearizationSpec & spec,
    const RHBMParameterVector & parameter_vector);

GaussianModel3DWithUncertainty DecodeParameterVector(
    const LinearizationSpec & spec,
    const RHBMParameterVector & parameter_vector,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix);

} // namespace rhbm_gem::linearization_service
