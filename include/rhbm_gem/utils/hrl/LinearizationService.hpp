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
    int basis_size{ 2 };
    GaussianModelKind model_kind{ GaussianModelKind::MODEL_3D };

    static LinearizationSpec DefaultDataset() { return LinearizationSpec{}; }
    static LinearizationSpec DefaultMetricModel() { return DefaultDataset(); }
    static LinearizationSpec AtomLocalDecode() { return DefaultDataset(); }
    static LinearizationSpec AtomGroupDecode() { return DefaultDataset(); }
    static LinearizationSpec BondGroupDecode();
};

SeriesPointList BuildDatasetSeries(
    const LinearizationSpec & spec,
    const LocalPotentialSampleList & sampling_entries,
    double x_min,
    double x_max);

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
