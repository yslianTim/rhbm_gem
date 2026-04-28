#pragma once

#include <optional>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem::linearization_service
{

enum class LinearizationKind
{
    LOG_QUADRATIC,
    TAYLOR_EXPANSION
};

enum class GaussianModelKind
{
    MODEL_2D,
    MODEL_3D
};

struct LinearizationSpec
{
    int basis_size{ 2 };
    LinearizationKind linearization_kind{ LinearizationKind::LOG_QUADRATIC };
    GaussianModelKind model_kind{ GaussianModelKind::MODEL_3D };
    bool requires_local_context{ false };

    static LinearizationSpec DefaultDataset() { return LinearizationSpec{}; }
    static LinearizationSpec DefaultMetricModel() { return LinearizationSpec{}; }
    static LinearizationSpec AtomLocalDecode() { return LinearizationSpec{}; }
    static LinearizationSpec AtomGroupDecode() { return LinearizationSpec{}; }
    static LinearizationSpec BondGroupDecode();
};

struct LinearizationContext
{
    std::optional<GaussianModel3D> model{};

    bool HasModelParameters() const { return model.has_value(); }
    static LinearizationContext FromModel(const GaussianModel3D & model);
};

SeriesPointList BuildDatasetSeries(
    const LinearizationSpec & spec,
    const LocalPotentialSampleList & sampling_entries,
    double x_min,
    double x_max,
    const LinearizationContext & context = {});

RHBMParameterVector EncodeGaussianToParameterVector(
    const LinearizationSpec & spec,
    const GaussianModel3D & gaussian_model);

GaussianModel3D DecodeParameterVector(
    const LinearizationSpec & spec,
    const RHBMParameterVector & parameter_vector,
    const LinearizationContext & context = {});

GaussianModel3DWithUncertainty DecodeParameterVector(
    const LinearizationSpec & spec,
    const RHBMParameterVector & parameter_vector,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix);

} // namespace rhbm_gem::linearization_service
