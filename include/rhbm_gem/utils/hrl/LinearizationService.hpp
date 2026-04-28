#pragma once

#include <optional>
#include <tuple>

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

    static LinearizationSpec DefaultDataset();
    static LinearizationSpec DefaultMetricModel();
    static LinearizationSpec AtomLocalDecode();
    static LinearizationSpec AtomGroupDecode();
    static LinearizationSpec BondGroupDecode();
};

struct LinearizationContext
{
    std::optional<GaussianModel3D> model{};

    bool HasModelParameters() const;

    static LinearizationContext FromModelParameters(
        const GaussianParameterVector & model_parameters);
    static LinearizationContext FromModel(const GaussianModel3D & model);
};

SeriesPointList BuildDatasetSeries(
    const LinearizationSpec & spec,
    const LocalPotentialSampleList & sampling_entries,
    double x_min,
    double x_max,
    const LinearizationContext & context = {});

SeriesPointList BuildLinearModelSeries(
    const LinearizationSpec & spec,
    const LocalPotentialSampleList & sampling_entries,
    const LinearizationContext & context = {});

RHBMBetaVector EncodeGaussianToBeta(
    const LinearizationSpec & spec,
    const GaussianParameterVector & gaussian_parameters);

RHBMBetaVector EncodeGaussianToBeta(
    const LinearizationSpec & spec,
    const GaussianEstimate & gaussian_estimate);

RHBMBetaVector EncodeGaussianToBeta(
    const LinearizationSpec & spec,
    const GaussianModel3D & gaussian_model);

GaussianParameterVector DecodeLocalBeta(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model,
    const LinearizationContext & context = {});

GaussianParameterVector DecodeGroupBeta(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model);

GaussianModel3D DecodeLocalModel3D(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model,
    const LinearizationContext & context = {});

GaussianModel3D DecodeGroupModel3D(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model);

std::tuple<GaussianParameterVector, GaussianParameterVector> DecodePosterior(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix);

GaussianEstimate DecodeLocalEstimate(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model,
    const LinearizationContext & context = {});

GaussianEstimate DecodeGroupEstimate(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model);

GaussianEstimateWithUncertainty DecodeGaussianEstimateWithUncertainty(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix);

} // namespace rhbm_gem::linearization_service
