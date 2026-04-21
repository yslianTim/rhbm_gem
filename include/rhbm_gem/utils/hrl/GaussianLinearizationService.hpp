#pragma once

#include <tuple>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/GaussianStatistics.hpp>
#include <rhbm_gem/utils/hrl/HRLModelTypes.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem {

enum class GaussianLinearizationKind
{
    LOG_QUADRATIC,
    TAYLOR_EXPANSION
};

enum class GaussianModelKind
{
    MODEL_2D,
    MODEL_3D
};

struct GaussianLinearizationSpec
{
    int basis_size{ 2 };
    GaussianLinearizationKind linearization_kind{ GaussianLinearizationKind::LOG_QUADRATIC };
    GaussianModelKind model_kind{ GaussianModelKind::MODEL_3D };
    bool requires_local_context{ false };

    static GaussianLinearizationSpec DefaultDataset();
    static GaussianLinearizationSpec DefaultMetricModel();
    static GaussianLinearizationSpec AtomLocalDecode();
    static GaussianLinearizationSpec AtomGroupDecode();
    static GaussianLinearizationSpec BondGroupDecode();
};

struct GaussianLinearizationContext
{
    GaussianParameterVector model_parameters{};

    bool HasModelParameters() const;

    static GaussianLinearizationContext FromModelParameters(
        const GaussianParameterVector & model_parameters);
};

class GaussianLinearizationService
{
    GaussianLinearizationSpec m_spec;
    
public:
    explicit GaussianLinearizationService(GaussianLinearizationSpec spec);

    const GaussianLinearizationSpec & Spec() const { return m_spec; }

    SeriesPointList BuildDatasetSeries(
        const LocalPotentialSampleList & sampling_entries,
        double x_min,
        double x_max,
        const GaussianLinearizationContext & context = {}) const;
    SeriesPointList BuildLinearModelSeries(
        const LocalPotentialSampleList & sampling_entries,
        const GaussianLinearizationContext & context = {}) const;
    HRLMemberDataset BuildDataset(
        const LocalPotentialSampleList & sampling_entries,
        double x_min,
        double x_max,
        const GaussianLinearizationContext & context = {}) const;

    HRLBetaVector EncodeGaussianToBeta(const GaussianParameterVector & gaussian_parameters) const;
    HRLBetaVector EncodeGaussianToBeta(const GaussianEstimate & gaussian_estimate) const;

    GaussianParameterVector DecodeLocalBeta(
        const HRLBetaVector & linear_model,
        const GaussianLinearizationContext & context = {}) const;
    GaussianParameterVector DecodeGroupBeta(const HRLBetaVector & linear_model) const;

    std::tuple<GaussianParameterVector, GaussianParameterVector> DecodePosterior(
        const HRLBetaVector & linear_model,
        const HRLPosteriorCovarianceMatrix & covariance_matrix) const;

    GaussianEstimate DecodeLocalEstimate(
        const HRLBetaVector & linear_model,
        const GaussianLinearizationContext & context = {}) const;
    GaussianEstimate DecodeGroupEstimate(const HRLBetaVector & linear_model) const;
    GaussianPosterior DecodePosteriorEstimate(
        const HRLBetaVector & linear_model,
        const HRLPosteriorCovarianceMatrix & covariance_matrix) const;

private:
    void ValidateSpec() const;
    void ValidateContextIfRequired(const GaussianLinearizationContext & context) const;
    GaussianParameterVector BuildGaussianVector(
        const HRLBetaVector & linear_model,
        const GaussianLinearizationContext * context) const;

};

} // namespace rhbm_gem
