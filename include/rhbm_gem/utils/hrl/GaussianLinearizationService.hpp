#pragma once

#include <tuple>

#include <Eigen/Dense>

#include <rhbm_gem/data/object/GaussianStatistics.hpp>
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
    Eigen::VectorXd model_parameters{};

    bool HasModelParameters() const;

    static GaussianLinearizationContext FromModelParameters(const Eigen::VectorXd & model_parameters);
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

    Eigen::VectorXd EncodeGaussianToBeta(const Eigen::VectorXd & gaussian_parameters) const;
    Eigen::VectorXd EncodeGaussianToBeta(const GaussianEstimate & gaussian_estimate) const;

    Eigen::VectorXd DecodeLocalBeta(
        const Eigen::VectorXd & linear_model,
        const GaussianLinearizationContext & context = {}) const;
    Eigen::VectorXd DecodeGroupBeta(const Eigen::VectorXd & linear_model) const;

    std::tuple<Eigen::VectorXd, Eigen::VectorXd> DecodePosterior(
        const Eigen::VectorXd & linear_model,
        const Eigen::MatrixXd & covariance_matrix) const;

    GaussianEstimate DecodeLocalEstimate(
        const Eigen::VectorXd & linear_model,
        const GaussianLinearizationContext & context = {}) const;
    GaussianEstimate DecodeGroupEstimate(const Eigen::VectorXd & linear_model) const;
    GaussianPosterior DecodePosteriorEstimate(
        const Eigen::VectorXd & linear_model,
        const Eigen::MatrixXd & covariance_matrix) const;

private:
    void ValidateSpec() const;
    void ValidateContextIfRequired(const GaussianLinearizationContext & context) const;
    Eigen::VectorXd BuildGaussianVector(
        const Eigen::VectorXd & linear_model,
        const GaussianLinearizationContext * context) const;

};

} // namespace rhbm_gem
