#pragma once

#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem {

struct LocalPotentialAnnotation
{
    GaussianPosterior posterior{};
    bool is_outlier{ false };
    double statistical_distance{ 0.0 };
};

class LocalPotentialEntry
{
    double m_alpha_r;
    LocalPotentialSampleList m_sampling_entries;
    GaussianEstimate m_gaus_estimate_ols;
    GaussianEstimate m_gaus_estimate_mdpde;
    std::unordered_map<std::string, LocalPotentialAnnotation> m_annotation_map;
    std::optional<RHBMMemberDataset> m_dataset;
    std::optional<RHBMBetaEstimateResult> m_fit_result;

public:
    LocalPotentialEntry();
    ~LocalPotentialEntry();

    void SetAlphaR(double value) { m_alpha_r = value; }
    void SetSamplingEntries(LocalPotentialSampleList value);
    void SetEstimateOLS(const GaussianEstimate & estimate) { m_gaus_estimate_ols = estimate; }
    void SetEstimateMDPDE(const GaussianEstimate & estimate) { m_gaus_estimate_mdpde = estimate; }
    void SetDataset(RHBMMemberDataset dataset);
    void SetFitResult(RHBMBetaEstimateResult value);
    void SetAnnotation(const std::string & key, LocalPotentialAnnotation annotation);
    void ClearTransientFitState();

    int GetSamplingEntryCount() const;
    double GetAlphaR() const { return m_alpha_r; }
    bool HasDataset() const { return m_dataset.has_value(); }
    bool HasFitResult() const { return m_fit_result.has_value(); }
    const GaussianEstimate & GetEstimateOLS() const { return m_gaus_estimate_ols; }
    const GaussianEstimate & GetEstimateMDPDE() const { return m_gaus_estimate_mdpde; }
    const RHBMMemberDataset & GetDataset() const;
    const RHBMBetaEstimateResult & GetFitResult() const;
    LocalPotentialAnnotation * FindAnnotation(const std::string & key);
    const LocalPotentialAnnotation * FindAnnotation(const std::string & key) const;
    const std::unordered_map<std::string, LocalPotentialAnnotation> & Annotations() const;
    const LocalPotentialSampleList & GetSamplingEntries() const;
    std::tuple<float, float> GetDistanceRange(double margin_rate = 0.0) const;
    std::tuple<float, float> GetResponseRange(double margin_rate = 0.0) const;
    SeriesPointList GetBinnedDistanceResponseSeries(
        int bin_size = 15,
        double x_min = 0.0,
        double x_max = 1.5) const;
    SeriesPointList GetLinearModelSeries() const;
    double GetMapValueNearCenter() const;
    double GetMomentZeroEstimate() const;
    double GetMomentTwoEstimate() const;
    double CalculateQScore(int par_choice) const;
};

} // namespace rhbm_gem
