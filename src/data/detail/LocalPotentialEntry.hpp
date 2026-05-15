#pragma once

#include <string>
#include <tuple>
#include <unordered_map>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem {

struct LocalPotentialAnnotation
{
    GaussianModel3DWithUncertainty gaussian{
        GaussianModel3D{ 0.0, 0.0 },
        GaussianModel3DUncertainty{}
    };
    bool is_outlier{ false };
    double statistical_distance{ 0.0 };
};

class LocalPotentialEntry
{
    LocalPotentialSampleList m_sampling_entries;
    LocalGaussianResult m_gaussian_result;
    std::unordered_map<std::string, LocalPotentialAnnotation> m_annotation_map;

public:
    LocalPotentialEntry();
    ~LocalPotentialEntry();

    void SetAlphaR(double value) { m_gaussian_result.alpha_r = value; }
    void SetSamplingEntries(LocalPotentialSampleList value);
    void SetEstimateOLS(const GaussianModel3D & estimate);
    void SetEstimateMDPDE(const GaussianModel3D & estimate);
    void SetGaussianResult(LocalGaussianResult value);
    void SetAnnotation(const std::string & key, LocalPotentialAnnotation annotation);
    void ClearTransientFitState();

    int GetSamplingEntryCount() const;
    double GetAlphaR() const { return m_gaussian_result.alpha_r; }
    const GaussianModel3D & GetEstimateOLS() const { return m_gaussian_result.ols.GetModel(); }
    const GaussianModel3D & GetEstimateMDPDE() const { return m_gaussian_result.mdpde.GetModel(); }
    const LocalGaussianResult & GetGaussianResult() const { return m_gaussian_result; }
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
    double GetMapValueNearCenter() const;
    double GetMomentZeroEstimate() const;
    double GetMomentTwoEstimate() const;
    double CalculateQScore(int par_choice) const;
};

} // namespace rhbm_gem
