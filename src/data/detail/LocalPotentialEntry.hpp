#pragma once

#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "data/detail/GaussianStatistics.hpp"

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
    std::vector<std::tuple<float, float>> m_distance_and_map_value_list;
    GaussianEstimate m_gaus_estimate_ols;
    GaussianEstimate m_gaus_estimate_mdpde;
    std::unordered_map<std::string, LocalPotentialAnnotation> m_annotation_map;

public:
    LocalPotentialEntry();
    ~LocalPotentialEntry();

    void SetAlphaR(double value) { m_alpha_r = value; }
    void SetDistanceAndMapValueList(std::vector<std::tuple<float, float>> value);
    void SetEstimateOLS(const GaussianEstimate & estimate) { m_gaus_estimate_ols = estimate; }
    void SetEstimateMDPDE(const GaussianEstimate & estimate) { m_gaus_estimate_mdpde = estimate; }
    void SetAnnotation(const std::string & key, LocalPotentialAnnotation annotation);

    int GetDistanceAndMapValueListSize() const;
    double GetAlphaR() const { return m_alpha_r; }
    const GaussianEstimate & GetEstimateOLS() const { return m_gaus_estimate_ols; }
    const GaussianEstimate & GetEstimateMDPDE() const { return m_gaus_estimate_mdpde; }
    LocalPotentialAnnotation * FindAnnotation(const std::string & key);
    const LocalPotentialAnnotation * FindAnnotation(const std::string & key) const;
    const std::unordered_map<std::string, LocalPotentialAnnotation> & Annotations() const;
    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList() const;
    double GetMapValueNearCenter() const;
    double GetMomentZeroEstimate() const;
    double GetMomentTwoEstimate() const;
    double CalculateQScore(int par_choice) const;
};

} // namespace rhbm_gem
