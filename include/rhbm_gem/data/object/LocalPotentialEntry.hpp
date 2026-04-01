#pragma once

#include <vector>
#include <tuple>
#include <string>
#include <unordered_map>

#include <rhbm_gem/data/object/GaussianStatistics.hpp>

namespace rhbm_gem {

class LocalPotentialEntry
{
    double m_alpha_r;
    std::vector<std::tuple<float, float>> m_distance_and_map_value_list;
    GaussianEstimate m_gaus_estimate_ols;
    GaussianEstimate m_gaus_estimate_mdpde;
    std::unordered_map<std::string, GaussianPosterior> m_gaus_posterior_map;
    std::unordered_map<std::string, bool> m_outlier_tag_map;
    std::unordered_map<std::string, double> m_statistical_distance_map;

public:
    LocalPotentialEntry();
    ~LocalPotentialEntry();

    void SetAlphaR(double value) { m_alpha_r = value; }
    void SetDistanceAndMapValueList(std::vector<std::tuple<float, float>> value);
    void SetEstimateOLS(const GaussianEstimate & estimate) { m_gaus_estimate_ols = estimate; }
    void SetEstimateMDPDE(const GaussianEstimate & estimate) { m_gaus_estimate_mdpde = estimate; }
    void SetPosterior(const std::string & key, const GaussianPosterior & posterior);
    void AddOutlierTag(const std::string & key, bool value);
    void AddStatisticalDistance(const std::string & key, double value);

    int GetDistanceAndMapValueListSize() const;
    double GetAlphaR() const { return m_alpha_r; }
    const GaussianEstimate & GetEstimateOLS() const { return m_gaus_estimate_ols; }
    const GaussianEstimate & GetEstimateMDPDE() const { return m_gaus_estimate_mdpde; }
    const GaussianPosterior & GetPosterior(const std::string & key) const;
    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList() const;
    double GetMapValueNearCenter() const;
    double GetMomentZeroEstimate() const;
    double GetMomentTwoEstimate() const;
    bool GetOutlierTag(const std::string & key) const;
    double GetStatisticalDistance(const std::string & key) const;
    double CalculateQScore(int par_choice) const;
};

} // namespace rhbm_gem
