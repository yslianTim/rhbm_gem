#pragma once

#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include <vector>
#include <tuple>
#include <string>
#include <unordered_map>
#include <Eigen/Dense>

class DataObjectBase;

class AtomicPotentialEntry
{
    std::vector<std::tuple<float, float>> m_distance_and_map_value_list;
    std::tuple<double, double> m_gaus_estimate_ols;
    std::tuple<double, double> m_gaus_estimate_mdpde;
    std::unordered_map<std::string, std::tuple<double, double>> m_gaus_estimate_posterior_map;
    std::unordered_map<std::string, std::tuple<double, double>> m_gaus_variance_posterior_map;
    std::unordered_map<std::string, bool> m_outlier_tag_map;
    std::unordered_map<std::string, double> m_statistical_distance_map;

public:
    AtomicPotentialEntry(void) = default;
    ~AtomicPotentialEntry() = default;

    void AddDistanceAndMapValueList(std::vector<std::tuple<float, float>> list) { m_distance_and_map_value_list = std::move(list); }
    void AddGausEstimateOLS(Eigen::VectorXd vec) { m_gaus_estimate_ols = std::make_tuple(vec(0), vec(1)); }
    void AddGausEstimateMDPDE(Eigen::VectorXd vec) { m_gaus_estimate_mdpde = std::make_tuple(vec(0), vec(1)); }
    void AddGausEstimatePosterior(const std::string & key, Eigen::VectorXd vec) { m_gaus_estimate_posterior_map[key] = std::make_tuple(vec(0), vec(1)); }
    void AddGausVariancePosterior(const std::string & key, Eigen::VectorXd vec) { m_gaus_variance_posterior_map[key] = std::make_tuple(vec(0), vec(1)); }
    void AddOutlierTag(const std::string & key, bool value) { m_outlier_tag_map[key] = value; }
    void AddStatisticalDistance(const std::string & key, double value) { m_statistical_distance_map[key] = value; }

    int GetDistanceAndMapValueListSize(void) const { return static_cast<int>(m_distance_and_map_value_list.size()); }
    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList(void) const { return m_distance_and_map_value_list; }

};