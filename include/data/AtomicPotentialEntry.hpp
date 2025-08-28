#pragma once

#include <vector>
#include <tuple>
#include <string>
#include <unordered_map>

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
    AtomicPotentialEntry(void);
    ~AtomicPotentialEntry();

    void AddDistanceAndMapValueList(std::vector<std::tuple<float, float>> && list);
    void ClearDistanceAndMapValueList(void);
    void AddGausEstimateOLS(double v0, double v1);
    void AddGausEstimateMDPDE(double v0, double v1);
    void AddGausEstimatePosterior(const std::string & key, double v0, double v1);
    void AddGausVariancePosterior(const std::string & key, double v0, double v1);
    void AddOutlierTag(const std::string & key, bool value);
    void AddStatisticalDistance(const std::string & key, double value);

    int GetDistanceAndMapValueListSize(void) const;
    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList(void) const;
    double GetGausEstimateOLS(int par_id) const;
    double GetAmplitudeEstimateOLS(void) const;
    double GetWidthEstimateOLS(void) const;
    double GetIntensityEstimateOLS(void) const;
    double GetGausEstimateMDPDE(int par_id) const;
    double GetAmplitudeEstimateMDPDE(void) const;
    double GetWidthEstimateMDPDE(void) const;
    double GetIntensityEstimateMDPDE(void) const;
    double GetGausEstimatePosterior(const std::string & key, int par_id) const;
    double GetAmplitudeEstimatePosterior(const std::string & key) const;
    double GetWidthEstimatePosterior(const std::string & key) const;
    double GetIntensityEstimatePosterior(const std::string & key) const;
    double GetGausVariancePosterior(const std::string & key, int par_id) const;
    double GetAmplitudeVariancePosterior(const std::string & key) const;
    double GetWidthVariancePosterior(const std::string & key) const;
    double GetIntensityVariancePosterior(const std::string & key) const;
    bool GetOutlierTag(const std::string & key) const;
    double GetStatisticalDistance(const std::string & key) const;

private:
    double CalculateIntensityEstimate(double amplitude, double width) const;
    double CalculateIntensityVariance(double amplitude, double sigma_amplitude, double width, double sigma_width) const;

};
