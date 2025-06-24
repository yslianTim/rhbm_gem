#pragma once

#include <cstddef>
#include <vector>
#include <tuple>
#include <string>
#include <unordered_map>

class DataObjectBase;

class AtomicChargeEntry
{
    std::vector<std::tuple<float, float, float>> m_regression_data_list;
    std::vector<std::tuple<float, float>> m_distance_and_map_value_list;
    std::vector<std::tuple<float, float>> m_distance_and_neutral_map_value_list;
    std::vector<std::tuple<float, float>> m_distance_and_positive_map_value_list;
    std::vector<std::tuple<float, float>> m_distance_and_negative_map_value_list;
    std::tuple<double, double, double> m_model_estimate_ols;
    std::tuple<double, double, double> m_model_estimate_mdpde;
    std::unordered_map<std::string, std::tuple<double, double, double>> m_model_estimate_posterior_map;
    std::unordered_map<std::string, std::tuple<double, double, double>> m_model_variance_posterior_map;
    std::unordered_map<std::string, bool> m_outlier_tag_map;
    std::unordered_map<std::string, double> m_statistical_distance_map;

public:
    AtomicChargeEntry(void);
    ~AtomicChargeEntry();

    bool IsDataListSizeValid(void) const { return CheckDataListSize(); }
    void AddRegressionDataList(std::vector<std::tuple<float, float, float>> list);
    void AddDistanceAndMapValueList(std::vector<std::tuple<float, float>> list);
    void AddDistanceAndNeutralMapValueList(std::vector<std::tuple<float, float>> list);
    void AddDistanceAndPositiveMapValueList(std::vector<std::tuple<float, float>> list);
    void AddDistanceAndNegativeMapValueList(std::vector<std::tuple<float, float>> list);
    void AddModelEstimateOLS(double v0, double v1, double v2);
    void AddModelEstimateMDPDE(double v0, double v1, double v2);
    void AddModelEstimatePosterior(const std::string & key, double v0, double v1, double v2);
    void AddModelVariancePosterior(const std::string & key, double v0, double v1, double v2);
    void AddOutlierTag(const std::string & key, bool value);
    void AddStatisticalDistance(const std::string & key, double value);

    size_t GetDistanceAndMapValueListSize(void) const;
    const std::vector<std::tuple<float, float, float>> & GetRegressionDataList(void) const;
    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList(void) const;
    const std::vector<std::tuple<float, float>> & GetDistanceAndNeutralMapValueList(void) const;
    const std::vector<std::tuple<float, float>> & GetDistanceAndPositiveMapValueList(void) const;
    const std::vector<std::tuple<float, float>> & GetDistanceAndNegativeMapValueList(void) const;
    double GetModelEstimateOLS(int par_id) const;
    double GetModelEstimateMDPDE(int par_id) const;
    double GetModelEstimatePosterior(const std::string & key, int par_id) const;
    double GetModelVariancePosterior(const std::string & key, int par_id) const;
    bool GetOutlierTag(const std::string & key) const;
    double GetStatisticalDistance(const std::string & key) const;

private:
    bool CheckDataListSize(void) const;

};
