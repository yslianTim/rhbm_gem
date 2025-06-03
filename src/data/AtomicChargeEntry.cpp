#include "AtomicChargeEntry.hpp"
#include "Constants.hpp"

#include <iostream>
#include <cmath>

AtomicChargeEntry::AtomicChargeEntry(void)
{

}

AtomicChargeEntry::~AtomicChargeEntry()
{

}

void AtomicChargeEntry::AddRegressionDataList(std::vector<std::tuple<float, float, float>> list)
{
    m_regression_data_list = std::move(list);
}

void AtomicChargeEntry::AddDistanceAndMapValueList(std::vector<std::tuple<float, float>> list)
{
    m_distance_and_map_value_list = std::move(list);
}

void AtomicChargeEntry::AddDistanceAndNeutralMapValueList(std::vector<std::tuple<float, float>> list)
{
    m_distance_and_neutral_map_value_list = std::move(list);
}

void AtomicChargeEntry::AddDistanceAndPositiveMapValueList(std::vector<std::tuple<float, float>> list)
{
    m_distance_and_positive_map_value_list = std::move(list);
}

void AtomicChargeEntry::AddDistanceAndNegativeMapValueList(std::vector<std::tuple<float, float>> list)
{
    m_distance_and_negative_map_value_list = std::move(list);
}

void AtomicChargeEntry::AddModelEstimateOLS(double v0, double v1, double v2)
{
    m_model_estimate_ols = std::make_tuple(v0, v1, v2);
}

void AtomicChargeEntry::AddModelEstimateMDPDE(double v0, double v1, double v2)
{
    m_model_estimate_mdpde = std::make_tuple(v0, v1, v2);
}

void AtomicChargeEntry::AddModelEstimatePosterior(const std::string & key, double v0, double v1, double v2)
{
    m_model_estimate_posterior_map[key] = std::make_tuple(v0, v1, v2);
}

void AtomicChargeEntry::AddModelVariancePosterior(const std::string & key, double v0, double v1, double v2)
{
    m_model_variance_posterior_map[key] = std::make_tuple(v0, v1, v2);
}

void AtomicChargeEntry::AddOutlierTag(const std::string & key, bool value)
{
    m_outlier_tag_map[key] = value;
}

void AtomicChargeEntry::AddStatisticalDistance(const std::string & key, double value)
{
    m_statistical_distance_map[key] = value;
}

size_t AtomicChargeEntry::GetDistanceAndMapValueListSize(void) const
{
    return m_distance_and_map_value_list.size();
}

const std::vector<std::tuple<float, float, float>> & AtomicChargeEntry::GetRegressionDataList(void) const
{
    return m_regression_data_list;
}

const std::vector<std::tuple<float, float>> & AtomicChargeEntry::GetDistanceAndMapValueList(void) const
{
    return m_distance_and_map_value_list;
}

const std::vector<std::tuple<float, float>> & AtomicChargeEntry::GetDistanceAndNeutralMapValueList(void) const
{
    return m_distance_and_neutral_map_value_list;
}

const std::vector<std::tuple<float, float>> & AtomicChargeEntry::GetDistanceAndPositiveMapValueList(void) const
{
    return m_distance_and_positive_map_value_list;
}

const std::vector<std::tuple<float, float>> & AtomicChargeEntry::GetDistanceAndNegativeMapValueList(void) const
{
    return m_distance_and_negative_map_value_list;
}

double AtomicChargeEntry::GetModelEstimateOLS(int par_id) const
{
    switch (par_id)
    {
        case 0: return std::get<0>(m_model_estimate_ols);
        case 1: return std::get<1>(m_model_estimate_ols);
        case 2: return std::get<2>(m_model_estimate_ols);
        default:
            std::cerr <<"Invalid parameter index: "<< par_id << std::endl;
            return 0.0;
    }
}

double AtomicChargeEntry::GetModelEstimateMDPDE(int par_id) const
{
    switch (par_id)
    {
        case 0: return std::get<0>(m_model_estimate_mdpde);
        case 1: return std::get<1>(m_model_estimate_mdpde);
        case 2: return std::get<2>(m_model_estimate_mdpde);
        default:
            std::cerr <<"Invalid parameter index: "<< par_id << std::endl;
            return 0.0;
    }
}

double AtomicChargeEntry::GetModelEstimatePosterior(const std::string & key, int par_id) const
{
    if (m_model_estimate_posterior_map.find(key) == m_model_estimate_posterior_map.end())
    {
        std::cerr << "Key not found in model estimate posterior map: " << key << std::endl;
        return 0.0;
    }
    switch (par_id)
    {
        case 0: return std::get<0>(m_model_estimate_posterior_map.at(key));
        case 1: return std::get<1>(m_model_estimate_posterior_map.at(key));
        case 2: return std::get<2>(m_model_estimate_posterior_map.at(key));
        default:
            std::cerr <<"Invalid parameter index: "<< par_id << std::endl;
            return 0.0;
    }
}

double AtomicChargeEntry::GetModelVariancePosterior(const std::string & key, int par_id) const
{
    if (m_model_variance_posterior_map.find(key) == m_model_variance_posterior_map.end())
    {
        std::cerr << "Key not found in model variance posterior map: " << key << std::endl;
        return 0.0;
    }
    switch (par_id)
    {
        case 0: return std::get<0>(m_model_variance_posterior_map.at(key));
        case 1: return std::get<1>(m_model_variance_posterior_map.at(key));
        case 2: return std::get<2>(m_model_variance_posterior_map.at(key));
        default:
            std::cerr <<"Invalid parameter index: "<< par_id << std::endl;
            return 0.0;
    }
}

bool AtomicChargeEntry::GetOutlierTag(const std::string & key) const
{
    return m_outlier_tag_map.at(key);
}

double AtomicChargeEntry::GetStatisticalDistance(const std::string & key) const
{
    return m_statistical_distance_map.at(key);
}

bool AtomicChargeEntry::CheckDataListSize(void) const
{
    if (m_distance_and_map_value_list.size() != m_distance_and_neutral_map_value_list.size() ||
        m_distance_and_map_value_list.size() != m_distance_and_positive_map_value_list.size() ||
        m_distance_and_map_value_list.size() != m_distance_and_negative_map_value_list.size())
    {
        std::cerr << "Data list sizes do not match!" << std::endl;
        return false;
    }
    return true;
}