#include "AtomicPotentialEntry.hpp"
#include "Constants.hpp"

#include <iostream>
#include <cmath>

AtomicPotentialEntry::AtomicPotentialEntry(void)
{

}

AtomicPotentialEntry::~AtomicPotentialEntry()
{

}

void AtomicPotentialEntry::AddDistanceAndMapValueList(std::vector<std::tuple<float, float>> list)
{
    m_distance_and_map_value_list = std::move(list);
}

void AtomicPotentialEntry::AddGausEstimateOLS(double v0, double v1)
{
    m_gaus_estimate_ols = std::make_tuple(v0, v1);
}

void AtomicPotentialEntry::AddGausEstimateMDPDE(double v0, double v1)
{
    m_gaus_estimate_mdpde = std::make_tuple(v0, v1);
}

void AtomicPotentialEntry::AddGausEstimatePosterior(const std::string & key, double v0, double v1)
{
    m_gaus_estimate_posterior_map[key] = std::make_tuple(v0, v1);
}

void AtomicPotentialEntry::AddGausVariancePosterior(const std::string & key, double v0, double v1)
{
    m_gaus_variance_posterior_map[key] = std::make_tuple(v0, v1);
}

void AtomicPotentialEntry::AddOutlierTag(const std::string & key, bool value)
{
    m_outlier_tag_map[key] = value;
}

void AtomicPotentialEntry::AddStatisticalDistance(const std::string & key, double value)
{
    m_statistical_distance_map[key] = value;
}

int AtomicPotentialEntry::GetDistanceAndMapValueListSize(void) const
{
    return static_cast<int>(m_distance_and_map_value_list.size());
}

const std::vector<std::tuple<float, float>> & AtomicPotentialEntry::GetDistanceAndMapValueList(void) const
{
    return m_distance_and_map_value_list;
}

double AtomicPotentialEntry::GetGausEstimateOLS(int par_id) const
{
    switch (par_id)
    {
        case 0: return GetAmplitudeEstimateOLS();
        case 1: return GetWidthEstimateOLS();
        case 2: return GetIntensityEstimateOLS();
        default:
            std::cerr <<"Invalid parameter index: "<< par_id << std::endl;
            return 0.0;
    }
}

double AtomicPotentialEntry::GetAmplitudeEstimateOLS(void) const
{
    return std::get<0>(m_gaus_estimate_ols);
}

double AtomicPotentialEntry::GetWidthEstimateOLS(void) const
{
    return std::get<1>(m_gaus_estimate_ols);
}

double AtomicPotentialEntry::GetIntensityEstimateOLS(void) const
{
    return CalculateIntensityEstimate(GetAmplitudeEstimateOLS(), GetWidthEstimateOLS());
}

double AtomicPotentialEntry::GetGausEstimateMDPDE(int par_id) const
{
    switch (par_id)
    {
        case 0: return GetAmplitudeEstimateMDPDE();
        case 1: return GetWidthEstimateMDPDE();
        case 2: return GetIntensityEstimateMDPDE();
        default:
            std::cerr <<"Invalid parameter index: "<< par_id << std::endl;
            return 0.0;
    }
}

double AtomicPotentialEntry::GetAmplitudeEstimateMDPDE(void) const
{
    return std::get<0>(m_gaus_estimate_mdpde);
}

double AtomicPotentialEntry::GetWidthEstimateMDPDE(void) const
{
    return std::get<1>(m_gaus_estimate_mdpde);
}

double AtomicPotentialEntry::GetIntensityEstimateMDPDE(void) const
{
    return CalculateIntensityEstimate(GetAmplitudeEstimateMDPDE(), GetWidthEstimateMDPDE());
}

double AtomicPotentialEntry::GetGausEstimatePosterior(const std::string & key, int par_id) const
{
    switch (par_id)
    {
        case 0: return GetAmplitudeEstimatePosterior(key);
        case 1: return GetWidthEstimatePosterior(key);
        case 2: return GetIntensityEstimatePosterior(key);
        default:
            std::cerr <<"Invalid parameter index: "<< par_id << std::endl;
            return 0.0;
    }
}

double AtomicPotentialEntry::GetAmplitudeEstimatePosterior(const std::string & key) const
{
    return std::get<0>(m_gaus_estimate_posterior_map.at(key));
}

double AtomicPotentialEntry::GetWidthEstimatePosterior(const std::string & key) const
{
    return std::get<1>(m_gaus_estimate_posterior_map.at(key));
}

double AtomicPotentialEntry::GetIntensityEstimatePosterior(const std::string & key) const
{
    return CalculateIntensityEstimate(
        GetAmplitudeEstimatePosterior(key), GetWidthEstimatePosterior(key));
}

double AtomicPotentialEntry::GetGausVariancePosterior(const std::string & key, int par_id) const
{
    switch (par_id)
    {
        case 0: return GetAmplitudeVariancePosterior(key);
        case 1: return GetWidthVariancePosterior(key);
        case 2: return GetIntensityVariancePosterior(key);
        default:
            std::cerr <<"Invalid parameter index: "<< par_id << std::endl;
            return 0.0;
    }
}

double AtomicPotentialEntry::GetAmplitudeVariancePosterior(const std::string & key) const
{
    return std::get<0>(m_gaus_variance_posterior_map.at(key));
}

double AtomicPotentialEntry::GetWidthVariancePosterior(const std::string & key) const
{
    return std::get<1>(m_gaus_variance_posterior_map.at(key));
}

double AtomicPotentialEntry::GetIntensityVariancePosterior(const std::string & key) const
{
    return CalculateIntensityVariance(
        GetAmplitudeEstimatePosterior(key), GetAmplitudeVariancePosterior(key),
        GetWidthEstimatePosterior(key), GetWidthVariancePosterior(key));
}

bool AtomicPotentialEntry::GetOutlierTag(const std::string & key) const
{
    return m_outlier_tag_map.at(key);
}

double AtomicPotentialEntry::GetStatisticalDistance(const std::string & key) const
{
    return m_statistical_distance_map.at(key);
}

double AtomicPotentialEntry::CalculateIntensityEstimate(double amplitude, double width) const
{
    if (width == 0.0) return 0.0;
    return amplitude * std::pow(Constants::two_pi * width * width, -1.5);
}

double AtomicPotentialEntry::CalculateIntensityVariance(
    double amplitude, double sigma_amplitude, double width, double sigma_width) const
{
    return std::sqrt(
        std::pow(std::pow(Constants::two_pi * width * width, -1.5) * sigma_amplitude, 2) +
        std::pow(-3.0 * amplitude * std::pow(Constants::two_pi, -1.5) * std::pow(width, -4) * sigma_width, 2)
    );
}