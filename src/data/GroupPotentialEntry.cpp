#include "GroupPotentialEntry.hpp"
#include "AtomObject.hpp"
#include "Constants.hpp"

#include <iostream>
#include <cmath>

GroupPotentialEntry::GroupPotentialEntry(void)
{

}

GroupPotentialEntry::~GroupPotentialEntry()
{

}

void GroupPotentialEntry::InsertGroupKey(uint64_t group_key)
{
    m_group_key_set.insert(group_key);
}

void GroupPotentialEntry::ReserveAtomObjectPtrList(uint64_t group_key, int size)
{
    m_group_atom_list_map[group_key].reserve(static_cast<size_t>(size));
}

void GroupPotentialEntry::AddAtomObjectPtr(uint64_t group_key, AtomObject * atom_object)
{
    m_group_atom_list_map[group_key].emplace_back(atom_object);
}

void GroupPotentialEntry::AddGausEstimateMean(uint64_t group_key, double v0, double v1)
{
    m_gaus_estimate_mean_map[group_key] = std::make_tuple(v0, v1);
}

void GroupPotentialEntry::AddGausEstimateMDPDE(uint64_t group_key, double v0, double v1)
{
    m_gaus_estimate_mdpde_map[group_key] = std::make_tuple(v0, v1);
}

void GroupPotentialEntry::AddGausEstimatePrior(uint64_t group_key, double v0, double v1)
{
    m_gaus_estimate_prior_map[group_key] = std::make_tuple(v0, v1);
}

void GroupPotentialEntry::AddGausVariancePrior(uint64_t group_key, double v0, double v1)
{
    m_gaus_variance_prior_map[group_key] = std::make_tuple(v0, v1);
}

int GroupPotentialEntry::GetAtomObjectPtrListSize(uint64_t group_key) const
{
    return static_cast<int>(m_group_atom_list_map.at(group_key).size());
}

const std::vector<AtomObject *> & GroupPotentialEntry::GetAtomObjectPtrList(uint64_t group_key) const
{
    return m_group_atom_list_map.at(group_key);
}

std::tuple<double, double> GroupPotentialEntry::GetGausEstimateMean(uint64_t group_key) const
{
    return m_gaus_estimate_mean_map.at(group_key);
}

std::tuple<double, double> GroupPotentialEntry::GetGausEstimateMDPDE(uint64_t group_key) const
{
    return m_gaus_estimate_mdpde_map.at(group_key);
}

std::tuple<double, double> GroupPotentialEntry::GetGausEstimatePrior(uint64_t group_key) const
{
    return m_gaus_estimate_prior_map.at(group_key);
}

std::tuple<double, double> GroupPotentialEntry::GetGausVariancePrior(uint64_t group_key) const
{
    return m_gaus_variance_prior_map.at(group_key);
}

const std::unordered_set<uint64_t> & GroupPotentialEntry::GetGroupKeySet(void) const
{
    return m_group_key_set;
}

double GroupPotentialEntry::GetGausEstimateMean(uint64_t group_key, int par_id) const
{
    return GetGausEstimateFromTuple(GetGausEstimateMean(group_key), par_id);
}

double GroupPotentialEntry::GetGausEstimateMDPDE(uint64_t group_key, int par_id) const
{
    return GetGausEstimateFromTuple(GetGausEstimateMDPDE(group_key), par_id);
}

double GroupPotentialEntry::GetGausEstimatePrior(uint64_t group_key, int par_id) const
{
    return GetGausEstimateFromTuple(GetGausEstimatePrior(group_key), par_id);
}

double GroupPotentialEntry::GetGausVariancePrior(uint64_t group_key, int par_id) const
{
    return GetGausVarianceFromTuple(
        GetGausEstimatePrior(group_key), GetGausVariancePrior(group_key), par_id);
}

double GroupPotentialEntry::GetGausEstimateFromTuple(
    const std::tuple<double, double> & estimate, int par_id) const
{
    switch (par_id)
    {
        case 0: return std::get<0>(estimate);
        case 1: return std::get<1>(estimate);
        case 2: return CalculateIntensityEstimate(estimate);
        default:
            std::cerr <<"Invalid parameter index: "<< par_id << std::endl;
            return 0.0;
    }
}

double GroupPotentialEntry::GetGausVarianceFromTuple(
    const std::tuple<double, double> & estimate,
    const std::tuple<double, double> & variance, int par_id) const
{
    switch (par_id)
    {
        case 0: return std::get<0>(variance);
        case 1: return std::get<1>(variance);
        case 2: return CalculateIntensityVariance(estimate, variance);
        default:
            std::cerr <<"Invalid parameter index: "<< par_id << std::endl;
            return 0.0;
    }
}

double GroupPotentialEntry::CalculateIntensityEstimate(
    const std::tuple<double, double> & estimate) const
{
    auto amplitude{ std::get<0>(estimate) };
    auto width{ std::get<1>(estimate) };
    if (width == 0.0) return 0.0;
    return amplitude * std::pow(Constants::two_pi * width * width, -1.5);
}

double GroupPotentialEntry::CalculateIntensityVariance(
    const std::tuple<double, double> & estimate,
    const std::tuple<double, double> & variance) const
{
    auto amplitude{ std::get<0>(estimate) };
    auto width{ std::get<1>(estimate) };
    auto sigma_amplitude{ std::get<0>(variance) };
    auto sigma_width{ std::get<1>(variance) };
    return std::sqrt(
        std::pow(std::pow(Constants::two_pi * width * width, -1.5) * sigma_amplitude, 2) +
        std::pow(-3.0 * amplitude * std::pow(Constants::two_pi, -1.5) * std::pow(width, -4) * sigma_width, 2)
    );
}
