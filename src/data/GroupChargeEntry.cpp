#include "GroupChargeEntry.hpp"
#include "AtomObject.hpp"
#include "Constants.hpp"

#include <iostream>
#include <cmath>

GroupChargeEntry::GroupChargeEntry(void)
{

}

GroupChargeEntry::~GroupChargeEntry()
{

}

void GroupChargeEntry::InsertGroupKey(uint64_t group_key)
{
    m_group_key_set.insert(group_key);
}

void GroupChargeEntry::ReserveAtomObjectPtrList(uint64_t group_key, int size)
{
    m_group_atom_list_map[group_key].reserve(static_cast<size_t>(size));
}

void GroupChargeEntry::AddAtomObjectPtr(uint64_t group_key, AtomObject * atom_object)
{
    m_group_atom_list_map[group_key].emplace_back(atom_object);
}

void GroupChargeEntry::AddModelEstimateMean(uint64_t group_key, double v0, double v1, double v2)
{
    m_model_estimate_mean_map[group_key] = std::make_tuple(v0, v1, v2);
}

void GroupChargeEntry::AddModelEstimateMDPDE(uint64_t group_key, double v0, double v1, double v2)
{
    m_model_estimate_mdpde_map[group_key] = std::make_tuple(v0, v1, v2);
}

void GroupChargeEntry::AddModelEstimatePrior(uint64_t group_key, double v0, double v1, double v2)
{
    m_model_estimate_prior_map[group_key] = std::make_tuple(v0, v1, v2);
}

void GroupChargeEntry::AddModelVariancePrior(uint64_t group_key, double v0, double v1, double v2)
{
    m_model_variance_prior_map[group_key] = std::make_tuple(v0, v1, v2);
}

int GroupChargeEntry::GetAtomObjectPtrListSize(uint64_t group_key) const
{
    return static_cast<int>(m_group_atom_list_map.at(group_key).size());
}

const std::vector<AtomObject *> & GroupChargeEntry::GetAtomObjectPtrList(uint64_t group_key) const
{
    return m_group_atom_list_map.at(group_key);
}

std::tuple<double, double, double> GroupChargeEntry::GetModelEstimateMean(uint64_t group_key) const
{
    return m_model_estimate_mean_map.at(group_key);
}

std::tuple<double, double, double> GroupChargeEntry::GetModelEstimateMDPDE(uint64_t group_key) const
{
    return m_model_estimate_mdpde_map.at(group_key);
}

std::tuple<double, double, double> GroupChargeEntry::GetModelEstimatePrior(uint64_t group_key) const
{
    return m_model_estimate_prior_map.at(group_key);
}

std::tuple<double, double, double> GroupChargeEntry::GetModelVariancePrior(uint64_t group_key) const
{
    return m_model_variance_prior_map.at(group_key);
}

const std::unordered_set<uint64_t> & GroupChargeEntry::GetGroupKeySet(void) const
{
    return m_group_key_set;
}

double GroupChargeEntry::GetModelEstimateMean(uint64_t group_key, int par_id) const
{
    return GetModelEstimateFromTuple(GetModelEstimateMean(group_key), par_id);
}

double GroupChargeEntry::GetModelEstimateMDPDE(uint64_t group_key, int par_id) const
{
    return GetModelEstimateFromTuple(GetModelEstimateMDPDE(group_key), par_id);
}

double GroupChargeEntry::GetModelEstimatePrior(uint64_t group_key, int par_id) const
{
    return GetModelEstimateFromTuple(GetModelEstimatePrior(group_key), par_id);
}

double GroupChargeEntry::GetModelVariancePrior(uint64_t group_key, int par_id) const
{
    return GetModelVarianceFromTuple(GetModelVariancePrior(group_key), par_id);
}

double GroupChargeEntry::GetModelEstimateFromTuple(
    const std::tuple<double, double, double> & estimate, int par_id) const
{
    switch (par_id)
    {
        case 0: return std::get<0>(estimate);
        case 1: return std::get<1>(estimate);
        case 2: return std::get<2>(estimate);
        default:
            std::cerr <<"Invalid parameter index: "<< par_id << std::endl;
            return 0.0;
    }
}

double GroupChargeEntry::GetModelVarianceFromTuple(
    const std::tuple<double, double, double> & variance, int par_id) const
{
    switch (par_id)
    {
        case 0: return std::get<0>(variance);
        case 1: return std::get<1>(variance);
        case 2: return std::get<2>(variance);
        default:
            std::cerr <<"Invalid parameter index: "<< par_id << std::endl;
            return 0.0;
    }
}