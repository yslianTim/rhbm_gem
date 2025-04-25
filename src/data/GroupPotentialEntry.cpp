#include "GroupPotentialEntry.hpp"
#include "AtomObject.hpp"

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

void GroupPotentialEntry::AddGausEstimateMean(uint64_t group_key, Eigen::VectorXd vec)
{
    m_gaus_estimate_mean_map[group_key] = std::make_tuple(vec(0), vec(1));
}

void GroupPotentialEntry::AddGausEstimateMean(uint64_t group_key, double v0, double v1)
{
    m_gaus_estimate_mean_map[group_key] = std::make_tuple(v0, v1);
}

void GroupPotentialEntry::AddGausEstimateMDPDE(uint64_t group_key, Eigen::VectorXd vec)
{
    m_gaus_estimate_mdpde_map[group_key] = std::make_tuple(vec(0), vec(1));
}

void GroupPotentialEntry::AddGausEstimateMDPDE(uint64_t group_key, double v0, double v1)
{
    m_gaus_estimate_mdpde_map[group_key] = std::make_tuple(v0, v1);
}

void GroupPotentialEntry::AddGausEstimatePrior(uint64_t group_key, Eigen::VectorXd vec)
{
    m_gaus_estimate_prior_map[group_key] = std::make_tuple(vec(0), vec(1));
}

void GroupPotentialEntry::AddGausEstimatePrior(uint64_t group_key, double v0, double v1)
{
    m_gaus_estimate_prior_map[group_key] = std::make_tuple(v0, v1);
}

void GroupPotentialEntry::AddGausVariancePrior(uint64_t group_key, Eigen::VectorXd vec)
{
    m_gaus_variance_prior_map[group_key] = std::make_tuple(vec(0), vec(1));
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