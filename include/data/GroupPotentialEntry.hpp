#pragma once

#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include <set>
#include <tuple>
#include <string>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>
#include "GroupPotentialEntryBase.hpp"
#include "TupleHashHelper.hpp"
#include "AtomObject.hpp"

template <typename GroupKey>
class GroupPotentialEntry : public GroupPotentialEntryBase
{
    std::set<GroupKey> m_group_key_set;
    std::unordered_map<GroupKey, std::vector<AtomObject *>, GenericTupleHash<GroupKey>> m_group_atom_list_map;
    std::unordered_map<GroupKey, std::tuple<double, double>, GenericTupleHash<GroupKey>> m_gaus_estimate_mean_map;
    std::unordered_map<GroupKey, std::tuple<double, double>, GenericTupleHash<GroupKey>> m_gaus_estimate_mdpde_map;
    std::unordered_map<GroupKey, std::tuple<double, double>, GenericTupleHash<GroupKey>> m_gaus_estimate_prior_map;
    std::unordered_map<GroupKey, std::tuple<double, double>, GenericTupleHash<GroupKey>> m_gaus_variance_prior_map;

public:
    GroupPotentialEntry(void) = default;
    ~GroupPotentialEntry() = default;
    void InsertGroupKey(const void * group_key) override { m_group_key_set.insert(GetKey(group_key)); }
    void ReserveAtomObjectPtrList(const void * group_key, int size) override { m_group_atom_list_map[GetKey(group_key)].reserve(size); }
    void AddAtomObjectPtr(const void * group_key, AtomObject * atom_object) override { m_group_atom_list_map[GetKey(group_key)].emplace_back(atom_object); }
    void AddGausEstimateMean(const void * group_key, Eigen::VectorXd vec) override { m_gaus_estimate_mean_map[GetKey(group_key)] = std::make_tuple(vec(0), vec(1)); }
    void AddGausEstimateMean(const void * group_key, double v0, double v1) override { m_gaus_estimate_mean_map[GetKey(group_key)] = std::make_tuple(v0, v1); }
    void AddGausEstimateMDPDE(const void * group_key, Eigen::VectorXd vec) override { m_gaus_estimate_mdpde_map[GetKey(group_key)] = std::make_tuple(vec(0), vec(1)); }
    void AddGausEstimateMDPDE(const void * group_key, double v0, double v1) override { m_gaus_estimate_mdpde_map[GetKey(group_key)] = std::make_tuple(v0, v1); }
    void AddGausEstimatePrior(const void * group_key, Eigen::VectorXd vec) override { m_gaus_estimate_prior_map[GetKey(group_key)] = std::make_tuple(vec(0), vec(1)); }
    void AddGausEstimatePrior(const void * group_key, double v0, double v1) override { m_gaus_estimate_prior_map[GetKey(group_key)] = std::make_tuple(v0, v1); }
    void AddGausVariancePrior(const void * group_key, Eigen::VectorXd vec) override { m_gaus_variance_prior_map[GetKey(group_key)] = std::make_tuple(vec(0), vec(1)); }
    void AddGausVariancePrior(const void * group_key, double v0, double v1) override { m_gaus_variance_prior_map[GetKey(group_key)] = std::make_tuple(v0, v1); }
    int GetAtomObjectPtrListSize(const void * group_key) const override { return static_cast<int>(m_group_atom_list_map.at(GetKey(group_key)).size()); }
    const std::vector<AtomObject *> & GetAtomObjectPtrList(const void * group_key) const override { return m_group_atom_list_map.at(GetKey(group_key)); }
    std::tuple<double, double> GetGausEstimateMean(const void * group_key) const override { return m_gaus_estimate_mean_map.at(GetKey(group_key)); }
    std::tuple<double, double> GetGausEstimateMDPDE(const void * group_key) const override { return m_gaus_estimate_mdpde_map.at(GetKey(group_key)); }
    std::tuple<double, double> GetGausEstimatePrior(const void * group_key) const override { return m_gaus_estimate_prior_map.at(GetKey(group_key)); }
    std::tuple<double, double> GetGausVariancePrior(const void * group_key) const override { return m_gaus_variance_prior_map.at(GetKey(group_key)); }

    const std::set<GroupKey> & GetGroupKeySet(void) const { return m_group_key_set; }

private:
    GroupKey GetKey(const void * group_key) const { return *reinterpret_cast<const GroupKey*>(group_key); }

};