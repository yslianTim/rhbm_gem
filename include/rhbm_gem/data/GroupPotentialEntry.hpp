#pragma once

#include <cstdint>
#include <tuple>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include <rhbm_gem/utils/GlobalEnumClass.hpp>

namespace rhbm_gem {

class AtomObject;
class BondObject;

class GroupPotentialEntry
{
    std::unordered_set<GroupKey> m_group_key_set;
    std::unordered_map<GroupKey, std::vector<AtomObject *>> m_group_atom_list_map;
    std::unordered_map<GroupKey, std::vector<BondObject *>> m_group_bond_list_map;
    std::unordered_map<GroupKey, std::tuple<double, double>> m_gaus_estimate_mean_map;
    std::unordered_map<GroupKey, std::tuple<double, double>> m_gaus_estimate_mdpde_map;
    std::unordered_map<GroupKey, std::tuple<double, double>> m_gaus_estimate_prior_map;
    std::unordered_map<GroupKey, std::tuple<double, double>> m_gaus_variance_prior_map;
    std::unordered_map<GroupKey, double> m_alpha_g_map;

public:
    GroupPotentialEntry();
    ~GroupPotentialEntry();
    void InsertGroupKey(GroupKey group_key);
    void ReserveAtomObjectPtrList(GroupKey group_key, int size);
    void ReserveBondObjectPtrList(GroupKey group_key, int size);
    void AddAtomObjectPtr(GroupKey group_key, AtomObject * atom_object);
    void AddBondObjectPtr(GroupKey group_key, BondObject * bond_object);
    void AddGausEstimateMean(GroupKey group_key, double v0, double v1);
    void AddGausEstimateMDPDE(GroupKey group_key, double v0, double v1);
    void AddGausEstimatePrior(GroupKey group_key, double v0, double v1);
    void AddGausVariancePrior(GroupKey group_key, double v0, double v1);
    void AddAlphaG(GroupKey group_key, double value);
    int GetAtomObjectPtrListSize(GroupKey group_key) const;
    int GetBondObjectPtrListSize(GroupKey group_key) const;
    double GetAlphaG(GroupKey group_key) const;
    const std::vector<AtomObject *> & GetAtomObjectPtrList(GroupKey group_key) const;
    const std::vector<BondObject *> & GetBondObjectPtrList(GroupKey group_key) const;
    std::tuple<double, double> GetGausEstimateMean(GroupKey group_key) const;
    std::tuple<double, double> GetGausEstimateMDPDE(GroupKey group_key) const;
    std::tuple<double, double> GetGausEstimatePrior(GroupKey group_key) const;
    std::tuple<double, double> GetGausVariancePrior(GroupKey group_key) const;
    double GetMuEstimateMean(GroupKey group_key, int par_id) const;
    double GetMuEstimateMDPDE(GroupKey group_key, int par_id) const;
    double GetMuEstimatePrior(GroupKey group_key, int par_id) const;
    double GetGausEstimateMean(GroupKey group_key, int par_id) const;
    double GetGausEstimateMDPDE(GroupKey group_key, int par_id) const;
    double GetGausEstimatePrior(GroupKey group_key, int par_id) const;
    double GetGausVariancePrior(GroupKey group_key, int par_id) const;
    const std::unordered_set<GroupKey> & GetGroupKeySet() const;

private:
    double GetGausEstimateFromTuple(const std::tuple<double, double> & estimate, int par_id) const;
    double GetGausVarianceFromTuple(const std::tuple<double, double> & estimate, const std::tuple<double, double> & variance, int par_id) const;
    double CalculateIntensityEstimate(const std::tuple<double, double> & estimate) const;
    double CalculateIntensityVariance(const std::tuple<double, double> & estimate, const std::tuple<double, double> & variance) const;

};

} // namespace rhbm_gem
