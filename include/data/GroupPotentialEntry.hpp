#pragma once

#include <cstdint>
#include <tuple>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

class AtomObject;
class BondObject;

class GroupPotentialEntry
{
    std::unordered_set<uint64_t> m_group_key_set;
    std::unordered_map<uint64_t, std::vector<AtomObject *>> m_group_atom_list_map;
    std::unordered_map<uint64_t, std::vector<BondObject *>> m_group_bond_list_map;
    std::unordered_map<uint64_t, std::tuple<double, double>> m_gaus_estimate_mean_map;
    std::unordered_map<uint64_t, std::tuple<double, double>> m_gaus_estimate_mdpde_map;
    std::unordered_map<uint64_t, std::tuple<double, double>> m_gaus_estimate_prior_map;
    std::unordered_map<uint64_t, std::tuple<double, double>> m_gaus_variance_prior_map;

public:
    GroupPotentialEntry(void);
    ~GroupPotentialEntry();
    void InsertGroupKey(uint64_t group_key);
    void ReserveAtomObjectPtrList(uint64_t group_key, int size);
    void ReserveBondObjectPtrList(uint64_t group_key, int size);
    void AddAtomObjectPtr(uint64_t group_key, AtomObject * atom_object);
    void AddBondObjectPtr(uint64_t group_key, BondObject * bond_object);
    void AddGausEstimateMean(uint64_t group_key, double v0, double v1);
    void AddGausEstimateMDPDE(uint64_t group_key, double v0, double v1);
    void AddGausEstimatePrior(uint64_t group_key, double v0, double v1);
    void AddGausVariancePrior(uint64_t group_key, double v0, double v1);
    int GetAtomObjectPtrListSize(uint64_t group_key) const;
    int GetBondObjectPtrListSize(uint64_t group_key) const;
    const std::vector<AtomObject *> & GetAtomObjectPtrList(uint64_t group_key) const;
    const std::vector<BondObject *> & GetBondObjectPtrList(uint64_t group_key) const;
    std::tuple<double, double> GetGausEstimateMean(uint64_t group_key) const;
    std::tuple<double, double> GetGausEstimateMDPDE(uint64_t group_key) const;
    std::tuple<double, double> GetGausEstimatePrior(uint64_t group_key) const;
    std::tuple<double, double> GetGausVariancePrior(uint64_t group_key) const;
    double GetGausEstimateMean(uint64_t group_key, int par_id) const;
    double GetGausEstimateMDPDE(uint64_t group_key, int par_id) const;
    double GetGausEstimatePrior(uint64_t group_key, int par_id) const;
    double GetGausVariancePrior(uint64_t group_key, int par_id) const;
    const std::unordered_set<uint64_t> & GetGroupKeySet(void) const;

private:
    double GetGausEstimateFromTuple(const std::tuple<double, double> & estimate, int par_id) const;
    double GetGausVarianceFromTuple(const std::tuple<double, double> & estimate, const std::tuple<double, double> & variance, int par_id) const;
    double CalculateIntensityEstimate(const std::tuple<double, double> & estimate) const;
    double CalculateIntensityVariance(const std::tuple<double, double> & estimate, const std::tuple<double, double> & variance) const;

};
