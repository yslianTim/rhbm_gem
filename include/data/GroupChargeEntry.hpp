#pragma once

#include <cstdint>
#include <tuple>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

class AtomObject;
class GroupChargeEntry
{
    std::unordered_set<uint64_t> m_group_key_set;
    std::unordered_map<uint64_t, std::vector<AtomObject *>> m_group_atom_list_map;
    std::unordered_map<uint64_t, std::tuple<double, double, double>> m_model_estimate_mean_map;
    std::unordered_map<uint64_t, std::tuple<double, double, double>> m_model_estimate_mdpde_map;
    std::unordered_map<uint64_t, std::tuple<double, double, double>> m_model_estimate_prior_map;
    std::unordered_map<uint64_t, std::tuple<double, double, double>> m_model_variance_prior_map;

public:
    GroupChargeEntry(void);
    ~GroupChargeEntry();
    void InsertGroupKey(uint64_t group_key);
    void ReserveAtomObjectPtrList(uint64_t group_key, int size);
    void AddAtomObjectPtr(uint64_t group_key, AtomObject * atom_object);
    void AddModelEstimateMean(uint64_t group_key, double v0, double v1, double v2);
    void AddModelEstimateMDPDE(uint64_t group_key, double v0, double v1, double v2);
    void AddModelEstimatePrior(uint64_t group_key, double v0, double v1, double v2);
    void AddModelVariancePrior(uint64_t group_key, double v0, double v1, double v2);
    int GetAtomObjectPtrListSize(uint64_t group_key) const;
    const std::vector<AtomObject *> & GetAtomObjectPtrList(uint64_t group_key) const;
    std::tuple<double, double, double> GetModelEstimateMean(uint64_t group_key) const;
    std::tuple<double, double, double> GetModelEstimateMDPDE(uint64_t group_key) const;
    std::tuple<double, double, double> GetModelEstimatePrior(uint64_t group_key) const;
    std::tuple<double, double, double> GetModelVariancePrior(uint64_t group_key) const;
    double GetModelEstimateMean(uint64_t group_key, int par_id) const;
    double GetModelEstimateMDPDE(uint64_t group_key, int par_id) const;
    double GetModelEstimatePrior(uint64_t group_key, int par_id) const;
    double GetModelVariancePrior(uint64_t group_key, int par_id) const;
    const std::unordered_set<uint64_t> & GetGroupKeySet(void) const;

private:
    double GetModelEstimateFromTuple(const std::tuple<double, double, double> & estimate, int par_id) const;
    double GetModelVarianceFromTuple(const std::tuple<double, double, double> & variance, int par_id) const;

};
