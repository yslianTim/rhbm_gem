#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

#include <rhbm_gem/data/object/GaussianStatistics.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>

namespace rhbm_gem {

class AtomObject;
class BondObject;

struct GroupPotentialBucket
{
    std::vector<AtomObject *> atom_members;
    std::vector<BondObject *> bond_members;
    GaussianEstimate mean{};
    GaussianEstimate mdpde{};
    GaussianEstimate prior{};
    GaussianEstimate prior_variance{};
    double alpha_g{ 0.0 };
};

class GroupPotentialEntry
{
    std::unordered_map<GroupKey, GroupPotentialBucket> m_group_map;

public:
    GroupPotentialEntry();
    ~GroupPotentialEntry();
    GroupPotentialBucket & EnsureGroup(GroupKey group_key);
    bool HasGroup(GroupKey group_key) const;
    const GroupPotentialBucket & GetGroup(GroupKey group_key) const;
    const std::unordered_map<GroupKey, GroupPotentialBucket> & GetGroups() const;
    std::vector<GroupKey> GetGroupKeys() const;
    void AddAtomMember(GroupKey group_key, AtomObject * atom_object);
    void AddBondMember(GroupKey group_key, BondObject * bond_object);
    void SetMeanEstimate(GroupKey group_key, const GaussianEstimate & estimate);
    void SetMDPDEEstimate(GroupKey group_key, const GaussianEstimate & estimate);
    void SetPriorEstimate(GroupKey group_key, const GaussianEstimate & estimate);
    void SetPriorVariance(GroupKey group_key, const GaussianEstimate & variance);
    void SetAlphaG(GroupKey group_key, double value);
    int GetAtomMemberCount(GroupKey group_key) const;
    int GetBondMemberCount(GroupKey group_key) const;
    double GetAlphaG(GroupKey group_key) const;
    const std::vector<AtomObject *> & GetAtomMembers(GroupKey group_key) const;
    const std::vector<BondObject *> & GetBondMembers(GroupKey group_key) const;
    const GaussianEstimate & GetMeanEstimate(GroupKey group_key) const;
    const GaussianEstimate & GetMDPDEEstimate(GroupKey group_key) const;
    const GaussianEstimate & GetPriorEstimate(GroupKey group_key) const;
    const GaussianEstimate & GetPriorVariance(GroupKey group_key) const;
    double GetMuEstimateMean(GroupKey group_key, int par_id) const;
    double GetMuEstimateMDPDE(GroupKey group_key, int par_id) const;
    double GetMuEstimatePrior(GroupKey group_key, int par_id) const;
    double GetGausEstimateMean(GroupKey group_key, int par_id) const;
    double GetGausEstimateMDPDE(GroupKey group_key, int par_id) const;
    double GetGausEstimatePrior(GroupKey group_key, int par_id) const;
    double GetGausVariancePrior(GroupKey group_key, int par_id) const;

private:
    double GetGausEstimate(const GaussianEstimate & estimate, int par_id) const;
    double GetGausVariance(
        const GaussianEstimate & estimate,
        const GaussianEstimate & variance,
        int par_id) const;
};

} // namespace rhbm_gem
