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
};

} // namespace rhbm_gem
