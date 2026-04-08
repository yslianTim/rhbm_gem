#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "data/detail/GaussianStatistics.hpp"
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>

namespace rhbm_gem {

class AtomObject;
class BondObject;

class GroupPotentialEntry
{
    enum class EntryKind
    {
        UNKNOWN,
        ATOM,
        BOND,
    };

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

    EntryKind m_entry_kind;
    std::unordered_map<GroupKey, GroupPotentialBucket> m_group_map;

public:
    GroupPotentialEntry();
    ~GroupPotentialEntry();
    bool HasGroup(GroupKey group_key) const;
    std::vector<GroupKey> CollectGroupKeys() const;
    size_t GroupCount() const;

    void AddAtomMember(GroupKey group_key, AtomObject & atom_object);
    void AddBondMember(GroupKey group_key, BondObject & bond_object);
    void ReserveAtomMembers(GroupKey group_key, size_t member_count);
    void ReserveBondMembers(GroupKey group_key, size_t member_count);
    const std::vector<AtomObject *> & GetAtomMembers(GroupKey group_key) const;
    const std::vector<BondObject *> & GetBondMembers(GroupKey group_key) const;
    size_t GetAtomMemberCount(GroupKey group_key) const;
    size_t GetBondMemberCount(GroupKey group_key) const;

    void SetGroupStatistics(
        GroupKey group_key,
        const GaussianEstimate & mean,
        const GaussianEstimate & mdpde,
        const GaussianEstimate & prior,
        const GaussianEstimate & prior_variance,
        double alpha_g);
    void SetAlphaG(GroupKey group_key, double alpha_g);
    const GaussianEstimate & GetMean(GroupKey group_key) const;
    const GaussianEstimate & GetMDPDE(GroupKey group_key) const;
    const GaussianEstimate & GetPrior(GroupKey group_key) const;
    const GaussianEstimate & GetPriorVariance(GroupKey group_key) const;
    GaussianPosterior BuildPriorPosterior(GroupKey group_key) const;
    double GetAlphaG(GroupKey group_key) const;

    void MarkAsAtomEntry();
    void MarkAsBondEntry();

private:
    void EnsureKind(EntryKind entry_kind);
    GroupPotentialBucket & EnsureGroup(GroupKey group_key);
    GroupPotentialBucket & RequireGroup(GroupKey group_key);
    const GroupPotentialBucket & RequireGroup(GroupKey group_key) const;
};

} // namespace rhbm_gem
