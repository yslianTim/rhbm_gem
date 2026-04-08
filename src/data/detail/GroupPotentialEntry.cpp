#include "data/detail/GroupPotentialEntry.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>

#include <stdexcept>

namespace rhbm_gem {

GroupPotentialEntry::GroupPotentialEntry()
    : m_entry_kind{ EntryKind::UNKNOWN }
{
}

GroupPotentialEntry::~GroupPotentialEntry()
{
}

bool GroupPotentialEntry::HasGroup(GroupKey group_key) const
{
    return m_group_map.find(group_key) != m_group_map.end();
}

std::vector<GroupKey> GroupPotentialEntry::CollectGroupKeys() const
{
    std::vector<GroupKey> group_keys;
    group_keys.reserve(m_group_map.size());
    for (const auto & [group_key, bucket] : m_group_map)
    {
        (void)bucket;
        group_keys.emplace_back(group_key);
    }
    return group_keys;
}

size_t GroupPotentialEntry::GroupCount() const
{
    return m_group_map.size();
}

void GroupPotentialEntry::AddAtomMember(GroupKey group_key, AtomObject & atom_object)
{
    EnsureKind(EntryKind::ATOM);
    EnsureGroup(group_key).atom_members.emplace_back(&atom_object);
}

void GroupPotentialEntry::AddBondMember(GroupKey group_key, BondObject & bond_object)
{
    EnsureKind(EntryKind::BOND);
    EnsureGroup(group_key).bond_members.emplace_back(&bond_object);
}

void GroupPotentialEntry::ReserveAtomMembers(GroupKey group_key, size_t member_count)
{
    EnsureKind(EntryKind::ATOM);
    EnsureGroup(group_key).atom_members.reserve(member_count);
}

void GroupPotentialEntry::ReserveBondMembers(GroupKey group_key, size_t member_count)
{
    EnsureKind(EntryKind::BOND);
    EnsureGroup(group_key).bond_members.reserve(member_count);
}

const std::vector<AtomObject *> & GroupPotentialEntry::GetAtomMembers(GroupKey group_key) const
{
    if (m_entry_kind != EntryKind::ATOM)
    {
        throw std::runtime_error("Atom members are not available for this group entry.");
    }
    return RequireGroup(group_key).atom_members;
}

const std::vector<BondObject *> & GroupPotentialEntry::GetBondMembers(GroupKey group_key) const
{
    if (m_entry_kind != EntryKind::BOND)
    {
        throw std::runtime_error("Bond members are not available for this group entry.");
    }
    return RequireGroup(group_key).bond_members;
}

size_t GroupPotentialEntry::GetAtomMemberCount(GroupKey group_key) const
{
    return GetAtomMembers(group_key).size();
}

size_t GroupPotentialEntry::GetBondMemberCount(GroupKey group_key) const
{
    return GetBondMembers(group_key).size();
}

void GroupPotentialEntry::SetGroupStatistics(
    GroupKey group_key,
    const GaussianEstimate & mean,
    const GaussianEstimate & mdpde,
    const GaussianEstimate & prior,
    const GaussianEstimate & prior_variance,
    double alpha_g)
{
    auto & group{ EnsureGroup(group_key) };
    group.mean = mean;
    group.mdpde = mdpde;
    group.prior = prior;
    group.prior_variance = prior_variance;
    group.alpha_g = alpha_g;
}

void GroupPotentialEntry::SetAlphaG(GroupKey group_key, double alpha_g)
{
    EnsureGroup(group_key).alpha_g = alpha_g;
}

const GaussianEstimate & GroupPotentialEntry::GetMean(GroupKey group_key) const
{
    return RequireGroup(group_key).mean;
}

const GaussianEstimate & GroupPotentialEntry::GetMDPDE(GroupKey group_key) const
{
    return RequireGroup(group_key).mdpde;
}

const GaussianEstimate & GroupPotentialEntry::GetPrior(GroupKey group_key) const
{
    return RequireGroup(group_key).prior;
}

const GaussianEstimate & GroupPotentialEntry::GetPriorVariance(GroupKey group_key) const
{
    return RequireGroup(group_key).prior_variance;
}

GaussianPosterior GroupPotentialEntry::BuildPriorPosterior(GroupKey group_key) const
{
    GaussianPosterior posterior;
    posterior.estimate = GetPrior(group_key);
    posterior.variance = GetPriorVariance(group_key);
    return posterior;
}

double GroupPotentialEntry::GetAlphaG(GroupKey group_key) const
{
    return RequireGroup(group_key).alpha_g;
}

void GroupPotentialEntry::MarkAsAtomEntry()
{
    EnsureKind(EntryKind::ATOM);
}

void GroupPotentialEntry::MarkAsBondEntry()
{
    EnsureKind(EntryKind::BOND);
}

void GroupPotentialEntry::EnsureKind(EntryKind entry_kind)
{
    if (m_entry_kind == EntryKind::UNKNOWN)
    {
        m_entry_kind = entry_kind;
        return;
    }
    if (m_entry_kind != entry_kind)
    {
        throw std::runtime_error("Group potential entry kind mismatch.");
    }
}

GroupPotentialEntry::GroupPotentialBucket & GroupPotentialEntry::EnsureGroup(GroupKey group_key)
{
    return m_group_map[group_key];
}

GroupPotentialEntry::GroupPotentialBucket & GroupPotentialEntry::RequireGroup(GroupKey group_key)
{
    const auto iter{ m_group_map.find(group_key) };
    if (iter == m_group_map.end())
    {
        throw std::runtime_error("Group key is not available.");
    }
    return iter->second;
}

const GroupPotentialEntry::GroupPotentialBucket & GroupPotentialEntry::RequireGroup(
    GroupKey group_key) const
{
    const auto iter{ m_group_map.find(group_key) };
    if (iter == m_group_map.end())
    {
        throw std::runtime_error("Group key is not available.");
    }
    return iter->second;
}

} // namespace rhbm_gem
