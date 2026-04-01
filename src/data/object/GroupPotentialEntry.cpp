#include <rhbm_gem/data/object/GroupPotentialEntry.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>

namespace rhbm_gem {

GroupPotentialEntry::GroupPotentialEntry()
{

}

GroupPotentialEntry::~GroupPotentialEntry()
{

}

GroupPotentialBucket & GroupPotentialEntry::EnsureGroup(GroupKey group_key)
{
    return m_group_map[group_key];
}

bool GroupPotentialEntry::HasGroup(GroupKey group_key) const
{
    return m_group_map.find(group_key) != m_group_map.end();
}

const GroupPotentialBucket & GroupPotentialEntry::GetGroup(GroupKey group_key) const
{
    return m_group_map.at(group_key);
}

const std::unordered_map<GroupKey, GroupPotentialBucket> & GroupPotentialEntry::GetGroups() const
{
    return m_group_map;
}

std::vector<GroupKey> GroupPotentialEntry::GetGroupKeys() const
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

void GroupPotentialEntry::AddAtomMember(GroupKey group_key, AtomObject * atom_object)
{
    EnsureGroup(group_key).atom_members.emplace_back(atom_object);
}

void GroupPotentialEntry::AddBondMember(GroupKey group_key, BondObject * bond_object)
{
    EnsureGroup(group_key).bond_members.emplace_back(bond_object);
}

void GroupPotentialEntry::SetMeanEstimate(
    GroupKey group_key,
    const GaussianEstimate & estimate)
{
    EnsureGroup(group_key).mean = estimate;
}

void GroupPotentialEntry::SetMDPDEEstimate(
    GroupKey group_key,
    const GaussianEstimate & estimate)
{
    EnsureGroup(group_key).mdpde = estimate;
}

void GroupPotentialEntry::SetPriorEstimate(
    GroupKey group_key,
    const GaussianEstimate & estimate)
{
    EnsureGroup(group_key).prior = estimate;
}

void GroupPotentialEntry::SetPriorVariance(
    GroupKey group_key,
    const GaussianEstimate & variance)
{
    EnsureGroup(group_key).prior_variance = variance;
}

void GroupPotentialEntry::SetAlphaG(GroupKey group_key, double value)
{
    EnsureGroup(group_key).alpha_g = value;
}

} // namespace rhbm_gem
