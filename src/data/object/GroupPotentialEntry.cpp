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

GroupPotentialBucket * GroupPotentialEntry::FindGroup(GroupKey group_key)
{
    const auto iter{ m_group_map.find(group_key) };
    return (iter == m_group_map.end()) ? nullptr : &iter->second;
}

const GroupPotentialBucket * GroupPotentialEntry::FindGroup(GroupKey group_key) const
{
    const auto iter{ m_group_map.find(group_key) };
    return (iter == m_group_map.end()) ? nullptr : &iter->second;
}

const std::unordered_map<GroupKey, GroupPotentialBucket> & GroupPotentialEntry::Entries() const
{
    return m_group_map;
}

} // namespace rhbm_gem
