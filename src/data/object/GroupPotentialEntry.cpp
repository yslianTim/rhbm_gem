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

} // namespace rhbm_gem
