#pragma once

#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>

namespace rhbm_gem {

class AtomObject;
class BondObject;

template <typename MemberT>
class GroupPotentialEntry
{
    struct GroupPotentialBucket
    {
        std::vector<MemberT *> members;
        GaussianModel3D mean{ 0.0, 0.0 };
        GaussianModel3D mdpde{ 0.0, 0.0 };
        GaussianModel3DWithUncertainty prior{
            GaussianModel3D{ 0.0, 0.0 },
            GaussianModel3DUncertainty{}
        };
        double alpha_g{ 0.0 };
    };

    std::unordered_map<GroupKey, GroupPotentialBucket> m_group_map;

public:
    GroupPotentialEntry() = default;
    ~GroupPotentialEntry() = default;

    bool HasGroup(GroupKey group_key) const
    {
        return m_group_map.find(group_key) != m_group_map.end();
    }

    std::vector<GroupKey> CollectGroupKeys() const
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

    size_t GroupCount() const
    {
        return m_group_map.size();
    }

    void AddMember(GroupKey group_key, MemberT & member)
    {
        EnsureGroup(group_key).members.emplace_back(&member);
    }

    void ReserveMembers(GroupKey group_key, size_t member_count)
    {
        EnsureGroup(group_key).members.reserve(member_count);
    }

    const std::vector<MemberT *> & GetMembers(GroupKey group_key) const
    {
        return RequireGroup(group_key).members;
    }

    size_t GetMemberCount(GroupKey group_key) const
    {
        return GetMembers(group_key).size();
    }

    void SetGaussianResult(GroupKey group_key, const GroupGaussianResult & result)
    {
        auto & group{ EnsureGroup(group_key) };
        group.mean = result.mean;
        group.mdpde = result.mdpde;
        group.prior = result.prior;
        group.alpha_g = result.alpha_g;
    }

    void SetAlphaG(GroupKey group_key, double alpha_g)
    {
        EnsureGroup(group_key).alpha_g = alpha_g;
    }

    const GaussianModel3D & GetMean(GroupKey group_key) const
    {
        return RequireGroup(group_key).mean;
    }

    const GaussianModel3D & GetMDPDE(GroupKey group_key) const
    {
        return RequireGroup(group_key).mdpde;
    }

    const GaussianModel3D & GetPrior(GroupKey group_key) const
    {
        return RequireGroup(group_key).prior.GetModel();
    }

    const GaussianModel3DUncertainty & GetPriorStandardDeviation(GroupKey group_key) const
    {
        return RequireGroup(group_key).prior.GetStandardDeviationModel();
    }

    GaussianModel3DWithUncertainty GetPriorWithUncertainty(GroupKey group_key) const
    {
        return RequireGroup(group_key).prior;
    }

    double GetAlphaG(GroupKey group_key) const
    {
        return RequireGroup(group_key).alpha_g;
    }

private:
    GroupPotentialBucket & EnsureGroup(GroupKey group_key)
    {
        return m_group_map[group_key];
    }

    GroupPotentialBucket & RequireGroup(GroupKey group_key)
    {
        const auto iter{ m_group_map.find(group_key) };
        if (iter == m_group_map.end())
        {
            throw std::runtime_error("Group key is not available.");
        }
        return iter->second;
    }

    const GroupPotentialBucket & RequireGroup(GroupKey group_key) const
    {
        const auto iter{ m_group_map.find(group_key) };
        if (iter == m_group_map.end())
        {
            throw std::runtime_error("Group key is not available.");
        }
        return iter->second;
    }
};

using AtomGroupPotentialEntry = GroupPotentialEntry<AtomObject>;
using BondGroupPotentialEntry = GroupPotentialEntry<BondObject>;

} // namespace rhbm_gem
