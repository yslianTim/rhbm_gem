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

    using GroupMap = std::unordered_map<GroupKey, GroupPotentialBucket>;

    GroupMap m_group_map_1st;
    GroupMap m_group_map_2nd;
    GroupMap m_group_map_3rd;

public:
    GroupPotentialEntry() = default;
    ~GroupPotentialEntry() = default;

    bool HasGroup(FittingStage stage, GroupKey group_key) const
    {
        const auto & group_map{ GetGroupMap(stage) };
        return group_map.find(group_key) != group_map.end();
    }

    std::vector<GroupKey> CollectGroupKeys(FittingStage stage) const
    {
        const auto & group_map{ GetGroupMap(stage) };
        std::vector<GroupKey> group_keys;
        group_keys.reserve(group_map.size());
        for (const auto & [group_key, bucket] : group_map)
        {
            (void)bucket;
            group_keys.emplace_back(group_key);
        }
        return group_keys;
    }

    size_t GroupCount(FittingStage stage) const
    {
        return GetGroupMap(stage).size();
    }

    void AddMember(FittingStage stage, GroupKey group_key, MemberT & member)
    {
        EnsureGroup(stage, group_key).members.emplace_back(&member);
    }

    void ReserveMembers(
        FittingStage stage,
        GroupKey group_key,
        size_t member_count)
    {
        EnsureGroup(stage, group_key).members.reserve(member_count);
    }

    const std::vector<MemberT *> & GetMembers(
        FittingStage stage,
        GroupKey group_key) const
    {
        return RequireGroup(stage, group_key).members;
    }

    size_t GetMemberCount(FittingStage stage, GroupKey group_key) const
    {
        return GetMembers(stage, group_key).size();
    }

    void SetGaussianResult(
        FittingStage stage,
        GroupKey group_key,
        const GroupGaussianResult & result)
    {
        auto & group{ EnsureGroup(stage, group_key) };
        group.mean = result.mean;
        group.mdpde = result.mdpde;
        group.prior = result.prior;
        group.alpha_g = result.alpha_g;
    }

    void SetAlphaG(FittingStage stage, GroupKey group_key, double alpha_g)
    {
        EnsureGroup(stage, group_key).alpha_g = alpha_g;
    }

    const GaussianModel3D & GetMean(FittingStage stage, GroupKey group_key) const
    {
        return RequireGroup(stage, group_key).mean;
    }

    const GaussianModel3D & GetMDPDE(FittingStage stage, GroupKey group_key) const
    {
        return RequireGroup(stage, group_key).mdpde;
    }

    const GaussianModel3D & GetPrior(FittingStage stage, GroupKey group_key) const
    {
        return RequireGroup(stage, group_key).prior.GetModel();
    }

    const GaussianModel3DUncertainty & GetPriorStandardDeviation(
        FittingStage stage,
        GroupKey group_key) const
    {
        return RequireGroup(stage, group_key).prior.GetStandardDeviationModel();
    }

    GaussianModel3DWithUncertainty GetPriorWithUncertainty(
        FittingStage stage,
        GroupKey group_key) const
    {
        return RequireGroup(stage, group_key).prior;
    }

    double GetAlphaG(FittingStage stage, GroupKey group_key) const
    {
        return RequireGroup(stage, group_key).alpha_g;
    }

    void CopyStage(FittingStage source_stage, FittingStage destination_stage)
    {
        GetGroupMap(destination_stage) = GetGroupMap(source_stage);
    }

private:
    GroupMap & GetGroupMap(FittingStage stage)
    {
        switch (stage)
        {
            case FittingStage::First: return m_group_map_1st;
            case FittingStage::Second: return m_group_map_2nd;
            case FittingStage::Third: return m_group_map_3rd;
        }
        throw std::invalid_argument("Unknown local fitting stage.");
    }

    const GroupMap & GetGroupMap(FittingStage stage) const
    {
        switch (stage)
        {
            case FittingStage::First: return m_group_map_1st;
            case FittingStage::Second: return m_group_map_2nd;
            case FittingStage::Third: return m_group_map_3rd;
        }
        throw std::invalid_argument("Unknown local fitting stage.");
    }

    GroupPotentialBucket & EnsureGroup(FittingStage stage, GroupKey group_key)
    {
        return GetGroupMap(stage)[group_key];
    }

    GroupPotentialBucket & RequireGroup(FittingStage stage, GroupKey group_key)
    {
        auto & group_map{ GetGroupMap(stage) };
        const auto iter{ group_map.find(group_key) };
        if (iter == group_map.end())
        {
            throw std::runtime_error("Group key is not available.");
        }
        return iter->second;
    }

    const GroupPotentialBucket & RequireGroup(
        FittingStage stage,
        GroupKey group_key) const
    {
        const auto & group_map{ GetGroupMap(stage) };
        const auto iter{ group_map.find(group_key) };
        if (iter == group_map.end())
        {
            throw std::runtime_error("Group key is not available.");
        }
        return iter->second;
    }
};

using AtomGroupPotentialEntry = GroupPotentialEntry<AtomObject>;
using BondGroupPotentialEntry = GroupPotentialEntry<BondObject>;

} // namespace rhbm_gem
