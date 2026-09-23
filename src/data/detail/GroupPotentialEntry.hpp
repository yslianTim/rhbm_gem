#pragma once

#include <cstddef>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>

namespace rhbm_gem {

class AtomObject;

class AtomGroupPotentialEntry
{
    struct GroupPotentialBucket
    {
        std::vector<AtomObject *> members;
        GaussianModel3D mean{ 0.0, 0.0 };
        GaussianModel3D mdpde{ 0.0, 0.0 };
        GaussianModel3DWithUncertainty prior{
            GaussianModel3D{ 0.0, 0.0 },
            GaussianModel3DUncertainty{}
        };
        double alpha_g{ 0.0 };
        std::optional<GroupParameterSummary> parameter_summary;
    };

    using GroupMap = std::unordered_map<GroupKey, GroupPotentialBucket>;

    GroupMap m_group_map{};

public:
    AtomGroupPotentialEntry() = default;
    ~AtomGroupPotentialEntry() = default;

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

    void ClearMembers()
    {
        for (auto & [key, bucket] : m_group_map) { (void)key; bucket.members.clear(); }
    }

    void AddMember(GroupKey group_key, AtomObject & member)
    {
        EnsureGroup(group_key).members.emplace_back(&member);
    }

    void ReserveMembers(
        GroupKey group_key,
        size_t member_count)
    {
        EnsureGroup(group_key).members.reserve(member_count);
    }

    const std::vector<AtomObject *> & GetMembers(GroupKey group_key) const
    {
        return RequireGroup(group_key).members;
    }

    size_t GetMemberCount(GroupKey group_key) const
    {
        return GetMembers(group_key).size();
    }

    void SetGaussianResult(
        GroupKey group_key,
        const GroupGaussianResult & result)
    {
        auto & group{ EnsureGroup(group_key) };
        group.parameter_summary.reset();
        group.mean = result.mean;
        group.mdpde = result.mdpde;
        group.prior = result.prior;
        group.alpha_g = result.alpha_g;
    }

    const std::optional<GroupParameterSummary> & GetParameterSummary(GroupKey key) const
    {
        return RequireGroup(key).parameter_summary;
    }
    void SetParameterSummary(GroupKey key, GroupParameterSummary value)
    {
        auto & group = EnsureGroup(key);
        if (value.inference) std::vector<GroupGaussianMemberResult>{}.swap(value.inference->member_results);
        group.mean = GaussianModel3D{0, 0}; group.mdpde = GaussianModel3D{0, 0};
        group.prior = GaussianModel3DWithUncertainty{GaussianModel3D{0, 0}, {}};
        group.parameter_summary = std::move(value);
    }
    void ClearResult(GroupKey key)
    {
        auto & group = RequireGroup(key);
        group.mean = GaussianModel3D{0, 0};
        group.mdpde = GaussianModel3D{0, 0};
        group.prior = GaussianModel3DWithUncertainty{GaussianModel3D{0, 0}, {}};
        group.parameter_summary.reset();
    }
    void SetAlphaG(GroupKey group_key, double alpha_g)
    {
        EnsureGroup(group_key).alpha_g = alpha_g;
    }

    const GaussianModel3D & GetMean(GroupKey group_key) const
    {
        const auto & group = RequireGroup(group_key);
        return group.parameter_summary && group.parameter_summary->descriptive_mean ?
            *group.parameter_summary->descriptive_mean : group.mean;
    }

    const GaussianModel3D & GetMDPDE(GroupKey group_key) const
    {
        const auto & group = RequireGroup(group_key);
        return group.parameter_summary && group.parameter_summary->inference ?
            group.parameter_summary->inference->mdpde : group.mdpde;
    }

    const GaussianModel3D & GetPrior(GroupKey group_key) const
    {
        const auto & group = RequireGroup(group_key);
        return group.parameter_summary && group.parameter_summary->inference ?
            group.parameter_summary->inference->prior.GetModel() : group.prior.GetModel();
    }

    const GaussianModel3DUncertainty & GetPriorStandardDeviation(GroupKey group_key) const
    {
        const auto & group = RequireGroup(group_key);
        return group.parameter_summary && group.parameter_summary->inference ?
            group.parameter_summary->inference->prior.GetStandardDeviationModel() : group.prior.GetStandardDeviationModel();
    }

    GaussianModel3DWithUncertainty GetPriorWithUncertainty(GroupKey group_key) const
    {
        const auto & group = RequireGroup(group_key);
        return group.parameter_summary && group.parameter_summary->inference ?
            group.parameter_summary->inference->prior : group.prior;
    }

    double GetAlphaG(GroupKey group_key) const
    {
        const auto & group = RequireGroup(group_key);
        return group.parameter_summary && group.parameter_summary->inference ?
            group.parameter_summary->inference->alpha_g : group.alpha_g;
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

} // namespace rhbm_gem
