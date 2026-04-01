#include <rhbm_gem/data/object/GroupPotentialEntry.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

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

int GroupPotentialEntry::GetAtomMemberCount(GroupKey group_key) const
{
    return static_cast<int>(GetGroup(group_key).atom_members.size());
}

int GroupPotentialEntry::GetBondMemberCount(GroupKey group_key) const
{
    return static_cast<int>(GetGroup(group_key).bond_members.size());
}

const std::vector<AtomObject *> & GroupPotentialEntry::GetAtomMembers(GroupKey group_key) const
{
    return GetGroup(group_key).atom_members;
}

const std::vector<BondObject *> & GroupPotentialEntry::GetBondMembers(GroupKey group_key) const
{
    return GetGroup(group_key).bond_members;
}

const GaussianEstimate & GroupPotentialEntry::GetMeanEstimate(GroupKey group_key) const
{
    return GetGroup(group_key).mean;
}

const GaussianEstimate & GroupPotentialEntry::GetMDPDEEstimate(GroupKey group_key) const
{
    return GetGroup(group_key).mdpde;
}

const GaussianEstimate & GroupPotentialEntry::GetPriorEstimate(GroupKey group_key) const
{
    return GetGroup(group_key).prior;
}

const GaussianEstimate & GroupPotentialEntry::GetPriorVariance(GroupKey group_key) const
{
    return GetGroup(group_key).prior_variance;
}

double GroupPotentialEntry::GetAlphaG(GroupKey group_key) const
{
    if (!HasGroup(group_key))
    {
        Logger::Log(
            LogLevel::Warning,
            "GroupPotentialEntry::GetAlphaG() - Group key not found: "
                + std::to_string(group_key));
        return 0.0;
    }
    return GetGroup(group_key).alpha_g;
}

double GroupPotentialEntry::GetMuEstimateMean(GroupKey group_key, int par_id) const
{
    return GetMeanEstimate(group_key).ToBeta()(par_id);
}

double GroupPotentialEntry::GetMuEstimateMDPDE(GroupKey group_key, int par_id) const
{
    return GetMDPDEEstimate(group_key).ToBeta()(par_id);
}

double GroupPotentialEntry::GetMuEstimatePrior(GroupKey group_key, int par_id) const
{
    return GetPriorEstimate(group_key).ToBeta()(par_id);
}

double GroupPotentialEntry::GetGausEstimateMean(GroupKey group_key, int par_id) const
{
    return GetGausEstimate(GetMeanEstimate(group_key), par_id);
}

double GroupPotentialEntry::GetGausEstimateMDPDE(GroupKey group_key, int par_id) const
{
    return GetGausEstimate(GetMDPDEEstimate(group_key), par_id);
}

double GroupPotentialEntry::GetGausEstimatePrior(GroupKey group_key, int par_id) const
{
    return GetGausEstimate(GetPriorEstimate(group_key), par_id);
}

double GroupPotentialEntry::GetGausVariancePrior(GroupKey group_key, int par_id) const
{
    return GetGausVariance(
        GetPriorEstimate(group_key),
        GetPriorVariance(group_key),
        par_id);
}

double GroupPotentialEntry::GetGausEstimate(
    const GaussianEstimate & estimate,
    int par_id) const
{
    return estimate.GetParameter(par_id);
}

double GroupPotentialEntry::GetGausVariance(
    const GaussianEstimate & estimate,
    const GaussianEstimate & variance,
    int par_id) const
{
    GaussianPosterior posterior;
    posterior.estimate = estimate;
    posterior.variance = variance;
    switch (par_id)
    {
    case 0:
        return variance.amplitude;
    case 1:
        return variance.width;
    case 2:
        return posterior.IntensityVariance();
    default:
        Logger::Log(LogLevel::Error, "Invalid parameter index: " + std::to_string(par_id));
        return 0.0;
    }
}

} // namespace rhbm_gem
