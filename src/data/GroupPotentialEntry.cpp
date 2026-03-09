#include <rhbm_gem/data/GroupPotentialEntry.hpp>
#include <rhbm_gem/utils/GausLinearTransformHelper.hpp>
#include <rhbm_gem/data/AtomObject.hpp>
#include <rhbm_gem/data/BondObject.hpp>
#include <rhbm_gem/utils/Constants.hpp>
#include <rhbm_gem/utils/Logger.hpp>

#include <cmath>

namespace rhbm_gem {

GroupPotentialEntry::GroupPotentialEntry()
{

}

GroupPotentialEntry::~GroupPotentialEntry()
{

}

void GroupPotentialEntry::InsertGroupKey(GroupKey group_key)
{
    m_group_key_set.insert(group_key);
}

void GroupPotentialEntry::ReserveAtomObjectPtrList(GroupKey group_key, int size)
{
    m_group_atom_list_map[group_key].reserve(static_cast<size_t>(size));
}

void GroupPotentialEntry::ReserveBondObjectPtrList(GroupKey group_key, int size)
{
    m_group_bond_list_map[group_key].reserve(static_cast<size_t>(size));
}

void GroupPotentialEntry::AddAtomObjectPtr(GroupKey group_key, AtomObject * atom_object)
{
    m_group_atom_list_map[group_key].emplace_back(atom_object);
}

void GroupPotentialEntry::AddBondObjectPtr(GroupKey group_key, BondObject * bond_object)
{
    m_group_bond_list_map[group_key].emplace_back(bond_object);
}

void GroupPotentialEntry::AddGausEstimateMean(GroupKey group_key, double v0, double v1)
{
    m_gaus_estimate_mean_map[group_key] = std::make_tuple(v0, v1);
}

void GroupPotentialEntry::AddGausEstimateMDPDE(GroupKey group_key, double v0, double v1)
{
    m_gaus_estimate_mdpde_map[group_key] = std::make_tuple(v0, v1);
}

void GroupPotentialEntry::AddGausEstimatePrior(GroupKey group_key, double v0, double v1)
{
    m_gaus_estimate_prior_map[group_key] = std::make_tuple(v0, v1);
}

void GroupPotentialEntry::AddGausVariancePrior(GroupKey group_key, double v0, double v1)
{
    m_gaus_variance_prior_map[group_key] = std::make_tuple(v0, v1);
}

void GroupPotentialEntry::AddAlphaG(GroupKey group_key, double value)
{
    m_alpha_g_map[group_key] = value;
}

int GroupPotentialEntry::GetAtomObjectPtrListSize(GroupKey group_key) const
{
    return static_cast<int>(m_group_atom_list_map.at(group_key).size());
}

int GroupPotentialEntry::GetBondObjectPtrListSize(GroupKey group_key) const
{
    return static_cast<int>(m_group_bond_list_map.at(group_key).size());
}

const std::vector<AtomObject *> & GroupPotentialEntry::GetAtomObjectPtrList(GroupKey group_key) const
{
    return m_group_atom_list_map.at(group_key);
}

const std::vector<BondObject *> & GroupPotentialEntry::GetBondObjectPtrList(GroupKey group_key) const
{
    return m_group_bond_list_map.at(group_key);
}

std::tuple<double, double> GroupPotentialEntry::GetGausEstimateMean(GroupKey group_key) const
{
    return m_gaus_estimate_mean_map.at(group_key);
}

std::tuple<double, double> GroupPotentialEntry::GetGausEstimateMDPDE(GroupKey group_key) const
{
    return m_gaus_estimate_mdpde_map.at(group_key);
}

std::tuple<double, double> GroupPotentialEntry::GetGausEstimatePrior(GroupKey group_key) const
{
    return m_gaus_estimate_prior_map.at(group_key);
}

std::tuple<double, double> GroupPotentialEntry::GetGausVariancePrior(GroupKey group_key) const
{
    return m_gaus_variance_prior_map.at(group_key);
}

double GroupPotentialEntry::GetAlphaG(GroupKey group_key) const
{
    if (m_alpha_g_map.find(group_key) == m_alpha_g_map.end())
    {
        Logger::Log(LogLevel::Warning,
            "GroupPotentialEntry::GetAlphaG() - Group key not found: " + std::to_string(group_key));
        return 0.0;
    }
    return m_alpha_g_map.at(group_key);
}

const std::unordered_set<GroupKey> & GroupPotentialEntry::GetGroupKeySet() const
{
    return m_group_key_set;
}

double GroupPotentialEntry::GetMuEstimateMean(GroupKey group_key, int par_id) const
{
    auto amplitude{ GetGausEstimateMean(group_key, 0) };
    auto width{ GetGausEstimateMean(group_key, 1) };
    Eigen::VectorXd coefficient_vector{
        GausLinearTransformHelper::BuildLinearModelCoefficentVector(amplitude, width)
    };
    return coefficient_vector(par_id);
}

double GroupPotentialEntry::GetMuEstimateMDPDE(GroupKey group_key, int par_id) const
{
    auto amplitude{ GetGausEstimateMDPDE(group_key, 0) };
    auto width{ GetGausEstimateMDPDE(group_key, 1) };
    Eigen::VectorXd coefficient_vector{
        GausLinearTransformHelper::BuildLinearModelCoefficentVector(amplitude, width)
    };
    return coefficient_vector(par_id);
}

double GroupPotentialEntry::GetMuEstimatePrior(GroupKey group_key, int par_id) const
{
    auto amplitude{ GetGausEstimatePrior(group_key, 0) };
    auto width{ GetGausEstimatePrior(group_key, 1) };
    Eigen::VectorXd coefficient_vector{
        GausLinearTransformHelper::BuildLinearModelCoefficentVector(amplitude, width)
    };
    return coefficient_vector(par_id);
}

double GroupPotentialEntry::GetGausEstimateMean(GroupKey group_key, int par_id) const
{
    return GetGausEstimateFromTuple(GetGausEstimateMean(group_key), par_id);
}

double GroupPotentialEntry::GetGausEstimateMDPDE(GroupKey group_key, int par_id) const
{
    return GetGausEstimateFromTuple(GetGausEstimateMDPDE(group_key), par_id);
}

double GroupPotentialEntry::GetGausEstimatePrior(GroupKey group_key, int par_id) const
{
    return GetGausEstimateFromTuple(GetGausEstimatePrior(group_key), par_id);
}

double GroupPotentialEntry::GetGausVariancePrior(GroupKey group_key, int par_id) const
{
    return GetGausVarianceFromTuple(
        GetGausEstimatePrior(group_key), GetGausVariancePrior(group_key), par_id);
}

double GroupPotentialEntry::GetGausEstimateFromTuple(
    const std::tuple<double, double> & estimate, int par_id) const
{
    switch (par_id)
    {
        case 0: return std::get<0>(estimate);
        case 1: return std::get<1>(estimate);
        case 2: return CalculateIntensityEstimate(estimate);
        default:
            Logger::Log(LogLevel::Error, "Invalid parameter index: " + std::to_string(par_id));
            return 0.0;
    }
}

double GroupPotentialEntry::GetGausVarianceFromTuple(
    const std::tuple<double, double> & estimate,
    const std::tuple<double, double> & variance, int par_id) const
{
    switch (par_id)
    {
        case 0: return std::get<0>(variance);
        case 1: return std::get<1>(variance);
        case 2: return CalculateIntensityVariance(estimate, variance);
        default:
            Logger::Log(LogLevel::Error, "Invalid parameter index: " + std::to_string(par_id));
            return 0.0;
    }
}

double GroupPotentialEntry::CalculateIntensityEstimate(
    const std::tuple<double, double> & estimate) const
{
    auto amplitude{ std::get<0>(estimate) };
    auto width{ std::get<1>(estimate) };
    if (width == 0.0) return 0.0;
    return amplitude * std::pow(Constants::two_pi * width * width, -1.5);
}

double GroupPotentialEntry::CalculateIntensityVariance(
    const std::tuple<double, double> & estimate,
    const std::tuple<double, double> & variance) const
{
    auto amplitude{ std::get<0>(estimate) };
    auto width{ std::get<1>(estimate) };
    auto sigma_amplitude{ std::get<0>(variance) };
    auto sigma_width{ std::get<1>(variance) };
    return std::sqrt(
        std::pow(std::pow(Constants::two_pi * width * width, -1.5) * sigma_amplitude, 2) +
        std::pow(-3.0 * amplitude * std::pow(Constants::two_pi, -1.5) * std::pow(width, -4) * sigma_width, 2)
    );
}

} // namespace rhbm_gem
