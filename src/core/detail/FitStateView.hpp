#pragma once

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <vector>

#include <rhbm_gem/core/GaussianEstimator.hpp>

namespace rhbm_gem::core::detail {

using FitState = std::vector<LocalGaussianResult>;
using PolishProvenance = std::vector<char>;
using ClusterKey = std::vector<std::size_t>;

inline const GaussianModel3D & GetFitModel(const FitState & state, std::size_t atom_index)
{
    return state.at(atom_index).mdpde.GetModel();
}

struct FitStatePatch
{
    ClusterKey atom_index_list{};
    std::vector<GaussianModel3DWithUncertainty> mdpde_list{};

    const GaussianModel3DWithUncertainty * Find(std::size_t atom_index) const
    {
        const auto iter{
            std::lower_bound(
                atom_index_list.begin(),
                atom_index_list.end(),
                atom_index)
        };
        if (iter == atom_index_list.end() || *iter != atom_index)
        {
            return nullptr;
        }
        return &mdpde_list.at(static_cast<std::size_t>(
            std::distance(atom_index_list.begin(), iter)));
    }

    void ApplyTo(FitState & state) const
    {
        if (atom_index_list.size() != mdpde_list.size())
        {
            throw std::invalid_argument(
                "Local fitting state patch sizes are inconsistent.");
        }
        for (std::size_t i = 0; i < atom_index_list.size(); i++)
        {
            state.at(atom_index_list.at(i)).mdpde = mdpde_list.at(i);
        }
    }
};

class FitStateView
{
public:
    explicit FitStateView(const FitState & base_state) : m_base_state{ &base_state }
    {
    }

    FitStateView(const FitState & base_state, const FitStatePatch & patch)
        : m_base_state{ &base_state }, m_patch{ &patch }
    {
    }

    const GaussianModel3DWithUncertainty & GetMdpde(std::size_t atom_index) const
    {
        if (m_patch != nullptr)
        {
            const auto * value{ m_patch->Find(atom_index) };
            if (value != nullptr) return *value;
        }
        return m_base_state->at(atom_index).mdpde;
    }

    const GaussianModel3D & GetModel(std::size_t atom_index) const
    {
        return GetMdpde(atom_index).GetModel();
    }

    std::size_t GetSize() const { return m_base_state->size(); }

private:
    const FitState * m_base_state{ nullptr };
    const FitStatePatch * m_patch{ nullptr };
};

inline const GaussianModel3D & GetFitModel(const FitStateView & state, std::size_t atom_index)
{
    return state.GetModel(atom_index);
}

} // namespace rhbm_gem::core::detail
