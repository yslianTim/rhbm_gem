#pragma once

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>

#include "core/detail/SecondStageContext.hpp"

namespace rhbm_gem::core::detail {

using FitState = std::vector<LocalGaussianResult>;

inline const GaussianModel3D & GetFitModel(const FitState & state, std::size_t atom_index)
{
    return state.at(atom_index).mdpde.GetModel();
}

struct FitStatePatch
{
    ClusterKey atom_index_list{};
    std::vector<GaussianModel3DWithUncertainty> mdpde_list{};

    static FitStatePatch FromState(const FitState & state, ClusterKey atom_index_list)
    {
        std::ranges::sort(atom_index_list);
        atom_index_list.erase(std::ranges::unique(atom_index_list).begin(), atom_index_list.end());

        FitStatePatch patch;
        patch.atom_index_list = std::move(atom_index_list);
        patch.mdpde_list.reserve(patch.atom_index_list.size());
        for (const auto atom_index : patch.atom_index_list)
        {
            patch.mdpde_list.emplace_back(state.at(atom_index).mdpde);
        }
        return patch;
    }

    const GaussianModel3DWithUncertainty * Find(std::size_t atom_index) const
    {
        const auto iter{ std::ranges::lower_bound(atom_index_list, atom_index) };
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
    const FitState & m_base_state;
    const FitStatePatch & m_patch;

public:
    FitStateView(const FitState & base_state, const FitStatePatch & patch)
        : m_base_state{ base_state }, m_patch{ patch }
    {
    }

    const GaussianModel3DWithUncertainty & GetMdpde(std::size_t atom_index) const
    {
        const auto * value{ FindOverride(atom_index) };
        if (value != nullptr) return *value;
        return m_base_state.at(atom_index).mdpde;
    }

    const GaussianModel3D & GetModel(std::size_t atom_index) const
    {
        return GetMdpde(atom_index).GetModel();
    }

    const GaussianModel3D & GetBaseModel(std::size_t atom_index) const
    {
        return m_base_state.at(atom_index).mdpde.GetModel();
    }

    const GaussianModel3DWithUncertainty * FindOverride(std::size_t atom_index) const
    {
        return m_patch.Find(atom_index);
    }

    const ClusterKey & GetOverrideAtomIndexList() const
    {
        return m_patch.atom_index_list;
    }

    FitState Materialize() const
    {
        auto state{ m_base_state };
        m_patch.ApplyTo(state);
        return state;
    }

    std::size_t size() const { return m_base_state.size(); }
};

inline const GaussianModel3D & GetFitModel(const FitStateView & state, std::size_t atom_index)
{
    return state.GetModel(atom_index);
}

inline const GaussianModel3D & GetFitModel(
    const std::vector<GaussianModel3D> & state,
    std::size_t atom_index)
{
    return state.at(atom_index);
}

} // namespace rhbm_gem::core::detail
