#pragma once

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <vector>

#include <rhbm_gem/core/GaussianEstimator.hpp>

namespace rhbm_gem::core::detail {

using LocalFittingState = std::vector<LocalGaussianResult>;
using LocalFittingPolishProvenance = std::vector<char>;
using LocalFittingClusterKey = std::vector<std::size_t>;

struct LocalFittingStatePatch
{
    LocalFittingClusterKey atom_index_list{};
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

    void ApplyTo(LocalFittingState & state) const
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

class LocalFittingStateView
{
public:
    explicit LocalFittingStateView(const LocalFittingState & base_state)
        : m_base_state{ &base_state }
    {
    }

    LocalFittingStateView(
        const LocalFittingState & base_state,
        const LocalFittingStatePatch & patch)
        : m_base_state{ &base_state },
          m_patch{ &patch }
    {
    }

    const GaussianModel3DWithUncertainty & GetMdpde(
        std::size_t atom_index) const
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
    const LocalFittingState * m_base_state{ nullptr };
    const LocalFittingStatePatch * m_patch{ nullptr };
};

} // namespace rhbm_gem::core::detail
