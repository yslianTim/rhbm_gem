#pragma once

#include <cstddef>
#include <optional>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>

#include "core/detail/CouplingGraph.hpp"

namespace rhbm_gem::core::detail {

struct NeighborAtomSample
{
    bool is_selected{ true };
    std::size_t atom_index{ 0 };
    double distance{ 0.0 };
};

struct UnselectedAtomContributor
{
    const AtomObject * atom{ nullptr };
    GroupKey group_key{};
    std::optional<std::size_t> selected_group_id{};
    std::optional<GaussianModel3DWithUncertainty> initial_seed{};
};

struct AtomContext
{
    const AtomObject * atom{ nullptr };
    GroupKey group_key{};
    std::size_t group_id{ 0 };
    GraphResidueKey residue_key{};
    LocalPotentialSampleList raw_sampling_entries{};
    LocalGaussianResult initial_result{};
    std::optional<GaussianModel3DWithUncertainty> group_prior{};
    std::vector<NeighborAtomSample> neighbor_atom_sample_list{};
    std::vector<std::size_t> neighbor_atom_sample_offset_list{};
    LocalGaussianDesignTemplate refit_design_template{};
    double alpha_r{ 0.0 };
    int neighbor_count_for_peeling{ 0 };

    auto NeighborBegin(std::size_t sample_index) const
    {
        return neighbor_atom_sample_list.begin() +
            static_cast<std::ptrdiff_t>(
                neighbor_atom_sample_offset_list.at(sample_index));
    }

    auto NeighborEnd(std::size_t sample_index) const
    {
        return neighbor_atom_sample_list.begin() +
            static_cast<std::ptrdiff_t>(
                neighbor_atom_sample_offset_list.at(sample_index + 1));
    }
};

struct SecondStageContext
{
    std::vector<AtomContext> selected_atom_list{};
    std::vector<UnselectedAtomContributor> unselected_atom_list{};
    std::unordered_map<GroupKey, std::size_t> selected_group_id_by_key{};
    std::vector<std::vector<std::size_t>> selected_atom_index_list_by_group{};

    std::size_t size() const { return selected_atom_list.size(); }
    AtomContext & at(std::size_t index) { return selected_atom_list.at(index); }
    const AtomContext & at(std::size_t index) const { return selected_atom_list.at(index); }
    auto begin() { return selected_atom_list.begin(); }
    auto end() { return selected_atom_list.end(); }
    auto begin() const { return selected_atom_list.begin(); }
    auto end() const { return selected_atom_list.end(); }
};

} // namespace rhbm_gem::core::detail
