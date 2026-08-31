#pragma once

#include "core/detail/GaussianModelOperations.hpp"
#include "core/detail/PreparedLocalGaussianFit.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace rhbm_gem {
class AtomObject;
}

namespace rhbm_gem::core::detail {

using ClusterKey = std::vector<std::size_t>;
using FitState = std::vector<LocalGaussianResult>;
using SecondStageAdjustedResponseCache = std::vector<std::vector<double>>;
using FittedGaussianSnapshot = std::vector<GaussianModel3D>;

struct SampleRef
{
    std::size_t atom_index{ 0 };
    std::size_t sample_index{ 0 };
    friend auto operator<=>(const SampleRef &, const SampleRef &) = default;
};

struct NeighborAtomSample
{
    bool is_selected{ true };
    std::size_t atom_index{ 0 };
    double distance{ 0.0 };
};

struct UnselectedAtomContributor
{
    std::optional<std::size_t> selected_group_id{};
    GaussianModel3D initial_seed{};
};

struct AtomContext
{
    const AtomObject * atom{ nullptr };
    std::size_t group_id{ 0 };
    LocalPotentialSampleList raw_sampling_entries{};
    std::vector<NeighborAtomSample> neighbor_atom_sample_list{};
    std::vector<std::size_t> neighbor_atom_sample_offset_list{};
    PreparedLocalGaussianDesign refit_design{};
    double alpha_r{ 0.0 };

    std::span<const NeighborAtomSample> Neighbors(std::size_t sample_index) const
    {
        const auto begin{ neighbor_atom_sample_offset_list.at(sample_index) };
        const auto end{ neighbor_atom_sample_offset_list.at(sample_index + 1) };
        return std::span<const NeighborAtomSample>{ neighbor_atom_sample_list }
            .subspan(begin, end - begin);
    }
};

struct SecondStageContext
{
    std::vector<AtomContext> selected_atom_list{};
    std::vector<UnselectedAtomContributor> unselected_atom_list{};
    std::vector<std::vector<std::size_t>> selected_atom_index_list_by_group{};

    std::size_t size() const { return selected_atom_list.size(); }
    AtomContext & at(std::size_t index) { return selected_atom_list.at(index); }
    const AtomContext & at(std::size_t index) const { return selected_atom_list.at(index); }
    auto begin() { return selected_atom_list.begin(); }
    auto end() { return selected_atom_list.end(); }
    auto begin() const { return selected_atom_list.begin(); }
    auto end() const { return selected_atom_list.end(); }
};

struct FitStatePatch
{
    ClusterKey atom_index_list{};
    std::vector<GaussianModel3DWithUncertainty> mdpde_list{};

    static FitStatePatch FromState(const FitState & state, ClusterKey atom_index_list);
    const GaussianModel3DWithUncertainty * Find(std::size_t atom_index) const;
    void ApplyTo(FitState & state) const;
};

struct FitStateProposal
{
    FitStatePatch patch{};
    double effective_damping{ 0.0 };
    double step_norm{ 0.0 };
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

    const GaussianModel3DWithUncertainty * FindOverride(
        std::size_t atom_index) const
    {
        return m_patch.Find(atom_index);
    }

    const ClusterKey & GetOverrideAtomIndexList() const
    {
        return m_patch.atom_index_list;
    }

    std::size_t size() const { return m_base_state.size(); }
};

const GaussianModel3D & GetFitModel(
    const FitState & state,
    std::size_t atom_index);

const GaussianModel3D & GetFitModel(
    const FitStateView & state,
    std::size_t atom_index);

const GaussianModel3D & GetFitModel(
    const FittedGaussianSnapshot & state,
    std::size_t atom_index);

FittedGaussianSnapshot BuildFittedGaussianSnapshot(const FitState & state);
FittedGaussianSnapshot BuildFittedGaussianSnapshot(const FitStateView & state);

struct SecondStageModelSnapshot
{
    FittedGaussianSnapshot selected{};
    FittedGaussianSnapshot unselected{};
};

const GaussianModel3D & ResolveNeighborAtomModel(
    const NeighborAtomSample & neighbor_atom_sample,
    const SecondStageModelSnapshot & model_snapshot);

SecondStageModelSnapshot BuildSecondStageModelSnapshot(
    const SecondStageContext & context,
    FittedGaussianSnapshot selected_snapshot);

SecondStageModelSnapshot BuildSecondStageModelSnapshot(
    const SecondStageContext & context,
    const FitState & state);

struct ResidualSample
{
    double adjusted_response{ 0.0 };
    double residual{ 0.0 };
};

struct ResidualBaseline
{
    SecondStageModelSnapshot model_snapshot{};
    std::vector<std::vector<std::optional<ResidualSample>>> sample_list{};

    std::optional<ResidualSample> operator()(const SampleRef & sample_ref) const
    {
        return sample_list.at(sample_ref.atom_index).at(sample_ref.sample_index);
    }

    const FittedGaussianSnapshot & GetState() const
    {
        return model_snapshot.selected;
    }
};

SecondStageAdjustedResponseCache BuildSecondStageAdjustedResponseCache(
    const SecondStageContext & context,
    const SecondStageModelSnapshot & model_snapshot);

LocalPotentialSampleList BuildSecondStageAdjustedSamples(
    const AtomContext & atom_context,
    const std::vector<double> & adjusted_response_list);

LocalPotentialSampleList BuildSecondStageAdjustedSamples(
    const AtomContext & atom_context,
    const SecondStageModelSnapshot & model_snapshot);

std::optional<ResidualSample> EvaluateResidualSample(
    const SecondStageContext & context,
    const SampleRef & sample_ref,
    const SecondStageModelSnapshot & model_snapshot);

struct SnapshotResidualEvaluator
{
    const SecondStageContext & context;
    const SecondStageModelSnapshot & model_snapshot;

    std::optional<ResidualSample> operator()(const SampleRef & sample_ref) const;
    const FittedGaussianSnapshot & GetState() const
    {
        return model_snapshot.selected;
    }
};

ResidualBaseline BuildResidualBaseline(
    const SecondStageContext & context,
    const FitState & state);

TransformedChangeSummary SummarizeTransformedChanges(
    const FitState & current_state,
    const FitState & previous_state,
    const std::vector<std::size_t> & index_list);

TransformedChangeSummary SummarizeTransformedChanges(
    const FitStateView & current_state,
    const FittedGaussianSnapshot & previous_state,
    const std::vector<std::size_t> & index_list);

TransformedChangeSummary SummarizeTransformedChangesByParameter(
    const FitState & current_state,
    const FitState & previous_state,
    const TransformedChangeIndexListByParameter & index_list_by_parameter);

} // namespace rhbm_gem::core::detail
