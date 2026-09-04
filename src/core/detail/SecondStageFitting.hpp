#pragma once

#include "core/detail/GaussianModelOperations.hpp"
#include "core/detail/PreparedLocalGaussianFit.hpp"

#include <cstddef>
#include <memory>
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
    std::size_t atom_index{ 0 };
    double distance{ 0.0 };
};

struct AtomContext
{
    const AtomObject * atom{ nullptr };
    LocalPotentialSampleList raw_sampling_entries{};
    std::vector<NeighborAtomSample> neighbor_atom_sample_list{};
    std::vector<std::size_t> neighbor_atom_sample_offset_list{};
    std::vector<std::vector<double>> unselected_distance_list_by_sample{};
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

struct FrozenBackground
{
    FittedGaussianSnapshot model_by_atom{};
    SecondStageAdjustedResponseCache response_by_atom{};
};

struct SecondStageContext
{
    std::vector<AtomContext> atom_list{};
    std::shared_ptr<const FrozenBackground> frozen_background{};

    std::size_t size() const { return atom_list.size(); }
    AtomContext & at(std::size_t index) { return atom_list.at(index); }
    const AtomContext & at(std::size_t index) const { return atom_list.at(index); }
    auto begin() const { return atom_list.begin(); }
    auto end() const { return atom_list.end(); }
};

std::shared_ptr<const FrozenBackground> BuildFrozenBackground(
    const SecondStageContext & context,
    const FitState & state,
    const std::vector<ClusterKey> & cluster_key_list);

double GetFrozenBackgroundResponse(const SecondStageContext & context, const SampleRef & sample_ref);

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
    FittedGaussianSnapshot node{};
    std::shared_ptr<const FrozenBackground> frozen_background{};
};

SecondStageModelSnapshot BuildSecondStageModelSnapshot(
    const SecondStageContext & context,
    FittedGaussianSnapshot node_snapshot);

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
        return model_snapshot.node;
    }
};

SecondStageAdjustedResponseCache BuildSecondStageAdjustedResponseCache(
    const SecondStageContext & context,
    const SecondStageModelSnapshot & model_snapshot);

LocalPotentialSampleList BuildSecondStageAdjustedSamples(
    const AtomContext & atom_context,
    const std::vector<double> & adjusted_response_list);

LocalPotentialSampleList BuildSecondStageAdjustedSamples(
    const SecondStageContext & context,
    std::size_t atom_index,
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
        return model_snapshot.node;
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

} // namespace rhbm_gem::core::detail
