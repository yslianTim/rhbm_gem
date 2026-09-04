#include "core/detail/SecondStageFitting.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <stdexcept>
#include <utility>

namespace rhbm_gem::core::detail {

namespace {

double CalculateSecondStageAdjustedResponse(
    const AtomContext & atom_context,
    const SampleRef & sample_ref,
    const SecondStageModelSnapshot & model_snapshot)
{
    const auto sample_index{ sample_ref.sample_index };
    auto response_value{ atom_context.raw_sampling_entries.at(sample_index).response };
    if (model_snapshot.frozen_background)
    {
        response_value -= model_snapshot.frozen_background->response_by_atom
            .at(sample_ref.atom_index).at(sample_index);
    }
    for (const auto & neighbor_atom_sample : atom_context.Neighbors(sample_index))
    {
        response_value -= GetFitModel(
            model_snapshot.node,
            neighbor_atom_sample.atom_index).ResponseAtDistance(neighbor_atom_sample.distance);
    }
    return response_value;
}

template<typename State>
FittedGaussianSnapshot BuildFittedGaussianSnapshotImpl(const State & state)
{
    FittedGaussianSnapshot snapshot;
    snapshot.reserve(state.size());
    for (std::size_t i = 0; i < state.size(); i++)
    {
        snapshot.emplace_back(GetFitModel(state, i));
    }
    return snapshot;
}

template<typename ResponseProvider>
LocalPotentialSampleList BuildSecondStageAdjustedSamplesImpl(
    const AtomContext & atom_context,
    ResponseProvider response_provider)
{
    auto adjusted_sampling_entries{ atom_context.raw_sampling_entries };
    for (std::size_t i = 0; i < adjusted_sampling_entries.size(); i++)
    {
        adjusted_sampling_entries.at(i).response = response_provider(i);
    }
    return adjusted_sampling_entries;
}

template<typename CurrentState, typename PreviousState>
TransformedChangeSummary SummarizeTransformedChangesImpl(
    const CurrentState & current_state,
    const PreviousState & previous_state,
    const std::vector<std::size_t> & index_list)
{
    std::vector<TransformedChange> change_list;
    change_list.reserve(index_list.size());
    for (const auto i : index_list)
    {
        change_list.emplace_back(CalculateTransformedChange(
            GetFitModel(current_state, i),
            GetFitModel(previous_state, i)));
    }
    return SummarizeTransformedChanges(change_list);
}

} // namespace

std::shared_ptr<const FrozenBackground> BuildFrozenBackground(
    const SecondStageContext & context,
    const FitState & state,
    const std::vector<ClusterKey> & cluster_key_list)
{
    if (state.size() != context.size()) return nullptr;
    auto background{ std::make_shared<FrozenBackground>() };
    background->model_by_atom.resize(context.size());
    background->response_by_atom.resize(context.size());
    std::vector<char> visited(context.size(), 0);
    for (const auto & key : cluster_key_list)
    {
        std::vector<GaussianModel3D> models;
        for (const auto atom_index : key)
        {
            if (atom_index >= state.size() || visited.at(atom_index) != 0) return nullptr;
            visited.at(atom_index) = 1;
            const auto & model{ state.at(atom_index).mdpde.GetModel() };
            if (!IsValidSecondStageGaussianModel(model)) return nullptr;
            models.emplace_back(model);
        }
        const auto median{ BuildGaussianParameterMedian(models) };
        if (!median.has_value()) return nullptr;
        for (const auto atom_index : key)
        {
            background->model_by_atom.at(atom_index) = *median;
            const auto & atom_context{ context.at(atom_index) };
            auto & responses{ background->response_by_atom.at(atom_index) };
            responses.assign(atom_context.raw_sampling_entries.size(), 0.0);
            if (atom_context.unselected_distance_list_by_sample.empty()) continue;
            if (atom_context.unselected_distance_list_by_sample.size() != responses.size()) return nullptr;
            for (std::size_t sample_index = 0; sample_index < responses.size(); sample_index++)
            {
                for (const auto distance : atom_context.unselected_distance_list_by_sample.at(sample_index))
                {
                    responses.at(sample_index) += median->ResponseAtDistance(distance);
                }
                if (!std::isfinite(responses.at(sample_index))) return nullptr;
            }
        }
    }
    if (std::ranges::find(visited, 0) != visited.end()) return nullptr;
    return background;
}

double GetFrozenBackgroundResponse(const SecondStageContext & context, const SampleRef & sample_ref)
{
    return context.frozen_background ?
        context.frozen_background->response_by_atom.at(sample_ref.atom_index).at(sample_ref.sample_index) :
        0.0;
}

FitStatePatch FitStatePatch::FromState(const FitState & state, ClusterKey atom_index_list)
{
    std::ranges::sort(atom_index_list);
    atom_index_list.erase(
        std::ranges::unique(atom_index_list).begin(),
        atom_index_list.end());

    FitStatePatch patch;
    patch.atom_index_list = std::move(atom_index_list);
    patch.mdpde_list.reserve(patch.atom_index_list.size());
    for (const auto atom_index : patch.atom_index_list)
    {
        patch.mdpde_list.emplace_back(state.at(atom_index).mdpde);
    }
    return patch;
}

const GaussianModel3DWithUncertainty * FitStatePatch::Find(std::size_t atom_index) const
{
    const auto iter{ std::ranges::lower_bound(atom_index_list, atom_index) };
    if (iter == atom_index_list.end() || *iter != atom_index) return nullptr;
    return &mdpde_list.at(static_cast<std::size_t>(std::distance(atom_index_list.begin(), iter)));
}

void FitStatePatch::ApplyTo(FitState & state) const
{
    if (atom_index_list.size() != mdpde_list.size())
    {
        throw std::invalid_argument("Local fitting state patch sizes are inconsistent.");
    }
    for (std::size_t i = 0; i < atom_index_list.size(); i++)
    {
        state.at(atom_index_list.at(i)).mdpde = mdpde_list.at(i);
    }
}

std::optional<ResidualSample> SnapshotResidualEvaluator::operator()(const SampleRef & sample_ref) const
{
    return EvaluateResidualSample(context, sample_ref, model_snapshot);
}

const GaussianModel3D & GetFitModel(const FitState & state, std::size_t atom_index)
{
    return state.at(atom_index).mdpde.GetModel();
}

const GaussianModel3D & GetFitModel(const FitStateView & state, std::size_t atom_index)
{
    return state.GetModel(atom_index);
}

const GaussianModel3D & GetFitModel(const FittedGaussianSnapshot & state, std::size_t atom_index)
{
    return state.at(atom_index);
}

FittedGaussianSnapshot BuildFittedGaussianSnapshot(const FitState & state)
{
    return BuildFittedGaussianSnapshotImpl(state);
}

FittedGaussianSnapshot BuildFittedGaussianSnapshot(const FitStateView & state)
{
    return BuildFittedGaussianSnapshotImpl(state);
}

SecondStageModelSnapshot BuildSecondStageModelSnapshot(
    const SecondStageContext & context,
    FittedGaussianSnapshot node_snapshot)
{
    if (node_snapshot.size() != context.size())
    {
        throw std::invalid_argument("Second-stage node snapshot size is inconsistent.");
    }
    return SecondStageModelSnapshot{ std::move(node_snapshot), context.frozen_background };
}

SecondStageModelSnapshot BuildSecondStageModelSnapshot(
    const SecondStageContext & context,
    const FitState & state)
{
    return BuildSecondStageModelSnapshot(context, BuildFittedGaussianSnapshot(state));
}

SecondStageAdjustedResponseCache BuildSecondStageAdjustedResponseCache(
    const SecondStageContext & context,
    const SecondStageModelSnapshot & model_snapshot)
{
    SecondStageAdjustedResponseCache cache(context.size());
    for (std::size_t i = 0; i < context.size(); i++)
    {
        const auto & atom_context{ context.at(i) };
        const auto sample_count{ atom_context.raw_sampling_entries.size() };
        auto & response_list{ cache.at(i) };
        response_list.reserve(sample_count);
        for (std::size_t j = 0; j < sample_count; j++)
        {
            response_list.emplace_back(
                CalculateSecondStageAdjustedResponse(atom_context, SampleRef{ i, j }, model_snapshot));
        }
    }
    return cache;
}

LocalPotentialSampleList BuildSecondStageAdjustedSamples(
    const AtomContext & atom_context,
    const std::vector<double> & adjusted_response_list)
{
    if (adjusted_response_list.size() != atom_context.raw_sampling_entries.size())
    {
        throw std::invalid_argument("Second-stage adjusted response count is inconsistent.");
    }
    return BuildSecondStageAdjustedSamplesImpl(
        atom_context,
        [&](std::size_t sample_index)
        {
            return adjusted_response_list.at(sample_index);
        });
}

LocalPotentialSampleList BuildSecondStageAdjustedSamples(
    const SecondStageContext & context,
    std::size_t atom_index,
    const SecondStageModelSnapshot & model_snapshot)
{
    const auto & atom_context{ context.at(atom_index) };
    return BuildSecondStageAdjustedSamplesImpl(
        atom_context,
        [&](std::size_t sample_index)
        {
            return CalculateSecondStageAdjustedResponse(atom_context, SampleRef{ atom_index, sample_index }, model_snapshot);
        });
}

std::optional<ResidualSample> EvaluateResidualSample(
    const SecondStageContext & context,
    const SampleRef & sample_ref,
    const SecondStageModelSnapshot & model_snapshot)
{
    const auto & atom_context{ context.at(sample_ref.atom_index) };
    const auto & sample{
        atom_context.raw_sampling_entries.at(sample_ref.sample_index)
    };
    const auto adjusted_response{
        CalculateSecondStageAdjustedResponse(
            atom_context,
            sample_ref,
            model_snapshot)
    };
    const auto expected_response{
        GetFitModel(
            model_snapshot.node,
            sample_ref.atom_index).ResponseAtDistance(sample.point.distance)
    };
    const auto residual{ adjusted_response - expected_response };
    if (!std::isfinite(adjusted_response) || !std::isfinite(residual)) return std::nullopt;
    return ResidualSample{ adjusted_response, residual };
}

ResidualBaseline BuildResidualBaseline(const SecondStageContext & context, const FitState & state)
{
    ResidualBaseline baseline{
        BuildSecondStageModelSnapshot(context, state),
        std::vector<std::vector<std::optional<ResidualSample>>>(context.size())
    };
    for (std::size_t i = 0; i < context.size(); i++)
    {
        const auto sample_count{ context.at(i).raw_sampling_entries.size() };
        baseline.sample_list.at(i).reserve(sample_count);
        for (std::size_t j = 0; j < sample_count; j++)
        {
            baseline.sample_list.at(i).emplace_back(
                EvaluateResidualSample(
                    context,
                    SampleRef{ i, j },
                    baseline.model_snapshot));
        }
    }
    return baseline;
}

TransformedChangeSummary SummarizeTransformedChanges(
    const FitState & current_state,
    const FitState & previous_state,
    const std::vector<std::size_t> & index_list)
{
    return SummarizeTransformedChangesImpl(current_state, previous_state, index_list);
}

TransformedChangeSummary SummarizeTransformedChanges(
    const FitStateView & current_state,
    const FittedGaussianSnapshot & previous_state,
    const std::vector<std::size_t> & index_list)
{
    return SummarizeTransformedChangesImpl(current_state, previous_state, index_list);
}

} // namespace rhbm_gem::core::detail
