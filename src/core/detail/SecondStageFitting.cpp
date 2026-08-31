#include "core/detail/SecondStageFitting.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <ranges>
#include <stdexcept>
#include <utility>

namespace rhbm_gem::core::detail {

namespace {

double CalculateSecondStageAdjustedResponse(
    const AtomContext & atom_context,
    std::size_t sample_index,
    const SecondStageModelSnapshot & model_snapshot)
{
    auto response_value{ atom_context.raw_sampling_entries.at(sample_index).response };
    for (const auto & neighbor_atom_sample : atom_context.Neighbors(sample_index))
    {
        response_value -= ResolveNeighborAtomModel(
            neighbor_atom_sample,
            model_snapshot).ResponseAtDistance(neighbor_atom_sample.distance);
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
    LocalPotentialSampleList adjusted_sampling_entries;
    adjusted_sampling_entries.reserve(atom_context.raw_sampling_entries.size());
    for (std::size_t i = 0; i < atom_context.raw_sampling_entries.size(); i++)
    {
        auto sample{ atom_context.raw_sampling_entries.at(i) };
        sample.response = response_provider(i);
        adjusted_sampling_entries.emplace_back(sample);
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

FittedGaussianSnapshot BuildUnselectedAtomContributorSnapshot(
    const SecondStageContext & context,
    const FittedGaussianSnapshot & selected_snapshot)
{
    if (selected_snapshot.size() != context.size())
    {
        throw std::invalid_argument(
            "Second-stage selected contributor snapshot size is inconsistent.");
    }

    std::vector<std::optional<GaussianModel3D>> median_model_by_group(
        context.selected_atom_index_list_by_group.size());
    std::vector<GaussianModel3D> model_list;
    for (std::size_t group_id = 0;
        group_id < context.selected_atom_index_list_by_group.size();
        group_id++)
    {
        const auto & atom_index_list{
            context.selected_atom_index_list_by_group.at(group_id)
        };
        model_list.clear();
        model_list.reserve(atom_index_list.size());
        for (const auto atom_index : atom_index_list)
        {
            model_list.emplace_back(GetFitModel(selected_snapshot, atom_index));
        }
        median_model_by_group.at(group_id) = BuildGaussianParameterMedian(model_list);
    }

    FittedGaussianSnapshot snapshot;
    snapshot.reserve(context.unselected_atom_list.size());
    for (const auto & contributor : context.unselected_atom_list)
    {
        if (contributor.selected_group_id.has_value() &&
            median_model_by_group.at(*contributor.selected_group_id).has_value())
        {
            snapshot.emplace_back(*median_model_by_group.at(*contributor.selected_group_id));
            continue;
        }
        if (!IsValidSecondStageGaussianModel(contributor.initial_seed))
        {
            throw std::logic_error("Second-stage unselected contributor seed is unavailable.");
        }
        snapshot.emplace_back(contributor.initial_seed);
    }
    return snapshot;
}

} // namespace

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

const GaussianModel3D & ResolveNeighborAtomModel(
    const NeighborAtomSample & neighbor_atom_sample,
    const SecondStageModelSnapshot & model_snapshot)
{
    return neighbor_atom_sample.is_selected ?
        GetFitModel(model_snapshot.selected, neighbor_atom_sample.atom_index) :
        GetFitModel(model_snapshot.unselected, neighbor_atom_sample.atom_index);
}

SecondStageModelSnapshot BuildSecondStageModelSnapshot(
    const SecondStageContext & context,
    FittedGaussianSnapshot selected_snapshot)
{
    auto unselected_snapshot{
        BuildUnselectedAtomContributorSnapshot(context, selected_snapshot)
    };
    return SecondStageModelSnapshot{
        std::move(selected_snapshot),
        std::move(unselected_snapshot)
    };
}

SecondStageModelSnapshot BuildSecondStageModelSnapshot(
    const SecondStageContext & context,
    const FitState & state)
{
    if (context.size() != state.size())
    {
        throw std::invalid_argument("Local fitting context and state sizes are inconsistent.");
    }
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
                CalculateSecondStageAdjustedResponse(atom_context, j, model_snapshot));
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
    const AtomContext & atom_context,
    const SecondStageModelSnapshot & model_snapshot)
{
    return BuildSecondStageAdjustedSamplesImpl(
        atom_context,
        [&](std::size_t sample_index)
        {
            return CalculateSecondStageAdjustedResponse(
                atom_context,
                sample_index,
                model_snapshot);
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
            sample_ref.sample_index,
            model_snapshot)
    };
    const auto expected_response{
        GetFitModel(
            model_snapshot.selected,
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
        std::vector<std::vector<std::optional<ResidualSample>>>(
            context.size())
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

TransformedChangeSummary SummarizeTransformedChangesByParameter(
    const FitState & current_state,
    const FitState & previous_state,
    const TransformedChangeIndexListByParameter & index_list_by_parameter)
{
    if (current_state.size() != previous_state.size())
    {
        throw std::invalid_argument(
            "Local fitting masked transformed change state sizes are inconsistent.");
    }
    std::vector<TransformedChange> change_list;
    change_list.reserve(current_state.size());
    for (std::size_t i = 0; i < current_state.size(); i++)
    {
        change_list.emplace_back(CalculateTransformedChange(
            GetFitModel(current_state, i),
            GetFitModel(previous_state, i)));
    }
    return SummarizeTransformedChangesByParameter(change_list, index_list_by_parameter);
}

} // namespace rhbm_gem::core::detail
