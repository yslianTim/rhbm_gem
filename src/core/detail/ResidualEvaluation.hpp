#pragma once

#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include "core/detail/TransformedGaussianModel.hpp"
#include "core/detail/FitStateView.hpp"
#include "core/detail/SecondStageContext.hpp"

namespace rhbm_gem::core::detail {

using SecondStageAdjustedResponseCache = std::vector<std::vector<double>>;

using FittedGaussianSnapshot = std::vector<GaussianModel3D>;

inline FittedGaussianSnapshot BuildFittedGaussianSnapshot(const FitState & state)
{
    FittedGaussianSnapshot snapshot;
    snapshot.reserve(state.size());
    for (const auto & result : state)
    {
        snapshot.emplace_back(result.mdpde.GetModel());
    }
    return snapshot;
}

inline FittedGaussianSnapshot BuildFittedGaussianSnapshot(const FitStateView & state)
{
    FittedGaussianSnapshot snapshot;
    snapshot.reserve(state.GetSize());
    for (std::size_t atom_index = 0; atom_index < state.GetSize(); atom_index++)
    {
        snapshot.emplace_back(state.GetModel(atom_index));
    }
    return snapshot;
}

struct SecondStageModelSnapshot
{
    FittedGaussianSnapshot selected{};
    FittedGaussianSnapshot unselected{};
};

inline const GaussianModel3D & ResolveNeighborAtomModel(
    const NeighborAtomSample & neighbor_atom_sample,
    const SecondStageModelSnapshot & model_snapshot)
{
    return neighbor_atom_sample.is_selected ?
        GetFitModel(model_snapshot.selected, neighbor_atom_sample.atom_index) :
        GetFitModel(model_snapshot.unselected, neighbor_atom_sample.atom_index);
}

inline FittedGaussianSnapshot BuildUnselectedAtomContributorSnapshot(
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
    for (const auto & unselected_atom_contributor : context.unselected_atom_list)
    {
        if (unselected_atom_contributor.selected_group_id.has_value() &&
            median_model_by_group.at(*unselected_atom_contributor.selected_group_id).has_value())
        {
            snapshot.emplace_back(
                *median_model_by_group.at(*unselected_atom_contributor.selected_group_id));
            continue;
        }
        if (!unselected_atom_contributor.initial_seed.has_value())
        {
            throw std::logic_error(
                "Second-stage unselected contributor seed is unavailable.");
        }
        snapshot.emplace_back(
            unselected_atom_contributor.initial_seed->GetModel());
    }
    return snapshot;
}

inline SecondStageModelSnapshot BuildSecondStageModelSnapshot(
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

inline SecondStageModelSnapshot BuildSecondStageModelSnapshot(
    const SecondStageContext & context,
    const FitState & state)
{
    if (context.size() != state.size())
    {
        throw std::invalid_argument(
            "Local fitting context and state sizes are inconsistent.");
    }
    return BuildSecondStageModelSnapshot(context, BuildFittedGaussianSnapshot(state));
}

struct ResidualSample
{
    double adjusted_response{ 0.0 };
    double residual{ 0.0 };
};

struct ResidualBaseline
{
    SecondStageModelSnapshot model_snapshot{};
    std::vector<std::vector<std::optional<ResidualSample>>> sample_list{};
};

inline double CalculateSecondStageAdjustedResponse(
    const AtomContext & atom_context,
    std::size_t sample_index,
    const SecondStageModelSnapshot & model_snapshot)
{
    auto response_value{
        static_cast<double>(atom_context.raw_sampling_entries.at(sample_index).response)
    };
    for (auto neighbor_iter = atom_context.NeighborBegin(sample_index);
        neighbor_iter != atom_context.NeighborEnd(sample_index);
        ++neighbor_iter)
    {
        const auto & neighbor_atom_sample{ *neighbor_iter };
        response_value -= ResolveNeighborAtomModel(
            neighbor_atom_sample,
            model_snapshot).ResponseAtDistance(neighbor_atom_sample.distance);
    }
    return response_value;
}

inline SecondStageAdjustedResponseCache BuildSecondStageAdjustedResponseCache(
    const SecondStageContext & context,
    const SecondStageModelSnapshot & model_snapshot)
{
    SecondStageAdjustedResponseCache cache(context.size());
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        const auto sample_count{
            context.at(atom_index).raw_sampling_entries.size()
        };
        auto & response_list{ cache.at(atom_index) };
        response_list.reserve(sample_count);
        for (std::size_t sample_index = 0; sample_index < sample_count; sample_index++)
        {
            response_list.emplace_back(static_cast<double>(static_cast<float>(
                CalculateSecondStageAdjustedResponse(
                    context.at(atom_index), sample_index, model_snapshot))));
        }
    }
    return cache;
}

inline LocalPotentialSampleList BuildSecondStageAdjustedSamples(
    const AtomContext & atom_context,
    const std::vector<double> & adjusted_response_list)
{
    if (adjusted_response_list.size() != atom_context.raw_sampling_entries.size())
    {
        throw std::invalid_argument(
            "Second-stage adjusted response count is inconsistent.");
    }
    LocalPotentialSampleList adjusted_sampling_entries;
    adjusted_sampling_entries.reserve(adjusted_response_list.size());
    for (std::size_t sample_index = 0; sample_index < adjusted_response_list.size(); sample_index++)
    {
        auto sample{ atom_context.raw_sampling_entries.at(sample_index) };
        sample.response = static_cast<float>(adjusted_response_list.at(sample_index));
        adjusted_sampling_entries.emplace_back(sample);
    }
    return adjusted_sampling_entries;
}

inline LocalPotentialSampleList BuildSecondStageAdjustedSamples(
    const AtomContext & atom_context,
    const SecondStageModelSnapshot & model_snapshot)
{
    LocalPotentialSampleList adjusted_sampling_entries;
    adjusted_sampling_entries.reserve(atom_context.raw_sampling_entries.size());
    for (std::size_t sample_index = 0;
        sample_index < atom_context.raw_sampling_entries.size();
        sample_index++)
    {
        auto sample{ atom_context.raw_sampling_entries.at(sample_index) };
        sample.response = static_cast<float>(
            CalculateSecondStageAdjustedResponse(atom_context, sample_index, model_snapshot));
        adjusted_sampling_entries.emplace_back(sample);
    }
    return adjusted_sampling_entries;
}

inline ResidualBaseline BuildResidualBaseline(const SecondStageContext & context, const FitState & state)
{
    ResidualBaseline baseline{
        BuildSecondStageModelSnapshot(context, state),
        std::vector<std::vector<std::optional<ResidualSample>>>(context.size())
    };
    for (std::size_t atom_index = 0; atom_index < context.size(); atom_index++)
    {
        const auto sample_count{
            context.at(atom_index).raw_sampling_entries.size()
        };
        baseline.sample_list.at(atom_index).reserve(sample_count);
        for (std::size_t sample_index = 0; sample_index < sample_count; sample_index++)
        {
            const auto adjusted_response{
                CalculateSecondStageAdjustedResponse(
                    context.at(atom_index),
                    sample_index,
                    baseline.model_snapshot)
            };
            const auto & sample{
                context.at(atom_index).raw_sampling_entries.at(sample_index)
            };
            const auto expected_response{
                GetFitModel(baseline.model_snapshot.selected, atom_index)
                .ResponseAtDistance(static_cast<double>(sample.point.distance))
            };
            const auto residual{ adjusted_response - expected_response };
            baseline.sample_list.at(atom_index).emplace_back(
                std::isfinite(adjusted_response) && std::isfinite(residual) ?
                    std::optional<ResidualSample>{
                        ResidualSample{
                            adjusted_response,
                            residual } } :
                    std::nullopt);
        }
    }
    return baseline;
}

template <typename State>
inline std::optional<ResidualSample> EvaluateResidualSample(
    const SecondStageContext & context,
    const State & state,
    const SampleRef & sample_ref,
    const SecondStageModelSnapshot & model_snapshot)
{
    const auto & atom_context{ context.at(sample_ref.atom_index) };
    const auto & sample{
        atom_context.raw_sampling_entries.at(sample_ref.sample_index)
    };
    const auto adjusted_response{
        CalculateSecondStageAdjustedResponse(atom_context, sample_ref.sample_index, model_snapshot)
    };
    const auto expected_response{
        GetFitModel(state, sample_ref.atom_index).ResponseAtDistance(static_cast<double>(sample.point.distance))
    };
    const auto residual{ adjusted_response - expected_response };
    if (!std::isfinite(adjusted_response) || !std::isfinite(residual))
    {
        return std::nullopt;
    }
    return ResidualSample{ adjusted_response, residual };
}

} // namespace rhbm_gem::core::detail
