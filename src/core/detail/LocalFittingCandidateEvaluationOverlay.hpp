#pragma once

#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include "core/detail/LocalFittingGroupMedian.hpp"
#include "core/detail/LocalFittingStateView.hpp"
#include "core/detail/SecondStageLocalFittingContext.hpp"

namespace rhbm_gem::core::detail {

using LocalFittingObjectiveSampleRef = GraphSampleId;
using FittedGaussianSnapshot = std::vector<GaussianModel3D>;

inline FittedGaussianSnapshot BuildFittedGaussianSnapshot(
    const LocalFittingState & state)
{
    FittedGaussianSnapshot snapshot;
    snapshot.reserve(state.size());
    for (const auto & result : state)
    {
        snapshot.emplace_back(result.mdpde.GetModel());
    }
    return snapshot;
}

inline FittedGaussianSnapshot BuildFittedGaussianSnapshot(
    const LocalFittingStateView & state)
{
    FittedGaussianSnapshot snapshot;
    snapshot.reserve(state.GetSize());
    for (std::size_t atom_index = 0;
        atom_index < state.GetSize();
        atom_index++)
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

inline const GaussianModel3D & ResolveSecondStageNeighborModel(
    const SecondStageNeighborSample & neighbor_sample,
    const FittedGaussianSnapshot & selected_snapshot,
    const FittedGaussianSnapshot & unselected_snapshot)
{
    return neighbor_sample.is_selected ?
        selected_snapshot.at(neighbor_sample.atom_index) :
        unselected_snapshot.at(neighbor_sample.atom_index);
}

inline FittedGaussianSnapshot BuildUnselectedContributorSnapshot(
    const SecondStageLocalFittingContext & context,
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
            model_list.emplace_back(selected_snapshot.at(atom_index));
        }
        median_model_by_group.at(group_id) =
            BuildLocalFittingGaussianParameterMedian(model_list);
    }

    FittedGaussianSnapshot snapshot;
    snapshot.reserve(context.unselected_atom_list.size());
    for (const auto & contributor : context.unselected_atom_list)
    {
        if (contributor.selected_group_id.has_value() &&
            median_model_by_group.at(*contributor.selected_group_id).has_value())
        {
            snapshot.emplace_back(
                *median_model_by_group.at(*contributor.selected_group_id));
            continue;
        }
        if (!contributor.initial_seed.has_value())
        {
            throw std::logic_error(
                "Second-stage unselected contributor seed is unavailable.");
        }
        snapshot.emplace_back(contributor.initial_seed->GetModel());
    }
    return snapshot;
}

inline SecondStageModelSnapshot BuildSecondStageModelSnapshot(
    const SecondStageLocalFittingContext & context,
    FittedGaussianSnapshot selected_snapshot)
{
    auto unselected_snapshot{
        BuildUnselectedContributorSnapshot(context, selected_snapshot)
    };
    return SecondStageModelSnapshot{
        std::move(selected_snapshot),
        std::move(unselected_snapshot)
    };
}

struct LocalFittingResidualSample
{
    double adjusted_response{ 0.0 };
    double residual{ 0.0 };
};

using LocalFittingResidualBaseline =
    std::vector<std::vector<std::optional<LocalFittingResidualSample>>>;

class LocalFittingCandidateEvaluationOverlay
{
public:
    LocalFittingCandidateEvaluationOverlay(
        const SecondStageLocalFittingContext & context,
        const SecondStageModelSnapshot & baseline_model_snapshot,
        const LocalFittingResidualBaseline & residual_baseline,
        const LocalFittingStateView & candidate_state,
        const LocalFittingStatePatch & patch)
        : m_context{ context },
          m_baseline_model_snapshot{ baseline_model_snapshot },
          m_residual_baseline{ residual_baseline },
          m_patch{ patch },
          m_changed_group_mask(
              context.selected_atom_index_list_by_group.size(),
              0),
          m_changed_group_median(
              context.selected_atom_index_list_by_group.size())
    {
        for (const auto atom_index : patch.atom_index_list)
        {
            m_changed_group_mask.at(m_context.at(atom_index).group_id) = 1;
        }
        std::vector<GaussianModel3D> model_list;
        for (std::size_t group_id = 0;
            group_id < m_changed_group_mask.size();
            group_id++)
        {
            if (m_changed_group_mask.at(group_id) == 0) continue;
            const auto & atom_index_list{
                m_context.selected_atom_index_list_by_group.at(group_id)
            };
            model_list.clear();
            model_list.reserve(atom_index_list.size());
            for (const auto atom_index : atom_index_list)
            {
                model_list.emplace_back(candidate_state.GetModel(atom_index));
            }
            m_changed_group_median.at(group_id) =
                BuildLocalFittingGaussianParameterMedian(model_list);
        }
    }

    std::optional<LocalFittingResidualSample> Evaluate(
        const LocalFittingStateView & candidate_state,
        const LocalFittingObjectiveSampleRef & sample_ref) const
    {
        const auto & baseline{
            m_residual_baseline.at(sample_ref.atom_index).at(
                sample_ref.sample_index)
        };
        if (!baseline.has_value()) return std::nullopt;
        const auto & atom_context{ m_context.at(sample_ref.atom_index) };
        const auto & sample{
            atom_context.raw_sampling_entries.at(sample_ref.sample_index)
        };
        auto adjusted_response{ baseline->adjusted_response };
        for (auto neighbor_iter =
                atom_context.NeighborBegin(sample_ref.sample_index);
            neighbor_iter != atom_context.NeighborEnd(sample_ref.sample_index);
            ++neighbor_iter)
        {
            const auto & neighbor_sample{ *neighbor_iter };
            const GaussianModel3D * candidate_model{ nullptr };
            const GaussianModel3D * baseline_model{ nullptr };
            if (neighbor_sample.is_selected)
            {
                if (m_patch.Find(neighbor_sample.atom_index) == nullptr)
                {
                    continue;
                }
                baseline_model = &m_baseline_model_snapshot.selected.at(
                    neighbor_sample.atom_index);
                candidate_model = &candidate_state.GetModel(
                    neighbor_sample.atom_index);
            }
            else
            {
                const auto & contributor{
                    m_context.unselected_atom_list.at(
                        neighbor_sample.atom_index)
                };
                if (!contributor.selected_group_id.has_value() ||
                    m_changed_group_mask.at(*contributor.selected_group_id) == 0)
                {
                    continue;
                }
                baseline_model = &m_baseline_model_snapshot.unselected.at(
                    neighbor_sample.atom_index);
                const auto & median{
                    m_changed_group_median.at(
                        *contributor.selected_group_id)
                };
                if (median.has_value())
                {
                    candidate_model = &*median;
                }
                else
                {
                    if (!contributor.initial_seed.has_value())
                    {
                        return std::nullopt;
                    }
                    candidate_model = &contributor.initial_seed->GetModel();
                }
            }
            adjusted_response +=
                baseline_model->ResponseAtDistance(neighbor_sample.distance) -
                candidate_model->ResponseAtDistance(neighbor_sample.distance);
        }
        const auto expected_response{
            m_patch.Find(sample_ref.atom_index) == nullptr ?
                baseline->adjusted_response - baseline->residual :
                candidate_state.GetModel(sample_ref.atom_index).ResponseAtDistance(
                    static_cast<double>(sample.point.distance))
        };
        const auto residual{ adjusted_response - expected_response };
        if (!std::isfinite(adjusted_response) || !std::isfinite(residual))
        {
            return std::nullopt;
        }
        return LocalFittingResidualSample{ adjusted_response, residual };
    }

private:
    const SecondStageLocalFittingContext & m_context;
    const SecondStageModelSnapshot & m_baseline_model_snapshot;
    const LocalFittingResidualBaseline & m_residual_baseline;
    const LocalFittingStatePatch & m_patch;
    std::vector<char> m_changed_group_mask{};
    std::vector<std::optional<GaussianModel3D>> m_changed_group_median{};
};

} // namespace rhbm_gem::core::detail
