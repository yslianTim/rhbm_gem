#pragma once

#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include "core/detail/LocalFittingGroupMedian.hpp"
#include "core/detail/ResidualEvaluation.hpp"
#include "core/detail/FitStateView.hpp"
#include "core/detail/SecondStageContext.hpp"

namespace rhbm_gem::core::detail {

class CandidateEvaluationOverlay
{
    const SecondStageContext & m_context;
    const SecondStageModelSnapshot & m_baseline_model_snapshot;
    const ResidualBaseline & m_residual_baseline;
    const FitStatePatch & m_patch;
    std::vector<char> m_changed_group_mask{};
    std::vector<std::optional<GaussianModel3D>> m_changed_group_median{};

public:
    CandidateEvaluationOverlay(
        const SecondStageContext & context,
        const SecondStageModelSnapshot & baseline_model_snapshot,
        const ResidualBaseline & residual_baseline,
        const FitStateView & candidate_state,
        const FitStatePatch & patch)
        : m_context{ context },
          m_baseline_model_snapshot{ baseline_model_snapshot },
          m_residual_baseline{ residual_baseline },
          m_patch{ patch },
          m_changed_group_mask(context.selected_atom_index_list_by_group.size(), 0),
          m_changed_group_median(context.selected_atom_index_list_by_group.size())
    {
        for (const auto atom_index : patch.atom_index_list)
        {
            m_changed_group_mask.at(m_context.at(atom_index).group_id) = 1;
        }
        std::vector<GaussianModel3D> model_list;
        for (std::size_t group_id = 0; group_id < m_changed_group_mask.size(); group_id++)
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
            m_changed_group_median.at(group_id) = BuildLocalFittingGaussianParameterMedian(model_list);
        }
    }

    std::optional<ResidualSample> Evaluate(
        const FitStateView & candidate_state,
        const ObjectiveSampleRef & sample_ref) const
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
        for (auto neighbor_iter = atom_context.NeighborBegin(sample_ref.sample_index);
            neighbor_iter != atom_context.NeighborEnd(sample_ref.sample_index);
            ++neighbor_iter)
        {
            const auto & neighbor_atom_sample{ *neighbor_iter };
            const GaussianModel3D * candidate_model{ nullptr };
            const GaussianModel3D * baseline_model{ nullptr };
            if (neighbor_atom_sample.is_selected)
            {
                if (m_patch.Find(neighbor_atom_sample.atom_index) == nullptr)
                {
                    continue;
                }
                baseline_model = &m_baseline_model_snapshot.selected.at(neighbor_atom_sample.atom_index);
                candidate_model = &candidate_state.GetModel(neighbor_atom_sample.atom_index);
            }
            else
            {
                const auto & unselected_atom_contributor{
                    m_context.unselected_atom_list.at(neighbor_atom_sample.atom_index)
                };
                if (!unselected_atom_contributor.selected_group_id.has_value() ||
                    m_changed_group_mask.at(*unselected_atom_contributor.selected_group_id) == 0)
                {
                    continue;
                }
                baseline_model = &m_baseline_model_snapshot.unselected.at(neighbor_atom_sample.atom_index);
                const auto & median{
                    m_changed_group_median.at(*unselected_atom_contributor.selected_group_id)
                };
                if (median.has_value())
                {
                    candidate_model = &*median;
                }
                else
                {
                    if (!unselected_atom_contributor.initial_seed.has_value())
                    {
                        return std::nullopt;
                    }
                    candidate_model = &unselected_atom_contributor.initial_seed->GetModel();
                }
            }
            adjusted_response +=
                baseline_model->ResponseAtDistance(neighbor_atom_sample.distance) -
                candidate_model->ResponseAtDistance(neighbor_atom_sample.distance);
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
        return ResidualSample{ adjusted_response, residual };
    }
};

} // namespace rhbm_gem::core::detail
