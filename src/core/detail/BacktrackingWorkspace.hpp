#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>

#include "core/detail/TransformedGaussianModel.hpp"
#include "core/detail/FitStateView.hpp"
#include "core/detail/TransformedChange.hpp"
#include "core/detail/SecondStageContext.hpp"

namespace rhbm_gem::core::detail {

enum class BacktrackingStepStatus
{
    CandidateReady,
    InvalidCandidate,
    Exhausted
};

struct BacktrackingStep
{
    BacktrackingStepStatus status{ BacktrackingStepStatus::Exhausted };
    double factor{ 0.0 };
    std::size_t trial_number{ 0 };
};

class BacktrackingWorkspace
{
    const FitState & m_previous_state;
    double m_minimum_transformed_change{ 0.0 };
    double m_next_factor{ 0.5 };
    std::size_t m_trial_number{ 1 };
    std::vector<GaussianModel3D> m_previous_model_list{};
    std::vector<GaussianModel3D> m_endpoint_model_list{};
    std::vector<double> m_previous_shared_offset_list{};
    std::vector<double> m_endpoint_shared_offset_list{};
    FitStatePatch m_candidate_patch{};

public:
    BacktrackingWorkspace(
        const SecondStageContext & context,
        const FitState & previous_state,
        const auto & endpoint_state,
        const std::vector<std::size_t> & active_index_list,
        double minimum_transformed_change)
        : m_previous_state{ previous_state },
          m_minimum_transformed_change{ minimum_transformed_change }
    {
        if (!std::isfinite(m_minimum_transformed_change) || m_minimum_transformed_change < 0.0)
        {
            throw std::invalid_argument(
                "Local fitting backtracking minimum transformed change is invalid.");
        }
        m_candidate_patch.atom_index_list = active_index_list;
        std::ranges::sort(m_candidate_patch.atom_index_list);
        m_candidate_patch.atom_index_list.erase(
            std::ranges::unique(m_candidate_patch.atom_index_list).begin(),
            m_candidate_patch.atom_index_list.end());
        std::vector<std::size_t> group_id_by_atom_position;
        group_id_by_atom_position.reserve(m_candidate_patch.atom_index_list.size());
        m_previous_model_list.reserve(m_candidate_patch.atom_index_list.size());
        m_endpoint_model_list.reserve(m_candidate_patch.atom_index_list.size());
        m_candidate_patch.mdpde_list.reserve(m_candidate_patch.atom_index_list.size());
        for (const auto atom_index : m_candidate_patch.atom_index_list)
        {
            group_id_by_atom_position.emplace_back(context.at(atom_index).group_id);
            m_previous_model_list.emplace_back(previous_state.at(atom_index).mdpde.GetModel());
            const auto & endpoint_mdpde{ GetEndpointMdpde(endpoint_state, atom_index) };
            m_endpoint_model_list.emplace_back(endpoint_mdpde.GetModel());
            m_candidate_patch.mdpde_list.emplace_back(endpoint_mdpde);
        }
        m_previous_shared_offset_list =
            BuildGroupMedianOffsetList(group_id_by_atom_position, m_previous_model_list);
        m_endpoint_shared_offset_list =
            BuildGroupMedianOffsetList(group_id_by_atom_position, m_endpoint_model_list);
    }

    BacktrackingWorkspace(const BacktrackingWorkspace &) = delete;
    BacktrackingWorkspace & operator=(const BacktrackingWorkspace &) = delete;

    BacktrackingStep BuildNextCandidate()
    {
        const auto factor{ m_next_factor };
        if (!std::isfinite(factor) || factor < std::numeric_limits<double>::epsilon())
        {
            return BacktrackingStep{
                BacktrackingStepStatus::Exhausted,
                factor,
                m_trial_number
            };
        }
        m_next_factor *= 0.5;

        if (!BuildCandidate(factor))
        {
            return BacktrackingStep{
                BacktrackingStepStatus::InvalidCandidate,
                factor,
                m_trial_number
            };
        }
        const auto maximum_transformed_change{ GetMaximumTransformedChange() };
        if (maximum_transformed_change < m_minimum_transformed_change)
        {
            return BacktrackingStep{
                BacktrackingStepStatus::Exhausted,
                factor,
                m_trial_number
            };
        }
        m_trial_number++;
        return BacktrackingStep{
            BacktrackingStepStatus::CandidateReady,
            factor,
            m_trial_number
        };
    }

    FitStatePatch TakeCandidatePatch() { return std::move(m_candidate_patch); }
    const FitStatePatch & GetCandidatePatch() const { return m_candidate_patch; }

    PolishProvenance BuildCandidatePolishProvenance(
        const PolishProvenance & previous_provenance,
        const PolishProvenance & endpoint_provenance) const
    {
        if (previous_provenance.size() != m_previous_state.size() ||
            endpoint_provenance.size() != m_previous_state.size())
        {
            throw std::invalid_argument(
                "Local fitting backtracking provenance sizes are inconsistent.");
        }
        auto provenance{ previous_provenance };
        for (std::size_t i = 0; i < m_candidate_patch.atom_index_list.size(); i++)
        {
            if (HasMaterialChange(i))
            {
                provenance.at(m_candidate_patch.atom_index_list.at(i)) =
                    endpoint_provenance.at(m_candidate_patch.atom_index_list.at(i));
            }
        }
        return provenance;
    }

    PolishProvenance BuildActiveCandidatePolishProvenance(
        const PolishProvenance & previous_provenance,
        const PolishProvenance & endpoint_provenance) const
    {
        if (previous_provenance.size() != m_candidate_patch.atom_index_list.size() ||
            endpoint_provenance.size() != m_candidate_patch.atom_index_list.size())
        {
            throw std::invalid_argument(
                "Local fitting active backtracking provenance sizes are inconsistent.");
        }
        auto provenance{ previous_provenance };
        for (std::size_t i = 0; i < m_candidate_patch.atom_index_list.size(); i++)
        {
            if (HasMaterialChange(i))
            {
                provenance.at(i) = endpoint_provenance.at(i);
            }
        }
        return provenance;
    }

private:
    bool BuildCandidate(double factor)
    {
        std::vector<GaussianModel3D> candidate_model_list;
        if (!TryBuildSharedOffsetDampedModelList(
                m_previous_model_list,
                m_endpoint_model_list,
                m_previous_shared_offset_list,
                m_endpoint_shared_offset_list,
                factor,
                candidate_model_list))
        {
            return false;
        }

        for (std::size_t i = 0; i < m_candidate_patch.atom_index_list.size(); i++)
        {
            const auto endpoint_uncertainty{
                m_candidate_patch.mdpde_list.at(i).GetStandardDeviationModel()
            };
            m_candidate_patch.mdpde_list.at(i) =
                GaussianModel3DWithUncertainty{
                    candidate_model_list.at(i),
                    endpoint_uncertainty
                };
        }
        return true;
    }

    bool HasMaterialChange(std::size_t atom_position) const
    {
        const auto change{
            CalculateTransformedChange(
                m_candidate_patch.mdpde_list.at(atom_position).GetModel(),
                m_previous_model_list.at(atom_position))
        };
        return IsTransformedChangeMaterial(change, m_minimum_transformed_change);
    }

    double GetMaximumTransformedChange() const
    {
        double maximum_change{ 0.0 };
        for (std::size_t i = 0; i < m_candidate_patch.atom_index_list.size(); i++)
        {
            maximum_change = std::max(
                maximum_change,
                detail::GetMaximumTransformedChange(
                    CalculateTransformedChange(
                        m_candidate_patch.mdpde_list.at(i).GetModel(),
                        m_previous_model_list.at(i)).value_list));
        }
        return maximum_change;
    }

    static const GaussianModel3DWithUncertainty & GetEndpointMdpde(const auto & state, std::size_t atom_index)
    {
        if constexpr (requires { state.GetMdpde(atom_index); })
        {
            return state.GetMdpde(atom_index);
        }
        else
        {
            return state.at(atom_index).mdpde;
        }
    }
};

} // namespace rhbm_gem::core::detail
