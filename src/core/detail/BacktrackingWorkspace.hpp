#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "core/detail/LocalFittingGroupMedian.hpp"
#include "core/detail/FitStateView.hpp"
#include "core/detail/LocalFittingTransformedChange.hpp"
#include "core/detail/SecondStageLocalFittingContext.hpp"

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
    double maximum_transformed_change{ std::numeric_limits<double>::infinity() };
    const FitStatePatch * candidate_patch{ nullptr };

    bool IsCandidateReady() const
    {
        return status == BacktrackingStepStatus::CandidateReady;
    }
};

struct BacktrackingWorkspace
{
    template <typename EndpointState>
    BacktrackingWorkspace(
        const SecondStageLocalFittingContext & context,
        const FitState & previous_state,
        const EndpointState & endpoint_state,
        const std::vector<std::size_t> & active_index_list,
        double minimum_transformed_change)
        : m_previous_state{ &previous_state },
          m_active_index_list{ active_index_list },
          m_minimum_transformed_change{ minimum_transformed_change }
    {
        if (!std::isfinite(m_minimum_transformed_change) || m_minimum_transformed_change < 0.0)
        {
            throw std::invalid_argument(
                "Local fitting backtracking minimum transformed change is invalid.");
        }
        std::sort(m_active_index_list.begin(), m_active_index_list.end());
        m_active_index_list.erase(
            std::unique(m_active_index_list.begin(), m_active_index_list.end()),
            m_active_index_list.end());
        m_group_key_by_atom_position.reserve(m_active_index_list.size());
        m_previous_model_list.reserve(m_active_index_list.size());
        m_endpoint_model_list.reserve(m_active_index_list.size());
        m_endpoint_uncertainty_list.reserve(m_active_index_list.size());
        for (const auto atom_index : m_active_index_list)
        {
            m_group_key_by_atom_position.emplace_back(context.at(atom_index).group_key);
            m_previous_model_list.emplace_back(previous_state.at(atom_index).mdpde.GetModel());
            m_endpoint_model_list.emplace_back(GetEndpointModel(endpoint_state, atom_index));
            m_endpoint_uncertainty_list.emplace_back(
                GetEndpointMdpde(endpoint_state, atom_index).GetStandardDeviationModel());
        }
        m_previous_shared_offset_model_list =
            BuildLocalFittingGroupMedianModelList(
                m_group_key_by_atom_position,
                m_previous_model_list);
        m_endpoint_shared_offset_model_list =
            BuildLocalFittingGroupMedianModelList(
                m_group_key_by_atom_position,
                m_endpoint_model_list);
        m_candidate_patch.atom_index_list = m_active_index_list;
        m_candidate_patch.mdpde_list.reserve(m_active_index_list.size());
        for (const auto atom_index : m_active_index_list)
        {
            m_candidate_patch.mdpde_list.emplace_back(GetEndpointMdpde(endpoint_state, atom_index));
        }
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
                m_trial_number,
                std::numeric_limits<double>::infinity(),
                nullptr
            };
        }
        m_next_factor *= 0.5;

        const auto * candidate_patch{ BuildCandidate(factor) };
        if (candidate_patch == nullptr)
        {
            return BacktrackingStep{
                BacktrackingStepStatus::InvalidCandidate,
                factor,
                m_trial_number,
                std::numeric_limits<double>::infinity(),
                nullptr
            };
        }
        const auto maximum_transformed_change{
            GetMaximumTransformedChange()
        };
        if (maximum_transformed_change < m_minimum_transformed_change)
        {
            return BacktrackingStep{
                BacktrackingStepStatus::Exhausted,
                factor,
                m_trial_number,
                maximum_transformed_change,
                nullptr
            };
        }
        m_trial_number++;
        return BacktrackingStep{
            BacktrackingStepStatus::CandidateReady,
            factor,
            m_trial_number,
            maximum_transformed_change,
            candidate_patch
        };
    }

    template <typename Evaluator>
    BacktrackingStep FindAcceptedCandidate(Evaluator && evaluator)
    {
        while (true)
        {
            const auto step{ BuildNextCandidate() };
            if (!step.IsCandidateReady()) return step;
            if (evaluator(step)) return step;
        }
    }

    FitStatePatch TakeCandidatePatch()
    {
        return std::move(m_candidate_patch);
    }

    FitState MaterializeCandidateState() const
    {
        auto state{ *m_previous_state };
        m_candidate_patch.ApplyTo(state);
        return state;
    }

    PolishProvenance BuildCandidatePolishProvenance(
        const PolishProvenance & previous_provenance,
        const PolishProvenance & endpoint_provenance) const
    {
        if (previous_provenance.size() != m_previous_state->size() ||
            endpoint_provenance.size() != m_previous_state->size())
        {
            throw std::invalid_argument(
                "Local fitting backtracking provenance sizes are inconsistent.");
        }
        auto provenance{ previous_provenance };
        for (std::size_t atom_position = 0; atom_position < m_active_index_list.size(); atom_position++)
        {
            if (HasMaterialChange(atom_position))
            {
                provenance.at(m_active_index_list.at(atom_position)) =
                    endpoint_provenance.at(m_active_index_list.at(atom_position));
            }
        }
        return provenance;
    }

    PolishProvenance BuildActiveCandidatePolishProvenance(
        const PolishProvenance & previous_provenance,
        const PolishProvenance & endpoint_provenance) const
    {
        if (previous_provenance.size() != m_active_index_list.size() ||
            endpoint_provenance.size() != m_active_index_list.size())
        {
            throw std::invalid_argument(
                "Local fitting active backtracking provenance sizes are inconsistent.");
        }
        auto provenance{ previous_provenance };
        for (std::size_t atom_position = 0; atom_position < m_active_index_list.size(); atom_position++)
        {
            if (HasMaterialChange(atom_position))
            {
                provenance.at(atom_position) = endpoint_provenance.at(atom_position);
            }
        }
        return provenance;
    }

private:
    const FitStatePatch * BuildCandidate(double factor)
    {
        const auto candidate_model_list{
            BuildLocalFittingSharedOffsetDampedModelList(
                m_previous_model_list,
                m_endpoint_model_list,
                m_previous_shared_offset_model_list,
                m_endpoint_shared_offset_model_list,
                factor)
        };
        if (!candidate_model_list.has_value()) return nullptr;

        for (std::size_t atom_position = 0; atom_position < m_active_index_list.size(); atom_position++)
        {
            m_candidate_patch.mdpde_list.at(atom_position) =
                GaussianModel3DWithUncertainty{
                    candidate_model_list->at(atom_position),
                    m_endpoint_uncertainty_list.at(atom_position)
                };
        }
        return &m_candidate_patch;
    }

    bool HasMaterialChange(std::size_t atom_position) const
    {
        const auto change{
            CalculateLocalFittingTransformedChange(
                m_candidate_patch.mdpde_list.at(atom_position).GetModel(),
                m_previous_model_list.at(atom_position))
        };
        return IsLocalFittingTransformedChangeMaterial(change, m_minimum_transformed_change);
    }

    double GetMaximumTransformedChange() const
    {
        std::vector<algorithm::ParameterChange> change_list;
        change_list.reserve(m_active_index_list.size());
        for (std::size_t atom_position = 0; atom_position < m_active_index_list.size(); atom_position++)
        {
            change_list.emplace_back(
                CalculateLocalFittingTransformedChange(
                    m_candidate_patch.mdpde_list.at(atom_position).GetModel(),
                    m_previous_model_list.at(atom_position)));
        }
        return GetMaximumLocalFittingTransformedChange(change_list);
    }

    template <typename EndpointState>
    static const GaussianModel3D & GetEndpointModel(const EndpointState & state, std::size_t atom_index)
    {
        if constexpr (std::is_same_v<std::decay_t<EndpointState>, FitStateView>)
        {
            return state.GetModel(atom_index);
        }
        else
        {
            return state.at(atom_index).mdpde.GetModel();
        }
    }

    template <typename EndpointState>
    static const GaussianModel3DWithUncertainty & GetEndpointMdpde(const EndpointState & state, std::size_t atom_index)
    {
        if constexpr (std::is_same_v<std::decay_t<EndpointState>, FitStateView>)
        {
            return state.GetMdpde(atom_index);
        }
        else
        {
            return state.at(atom_index).mdpde;
        }
    }

    const FitState * m_previous_state{ nullptr };
    std::vector<std::size_t> m_active_index_list{};
    double m_minimum_transformed_change{ 0.0 };
    double m_next_factor{ 0.5 };
    std::size_t m_trial_number{ 1 };
    std::vector<GroupKey> m_group_key_by_atom_position{};
    std::vector<GaussianModel3D> m_previous_model_list{};
    std::vector<GaussianModel3D> m_endpoint_model_list{};
    std::vector<GaussianModel3DUncertainty> m_endpoint_uncertainty_list{};
    std::vector<GaussianModel3D> m_previous_shared_offset_model_list{};
    std::vector<GaussianModel3D> m_endpoint_shared_offset_model_list{};
    FitStatePatch m_candidate_patch{};
};

} // namespace rhbm_gem::core::detail
