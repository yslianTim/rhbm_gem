#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "core/detail/LocalFittingGroupMedian.hpp"
#include "core/detail/LocalFittingStateView.hpp"
#include "core/detail/LocalFittingTransformedChange.hpp"
#include "core/detail/SecondStageLocalFittingContext.hpp"

namespace rhbm_gem::core::detail {

enum class LocalFittingBacktrackingStepStatus
{
    CandidateReady,
    InvalidCandidate,
    Exhausted
};

struct LocalFittingBacktrackingStep
{
    LocalFittingBacktrackingStepStatus status{
        LocalFittingBacktrackingStepStatus::Exhausted
    };
    double factor{ 0.0 };
    std::size_t trial_number{ 0 };
    double maximum_transformed_change{
        std::numeric_limits<double>::infinity()
    };
    const LocalFittingStatePatch * candidate_patch{ nullptr };

    bool IsCandidateReady() const
    {
        return status == LocalFittingBacktrackingStepStatus::CandidateReady;
    }
};

struct LocalFittingBacktrackingWorkspace
{
    template <typename EndpointState>
    LocalFittingBacktrackingWorkspace(
        const SecondStageLocalFittingContext & context,
        const LocalFittingState & previous_state,
        const EndpointState & endpoint_state,
        const std::vector<std::size_t> & active_index_list,
        double minimum_transformed_change)
        : m_previous_state{ &previous_state },
          m_active_index_list{ active_index_list },
          m_minimum_transformed_change{ minimum_transformed_change }
    {
        if (!std::isfinite(m_minimum_transformed_change) ||
            m_minimum_transformed_change < 0.0)
        {
            throw std::invalid_argument(
                "Local fitting backtracking minimum transformed change is invalid.");
        }
        std::sort(
            m_active_index_list.begin(),
            m_active_index_list.end());
        m_active_index_list.erase(
            std::unique(
                m_active_index_list.begin(),
                m_active_index_list.end()),
            m_active_index_list.end());
        m_group_key_by_atom_position.reserve(m_active_index_list.size());
        m_previous_model_list.reserve(m_active_index_list.size());
        m_endpoint_model_list.reserve(m_active_index_list.size());
        m_endpoint_uncertainty_list.reserve(m_active_index_list.size());
        m_previous_transformed_estimation_list.reserve(
            m_active_index_list.size());
        for (const auto atom_index : m_active_index_list)
        {
            m_group_key_by_atom_position.emplace_back(
                context.at(atom_index).group_key);
            m_previous_model_list.emplace_back(
                previous_state.at(atom_index).mdpde.GetModel());
            m_endpoint_model_list.emplace_back(
                GetEndpointModel(endpoint_state, atom_index));
            m_endpoint_uncertainty_list.emplace_back(
                GetEndpointMdpde(endpoint_state, atom_index)
                    .GetStandardDeviationModel());
            m_previous_transformed_estimation_list.emplace_back(
                EncodeLocalFittingTransformedCoordinates(
                    previous_state.at(atom_index).mdpde.GetModel()));
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
            m_candidate_patch.mdpde_list.emplace_back(
                GetEndpointMdpde(endpoint_state, atom_index));
        }
    }

    LocalFittingBacktrackingWorkspace(
        const LocalFittingBacktrackingWorkspace &) = delete;
    LocalFittingBacktrackingWorkspace & operator=(
        const LocalFittingBacktrackingWorkspace &) = delete;

    LocalFittingBacktrackingStep BuildNextCandidate()
    {
        const auto factor{ m_next_factor };
        if (!std::isfinite(factor) ||
            factor < std::numeric_limits<double>::epsilon())
        {
            return LocalFittingBacktrackingStep{
                LocalFittingBacktrackingStepStatus::Exhausted,
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
            return LocalFittingBacktrackingStep{
                LocalFittingBacktrackingStepStatus::InvalidCandidate,
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
            return LocalFittingBacktrackingStep{
                LocalFittingBacktrackingStepStatus::Exhausted,
                factor,
                m_trial_number,
                maximum_transformed_change,
                nullptr
            };
        }
        m_trial_number++;
        return LocalFittingBacktrackingStep{
            LocalFittingBacktrackingStepStatus::CandidateReady,
            factor,
            m_trial_number,
            maximum_transformed_change,
            candidate_patch
        };
    }

    template <typename Evaluator>
    LocalFittingBacktrackingStep FindAcceptedCandidate(Evaluator && evaluator)
    {
        while (true)
        {
            const auto step{ BuildNextCandidate() };
            if (!step.IsCandidateReady()) return step;
            if (evaluator(step)) return step;
        }
    }

    LocalFittingStatePatch TakeCandidatePatch()
    {
        return std::move(m_candidate_patch);
    }

    LocalFittingState MaterializeCandidateState() const
    {
        auto state{ *m_previous_state };
        m_candidate_patch.ApplyTo(state);
        return state;
    }

    LocalFittingPolishProvenance BuildCandidatePolishProvenance(
        const LocalFittingPolishProvenance & previous_provenance,
        const LocalFittingPolishProvenance & endpoint_provenance) const
    {
        if (previous_provenance.size() != m_previous_state->size() ||
            endpoint_provenance.size() != m_previous_state->size())
        {
            throw std::invalid_argument(
                "Local fitting backtracking provenance sizes are inconsistent.");
        }
        auto provenance{ previous_provenance };
        for (std::size_t atom_position = 0;
            atom_position < m_active_index_list.size();
            atom_position++)
        {
            if (HasMaterialChange(atom_position))
            {
                provenance.at(m_active_index_list.at(atom_position)) =
                    endpoint_provenance.at(
                        m_active_index_list.at(atom_position));
            }
        }
        return provenance;
    }

    LocalFittingPolishProvenance BuildActiveCandidatePolishProvenance(
        const LocalFittingPolishProvenance & previous_provenance,
        const LocalFittingPolishProvenance & endpoint_provenance) const
    {
        if (previous_provenance.size() != m_active_index_list.size() ||
            endpoint_provenance.size() != m_active_index_list.size())
        {
            throw std::invalid_argument(
                "Local fitting active backtracking provenance sizes are inconsistent.");
        }
        auto provenance{ previous_provenance };
        for (std::size_t atom_position = 0;
            atom_position < m_active_index_list.size();
            atom_position++)
        {
            if (HasMaterialChange(atom_position))
            {
                provenance.at(atom_position) = endpoint_provenance.at(
                    atom_position);
            }
        }
        return provenance;
    }

private:
    const LocalFittingStatePatch * BuildCandidate(double factor)
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

        for (std::size_t atom_position = 0;
            atom_position < m_active_index_list.size();
            atom_position++)
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
        const auto current_coordinates{
            EncodeLocalFittingTransformedCoordinates(
                m_candidate_patch.mdpde_list.at(atom_position).GetModel())
        };
        const auto & previous_coordinates{
            m_previous_transformed_estimation_list.at(atom_position)
        };
        if (!current_coordinates.has_value() ||
            !previous_coordinates.has_value())
        {
            return false;
        }
        return ((*current_coordinates - *previous_coordinates).array().abs() >=
            m_minimum_transformed_change).any();
    }

    double GetMaximumTransformedChange() const
    {
        std::array<double, kTransformedChangeSize> maximum_list{};
        for (std::size_t atom_position = 0;
            atom_position < m_active_index_list.size();
            atom_position++)
        {
            const auto current_coordinates{
                EncodeLocalFittingTransformedCoordinates(
                    m_candidate_patch.mdpde_list.at(atom_position).GetModel())
            };
            const auto & previous_coordinates{
                m_previous_transformed_estimation_list.at(atom_position)
            };
            if (!current_coordinates.has_value() ||
                !previous_coordinates.has_value())
            {
                return std::numeric_limits<double>::infinity();
            }
            for (std::size_t parameter_index = 0;
                parameter_index < kTransformedChangeSize;
                parameter_index++)
            {
                const auto parameter_value{
                    std::abs(
                        (*current_coordinates)(
                            static_cast<Eigen::Index>(parameter_index)) -
                        (*previous_coordinates)(
                            static_cast<Eigen::Index>(parameter_index)))
                };
                if (!std::isfinite(parameter_value))
                {
                    return std::numeric_limits<double>::infinity();
                }
                maximum_list.at(parameter_index) = std::max(
                    maximum_list.at(parameter_index),
                    parameter_value);
            }
        }
        return *std::max_element(maximum_list.begin(), maximum_list.end());
    }

    template <typename EndpointState>
    static const GaussianModel3D & GetEndpointModel(
        const EndpointState & state,
        std::size_t atom_index)
    {
        if constexpr (std::is_same_v<
                          std::decay_t<EndpointState>,
                          LocalFittingStateView>)
        {
            return state.GetModel(atom_index);
        }
        else
        {
            return state.at(atom_index).mdpde.GetModel();
        }
    }

    template <typename EndpointState>
    static const GaussianModel3DWithUncertainty & GetEndpointMdpde(
        const EndpointState & state,
        std::size_t atom_index)
    {
        if constexpr (std::is_same_v<
                          std::decay_t<EndpointState>,
                          LocalFittingStateView>)
        {
            return state.GetMdpde(atom_index);
        }
        else
        {
            return state.at(atom_index).mdpde;
        }
    }

    const LocalFittingState * m_previous_state{ nullptr };
    std::vector<std::size_t> m_active_index_list{};
    double m_minimum_transformed_change{ 0.0 };
    double m_next_factor{ 0.5 };
    std::size_t m_trial_number{ 1 };
    std::vector<GroupKey> m_group_key_by_atom_position{};
    std::vector<GaussianModel3D> m_previous_model_list{};
    std::vector<GaussianModel3D> m_endpoint_model_list{};
    std::vector<GaussianModel3DUncertainty> m_endpoint_uncertainty_list{};
    std::vector<std::optional<Eigen::Vector3d>>
        m_previous_transformed_estimation_list{};
    std::vector<GaussianModel3D> m_previous_shared_offset_model_list{};
    std::vector<GaussianModel3D> m_endpoint_shared_offset_model_list{};
    LocalFittingStatePatch m_candidate_patch{};
};

} // namespace rhbm_gem::core::detail
