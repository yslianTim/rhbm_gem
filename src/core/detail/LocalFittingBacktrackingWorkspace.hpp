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

struct LocalFittingBacktrackingWorkspace
{
    template <typename EndpointState>
    LocalFittingBacktrackingWorkspace(
        const SecondStageLocalFittingContext & context,
        const LocalFittingState & previous_state,
        const EndpointState & endpoint_state,
        const std::vector<std::size_t> & active_index_list)
        : m_previous_state{ &previous_state },
          m_active_index_list{ active_index_list }
    {
        m_group_key_by_atom_position.reserve(m_active_index_list.size());
        m_previous_model_list.reserve(m_active_index_list.size());
        m_endpoint_model_list.reserve(m_active_index_list.size());
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
        m_candidate_patch.mdpde_list.resize(m_active_index_list.size());
    }

    LocalFittingBacktrackingWorkspace(
        const LocalFittingBacktrackingWorkspace &) = delete;
    LocalFittingBacktrackingWorkspace & operator=(
        const LocalFittingBacktrackingWorkspace &) = delete;

    const LocalFittingStatePatch * BuildCandidate(double factor)
    {
        if (!std::isfinite(factor) || factor < 0.0 || factor > 1.0)
        {
            throw std::invalid_argument(
                "Local fitting objective backtracking factor must be in [0, 1].");
        }
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

    LocalFittingStatePatch TakeCandidatePatch()
    {
        return std::move(m_candidate_patch);
    }

    LocalFittingState TakeCandidateState()
    {
        auto state{ *m_previous_state };
        m_candidate_patch.ApplyTo(state);
        return state;
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

private:
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
