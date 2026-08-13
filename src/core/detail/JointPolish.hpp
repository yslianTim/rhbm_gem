#pragma once

#include "core/detail/JointOffset.hpp"
#include "core/detail/LocalFittingCandidateEvaluationOverlay.hpp"
#include "core/detail/LocalFittingRobustScale.hpp"
#include "core/detail/LocalFittingTransformedChange.hpp"
#include "core/detail/LocalFittingTrustRegion.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>
#include <rhbm_gem/utils/algorithm/WeightedRidgeSolver.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem::core::detail {

constexpr std::size_t kJointPolishShapeParameterSize{ 2 };

struct JointPolishParameterization
{
    std::vector<std::size_t> group_position_by_atom{};
    std::vector<std::vector<std::size_t>> atom_position_list_by_group{};
    Eigen::VectorXd seed_parameter{};

    std::size_t AtomCount() const
    {
        return group_position_by_atom.size();
    }

    std::size_t GroupCount() const
    {
        return atom_position_list_by_group.size();
    }

    Eigen::Index ParameterCount() const
    {
        return seed_parameter.size();
    }

    Eigen::Index ShapeColumn(
        std::size_t atom_position,
        std::size_t shape_parameter_index) const
    {
        return static_cast<Eigen::Index>(
            atom_position * kJointPolishShapeParameterSize +
            shape_parameter_index);
    }

    Eigen::Index OffsetColumn(std::size_t atom_position) const
    {
        return static_cast<Eigen::Index>(
            AtomCount() * kJointPolishShapeParameterSize +
            group_position_by_atom.at(atom_position));
    }

    std::optional<std::vector<GaussianModel3D>> DecodeModels(
        const Eigen::VectorXd & direction,
        double damping) const
    {
        if (direction.size() != ParameterCount() ||
            !std::isfinite(damping) || damping < 0.0 || damping > 1.0)
        {
            return std::nullopt;
        }
        const Eigen::VectorXd parameter{
            seed_parameter + damping * direction
        };
        if (!parameter.allFinite()) return std::nullopt;

        std::vector<GaussianModel3D> model_list;
        model_list.reserve(AtomCount());
        for (std::size_t atom_position = 0;
            atom_position < AtomCount();
            atom_position++)
        {
            const Eigen::Vector3d shape_coordinates{
                parameter(ShapeColumn(atom_position, 0)),
                parameter(ShapeColumn(atom_position, 1)),
                0.0
            };
            const auto shape_model{
                DecodeLocalFittingTransformedCoordinates(shape_coordinates)
            };
            if (!shape_model.has_value()) return std::nullopt;
            const auto model{
                shape_model->WithOffset(parameter(OffsetColumn(atom_position)))
            };
            if (!EncodeLocalFittingTransformedCoordinates(model).has_value())
            {
                return std::nullopt;
            }
            model_list.emplace_back(model);
        }
        return model_list;
    }
};

inline std::optional<JointPolishParameterization>
BuildJointPolishParameterization(
    const std::vector<GroupKey> & group_key_by_atom_position,
    const std::vector<GaussianModel3D> & base_model_list)
{
    if (group_key_by_atom_position.empty() ||
        group_key_by_atom_position.size() != base_model_list.size())
    {
        return std::nullopt;
    }

    std::map<GroupKey, std::size_t> group_position_by_key;
    for (const auto group_key : group_key_by_atom_position)
    {
        group_position_by_key.emplace(group_key, 0);
    }
    std::size_t group_position{ 0 };
    for (auto & [group_key, position] : group_position_by_key)
    {
        static_cast<void>(group_key);
        position = group_position++;
    }

    JointPolishParameterization parameterization;
    parameterization.group_position_by_atom.resize(base_model_list.size());
    parameterization.atom_position_list_by_group.resize(
        group_position_by_key.size());
    parameterization.seed_parameter = Eigen::VectorXd::Zero(
        static_cast<Eigen::Index>(
            base_model_list.size() *
                kJointPolishShapeParameterSize +
            group_position_by_key.size()));
    std::vector<std::vector<double>> offset_list_by_group(
        group_position_by_key.size());

    for (std::size_t atom_position = 0;
        atom_position < base_model_list.size();
        atom_position++)
    {
        const auto transformed{
            EncodeLocalFittingTransformedCoordinates(
                base_model_list.at(atom_position))
        };
        if (!transformed.has_value()) return std::nullopt;
        const auto atom_group_position{
            group_position_by_key.at(group_key_by_atom_position.at(atom_position))
        };
        parameterization.group_position_by_atom.at(atom_position) =
            atom_group_position;
        parameterization.atom_position_list_by_group.at(atom_group_position)
            .emplace_back(atom_position);
        parameterization.seed_parameter(
            parameterization.ShapeColumn(atom_position, 0)) =
            (*transformed)(static_cast<Eigen::Index>(
                kLogPeakHeightChangeIndex));
        parameterization.seed_parameter(
            parameterization.ShapeColumn(atom_position, 1)) =
            (*transformed)(static_cast<Eigen::Index>(
                kLogWidthChangeIndex));
        offset_list_by_group.at(atom_group_position).emplace_back(
            base_model_list.at(atom_position).GetOffset());
    }

    for (std::size_t current_group_position = 0;
        current_group_position < offset_list_by_group.size();
        current_group_position++)
    {
        auto & offset_list{ offset_list_by_group.at(current_group_position) };
        if (offset_list.empty()) return std::nullopt;
        std::sort(offset_list.begin(), offset_list.end());
        const auto middle{ offset_list.size() / 2 };
        const auto median{
            offset_list.size() % 2 == 0 ?
                0.5 * offset_list.at(middle - 1) +
                    0.5 * offset_list.at(middle) :
                offset_list.at(middle)
        };
        if (!std::isfinite(median)) return std::nullopt;
        const auto representative_atom_position{
            parameterization.atom_position_list_by_group
                .at(current_group_position).front()
        };
        parameterization.seed_parameter(
            parameterization.OffsetColumn(representative_atom_position)) = median;
    }
    if (!parameterization.seed_parameter.allFinite()) return std::nullopt;
    return parameterization;
}

struct SharedOffsetResponse
{
    double response{ 0.0 };
    Eigen::Vector2d shape_jacobian{ Eigen::Vector2d::Zero() };
    double offset_jacobian{ 0.0 };
};

inline std::optional<SharedOffsetResponse>
EvaluateSharedOffsetResponse(
    const GaussianModel3D & model,
    double distance)
{
    if (!EncodeLocalFittingTransformedCoordinates(model).has_value() ||
        !std::isfinite(distance) || distance < 0.0)
    {
        return std::nullopt;
    }

    const auto width{ model.GetWidth() };
    const auto signal{ model.SignalAtDistance(distance) };
    const auto offset_basis{ model.OffsetBasisAtDistance(distance) };
    const auto response{ signal + model.GetOffset() * offset_basis };
    if (!std::isfinite(signal) || !std::isfinite(offset_basis) ||
        !std::isfinite(response))
    {
        return std::nullopt;
    }

    const double center_offset_basis_scale{ std::sqrt(2.0 / M_PI) };
    const auto normalized_distance{ distance / width };
    auto log_width_derivative{
        signal * normalized_distance * normalized_distance
    };
    if (distance < 1.0e-5)
    {
        log_width_derivative -=
            model.GetOffset() * center_offset_basis_scale / width;
    }
    else
    {
        const auto exponent{
            -0.5 * normalized_distance * normalized_distance
        };
        log_width_derivative -= model.GetOffset() *
            center_offset_basis_scale / width * std::exp(exponent);
    }

    Eigen::Vector2d shape_jacobian{
        signal,
        log_width_derivative
    };
    if (!shape_jacobian.allFinite()) return std::nullopt;

    return SharedOffsetResponse{
        response,
        shape_jacobian,
        offset_basis
    };
}

struct TransformedResponse
{
    double response{ 0.0 };
    Eigen::Vector3d jacobian{ Eigen::Vector3d::Zero() };
};

struct TransformedModelInvariants
{
    GaussianModel3D model{};
    Eigen::Vector3d transformed{ Eigen::Vector3d::Zero() };
    double peak_height{ 0.0 };
};

inline std::optional<TransformedModelInvariants>
BuildTransformedModelInvariants(
    const GaussianModel3D & model)
{
    const auto transformed{ EncodeLocalFittingTransformedCoordinates(model) };
    if (!transformed.has_value()) return std::nullopt;

    const auto peak_height{
        std::exp((*transformed)(static_cast<Eigen::Index>(
            kLogPeakHeightChangeIndex)))
    };
    if (!std::isfinite(peak_height)) return std::nullopt;

    return TransformedModelInvariants{
        model,
        *transformed,
        peak_height
    };
}

inline std::optional<SharedOffsetResponse>
EvaluateSharedOffsetResponse(
    const TransformedModelInvariants & invariants,
    double distance)
{
    if (!std::isfinite(distance) || distance < 0.0) return std::nullopt;

    const auto evaluation{ invariants.model.EvaluateAtDistance(distance) };
    if (!std::isfinite(evaluation.signal) ||
        !std::isfinite(evaluation.offset_basis) ||
        !std::isfinite(evaluation.response))
    {
        return std::nullopt;
    }

    const auto width{ invariants.model.GetWidth() };
    const double center_offset_basis_scale{ std::sqrt(2.0 / M_PI) };
    const auto normalized_distance{ distance / width };
    auto log_width_derivative{
        evaluation.signal * normalized_distance * normalized_distance
    };
    if (distance < 1.0e-5)
    {
        log_width_derivative -=
            invariants.model.GetOffset() * center_offset_basis_scale / width;
    }
    else
    {
        const auto exponent{
            -0.5 * normalized_distance * normalized_distance
        };
        log_width_derivative -= invariants.model.GetOffset() *
            center_offset_basis_scale / width * std::exp(exponent);
    }

    const Eigen::Vector2d shape_jacobian{
        evaluation.signal,
        log_width_derivative
    };
    if (!shape_jacobian.allFinite()) return std::nullopt;

    return SharedOffsetResponse{
        evaluation.response,
        shape_jacobian,
        evaluation.offset_basis
    };
}

inline std::optional<TransformedResponse>
EvaluateTransformedResponse(
    const TransformedModelInvariants & invariants,
    double distance)
{
    const auto shared_offset_evaluation{
        EvaluateSharedOffsetResponse(invariants, distance)
    };
    if (!shared_offset_evaluation.has_value())
    {
        return std::nullopt;
    }

    const auto & model{ invariants.model };
    const auto width{ model.GetWidth() };
    const double center_offset_basis_scale{ std::sqrt(2.0 / M_PI) };
    Eigen::Vector3d jacobian{ Eigen::Vector3d::Zero() };
    jacobian(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex)) =
        shared_offset_evaluation->shape_jacobian(0) +
        model.GetOffset() * shared_offset_evaluation->offset_jacobian;
    jacobian(static_cast<Eigen::Index>(kLogWidthChangeIndex)) =
        shared_offset_evaluation->shape_jacobian(1) +
        model.GetOffset() * shared_offset_evaluation->offset_jacobian;
    jacobian(static_cast<Eigen::Index>(kOffsetToPeakRatioChangeIndex)) =
        invariants.peak_height * width *
        shared_offset_evaluation->offset_jacobian /
        center_offset_basis_scale;
    if (!jacobian.allFinite())
    {
        return std::nullopt;
    }

    return TransformedResponse{
        shared_offset_evaluation->response,
        jacobian
    };
}

inline std::optional<TransformedResponse>
EvaluateTransformedResponse(
    const GaussianModel3D & model,
    double distance)
{
    const auto invariants{ BuildTransformedModelInvariants(model) };
    if (!invariants.has_value()) return std::nullopt;
    return EvaluateTransformedResponse(*invariants, distance);
}

constexpr double kJointPolishTransformedChangeTolerance{ 1.0e-4 };

inline std::optional<Eigen::VectorXd> BuildJointPolishDirection(
    const SecondStageLocalFittingContext & context,
    const LocalFittingStateView & base_state,
    const std::vector<GaussianModel3D> & seed_model_list,
    const LocalFittingClusterKey & key,
    const std::vector<LocalFittingObjectiveSampleRef> & sample_ref_list,
    const std::vector<double> & ridge_multiplier_list,
    const JointPolishParameterization & parameterization,
    ReusableWeightedRidgeSolver & reusable_solver)
{
    if (key.empty() || sample_ref_list.empty() ||
        seed_model_list.size() != key.size() ||
        parameterization.AtomCount() != key.size())
    {
        return std::nullopt;
    }

    const auto column_count{ parameterization.ParameterCount() };
    std::unordered_map<std::size_t, std::size_t>
        local_position_by_atom_index;
    local_position_by_atom_index.reserve(key.size());
    std::unordered_map<GroupKey, std::size_t> local_position_by_group_key;
    for (std::size_t local_position = 0;
        local_position < key.size();
        local_position++)
    {
        const auto atom_index{ key.at(local_position) };
        local_position_by_atom_index.emplace(atom_index, local_position);
        local_position_by_group_key.emplace(
            context.at(atom_index).group_key,
            local_position);
    }
    auto selected_snapshot{ BuildFittedGaussianSnapshot(base_state) };
    for (std::size_t local_position = 0;
        local_position < key.size();
        local_position++)
    {
        selected_snapshot.at(key.at(local_position)) =
            seed_model_list.at(local_position);
    }
    const auto model_snapshot{
        BuildSecondStageModelSnapshot(context, std::move(selected_snapshot))
    };

    std::vector<Eigen::Triplet<double>> triplet_list;
    std::vector<double> residual_list;
    residual_list.reserve(sample_ref_list.size());
    for (const auto & sample_ref : sample_ref_list)
    {
        const auto & atom_context{ context.at(sample_ref.atom_index) };
        const auto & sample{
            atom_context.raw_sampling_entries.at(sample_ref.sample_index)
        };
        if (!std::isfinite(static_cast<double>(sample.response)))
        {
            return std::nullopt;
        }

        const auto row_index{ static_cast<Eigen::Index>(residual_list.size()) };
        double predicted_response{ 0.0 };
        const auto append_model = [&](std::size_t atom_index, double distance) -> bool
        {
            const auto local_position_iter{
                local_position_by_atom_index.find(atom_index)
            };
            const auto local_position_value{
                local_position_iter == local_position_by_atom_index.end() ?
                    -1 : static_cast<int>(local_position_iter->second)
            };
            const auto & model{
                local_position_value >= 0 ?
                    seed_model_list.at(static_cast<std::size_t>(
                        local_position_value)) :
                    base_state.GetModel(atom_index)
            };
            const auto evaluation{
                EvaluateSharedOffsetResponse(model, distance)
            };
            if (!evaluation.has_value()) return false;
            predicted_response += evaluation->response;

            if (local_position_value < 0) return true;
            const auto local_position{
                static_cast<std::size_t>(local_position_value)
            };
            for (std::size_t parameter_index = 0;
                parameter_index < kJointPolishShapeParameterSize;
                parameter_index++)
            {
                const auto column_index{
                    parameterization.ShapeColumn(
                        local_position,
                        parameter_index)
                };
                const auto derivative{
                    evaluation->shape_jacobian(
                        static_cast<Eigen::Index>(parameter_index))
                };
                if (std::abs(derivative) <=
                    std::numeric_limits<double>::epsilon())
                {
                    continue;
                }
                triplet_list.emplace_back(row_index, column_index, derivative);
            }
            if (std::abs(evaluation->offset_jacobian) >
                std::numeric_limits<double>::epsilon())
            {
                triplet_list.emplace_back(
                    row_index,
                    parameterization.OffsetColumn(local_position),
                    evaluation->offset_jacobian);
            }
            return true;
        };
        const auto append_unselected_model = [&](
            std::size_t contributor_index,
            double distance) -> bool
        {
            const auto & contributor{
                context.unselected_atom_list.at(contributor_index)
            };
            const auto & model{
                model_snapshot.unselected.at(contributor_index)
            };
            const auto evaluation{
                EvaluateSharedOffsetResponse(model, distance)
            };
            if (!evaluation.has_value()) return false;
            predicted_response += evaluation->response;

            const auto local_position_iter{
                local_position_by_group_key.find(contributor.group_key)
            };
            if (local_position_iter == local_position_by_group_key.end() ||
                std::abs(evaluation->offset_jacobian) <=
                    std::numeric_limits<double>::epsilon())
            {
                return true;
            }
            triplet_list.emplace_back(
                row_index,
                parameterization.OffsetColumn(local_position_iter->second),
                evaluation->offset_jacobian);
            return true;
        };

        if (!append_model(
                sample_ref.atom_index,
                static_cast<double>(sample.point.distance)))
        {
            return std::nullopt;
        }
        for (auto neighbor_iter =
                atom_context.NeighborBegin(sample_ref.sample_index);
            neighbor_iter != atom_context.NeighborEnd(sample_ref.sample_index);
            ++neighbor_iter)
        {
            const auto & neighbor_sample{ *neighbor_iter };
            const auto appended{
                neighbor_sample.is_selected ?
                    append_model(
                        neighbor_sample.atom_index,
                        neighbor_sample.distance) :
                    append_unselected_model(
                        neighbor_sample.atom_index,
                        neighbor_sample.distance)
            };
            if (!appended) return std::nullopt;
        }

        const auto residual{
            static_cast<double>(sample.response) - predicted_response
        };
        if (!std::isfinite(residual)) return std::nullopt;
        residual_list.emplace_back(residual);
    }

    const auto row_count{ static_cast<Eigen::Index>(residual_list.size()) };
    algorithm::WeightedRidgeSystem system;
    system.design_matrix.resize(row_count, column_count);
    system.design_matrix.setFromTriplets(
        triplet_list.begin(),
        triplet_list.end());
    Eigen::VectorXd column_square_sum{
        Eigen::VectorXd::Zero(column_count)
    };
    for (Eigen::Index column_index = 0;
        column_index < system.design_matrix.outerSize();
        column_index++)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator iter(
                system.design_matrix,
                column_index);
            iter;
            ++iter)
        {
            column_square_sum(column_index) += iter.value() * iter.value();
        }
    }
    system.response = Eigen::VectorXd::Zero(row_count);
    for (Eigen::Index row_index = 0;
        row_index < row_count;
        row_index++)
    {
        system.response(row_index) = residual_list.at(
            static_cast<std::size_t>(row_index));
    }
    system.previous_parameter = Eigen::VectorXd::Zero(column_count);
    system.ridge_diagonal = Eigen::VectorXd::Zero(column_count);

    const auto conditioning{
        EvaluateJointOffsetConditioning(
            system.design_matrix,
            kJointOffsetConditioningPivotRatioThreshold)
    };
    const auto conditioning_multiplier{
        conditioning.guard_required ?
            kCollinearJointOffsetRidgeMultiplier : 1.0
    };
    for (Eigen::Index column_index = 0;
        column_index < column_count;
        column_index++)
    {
        double parameter_multiplier{ 1.0 };
        const auto offset_column_base{
            static_cast<Eigen::Index>(
                key.size() * kJointPolishShapeParameterSize)
        };
        if (column_index < offset_column_base)
        {
            const auto local_position{
                static_cast<std::size_t>(column_index) /
                    kJointPolishShapeParameterSize
            };
            parameter_multiplier = ridge_multiplier_list.at(
                key.at(local_position));
        }
        else
        {
            const auto group_position{
                static_cast<std::size_t>(column_index - offset_column_base)
            };
            for (const auto local_position :
                parameterization.atom_position_list_by_group.at(
                    group_position))
            {
                parameter_multiplier = std::max(
                    parameter_multiplier,
                    ridge_multiplier_list.at(key.at(local_position)));
            }
        }
        const auto square_sum{ column_square_sum(column_index) };
        system.ridge_diagonal(column_index) =
            CalculateJointOffsetRidgeDiagonal(
                square_sum,
                std::max(parameter_multiplier, conditioning_multiplier));
    }

    const auto residual_scale{
        std::max(
            CalculateLocalFittingMedianAbsoluteDeviationScale(residual_list),
            kLocalFittingRobustScaleMin)
    };
    if (!std::isfinite(residual_scale)) return std::nullopt;
    Eigen::VectorXd weight{ Eigen::VectorXd::Ones(row_count) };
    for (Eigen::Index row_index = 0;
        row_index < row_count;
        row_index++)
    {
        weight(row_index) = algorithm::CalculateCauchyWeight(
            system.response(row_index),
            residual_scale,
            kRobustLossCutoffMultiplier);
    }

    Eigen::VectorXd direction;
    if (!reusable_solver.Solve(system, weight, direction) ||
        !direction.allFinite())
    {
        return std::nullopt;
    }

    return direction;
}

inline bool HasMaterialJointPolishChange(
    const std::vector<GaussianModel3D> & candidate_model_list,
    const std::vector<GaussianModel3D> & seed_model_list)
{
    if (candidate_model_list.size() != seed_model_list.size()) return false;
    for (std::size_t atom_position = 0;
        atom_position < candidate_model_list.size();
        atom_position++)
    {
        const auto change{
            CalculateLocalFittingTransformedChange(
                candidate_model_list.at(atom_position),
                seed_model_list.at(atom_position))
        };
        if (std::any_of(
                change.value_list.begin(),
                change.value_list.end(),
                [](double value)
                {
                    return std::isfinite(value) &&
                        value >=
                            kJointPolishTransformedChangeTolerance;
                }))
        {
            return true;
        }
    }
    return false;
}

struct JointPolishProposal
{
    LocalFittingStatePatch patch{};
    double effective_damping{ 0.0 };
    double step_norm{ 0.0 };
    std::vector<std::size_t> changed_atom_index_list{};
};

inline std::optional<JointPolishProposal>
BuildJointPolishProposal(
    const SecondStageLocalFittingContext & context,
    const LocalFittingState & outer_previous_state,
    const LocalFittingStateView & base_state,
    const LocalFittingClusterKey & key,
    const std::vector<LocalFittingObjectiveSampleRef> & sample_ref_list,
    const std::vector<double> & ridge_multiplier_list,
    ReusableWeightedRidgeSolver & reusable_solver,
    double trust_region_radius)
{
    constexpr double trust_region_tolerance{ 1.0e-12 };
    std::vector<GroupKey> group_key_by_atom_position;
    std::vector<GaussianModel3D> base_model_list;
    group_key_by_atom_position.reserve(key.size());
    base_model_list.reserve(key.size());
    for (const auto atom_index : key)
    {
        group_key_by_atom_position.emplace_back(
            context.at(atom_index).group_key);
        base_model_list.emplace_back(base_state.GetModel(atom_index));
    }
    const auto parameterization{
        BuildJointPolishParameterization(
            group_key_by_atom_position,
            base_model_list)
    };
    if (!parameterization.has_value()) return std::nullopt;

    const Eigen::VectorXd zero_direction{
        Eigen::VectorXd::Zero(parameterization->ParameterCount())
    };
    const auto seed_model_list{
        parameterization->DecodeModels(zero_direction, 0.0)
    };
    if (!seed_model_list.has_value() ||
        std::any_of(
            seed_model_list->begin(),
            seed_model_list->end(),
            [](const GaussianModel3D & model)
            {
                return !IsValidSecondStageGaussianModel(model);
            }))
    {
        return std::nullopt;
    }
    const auto seed_step_norm{
        CalculateLocalFittingClusterModelTrustRegionStepNorm(
            outer_previous_state,
            key,
            *seed_model_list)
    };
    if (!seed_step_norm.has_value() ||
        *seed_step_norm > trust_region_radius + trust_region_tolerance)
    {
        return std::nullopt;
    }
    const auto direction{
        BuildJointPolishDirection(
            context,
            base_state,
            *seed_model_list,
            key,
            sample_ref_list,
            ridge_multiplier_list,
            *parameterization,
            reusable_solver)
    };
    if (!direction.has_value()) return std::nullopt;

    double damping{ 1.0 };
    while (damping >= std::numeric_limits<double>::epsilon())
    {
        auto candidate_model_list{
            parameterization->DecodeModels(*direction, damping)
        };
        if (candidate_model_list.has_value() &&
            std::none_of(
                candidate_model_list->begin(),
                candidate_model_list->end(),
                [](const GaussianModel3D & model)
                {
                    return !IsValidSecondStageGaussianModel(model);
                }))
        {
            const auto step_norm{
                CalculateLocalFittingClusterModelTrustRegionStepNorm(
                    outer_previous_state,
                    key,
                    *candidate_model_list)
            };
            if (step_norm.has_value() &&
                *step_norm <= trust_region_radius + trust_region_tolerance)
            {
                if (!HasMaterialJointPolishChange(
                        *candidate_model_list,
                        *seed_model_list))
                {
                    return std::nullopt;
                }

                JointPolishProposal proposal;
                proposal.patch.atom_index_list = key;
                proposal.patch.mdpde_list.reserve(key.size());
                proposal.effective_damping = damping;
                proposal.step_norm = *step_norm;
                for (std::size_t atom_position = 0;
                    atom_position < key.size();
                    atom_position++)
                {
                    const auto atom_index{ key.at(atom_position) };
                    const auto & base_mdpde{ base_state.GetMdpde(atom_index) };
                    const auto & candidate_model{
                        candidate_model_list->at(atom_position)
                    };
                    proposal.patch.mdpde_list.emplace_back(
                        GaussianModel3DWithUncertainty{
                            candidate_model,
                            base_mdpde.GetStandardDeviationModel()
                        });
                    const auto base_coordinates{
                        EncodeLocalFittingTransformedCoordinates(
                            base_mdpde.GetModel())
                    };
                    const auto candidate_coordinates{
                        EncodeLocalFittingTransformedCoordinates(
                            candidate_model)
                    };
                    if (!base_coordinates.has_value() ||
                        !candidate_coordinates.has_value())
                    {
                        return std::nullopt;
                    }
                    if ((base_coordinates->array() !=
                        candidate_coordinates->array()).any())
                    {
                        proposal.changed_atom_index_list.emplace_back(atom_index);
                    }
                }
                return proposal;
            }
        }
        damping *= 0.5;
    }
    return std::nullopt;
}

} // namespace rhbm_gem::core::detail
