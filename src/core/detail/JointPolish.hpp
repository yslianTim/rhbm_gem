#pragma once

#include "core/detail/FitStateView.hpp"
#include "core/detail/JointFittingSolverPolicy.hpp"
#include "core/detail/ResidualEvaluation.hpp"
#include "core/detail/ClusterSolverWorkspace.hpp"
#include "core/detail/SecondStageContext.hpp"
#include "core/detail/TransformedChange.hpp"
#include "core/detail/TransformedGaussianModel.hpp"
#include "core/detail/TrustRegion.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <optional>
#include <ranges>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>
#include <rhbm_gem/utils/algorithm/WeightedRidgeSolver.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem::core::detail {

constexpr std::size_t kJointPolishShapeParameterSize{ 2 };
constexpr double kJointPolishResidualScaleMin{ 1.0e-12 };
constexpr double kJointPolishRobustLossCutoffMultiplier{ 1.345 };
constexpr double kJointPolishRidgeRatio{ 1.0e-3 };
constexpr double kJointPolishConditioningRidgeMultiplier{ 10.0 };
constexpr double kJointPolishConditioningPivotRatioThreshold{ 1.0e-8 };
constexpr double kJointPolishTransformedChangeTolerance{ 1.0e-4 };

struct JointPolishParameterization
{
    std::vector<std::size_t> group_position_by_atom{};
    Eigen::VectorXd seed_parameter{};

private:
    std::optional<std::vector<GaussianModel3D>> DecodeParameter(
        const Eigen::VectorXd & parameter) const
    {
        if (parameter.size() != seed_parameter.size() || !parameter.allFinite())
        {
            return std::nullopt;
        }

        std::vector<GaussianModel3D> model_list;
        model_list.reserve(group_position_by_atom.size());
        for (std::size_t atom_position = 0;
            atom_position < group_position_by_atom.size();
            atom_position++)
        {
            const Eigen::Vector3d shape_coordinates{
                parameter(ShapeColumn(atom_position, 0)),
                parameter(ShapeColumn(atom_position, 1)),
                0.0
            };
            const auto shape_model{
                DecodeTransformedCoordinates(shape_coordinates)
            };
            if (!shape_model.has_value()) return std::nullopt;
            const auto model{
                shape_model->WithOffset(parameter(OffsetColumn(atom_position)))
            };
            if (!EncodeTransformedCoordinates(model).has_value())
            {
                return std::nullopt;
            }
            model_list.emplace_back(model);
        }
        return model_list;
    }

public:
    Eigen::Index ShapeColumn(std::size_t atom_position, std::size_t shape_parameter_index) const
    {
        return static_cast<Eigen::Index>(
            atom_position * kJointPolishShapeParameterSize + shape_parameter_index);
    }

    Eigen::Index OffsetColumn(std::size_t atom_position) const
    {
        return static_cast<Eigen::Index>(
            group_position_by_atom.size() * kJointPolishShapeParameterSize +
            group_position_by_atom.at(atom_position));
    }

    std::optional<std::vector<GaussianModel3D>> DecodeModels(
        const Eigen::VectorXd & direction,
        double damping) const
    {
        if (direction.size() != seed_parameter.size() ||
            !std::isfinite(damping) || damping < 0.0 || damping > 1.0)
        {
            return std::nullopt;
        }
        const Eigen::VectorXd parameter{ seed_parameter + damping * direction };
        return DecodeParameter(parameter);
    }

    std::optional<std::vector<GaussianModel3D>> DecodeSeedModels() const
    {
        return DecodeParameter(seed_parameter);
    }
};

inline std::optional<JointPolishParameterization> BuildJointPolishParameterization(
    const std::vector<std::size_t> & group_id_by_atom_position,
    const std::vector<GaussianModel3D> & base_model_list)
{
    if (group_id_by_atom_position.empty() || group_id_by_atom_position.size() != base_model_list.size())
    {
        return std::nullopt;
    }

    std::map<std::size_t, std::size_t> group_position_by_id;
    for (const auto group_id : group_id_by_atom_position)
    {
        group_position_by_id.emplace(group_id, 0);
    }
    std::size_t group_position{ 0 };
    for (auto & group_entry : group_position_by_id)
    {
        group_entry.second = group_position++;
    }

    JointPolishParameterization parameterization;
    parameterization.group_position_by_atom.resize(base_model_list.size());
    parameterization.seed_parameter = Eigen::VectorXd::Zero(
        static_cast<Eigen::Index>(
            base_model_list.size() * kJointPolishShapeParameterSize + group_position_by_id.size()));
    std::vector<std::vector<double>> offset_list_by_group(group_position_by_id.size());

    for (std::size_t atom_position = 0; atom_position < base_model_list.size(); atom_position++)
    {
        const auto transformed{
            EncodeTransformedCoordinates(base_model_list.at(atom_position))
        };
        if (!transformed.has_value()) return std::nullopt;
        const auto atom_group_position{
            group_position_by_id.at(group_id_by_atom_position.at(atom_position))
        };
        parameterization.group_position_by_atom.at(atom_position) = atom_group_position;
        parameterization.seed_parameter(
            parameterization.ShapeColumn(atom_position, 0)) =
            (*transformed)(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex));
        parameterization.seed_parameter(
            parameterization.ShapeColumn(atom_position, 1)) =
            (*transformed)(static_cast<Eigen::Index>(kLogWidthChangeIndex));
        offset_list_by_group.at(atom_group_position).emplace_back(
            base_model_list.at(atom_position).GetOffset());
    }

    for (std::size_t current_group_position = 0;
        current_group_position < offset_list_by_group.size();
        current_group_position++)
    {
        auto & offset_list{ offset_list_by_group.at(current_group_position) };
        std::ranges::sort(offset_list);
        const auto middle{ offset_list.size() / 2 };
        const auto median{
            offset_list.size() % 2 == 0 ?
                0.5 * offset_list.at(middle - 1) + 0.5 * offset_list.at(middle) :
                offset_list.at(middle)
        };
        if (!std::isfinite(median)) return std::nullopt;
        parameterization.seed_parameter(
            static_cast<Eigen::Index>(
                parameterization.group_position_by_atom.size() * kJointPolishShapeParameterSize +
                current_group_position)) = median;
    }
    return parameterization;
}

inline std::optional<Eigen::VectorXd> BuildJointPolishDirection(
    const SecondStageContext & context,
    const FitStateView & base_state,
    const ClusterKey & key,
    const std::vector<SampleRef> & sample_ref_list,
    const std::vector<double> & ridge_multiplier_list,
    const JointPolishParameterization & parameterization,
    ReusableWeightedRidgeSolver & reusable_solver)
{
    if (key.empty() || sample_ref_list.empty() || parameterization.group_position_by_atom.size() != key.size())
    {
        return std::nullopt;
    }
    const auto seed_model_list{ parameterization.DecodeSeedModels() };
    if (!seed_model_list.has_value()) return std::nullopt;

    const auto column_count{ parameterization.seed_parameter.size() };
    std::unordered_map<std::size_t, std::size_t> local_position_by_atom_index;
    local_position_by_atom_index.reserve(key.size());
    std::unordered_map<std::size_t, Eigen::Index> offset_column_by_group_id;
    Eigen::VectorXd ridge_multiplier_by_column{ Eigen::VectorXd::Ones(column_count) };
    for (std::size_t local_position = 0; local_position < key.size(); local_position++)
    {
        const auto atom_index{ key.at(local_position) };
        local_position_by_atom_index.emplace(atom_index, local_position);
        const auto ridge_multiplier{ ridge_multiplier_list.at(atom_index) };
        for (std::size_t parameter_index = 0;
            parameter_index < kJointPolishShapeParameterSize;
            parameter_index++)
        {
            ridge_multiplier_by_column(
                parameterization.ShapeColumn(local_position, parameter_index)) =
                ridge_multiplier;
        }
        const auto offset_column{ parameterization.OffsetColumn(local_position) };
        offset_column_by_group_id.emplace(
            context.at(atom_index).group_id,
            offset_column);
        ridge_multiplier_by_column(offset_column) = std::max(
            ridge_multiplier_by_column(offset_column),
            ridge_multiplier);
    }
    auto selected_snapshot{ BuildFittedGaussianSnapshot(base_state) };
    for (std::size_t local_position = 0; local_position < key.size(); local_position++)
    {
        selected_snapshot.at(key.at(local_position)) = seed_model_list->at(local_position);
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
            const auto & model{
                local_position_iter != local_position_by_atom_index.end() ?
                    seed_model_list->at(local_position_iter->second) :
                    base_state.GetModel(atom_index)
            };
            const auto evaluation{ EvaluateSharedOffsetResponse(model, distance) };
            if (!evaluation.has_value()) return false;
            predicted_response += evaluation->response;

            if (local_position_iter == local_position_by_atom_index.end()) return true;
            const auto local_position{ local_position_iter->second };
            for (std::size_t parameter_index = 0;
                parameter_index < kJointPolishShapeParameterSize;
                parameter_index++)
            {
                const auto column_index{
                    parameterization.ShapeColumn(local_position, parameter_index)
                };
                const auto derivative{
                    evaluation->shape_jacobian(static_cast<Eigen::Index>(parameter_index))
                };
                if (std::abs(derivative) <= std::numeric_limits<double>::epsilon()) continue;
                triplet_list.emplace_back(row_index, column_index, derivative);
            }
            if (std::abs(evaluation->offset_jacobian) > std::numeric_limits<double>::epsilon())
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
            const auto & unselected_atom_contributor{
                context.unselected_atom_list.at(contributor_index)
            };
            const auto & model{
                GetFitModel(model_snapshot.unselected, contributor_index)
            };
            const auto evaluation{ EvaluateSharedOffsetResponse(model, distance) };
            if (!evaluation.has_value()) return false;
            predicted_response += evaluation->response;

            const auto selected_group_id{ unselected_atom_contributor.selected_group_id };
            const auto local_position_iter{
                selected_group_id.has_value() ?
                    offset_column_by_group_id.find(*selected_group_id) :
                    offset_column_by_group_id.end()
            };
            if (local_position_iter == offset_column_by_group_id.end() ||
                std::abs(evaluation->offset_jacobian) <= std::numeric_limits<double>::epsilon())
            {
                return true;
            }
            triplet_list.emplace_back(row_index, local_position_iter->second, evaluation->offset_jacobian);
            return true;
        };

        if (!append_model(sample_ref.atom_index, static_cast<double>(sample.point.distance)))
        {
            return std::nullopt;
        }
        for (const auto & neighbor_atom_sample : atom_context.Neighbors(sample_ref.sample_index))
        {
            const auto appended{
                neighbor_atom_sample.is_selected ?
                    append_model(
                        neighbor_atom_sample.atom_index,
                        neighbor_atom_sample.distance) :
                    append_unselected_model(
                        neighbor_atom_sample.atom_index,
                        neighbor_atom_sample.distance)
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
    system.design_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
    Eigen::VectorXd column_square_sum{ Eigen::VectorXd::Zero(column_count) };
    for (Eigen::Index column_index = 0; column_index < system.design_matrix.outerSize(); column_index++)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator iter(system.design_matrix, column_index);
            iter;
            ++iter)
        {
            column_square_sum(column_index) += iter.value() * iter.value();
        }
    }
    system.response = Eigen::VectorXd::Zero(row_count);
    for (Eigen::Index row_index = 0; row_index < row_count; row_index++)
    {
        system.response(row_index) = residual_list.at(static_cast<std::size_t>(row_index));
    }
    system.previous_parameter = Eigen::VectorXd::Zero(column_count);
    system.ridge_diagonal = Eigen::VectorXd::Zero(column_count);

    const auto conditioning{
        EvaluateJointFittingConditioning(
            system.design_matrix,
            kJointPolishConditioningPivotRatioThreshold)
    };
    if (conditioning.guard_required)
    {
        ridge_multiplier_by_column.array() = ridge_multiplier_by_column.array().max(kJointPolishConditioningRidgeMultiplier);
    }
    for (Eigen::Index column_index = 0; column_index < column_count; column_index++)
    {
        system.ridge_diagonal(column_index) =
            CalculateJointFittingRidgeDiagonal(
                column_square_sum(column_index),
                kJointPolishRidgeRatio,
                ridge_multiplier_by_column(column_index));
    }

    const auto residual_scale{
        std::max(
            array_helper::ComputeMedianAbsoluteDeviationScale(residual_list),
            kJointPolishResidualScaleMin)
    };
    if (!std::isfinite(residual_scale)) return std::nullopt;
    Eigen::VectorXd weight{ Eigen::VectorXd::Ones(row_count) };
    for (Eigen::Index row_index = 0; row_index < row_count; row_index++)
    {
        weight(row_index) = algorithm::CalculateCauchyWeight(
            system.response(row_index),
            residual_scale,
            kJointPolishRobustLossCutoffMultiplier);
    }

    Eigen::VectorXd direction;
    if (!reusable_solver.Solve(system, weight, direction))
    {
        return std::nullopt;
    }

    return direction;
}

struct JointPolishProposal
{
    FitStatePatch patch{};
    double effective_damping{ 0.0 };
    double step_norm{ 0.0 };
};

inline std::optional<JointPolishProposal> BuildJointPolishProposal(
    const SecondStageContext & context,
    const FitStateView & base_state,
    const ClusterKey & key,
    const std::vector<SampleRef> & sample_ref_list,
    const std::vector<double> & ridge_multiplier_list,
    ReusableWeightedRidgeSolver & reusable_solver,
    double trust_region_radius)
{
    std::vector<std::size_t> group_id_by_atom_position;
    std::vector<GaussianModel3D> outer_previous_model_list;
    std::vector<GaussianModel3D> base_model_list;
    group_id_by_atom_position.reserve(key.size());
    outer_previous_model_list.reserve(key.size());
    base_model_list.reserve(key.size());
    for (const auto atom_index : key)
    {
        group_id_by_atom_position.emplace_back(context.at(atom_index).group_id);
        outer_previous_model_list.emplace_back(base_state.GetBaseModel(atom_index));
        base_model_list.emplace_back(base_state.GetModel(atom_index));
    }
    const auto parameterization{
        BuildJointPolishParameterization(group_id_by_atom_position, base_model_list)
    };
    if (!parameterization.has_value()) return std::nullopt;

    const auto seed_model_list{ parameterization->DecodeSeedModels() };
    if (!seed_model_list.has_value()) return std::nullopt;
    const auto seed_step_norm{
        CalculateModelTrustRegionStepNorm(outer_previous_model_list, *seed_model_list)
    };
    if (!seed_step_norm.has_value() ||
        !IsTrustRegionStepWithinRadius(*seed_step_norm, trust_region_radius))
    {
        return std::nullopt;
    }
    const auto direction{
        BuildJointPolishDirection(
            context,
            base_state,
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
        if (candidate_model_list.has_value())
        {
            const auto step_norm{
                CalculateModelTrustRegionStepNorm(
                    outer_previous_model_list,
                    *candidate_model_list)
            };
            if (step_norm.has_value() &&
                IsTrustRegionStepWithinRadius(*step_norm, trust_region_radius))
            {
                bool has_material_change{ false };
                for (std::size_t atom_position = 0; atom_position < candidate_model_list->size(); atom_position++)
                {
                    const auto change{
                        CalculateTransformedChange(
                            candidate_model_list->at(atom_position),
                            seed_model_list->at(atom_position))
                    };
                    if (IsTransformedChangeMaterial(change, kJointPolishTransformedChangeTolerance))
                    {
                        has_material_change = true;
                        break;
                    }
                }
                if (!has_material_change)
                {
                    return std::nullopt;
                }

                JointPolishProposal proposal{
                    .patch{ .atom_index_list = key },
                    .effective_damping = damping,
                    .step_norm = *step_norm
                };
                proposal.patch.mdpde_list.reserve(key.size());
                for (std::size_t atom_position = 0; atom_position < key.size(); atom_position++)
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
                }
                return proposal;
            }
        }
        damping *= 0.5;
    }
    return std::nullopt;
}

} // namespace rhbm_gem::core::detail
