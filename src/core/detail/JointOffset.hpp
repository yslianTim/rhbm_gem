#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>
#include <rhbm_gem/utils/algorithm/Convergence.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include "core/detail/ResidualEvaluation.hpp"
#include "core/detail/JointFittingSolverPolicy.hpp"
#include "core/detail/JointOffsetSolveStatus.hpp"
#include "core/detail/ClusterSolverWorkspace.hpp"
#include "core/detail/SecondStageContext.hpp"

namespace rhbm_gem::core::detail {

constexpr int kRobustLossMaximumIterations{ 50 };
constexpr double kJointOffsetRobustLossCutoffMultiplier{ 1.345 };
constexpr double kJointOffsetResidualScaleMin{ 1.0e-12 };
constexpr double kJointOffsetRidgeRatio{ 1.0e-3 };
constexpr double kJointOffsetCollinearityOverlapThreshold{ 0.98 };
constexpr double kJointOffsetConditioningRidgeMultiplier{ 10.0 };
constexpr double kJointOffsetConditioningPivotRatioThreshold{ 1.0e-8 };
constexpr double kJointOffsetIrlsScaleFloor{ 1.0e-2 };
constexpr double kJointOffsetIrlsNormalizedChangeTolerance{ 1.0e-6 };
constexpr double kJointOffsetIrlsObjectiveRelativeTolerance{ 1.0e-10 };

struct JointOffsetParameterization
{
    std::vector<std::size_t> group_position_by_atom{};
    Eigen::VectorXd seed_offset{};

    std::size_t AtomCount() const
    {
        return group_position_by_atom.size();
    }

    Eigen::Index ParameterCount() const
    {
        return seed_offset.size();
    }

    Eigen::Index OffsetColumn(std::size_t atom_position) const
    {
        return static_cast<Eigen::Index>(group_position_by_atom.at(atom_position));
    }

    Eigen::VectorXd ExpandOffsets(const Eigen::VectorXd & group_offset) const
    {
        Eigen::VectorXd atom_offset{
            Eigen::VectorXd::Zero(static_cast<Eigen::Index>(AtomCount()))
        };
        for (std::size_t atom_position = 0; atom_position < AtomCount(); atom_position++)
        {
            atom_offset(static_cast<Eigen::Index>(atom_position)) =
                group_offset(OffsetColumn(atom_position));
        }
        return atom_offset;
    }
};

inline std::optional<JointOffsetParameterization> BuildJointOffsetParameterization(
    const std::vector<GroupKey> & group_key_by_atom_position,
    const Eigen::VectorXd & atom_offset)
{
    const auto atom_count{ static_cast<std::size_t>(atom_offset.size()) };
    if (group_key_by_atom_position.empty() || group_key_by_atom_position.size() != atom_count)
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

    JointOffsetParameterization parameterization;
    parameterization.group_position_by_atom.resize(atom_count);
    parameterization.seed_offset = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(group_position_by_key.size()));
    std::vector<std::vector<double>> offset_list_by_group(group_position_by_key.size());

    for (std::size_t atom_position = 0; atom_position < atom_count; atom_position++)
    {
        const auto offset{ atom_offset(static_cast<Eigen::Index>(atom_position)) };
        if (!std::isfinite(offset)) return std::nullopt;

        const auto atom_group_position{
            group_position_by_key.at(group_key_by_atom_position.at(atom_position))
        };
        parameterization.group_position_by_atom.at(atom_position) = atom_group_position;
        offset_list_by_group.at(atom_group_position).emplace_back(offset);
    }

    for (std::size_t current_group_position = 0;
        current_group_position < offset_list_by_group.size();
        current_group_position++)
    {
        auto & offset_list{ offset_list_by_group.at(current_group_position) };
        if (offset_list.empty()) return std::nullopt;
        std::sort(offset_list.begin(), offset_list.end());
        const auto middle{ offset_list.size() / 2 };
        const auto median{ offset_list.size() % 2 == 0 ?
            0.5 * offset_list.at(middle - 1) + 0.5 * offset_list.at(middle) :
            offset_list.at(middle)
        };
        if (!std::isfinite(median)) return std::nullopt;
        parameterization.seed_offset(static_cast<Eigen::Index>(current_group_position)) = median;
    }
    return parameterization;
}

struct JointOffsetSolveResult
{
    JointOffsetSolveStatus status{ JointOffsetSolveStatus::SystemBuildFailed };
    Eigen::VectorXd offset{};
};

inline algorithm::WeightedRidgeSystem BuildJointOffsetSystem(
    const SecondStageContext & context,
    const std::vector<std::size_t> & active_index_list,
    const SecondStageModelSnapshot & model_snapshot,
    const std::vector<double> & ridge_multiplier_list,
    const JointOffsetParameterization & parameterization,
    bool log_debug_diagnostics)
{
    std::unordered_map<std::size_t, Eigen::Index> active_offset_column_by_atom_index;
    active_offset_column_by_atom_index.reserve(active_index_list.size());
    std::unordered_map<GroupKey, Eigen::Index> active_offset_column_by_group_key;
    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto atom_index{ active_index_list.at(i) };
        const auto offset_column{ parameterization.OffsetColumn(i) };
        active_offset_column_by_atom_index.emplace(atom_index, offset_column);
        active_offset_column_by_group_key.emplace(
            context.at(atom_index).group_key,
            offset_column);
    }

    const auto column_count{ parameterization.ParameterCount() };
    std::vector<Eigen::Triplet<double>> triplet_list;
    std::vector<double> response_list;
    Eigen::VectorXd group_column_square_sum{ Eigen::VectorXd::Zero(column_count) };
    std::map<std::pair<Eigen::Index, Eigen::Index>, double> group_column_cross_sum_map;
    std::vector<std::pair<Eigen::Index, double>> group_row_basis_entries;
    group_row_basis_entries.reserve(static_cast<std::size_t>(column_count));
    for (const auto active_index : active_index_list)
    {
        const auto & atom_context{ context.at(active_index) };
        const auto & target_model{
            GetFitModel(model_snapshot.selected, active_index)
        };
        for (std::size_t sample_index = 0; sample_index < atom_context.raw_sampling_entries.size(); sample_index++)
        {
            const auto & sample{ atom_context.raw_sampling_entries.at(sample_index) };
            if (!std::isfinite(static_cast<double>(sample.response)))
            {
                throw std::runtime_error("Joint offset sample response is not finite.");
            }
            const auto target_distance{ static_cast<double>(sample.point.distance) };
            const auto target_signal{ target_model.SignalAtDistance(target_distance) };
            const auto target_basis{ target_model.OffsetBasisAtDistance(target_distance) };
            if (!std::isfinite(target_signal) || !std::isfinite(target_basis))
            {
                throw std::runtime_error("Joint offset target model evaluation is not finite.");
            }
            auto residual{ static_cast<double>(sample.response) - target_signal };
            Eigen::VectorXd group_basis{ Eigen::VectorXd::Zero(column_count) };
            if (std::abs(target_basis) > std::numeric_limits<double>::epsilon())
            {
                group_basis(active_offset_column_by_atom_index.at(active_index)) += target_basis;
            }

            for (auto neighbor_iter = atom_context.NeighborBegin(sample_index);
                neighbor_iter != atom_context.NeighborEnd(sample_index);
                ++neighbor_iter)
            {
                const auto & neighbor_atom_sample{ *neighbor_iter };
                const auto & neighbor_model{
                    ResolveNeighborAtomModel(neighbor_atom_sample, model_snapshot)
                };
                Eigen::Index neighbor_offset_column{ -1 };
                if (neighbor_atom_sample.is_selected)
                {
                    const auto neighbor_offset_column_iter{
                        active_offset_column_by_atom_index.find(neighbor_atom_sample.atom_index)
                    };
                    if (neighbor_offset_column_iter != active_offset_column_by_atom_index.end())
                    {
                        neighbor_offset_column = neighbor_offset_column_iter->second;
                    }
                }
                else
                {
                    const auto group_key{
                        context.unselected_atom_list.at(neighbor_atom_sample.atom_index).group_key
                    };
                    const auto offset_column_iter{
                        active_offset_column_by_group_key.find(group_key)
                    };
                    if (offset_column_iter != active_offset_column_by_group_key.end())
                    {
                        neighbor_offset_column = offset_column_iter->second;
                    }
                }
                if (neighbor_offset_column < 0)
                {
                    const auto response{
                        neighbor_model.ResponseAtDistance(neighbor_atom_sample.distance)
                    };
                    if (!std::isfinite(response))
                    {
                        throw std::runtime_error(
                            "Joint offset fixed neighbor model evaluation is not finite.");
                    }
                    residual -= response;
                    continue;
                }

                const auto signal{
                    neighbor_model.SignalAtDistance(neighbor_atom_sample.distance)
                };
                const auto basis{
                    neighbor_model.OffsetBasisAtDistance(neighbor_atom_sample.distance)
                };
                if (!std::isfinite(signal) || !std::isfinite(basis))
                {
                    throw std::runtime_error(
                        "Joint offset active neighbor model evaluation is not finite.");
                }
                residual -= signal;
                if (std::abs(basis) > std::numeric_limits<double>::epsilon())
                {
                    group_basis(neighbor_offset_column) += basis;
                }
            }
            if (!std::isfinite(residual))
            {
                throw std::runtime_error("Joint offset residual is not finite.");
            }
            if (!group_basis.allFinite())
            {
                throw std::runtime_error("Joint offset group basis is invalid.");
            }
            group_row_basis_entries.clear();
            for (Eigen::Index column_index = 0; column_index < group_basis.size(); column_index++)
            {
                const auto basis{ group_basis(column_index) };
                if (std::abs(basis) <= std::numeric_limits<double>::epsilon())
                {
                    continue;
                }
                group_row_basis_entries.emplace_back(column_index, basis);
            }
            if (group_row_basis_entries.empty()) continue;

            const auto row_index{
                static_cast<Eigen::Index>(response_list.size())
            };
            response_list.emplace_back(residual);
            for (const auto & [column_index, basis] : group_row_basis_entries)
            {
                triplet_list.emplace_back(row_index, column_index, basis);
                group_column_square_sum(column_index) += basis * basis;
            }
            for (std::size_t i = 0; i < group_row_basis_entries.size(); i++)
            {
                const auto [left_column, left_basis]{
                    group_row_basis_entries.at(i)
                };
                for (std::size_t j = i + 1; j < group_row_basis_entries.size(); j++)
                {
                    const auto [right_column, right_basis]{
                        group_row_basis_entries.at(j)
                    };
                    const auto column_pair{ std::minmax(left_column, right_column) };
                    group_column_cross_sum_map[column_pair] += left_basis * right_basis;
                }
            }
        }
    }

    Eigen::VectorXd proactive_ridge_multiplier{ Eigen::VectorXd::Ones(column_count) };
    for (const auto & [column_pair, cross_sum] : group_column_cross_sum_map)
    {
        const auto left_column{ column_pair.first };
        const auto right_column{ column_pair.second };
        const auto left_square_sum{ group_column_square_sum(left_column) };
        const auto right_square_sum{ group_column_square_sum(right_column) };
        if (left_square_sum <= std::numeric_limits<double>::epsilon() ||
            right_square_sum <= std::numeric_limits<double>::epsilon())
        {
            continue;
        }
        const auto overlap{
            std::abs(cross_sum) / std::sqrt(left_square_sum * right_square_sum)
        };
        if (!std::isfinite(overlap) || overlap < kJointOffsetCollinearityOverlapThreshold)
        {
            continue;
        }

        proactive_ridge_multiplier(left_column) = std::max(
            proactive_ridge_multiplier(left_column),
            kJointOffsetConditioningRidgeMultiplier);
        proactive_ridge_multiplier(right_column) = std::max(
            proactive_ridge_multiplier(right_column),
            kJointOffsetConditioningRidgeMultiplier);
    }

    const auto row_count{ static_cast<Eigen::Index>(response_list.size()) };
    algorithm::WeightedRidgeSystem system;
    system.design_matrix.resize(row_count, column_count);
    system.design_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
    system.response = Eigen::VectorXd::Zero(row_count);
    for (Eigen::Index row_index = 0; row_index < row_count; row_index++)
    {
        system.response(row_index) = response_list.at(static_cast<std::size_t>(row_index));
    }
    const auto conditioning{
        EvaluateJointFittingConditioning(
            system.design_matrix,
            kJointOffsetConditioningPivotRatioThreshold)
    };
    if (conditioning.guard_required)
    {
        proactive_ridge_multiplier.array() = proactive_ridge_multiplier.array().max(
            kJointOffsetConditioningRidgeMultiplier);
        if (log_debug_diagnostics)
        {
            std::ostringstream message;
            message
                << std::scientific << std::setprecision(2)
                << "Joint offset conditioning guard: columns = "
                << column_count
                << ", normalized LDLT pivot ratio = "
                << conditioning.pivot_ratio
                << ", proactive ridge multiplier = "
                << kJointOffsetConditioningRidgeMultiplier << ".";
            Logger::Log(LogLevel::Debug, message.str());
        }
    }
    system.previous_parameter = parameterization.seed_offset;
    system.ridge_diagonal = Eigen::VectorXd::Zero(column_count);
    Eigen::VectorXd group_ridge_multiplier{ Eigen::VectorXd::Ones(column_count) };
    for (std::size_t atom_position = 0; atom_position < active_index_list.size(); atom_position++)
    {
        const auto atom_index{ active_index_list.at(atom_position) };
        const auto offset_column{ parameterization.OffsetColumn(atom_position) };
        group_ridge_multiplier(offset_column) = std::max(
            group_ridge_multiplier(offset_column),
            ridge_multiplier_list.at(atom_index));
    }
    for (Eigen::Index column_index = 0; column_index < column_count; column_index++)
    {
        const auto square_sum{ group_column_square_sum(column_index) };
        const auto combined_multiplier{
            std::max(
                group_ridge_multiplier(column_index),
                proactive_ridge_multiplier(column_index))
        };
        system.ridge_diagonal(column_index) =
            CalculateJointFittingRidgeDiagonal(
                square_sum,
                kJointOffsetRidgeRatio,
                combined_multiplier);
    }
    return system;
}

inline double CalculateWeightedRidgeSurrogateObjective(
    const algorithm::WeightedRidgeSystem & system,
    const Eigen::VectorXd & weight,
    const Eigen::VectorXd & offset)
{
    if (system.response.size() != weight.size())
    {
        throw std::invalid_argument(
            "Weighted ridge objective input sizes are inconsistent.");
    }
    if (system.previous_parameter.size() != offset.size() ||
        system.ridge_diagonal.size() != offset.size())
    {
        throw std::invalid_argument(
            "Weighted ridge objective parameter sizes are inconsistent.");
    }
    if (system.response.size() == 0)
    {
        return std::numeric_limits<double>::infinity();
    }

    const Eigen::VectorXd residual{
        system.response - system.design_matrix * offset
    };
    const auto weighted_residual_loss{
        weight.cwiseProduct(residual.cwiseAbs2()).sum()
    };
    const Eigen::VectorXd offset_delta{
        offset - system.previous_parameter
    };
    const auto ridge_loss{
        system.ridge_diagonal.cwiseProduct(offset_delta.cwiseAbs2()).sum()
    };
    const auto objective{
        (weighted_residual_loss + ridge_loss) / static_cast<double>(system.response.size())
    };
    return std::isfinite(objective) ? objective : std::numeric_limits<double>::infinity();
}

inline bool IsJointOffsetObjectiveDeteriorated(double updated_objective, double current_objective)
{
    if (!std::isfinite(updated_objective)) return true;
    if (!std::isfinite(current_objective)) return false;
    const auto scale{
        std::max({
            std::abs(updated_objective),
            std::abs(current_objective),
            1.0
        })
    };
    return updated_objective > current_objective + kJointOffsetIrlsObjectiveRelativeTolerance * scale;
}

inline JointOffsetSolveResult EstimateJointOffsets(
    const SecondStageContext & context,
    const std::vector<std::size_t> & active_index_list,
    const SecondStageModelSnapshot & model_snapshot,
    const std::vector<double> & ridge_multiplier_list,
    ReusableWeightedRidgeSolver & reusable_solver,
    bool log_debug_diagnostics)
{
    Eigen::VectorXd previous_offset{
        Eigen::VectorXd::Zero(static_cast<Eigen::Index>(active_index_list.size()))
    };
    std::vector<GroupKey> group_key_by_atom_position;
    group_key_by_atom_position.reserve(active_index_list.size());
    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto atom_index{ active_index_list.at(i) };
        group_key_by_atom_position.emplace_back(context.at(atom_index).group_key);
        previous_offset(static_cast<Eigen::Index>(i)) =
            GetFitModel(model_snapshot.selected, atom_index).GetOffset();
    }
    auto parameterization{
        BuildJointOffsetParameterization(group_key_by_atom_position, previous_offset)
    };
    if (!parameterization.has_value())
    {
        return JointOffsetSolveResult{
            JointOffsetSolveStatus::SystemBuildFailed,
            previous_offset
        };
    }
    algorithm::WeightedRidgeSystem system;
    try
    {
        system = BuildJointOffsetSystem(
            context,
            active_index_list,
            model_snapshot,
            ridge_multiplier_list,
            *parameterization,
            log_debug_diagnostics);
    }
    catch (const std::runtime_error &)
    {
        return JointOffsetSolveResult{
            JointOffsetSolveStatus::SystemBuildFailed,
            previous_offset
        };
    }
    const auto make_progress_result = [&](
        JointOffsetSolveStatus status,
        const Eigen::VectorXd & group_offset)
    {
        auto atom_offset{ parameterization->ExpandOffsets(group_offset) };
        return JointOffsetSolveResult{
            status,
            std::move(atom_offset)
        };
    };
    if (system.response.size() == 0 || system.previous_parameter.size() == 0)
    {
        return JointOffsetSolveResult{
            JointOffsetSolveStatus::EmptySystem,
            previous_offset
        };
    }

    Eigen::VectorXd weight{ Eigen::VectorXd::Ones(system.response.size()) };
    Eigen::VectorXd offset;
    if (!reusable_solver.Solve(system, weight, offset))
    {
        return JointOffsetSolveResult{
            JointOffsetSolveStatus::InitialSolveFailed,
            previous_offset
        };
    }

    for (int iteration = 0; iteration < kRobustLossMaximumIterations; iteration++)
    {
        const Eigen::VectorXd residual{ system.response - system.design_matrix * offset };
        std::vector<double> residual_list(residual.data(), residual.data() + residual.size());
        const auto residual_scale{
            std::max(
                array_helper::ComputeMedianAbsoluteDeviationScale(residual_list),
                kJointOffsetResidualScaleMin)
        };
        for (Eigen::Index i = 0; i < residual.size(); i++)
        {
            weight(i) = algorithm::CalculateCauchyWeight(
                residual(i),
                residual_scale,
                kJointOffsetRobustLossCutoffMultiplier);
        }

        Eigen::VectorXd updated_offset;
        if (!reusable_solver.Solve(system, weight, updated_offset))
        {
            return JointOffsetSolveResult{
                JointOffsetSolveStatus::IrlsSolveFailed,
                previous_offset
            };
        }
        const auto current_objective{
            CalculateWeightedRidgeSurrogateObjective(system, weight, offset)
        };
        const auto updated_objective{
            CalculateWeightedRidgeSurrogateObjective(system, weight, updated_offset)
        };
        if (IsJointOffsetObjectiveDeteriorated(updated_objective, current_objective))
        {
            return make_progress_result(
                JointOffsetSolveStatus::IrlsObjectiveDeteriorated,
                offset);
        }
        const auto maximum_change{
            algorithm::CalculateMaximumNormalizedVectorChange(
                updated_offset,
                offset,
                kJointOffsetIrlsScaleFloor)
        };
        offset = std::move(updated_offset);
        if (maximum_change < kJointOffsetIrlsNormalizedChangeTolerance)
        {
            return make_progress_result(JointOffsetSolveStatus::Converged, offset);
        }
    }

    return make_progress_result(JointOffsetSolveStatus::IrlsMaximumIterationsReached, offset);
}

} // namespace rhbm_gem::core::detail
