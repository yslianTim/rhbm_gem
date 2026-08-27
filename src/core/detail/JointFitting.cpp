#include "core/detail/JointFitting.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ranges>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include <rhbm_gem/utils/algorithm/Convergence.hpp>
#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

namespace rhbm_gem::core::detail {

namespace {

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

constexpr double kJointPolishResidualScaleMin{ 1.0e-12 };
constexpr double kJointPolishRobustLossCutoffMultiplier{ 1.345 };
constexpr double kJointPolishRidgeRatio{ 1.0e-3 };
constexpr double kJointPolishConditioningRidgeMultiplier{ 10.0 };
constexpr double kJointPolishConditioningPivotRatioThreshold{ 1.0e-8 };
constexpr double kJointPolishTransformedChangeTolerance{ 1.0e-4 };

} // namespace

bool ReusableWeightedRidgeSolver::Solve(
    const algorithm::WeightedRidgeSystem & system,
    const Eigen::VectorXd & weight,
    Eigen::VectorXd & parameter)
{
    const auto pattern{ BuildPattern(system.design_matrix) };
    if (m_row_count != system.design_matrix.rows() ||
        m_column_count != system.design_matrix.cols() ||
        m_pattern != pattern)
    {
        if (!m_solver.AnalyzePattern(system)) return false;
        m_row_count = system.design_matrix.rows();
        m_column_count = system.design_matrix.cols();
        m_pattern = pattern;
    }
    return m_solver.SolveNumeric(system, weight, parameter);
}

std::vector<std::pair<Eigen::Index, Eigen::Index>>
ReusableWeightedRidgeSolver::BuildPattern(
    const Eigen::SparseMatrix<double> & matrix)
{
    std::vector<std::pair<Eigen::Index, Eigen::Index>> pattern;
    pattern.reserve(static_cast<std::size_t>(matrix.nonZeros()));
    for (Eigen::Index outer = 0; outer < matrix.outerSize(); outer++)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator iter(matrix, outer);
            iter;
            ++iter)
        {
            pattern.emplace_back(iter.row(), iter.col());
        }
    }
    return pattern;
}

Eigen::VectorXd JointOffsetParameterization::ExpandOffsets(
    const Eigen::VectorXd & group_offset) const
{
    Eigen::VectorXd atom_offset{
        Eigen::VectorXd::Zero(
            static_cast<Eigen::Index>(group_position_by_atom.size()))
    };
    for (std::size_t atom_position = 0;
        atom_position < group_position_by_atom.size();
        atom_position++)
    {
        atom_offset(static_cast<Eigen::Index>(atom_position)) =
            group_offset(OffsetColumn(atom_position));
    }
    return atom_offset;
}

std::optional<std::vector<GaussianModel3D>>
JointPolishParameterization::DecodeParameter(
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
        auto active_shape_coordinates{
            base_shape_coordinate_by_atom.at(atom_position)
        };
        if (HasShapeColumn(atom_position))
        {
            active_shape_coordinates(0) =
                parameter(ShapeColumn(atom_position, 0));
            active_shape_coordinates(1) =
                parameter(ShapeColumn(atom_position, 1));
        }
        const Eigen::Vector3d shape_coordinates{
            active_shape_coordinates(0),
            active_shape_coordinates(1),
            0.0
        };
        const auto shape_model{
            DecodeTransformedCoordinates(shape_coordinates)
        };
        if (!shape_model.has_value()) return std::nullopt;
        const auto group_position{ group_position_by_atom.at(atom_position) };
        const auto offset{
            HasOffsetColumn(atom_position) ?
                parameter(OffsetColumn(atom_position)) :
                base_offset_by_group.at(group_position)
        };
        const auto model{ shape_model->WithOffset(offset) };
        if (!EncodeTransformedCoordinates(model).has_value())
        {
            return std::nullopt;
        }
        model_list.emplace_back(model);
    }
    return model_list;
}

std::optional<std::vector<GaussianModel3D>>
JointPolishParameterization::DecodeModels(
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

std::optional<std::vector<GaussianModel3D>>
JointPolishParameterization::DecodeSeedModels() const
{
    return DecodeParameter(seed_parameter);
}

namespace {

struct JointFittingGroupLayout
{
    std::vector<std::size_t> group_position_by_atom{};
    std::size_t group_count{ 0 };
};

std::optional<JointFittingGroupLayout> BuildJointFittingGroupLayout(
    const std::vector<std::size_t> & group_id_by_atom_position,
    std::size_t atom_count)
{
    if (group_id_by_atom_position.empty() ||
        group_id_by_atom_position.size() != atom_count)
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

    JointFittingGroupLayout layout;
    layout.group_position_by_atom.reserve(atom_count);
    for (const auto group_id : group_id_by_atom_position)
    {
        layout.group_position_by_atom.emplace_back(
            group_position_by_id.at(group_id));
    }
    layout.group_count = group_position_by_id.size();
    return layout;
}

std::optional<double> CalculateJointFittingGroupMedian(
    std::vector<double> & value_list)
{
    if (value_list.empty()) return std::nullopt;

    std::ranges::sort(value_list);
    const auto middle{ value_list.size() / 2 };
    const auto median{
        value_list.size() % 2 == 0 ?
            0.5 * value_list.at(middle - 1) + 0.5 * value_list.at(middle) :
            value_list.at(middle)
    };
    return std::isfinite(median) ?
        std::optional<double>{ median } : std::nullopt;
}

} // namespace

JointFittingConditioning EvaluateJointFittingConditioning(
    const Eigen::SparseMatrix<double> & design_matrix,
    double pivot_ratio_threshold)
{
    if (design_matrix.cols() == 0)
    {
        return JointFittingConditioning{ true, 0.0 };
    }

    Eigen::SparseMatrix<double> normalized_design{ design_matrix };
    for (Eigen::Index column = 0; column < normalized_design.outerSize(); column++)
    {
        double square_sum{ 0.0 };
        for (Eigen::SparseMatrix<double>::InnerIterator entry(normalized_design, column);
            entry;
            ++entry)
        {
            square_sum += entry.value() * entry.value();
        }
        if (!std::isfinite(square_sum) ||
            square_sum <= std::numeric_limits<double>::epsilon())
        {
            return JointFittingConditioning{ true, 0.0 };
        }
        const auto scale{ std::sqrt(square_sum) };
        for (Eigen::SparseMatrix<double>::InnerIterator entry(normalized_design, column);
            entry;
            ++entry)
        {
            entry.valueRef() /= scale;
        }
    }

    Eigen::SparseMatrix<double> normalized_gram{
        normalized_design.transpose() * normalized_design
    };
    normalized_gram.makeCompressed();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(normalized_gram);
    if (solver.info() != Eigen::Success)
    {
        return JointFittingConditioning{ true, 0.0 };
    }

    const auto diagonal{ solver.vectorD().eval() };
    if (diagonal.size() == 0 || !diagonal.allFinite())
    {
        return JointFittingConditioning{ true, 0.0 };
    }
    if (diagonal.minCoeff() <= 0.0)
    {
        return JointFittingConditioning{ true, 0.0 };
    }
    const auto maximum_pivot{ diagonal.maxCoeff() };
    const auto minimum_pivot{ diagonal.minCoeff() };
    if (!std::isfinite(maximum_pivot) ||
        maximum_pivot <= std::numeric_limits<double>::epsilon())
    {
        return JointFittingConditioning{ true, 0.0 };
    }

    const auto pivot_ratio{ minimum_pivot / maximum_pivot };
    return JointFittingConditioning{
        !std::isfinite(pivot_ratio) || pivot_ratio <= pivot_ratio_threshold,
        std::isfinite(pivot_ratio) ? pivot_ratio : 0.0
    };
}

static double CalculateJointFittingRidgeDiagonal(
    double column_square_sum,
    double ridge_ratio,
    double multiplier)
{
    if (!std::isfinite(multiplier) || multiplier <= 0.0)
    {
        throw std::invalid_argument(
            "Local fitting ridge multiplier must be positive and finite.");
    }
    const auto base_ridge{
        column_square_sum > std::numeric_limits<double>::epsilon() ?
            ridge_ratio * column_square_sum : 1.0
    };
    return multiplier * base_ridge;
}

bool IsJointOffsetSolveHardFailure(JointOffsetSolveStatus status)
{
    switch (status)
    {
    case JointOffsetSolveStatus::Converged:
    case JointOffsetSolveStatus::IrlsObjectiveDeteriorated:
    case JointOffsetSolveStatus::IrlsMaximumIterationsReached:
        return false;
    case JointOffsetSolveStatus::SystemBuildFailed:
    case JointOffsetSolveStatus::EmptySystem:
    case JointOffsetSolveStatus::InitialSolveFailed:
    case JointOffsetSolveStatus::IrlsSolveFailed:
        return true;
    }
    throw std::logic_error("Joint offset solve status is invalid.");
}

const char * GetJointOffsetSolveStatusText(JointOffsetSolveStatus status)
{
    switch (status)
    {
    case JointOffsetSolveStatus::Converged:
        return "converged";
    case JointOffsetSolveStatus::SystemBuildFailed:
        return "system-build-failed";
    case JointOffsetSolveStatus::EmptySystem:
        return "empty-system";
    case JointOffsetSolveStatus::InitialSolveFailed:
        return "initial-solve-failed";
    case JointOffsetSolveStatus::IrlsSolveFailed:
        return "irls-solve-failed";
    case JointOffsetSolveStatus::IrlsObjectiveDeteriorated:
        return "irls-objective-deteriorated";
    case JointOffsetSolveStatus::IrlsMaximumIterationsReached:
        return "irls-maximum-iterations-reached";
    }
    throw std::logic_error("Joint offset solve status is invalid.");
}

void ResetClusterSolverWorkspace(
    const std::vector<ClusterKey> & cluster_key_list,
    ClusterSolverWorkspaceMap & workspace_by_key)
{
    workspace_by_key.clear();
    for (const auto & key : cluster_key_list)
    {
        workspace_by_key.try_emplace(key);
    }
}

bool AreClustersSolverQualified(const ClusterHealthMap & health_by_key)
{
    return std::ranges::all_of(
        health_by_key | std::views::values,
        &ClusterHealth::IsSolverQualified);
}

bool IsLocalRefitStatusSolverQualified(RHBMEstimationStatus status)
{
    switch (status)
    {
    case RHBMEstimationStatus::SUCCESS:
        return true;
    case RHBMEstimationStatus::MAX_ITERATIONS_REACHED:
    case RHBMEstimationStatus::SINGLE_MEMBER:
    case RHBMEstimationStatus::INSUFFICIENT_DATA:
    case RHBMEstimationStatus::NUMERICAL_FALLBACK:
        return false;
    }
    throw std::logic_error("Local Gaussian refit status is invalid.");
}

std::optional<JointOffsetParameterization> BuildJointOffsetParameterization(
    const std::vector<std::size_t> & group_id_by_atom_position,
    const Eigen::VectorXd & atom_offset)
{
    const auto atom_count{ static_cast<std::size_t>(atom_offset.size()) };
    auto group_layout{ BuildJointFittingGroupLayout(group_id_by_atom_position, atom_count) };
    if (!group_layout.has_value()) return std::nullopt;

    JointOffsetParameterization parameterization;
    parameterization.group_position_by_atom = std::move(group_layout->group_position_by_atom);
    parameterization.seed_offset = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(group_layout->group_count));
    std::vector<std::vector<double>> offset_list_by_group(group_layout->group_count);

    for (std::size_t atom_position = 0; atom_position < atom_count; atom_position++)
    {
        const auto offset{ atom_offset(static_cast<Eigen::Index>(atom_position)) };
        if (!std::isfinite(offset)) return std::nullopt;

        const auto atom_group_position{
            parameterization.group_position_by_atom.at(atom_position)
        };
        offset_list_by_group.at(atom_group_position).emplace_back(offset);
    }

    for (std::size_t current_group_position = 0;
        current_group_position < offset_list_by_group.size();
        current_group_position++)
    {
        auto & offset_list{ offset_list_by_group.at(current_group_position) };
        const auto median{ CalculateJointFittingGroupMedian(offset_list) };
        if (!median.has_value()) return std::nullopt;
        parameterization.seed_offset(static_cast<Eigen::Index>(current_group_position)) = *median;
    }
    return parameterization;
}

static algorithm::WeightedRidgeSystem BuildJointOffsetSystem(
    const SecondStageContext & context,
    const std::vector<std::size_t> & active_index_list,
    const SecondStageModelSnapshot & model_snapshot,
    const std::vector<double> & ridge_multiplier_list,
    const JointOffsetParameterization & parameterization,
    bool log_debug_diagnostics)
{
    const auto column_count{ parameterization.seed_offset.size() };
    std::unordered_map<std::size_t, Eigen::Index> active_offset_column_by_atom_index;
    active_offset_column_by_atom_index.reserve(active_index_list.size());
    std::unordered_map<std::size_t, Eigen::Index> active_offset_column_by_group_id;
    Eigen::VectorXd ridge_multiplier_by_group{ Eigen::VectorXd::Ones(column_count) };
    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto atom_index{ active_index_list.at(i) };
        const auto offset_column{ parameterization.OffsetColumn(i) };
        active_offset_column_by_atom_index.emplace(atom_index, offset_column);
        active_offset_column_by_group_id.emplace(context.at(atom_index).group_id, offset_column);
        ridge_multiplier_by_group(offset_column) = std::max(
            ridge_multiplier_by_group(offset_column),
            ridge_multiplier_list.at(atom_index));
    }

    std::vector<Eigen::Triplet<double>> triplet_list;
    std::vector<double> response_list;
    Eigen::VectorXd group_column_square_sum{ Eigen::VectorXd::Zero(column_count) };
    std::map<std::pair<Eigen::Index, Eigen::Index>, double> group_column_cross_sum_map;
    std::vector<std::pair<Eigen::Index, double>> group_row_basis_entries;
    group_row_basis_entries.reserve(static_cast<std::size_t>(column_count));
    for (std::size_t atom_position = 0; atom_position < active_index_list.size(); atom_position++)
    {
        const auto active_index{ active_index_list.at(atom_position) };
        const auto target_offset_column{ parameterization.OffsetColumn(atom_position) };
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
                group_basis(target_offset_column) += target_basis;
            }

            for (const auto & neighbor_atom_sample : atom_context.Neighbors(sample_index))
            {
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
                    const auto selected_group_id{
                        context.unselected_atom_list.at(neighbor_atom_sample.atom_index).selected_group_id
                    };
                    const auto offset_column_iter{
                        selected_group_id.has_value() ?
                            active_offset_column_by_group_id.find(*selected_group_id) :
                            active_offset_column_by_group_id.end()
                    };
                    if (offset_column_iter != active_offset_column_by_group_id.end())
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

        ridge_multiplier_by_group(left_column) = std::max(
            ridge_multiplier_by_group(left_column),
            kJointOffsetConditioningRidgeMultiplier);
        ridge_multiplier_by_group(right_column) = std::max(
            ridge_multiplier_by_group(right_column),
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
        ridge_multiplier_by_group.array() = ridge_multiplier_by_group.array().max(
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
    for (Eigen::Index column_index = 0; column_index < column_count; column_index++)
    {
        system.ridge_diagonal(column_index) =
            CalculateJointFittingRidgeDiagonal(
                group_column_square_sum(column_index),
                kJointOffsetRidgeRatio,
                ridge_multiplier_by_group(column_index));
    }
    return system;
}

static double CalculateWeightedRidgeSurrogateObjective(
    const algorithm::WeightedRidgeSystem & system,
    const Eigen::VectorXd & weight,
    const Eigen::VectorXd & offset)
{
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

static bool IsJointOffsetObjectiveDeteriorated(
    double updated_objective,
    double current_objective)
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

JointOffsetSolveResult EstimateJointOffsets(
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
    std::vector<std::size_t> group_id_by_atom_position;
    group_id_by_atom_position.reserve(active_index_list.size());
    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto atom_index{ active_index_list.at(i) };
        group_id_by_atom_position.emplace_back(context.at(atom_index).group_id);
        previous_offset(static_cast<Eigen::Index>(i)) =
            GetFitModel(model_snapshot.selected, atom_index).GetOffset();
    }
    auto parameterization{
        BuildJointOffsetParameterization(group_id_by_atom_position, previous_offset)
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
    if (system.response.size() == 0)
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

std::optional<JointPolishParameterization> BuildJointPolishParameterization(
    const std::vector<std::size_t> & group_id_by_atom_position,
    const std::vector<GaussianModel3D> & base_model_list)
{
    return BuildActiveSetJointPolishParameterization(
        group_id_by_atom_position,
        base_model_list,
        std::vector<char>(base_model_list.size(), 1),
        std::vector<char>(base_model_list.size(), 1));
}

std::optional<JointPolishParameterization> BuildActiveSetJointPolishParameterization(
    const std::vector<std::size_t> & group_id_by_atom_position,
    const std::vector<GaussianModel3D> & base_model_list,
    const std::vector<char> & shape_active_mask,
    const std::vector<char> & offset_active_mask)
{
    if (shape_active_mask.size() != base_model_list.size() ||
        offset_active_mask.size() != base_model_list.size())
    {
        return std::nullopt;
    }
    auto group_layout{
        BuildJointFittingGroupLayout(group_id_by_atom_position, base_model_list.size())
    };
    if (!group_layout.has_value()) return std::nullopt;

    JointPolishParameterization parameterization;
    parameterization.group_position_by_atom = std::move(group_layout->group_position_by_atom);
    parameterization.shape_position_by_atom.resize(base_model_list.size());
    parameterization.base_shape_coordinate_by_atom.reserve(base_model_list.size());
    parameterization.offset_position_by_group.resize(group_layout->group_count);
    parameterization.base_offset_by_group.resize(group_layout->group_count);
    for (const auto is_shape_active : shape_active_mask)
    {
        if (is_shape_active != 0) parameterization.shape_atom_count++;
    }
    std::vector<int> offset_activity_by_group(group_layout->group_count, -1);
    for (std::size_t atom_position = 0; atom_position < base_model_list.size(); atom_position++)
    {
        const auto group_position{
            parameterization.group_position_by_atom.at(atom_position)
        };
        const auto is_active{ offset_active_mask.at(atom_position) != 0 ? 1 : 0 };
        if (offset_activity_by_group.at(group_position) >= 0 &&
            offset_activity_by_group.at(group_position) != is_active)
        {
            return std::nullopt;
        }
        offset_activity_by_group.at(group_position) = is_active;
    }
    for (std::size_t group_position = 0;
        group_position < offset_activity_by_group.size();
        group_position++)
    {
        if (offset_activity_by_group.at(group_position) != 0)
        {
            parameterization.offset_position_by_group.at(group_position) =
                parameterization.offset_group_count++;
        }
        else
        {
            parameterization.offset_position_by_group.at(group_position) =
                group_layout->group_count;
        }
    }
    parameterization.seed_parameter = Eigen::VectorXd::Zero(
        static_cast<Eigen::Index>(
            parameterization.shape_atom_count * kJointPolishShapeParameterSize +
                parameterization.offset_group_count));
    std::vector<std::vector<double>> offset_list_by_group(group_layout->group_count);
    std::size_t shape_position{ 0 };

    for (std::size_t atom_position = 0; atom_position < base_model_list.size(); atom_position++)
    {
        const auto transformed{
            EncodeTransformedCoordinates(base_model_list.at(atom_position))
        };
        if (!transformed.has_value()) return std::nullopt;
        const auto atom_group_position{
            parameterization.group_position_by_atom.at(atom_position)
        };
        parameterization.base_shape_coordinate_by_atom.emplace_back(
            (*transformed)(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex)),
            (*transformed)(static_cast<Eigen::Index>(kLogWidthChangeIndex)));
        parameterization.shape_position_by_atom.at(atom_position) =
            shape_active_mask.at(atom_position) != 0 ?
                shape_position++ : parameterization.shape_atom_count;
        if (parameterization.HasShapeColumn(atom_position))
        {
            parameterization.seed_parameter(
                parameterization.ShapeColumn(atom_position, 0)) =
                (*transformed)(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex));
            parameterization.seed_parameter(
                parameterization.ShapeColumn(atom_position, 1)) =
                (*transformed)(static_cast<Eigen::Index>(kLogWidthChangeIndex));
        }
        offset_list_by_group.at(atom_group_position).emplace_back(
            base_model_list.at(atom_position).GetOffset());
    }

    for (std::size_t current_group_position = 0;
        current_group_position < offset_list_by_group.size();
        current_group_position++)
    {
        auto & offset_list{ offset_list_by_group.at(current_group_position) };
        const auto median{ CalculateJointFittingGroupMedian(offset_list) };
        if (!median.has_value()) return std::nullopt;
        parameterization.base_offset_by_group.at(current_group_position) = *median;
        if (parameterization.offset_position_by_group.at(current_group_position) <
            parameterization.offset_group_count)
        {
            parameterization.seed_parameter(
                static_cast<Eigen::Index>(
                    parameterization.shape_atom_count * kJointPolishShapeParameterSize +
                    parameterization.offset_position_by_group.at(current_group_position))) = *median;
        }
    }
    return parameterization;
}

std::optional<Eigen::VectorXd> BuildJointPolishDirection(
    const SecondStageContext & context,
    const FitStateView & base_state,
    const ClusterKey & key,
    const std::vector<SampleRef> & sample_ref_list,
    const std::vector<double> & ridge_multiplier_list,
    const JointPolishParameterization & parameterization,
    ReusableWeightedRidgeSolver & reusable_solver)
{
    if (key.empty() || sample_ref_list.empty() ||
        parameterization.group_position_by_atom.size() != key.size() ||
        parameterization.shape_position_by_atom.size() != key.size() ||
        parameterization.base_shape_coordinate_by_atom.size() != key.size() ||
        parameterization.offset_position_by_group.size() !=
            parameterization.base_offset_by_group.size())
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
        if (parameterization.HasShapeColumn(local_position))
        {
            for (std::size_t parameter_index = 0;
                parameter_index < kJointPolishShapeParameterSize;
                parameter_index++)
            {
                ridge_multiplier_by_column(
                    parameterization.ShapeColumn(
                        local_position,
                        parameter_index)) = ridge_multiplier;
            }
        }
        if (parameterization.HasOffsetColumn(local_position))
        {
            const auto offset_column{ parameterization.OffsetColumn(local_position) };
            offset_column_by_group_id.emplace(
                context.at(atom_index).group_id,
                offset_column);
            ridge_multiplier_by_column(offset_column) = std::max(
                ridge_multiplier_by_column(offset_column),
                ridge_multiplier);
        }
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
            if (parameterization.HasShapeColumn(local_position))
            {
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
                        evaluation->shape_jacobian(static_cast<Eigen::Index>(parameter_index))
                    };
                    if (std::abs(derivative) <= std::numeric_limits<double>::epsilon())
                    {
                        continue;
                    }
                    triplet_list.emplace_back(row_index, column_index, derivative);
                }
            }
            if (parameterization.HasOffsetColumn(local_position) &&
                std::abs(evaluation->offset_jacobian) > std::numeric_limits<double>::epsilon())
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

std::optional<FitStateProposal> BuildJointPolishProposal(
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

                FitStateProposal proposal{
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

const char * GetBoundaryJointCorrectionStatusText(
    BoundaryJointCorrectionStatus status)
{
    switch (status)
    {
    case BoundaryJointCorrectionStatus::CandidateReady:
        return "candidate-ready";
    case BoundaryJointCorrectionStatus::InvalidInput:
        return "invalid-input";
    case BoundaryJointCorrectionStatus::InvalidSeed:
        return "invalid-seed";
    case BoundaryJointCorrectionStatus::SystemBuildFailed:
        return "system-build-failed";
    case BoundaryJointCorrectionStatus::TrustRegionUnavailable:
        return "trust-region-unavailable";
    case BoundaryJointCorrectionStatus::NoMaterialChange:
        return "no-material-change";
    }
    return "invalid-input";
}

BoundaryJointCorrectionResult BuildBoundaryJointCorrection(
    const SecondStageContext & context,
    const FitState & previous_state,
    const FitStateView & endpoint_state,
    const std::vector<std::size_t> & shape_active_atom_index_list,
    const std::vector<std::size_t> & offset_active_atom_index_list,
    const std::vector<std::size_t> & offset_closure_atom_index_list,
    const std::vector<SampleRef> & sample_ref_list,
    const std::vector<double> & ridge_multiplier_list,
    const std::vector<BoundaryJointTrustRegion> & trust_region_list,
    ReusableWeightedRidgeSolver & reusable_solver)
{
    BoundaryJointCorrectionResult result;
    if ((shape_active_atom_index_list.empty() && offset_active_atom_index_list.empty()) ||
        offset_closure_atom_index_list.empty() ||
        sample_ref_list.empty() ||
        previous_state.size() != context.size() ||
        endpoint_state.size() != context.size() ||
        ridge_multiplier_list.size() != context.size() ||
        trust_region_list.empty() ||
        !std::ranges::is_sorted(shape_active_atom_index_list) ||
        !std::ranges::is_sorted(offset_active_atom_index_list) ||
        !std::ranges::is_sorted(offset_closure_atom_index_list))
    {
        return result;
    }
    for (const auto atom_index : shape_active_atom_index_list)
    {
        if (atom_index >= context.size() ||
            !std::ranges::binary_search(offset_closure_atom_index_list, atom_index))
        {
            return result;
        }
    }
    for (const auto atom_index : offset_active_atom_index_list)
    {
        if (atom_index >= context.size() ||
            !std::ranges::binary_search(offset_closure_atom_index_list, atom_index))
        {
            return result;
        }
    }

    std::vector<std::size_t> group_id_by_atom_position;
    std::vector<GaussianModel3D> endpoint_model_list;
    std::vector<char> shape_active_mask;
    std::vector<char> offset_active_mask;
    group_id_by_atom_position.reserve(offset_closure_atom_index_list.size());
    endpoint_model_list.reserve(offset_closure_atom_index_list.size());
    shape_active_mask.reserve(offset_closure_atom_index_list.size());
    offset_active_mask.reserve(offset_closure_atom_index_list.size());
    for (const auto atom_index : offset_closure_atom_index_list)
    {
        if (atom_index >= context.size()) return result;
        group_id_by_atom_position.emplace_back(context.at(atom_index).group_id);
        endpoint_model_list.emplace_back(endpoint_state.GetModel(atom_index));
        shape_active_mask.emplace_back(
            std::ranges::binary_search(
                shape_active_atom_index_list,
                atom_index) ? 1 : 0);
        offset_active_mask.emplace_back(
            std::ranges::binary_search(
                offset_active_atom_index_list,
                atom_index) ? 1 : 0);
    }
    const auto parameterization{
        BuildActiveSetJointPolishParameterization(
            group_id_by_atom_position,
            endpoint_model_list,
            shape_active_mask,
            offset_active_mask)
    };
    if (!parameterization.has_value())
    {
        result.status = BoundaryJointCorrectionStatus::InvalidSeed;
        return result;
    }
    result.parameter_count = static_cast<std::size_t>(parameterization->seed_parameter.size());
    const auto seed_model_list{ parameterization->DecodeSeedModels() };
    if (!seed_model_list.has_value())
    {
        result.status = BoundaryJointCorrectionStatus::InvalidSeed;
        return result;
    }
    const auto direction{
        BuildJointPolishDirection(
            context,
            endpoint_state,
            offset_closure_atom_index_list,
            sample_ref_list,
            ridge_multiplier_list,
            *parameterization,
            reusable_solver)
    };
    if (!direction.has_value())
    {
        result.status = BoundaryJointCorrectionStatus::SystemBuildFailed;
        return result;
    }

    double damping{ 1.0 };
    bool found_trust_region_candidate{ false };
    while (damping >= std::numeric_limits<double>::epsilon())
    {
        const auto candidate_model_list{
            parameterization->DecodeModels(*direction, damping)
        };
        if (!candidate_model_list.has_value())
        {
            damping *= 0.5;
            continue;
        }
        double maximum_normalized_trust_step{ 0.0 };
        bool trust_region_accepted{ true };
        for (const auto & trust_region : trust_region_list)
        {
            if (!std::isfinite(trust_region.radius) ||
                trust_region.radius <= 0.0 || trust_region.key.empty())
            {
                return result;
            }
            std::vector<GaussianModel3D> previous_model_list;
            std::vector<GaussianModel3D> candidate_cluster_model_list;
            previous_model_list.reserve(trust_region.key.size());
            candidate_cluster_model_list.reserve(trust_region.key.size());
            for (const auto atom_index : trust_region.key)
            {
                if (atom_index >= context.size()) return result;
                previous_model_list.emplace_back(previous_state.at(atom_index).mdpde.GetModel());
                const auto closure_iter{
                    std::ranges::lower_bound(
                        offset_closure_atom_index_list,
                        atom_index)
                };
                if (closure_iter != offset_closure_atom_index_list.end() && *closure_iter == atom_index)
                {
                    const auto position{ static_cast<std::size_t>(
                        std::distance(
                            offset_closure_atom_index_list.begin(),
                            closure_iter)) };
                    candidate_cluster_model_list.emplace_back(candidate_model_list->at(position));
                }
                else
                {
                    candidate_cluster_model_list.emplace_back(endpoint_state.GetModel(atom_index));
                }
            }
            const auto step_norm{
                CalculateModelTrustRegionStepNorm(
                    previous_model_list,
                    candidate_cluster_model_list)
            };
            if (!step_norm.has_value())
            {
                trust_region_accepted = false;
                break;
            }
            maximum_normalized_trust_step = std::max(
                maximum_normalized_trust_step,
                *step_norm / trust_region.radius);
            if (!IsTrustRegionStepWithinRadius(*step_norm, trust_region.radius))
            {
                trust_region_accepted = false;
                break;
            }
        }
        if (!trust_region_accepted)
        {
            damping *= 0.5;
            continue;
        }
        found_trust_region_candidate = true;
        bool has_material_change{ false };
        for (std::size_t position = 0; position < candidate_model_list->size(); position++)
        {
            const auto change{
                CalculateTransformedChange(
                    candidate_model_list->at(position),
                    endpoint_model_list.at(position))
            };
            if (IsTransformedChangeMaterial(change, kJointPolishTransformedChangeTolerance))
            {
                has_material_change = true;
                break;
            }
        }
        if (!has_material_change)
        {
            result.status = BoundaryJointCorrectionStatus::NoMaterialChange;
            return result;
        }

        FitStatePatch patch;
        patch.atom_index_list = offset_closure_atom_index_list;
        patch.mdpde_list.reserve(offset_closure_atom_index_list.size());
        for (std::size_t position = 0; position < offset_closure_atom_index_list.size(); position++)
        {
            const auto atom_index{
                offset_closure_atom_index_list.at(position)
            };
            patch.mdpde_list.emplace_back(GaussianModel3DWithUncertainty{
                candidate_model_list->at(position),
                endpoint_state.GetMdpde(atom_index).GetStandardDeviationModel()
            });
        }
        result.status = BoundaryJointCorrectionStatus::CandidateReady;
        result.patch = std::move(patch);
        result.damping = damping;
        result.maximum_normalized_trust_step = maximum_normalized_trust_step;
        return result;
    }
    result.status = found_trust_region_candidate ?
        BoundaryJointCorrectionStatus::NoMaterialChange :
        BoundaryJointCorrectionStatus::TrustRegionUnavailable;
    return result;
}

} // namespace rhbm_gem::core::detail
