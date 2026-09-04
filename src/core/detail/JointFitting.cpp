#include "core/detail/JointFitting.hpp"

#include "core/detail/GaussianModelOperations.hpp"

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
constexpr double kJointFittingRobustLossCutoffMultiplier{ 1.345 };
constexpr double kJointFittingResidualScaleMin{ 1.0e-12 };
constexpr double kJointFittingRidgeRatio{ 1.0e-3 };
constexpr double kJointFittingConditioningRidgeMultiplier{ 10.0 };
constexpr double kJointFittingConditioningPivotRatioThreshold{ 1.0e-8 };
constexpr double kJointOffsetCollinearityOverlapThreshold{ 0.98 };
constexpr double kJointOffsetIrlsScaleFloor{ 1.0e-2 };
constexpr double kJointOffsetIrlsNormalizedChangeTolerance{ 1.0e-6 };
constexpr double kJointOffsetIrlsObjectiveRelativeTolerance{ 1.0e-10 };
constexpr double kJointPolishTransformedChangeTolerance{ 1.0e-4 };

bool HasMaterialJointFittingChange(
    const std::vector<GaussianModel3D> & candidate_model_list,
    const std::vector<GaussianModel3D> & reference_model_list)
{
    for (std::size_t position = 0; position < candidate_model_list.size(); position++)
    {
        const auto change{
            CalculateTransformedChange(
                candidate_model_list.at(position),
                reference_model_list.at(position))
        };
        if (IsTransformedChangeMaterial(change, kJointPolishTransformedChangeTolerance))
        {
            return true;
        }
    }
    return false;
}

} // namespace

std::optional<std::vector<GaussianModel3D>>
JointPolishParameterization::DecodeParameter(const Eigen::VectorXd & parameter) const
{
    if (parameter.size() != m_seed_parameter.size() || !parameter.allFinite()) return std::nullopt;

    std::vector<GaussianModel3D> model_list;
    model_list.reserve(m_base_offset_by_atom.size());
    for (std::size_t atom_position = 0; atom_position < m_base_offset_by_atom.size(); atom_position++)
    {
        auto active_shape_coordinates{
            m_base_shape_coordinate_by_atom.at(atom_position)
        };
        if (HasShapeColumn(atom_position))
        {
            active_shape_coordinates(0) =
                parameter(ShapeColumn(atom_position, 0));
            active_shape_coordinates(1) =
                parameter(ShapeColumn(atom_position, 1));
        }
        const GaussianModel3D::TransformedCoordinates shape_coordinates{
            active_shape_coordinates(0),
            active_shape_coordinates(1),
            0.0
        };
        const auto shape_model{
            GaussianModel3D::FromTransformedCoordinates(shape_coordinates)
        };
        if (!shape_model.has_value()) return std::nullopt;
        const auto offset{
            HasOffsetColumn(atom_position) ?
                parameter(OffsetColumn(atom_position)) :
                m_base_offset_by_atom.at(atom_position)
        };
        const auto model{ shape_model->WithOffset(offset) };
        if (!model.ToTransformedCoordinates().has_value()) return std::nullopt;
        model_list.emplace_back(model);
    }
    return model_list;
}

std::optional<std::vector<GaussianModel3D>>
JointPolishParameterization::DecodeModels(const Eigen::VectorXd & direction, double damping) const
{
    if (direction.size() != m_seed_parameter.size() ||
        !std::isfinite(damping) || damping < 0.0 || damping > 1.0)
    {
        return std::nullopt;
    }
    const Eigen::VectorXd parameter{ m_seed_parameter + damping * direction };
    return DecodeParameter(parameter);
}

std::optional<std::vector<GaussianModel3D>>
JointPolishParameterization::DecodeSeedModels() const
{
    return DecodeParameter(m_seed_parameter);
}

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

static algorithm::WeightedRidgeSystem BuildJointOffsetSystem(
    const SecondStageContext & context,
    const std::vector<std::size_t> & active_index_list,
    const SecondStageModelSnapshot & model_snapshot,
    const std::vector<double> & ridge_multiplier_list,
    const Eigen::VectorXd & previous_offset,
    bool log_debug_diagnostics)
{
    const auto column_count{ previous_offset.size() };
    std::unordered_map<std::size_t, Eigen::Index> active_offset_column_by_atom_index;
    active_offset_column_by_atom_index.reserve(active_index_list.size());
    Eigen::VectorXd ridge_multiplier_by_column{ Eigen::VectorXd::Ones(column_count) };
    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto atom_index{ active_index_list.at(i) };
        const auto offset_column{ static_cast<Eigen::Index>(i) };
        active_offset_column_by_atom_index.emplace(atom_index, offset_column);
        ridge_multiplier_by_column(offset_column) = std::max(
            ridge_multiplier_by_column(offset_column),
            ridge_multiplier_list.at(atom_index));
    }

    std::vector<Eigen::Triplet<double>> triplet_list;
    std::vector<double> response_list;
    Eigen::VectorXd column_square_sum{ Eigen::VectorXd::Zero(column_count) };
    std::map<std::pair<Eigen::Index, Eigen::Index>, double> column_cross_sum_map;
    std::map<Eigen::Index, double> basis_by_column;
    for (std::size_t atom_position = 0; atom_position < active_index_list.size(); atom_position++)
    {
        const auto active_index{ active_index_list.at(atom_position) };
        const auto target_offset_column{ static_cast<Eigen::Index>(atom_position) };
        const auto & atom_context{ context.at(active_index) };
        const auto & target_model{
            GetFitModel(model_snapshot.node, active_index)
        };
        for (std::size_t sample_index = 0; sample_index < atom_context.raw_sampling_entries.size(); sample_index++)
        {
            const auto & sample{ atom_context.raw_sampling_entries.at(sample_index) };
            if (!std::isfinite(sample.response))
            {
                throw std::runtime_error("Joint offset sample response is not finite.");
            }
            const auto target_distance{ sample.point.distance };
            const auto target_signal{ target_model.SignalAtDistance(target_distance) };
            const auto target_basis{ target_model.OffsetBasisAtDistance(target_distance) };
            if (!std::isfinite(target_signal) || !std::isfinite(target_basis))
            {
                throw std::runtime_error("Joint offset target model evaluation is not finite.");
            }
            auto residual{ sample.response - target_signal -
                (model_snapshot.frozen_background ?
                    model_snapshot.frozen_background->response_by_atom.at(active_index).at(sample_index) : 0.0) };
            basis_by_column.clear();
            if (std::abs(target_basis) > std::numeric_limits<double>::epsilon())
            {
                basis_by_column[target_offset_column] += target_basis;
            }

            for (const auto & neighbor_atom_sample : atom_context.Neighbors(sample_index))
            {
                const auto & neighbor_model{
                    GetFitModel(model_snapshot.node, neighbor_atom_sample.atom_index)
                };
                Eigen::Index neighbor_offset_column{ -1 };
                const auto neighbor_offset_column_iter{
                    active_offset_column_by_atom_index.find(
                        neighbor_atom_sample.atom_index)
                };
                if (neighbor_offset_column_iter !=
                    active_offset_column_by_atom_index.end())
                {
                    neighbor_offset_column = neighbor_offset_column_iter->second;
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
                    basis_by_column[neighbor_offset_column] += basis;
                }
            }
            if (!std::isfinite(residual))
            {
                throw std::runtime_error("Joint offset residual is not finite.");
            }
            for (auto basis_iter = basis_by_column.begin();
                basis_iter != basis_by_column.end();)
            {
                if (!std::isfinite(basis_iter->second))
                {
                    throw std::runtime_error("Joint offset atom basis is invalid.");
                }
                if (std::abs(basis_iter->second) <= std::numeric_limits<double>::epsilon())
                {
                    basis_iter = basis_by_column.erase(basis_iter);
                    continue;
                }
                ++basis_iter;
            }
            if (basis_by_column.empty()) continue;

            const auto row_index{ static_cast<Eigen::Index>(response_list.size()) };
            response_list.emplace_back(residual);
            for (const auto & [column_index, basis] : basis_by_column)
            {
                triplet_list.emplace_back(row_index, column_index, basis);
                column_square_sum(column_index) += basis * basis;
            }
            for (auto left_iter = basis_by_column.begin();
                left_iter != basis_by_column.end();
                ++left_iter)
            {
                auto right_iter{ left_iter };
                for (++right_iter; right_iter != basis_by_column.end(); ++right_iter)
                {
                    column_cross_sum_map[
                        { left_iter->first, right_iter->first }] += left_iter->second * right_iter->second;
                }
            }
        }
    }

    for (const auto & [column_pair, cross_sum] : column_cross_sum_map)
    {
        const auto left_column{ column_pair.first };
        const auto right_column{ column_pair.second };
        const auto left_square_sum{ column_square_sum(left_column) };
        const auto right_square_sum{ column_square_sum(right_column) };
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

        ridge_multiplier_by_column(left_column) = std::max(
            ridge_multiplier_by_column(left_column),
            kJointFittingConditioningRidgeMultiplier);
        ridge_multiplier_by_column(right_column) = std::max(
            ridge_multiplier_by_column(right_column),
            kJointFittingConditioningRidgeMultiplier);
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
            kJointFittingConditioningPivotRatioThreshold)
    };
    if (conditioning.guard_required)
    {
        ridge_multiplier_by_column.array() = ridge_multiplier_by_column.array().max(
            kJointFittingConditioningRidgeMultiplier);
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
                << kJointFittingConditioningRidgeMultiplier << ".";
            Logger::Log(LogLevel::Debug, message.str());
        }
    }
    system.previous_parameter = previous_offset;
    system.ridge_diagonal = Eigen::VectorXd::Zero(column_count);
    for (Eigen::Index column_index = 0; column_index < column_count; column_index++)
    {
        system.ridge_diagonal(column_index) =
            CalculateJointFittingRidgeDiagonal(
                column_square_sum(column_index),
                kJointFittingRidgeRatio,
                ridge_multiplier_by_column(column_index));
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
    algorithm::WeightedRidgeSolver & reusable_solver,
    bool log_debug_diagnostics)
{
    Eigen::VectorXd previous_offset{
        Eigen::VectorXd::Zero(static_cast<Eigen::Index>(active_index_list.size()))
    };
    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto atom_index{ active_index_list.at(i) };
        previous_offset(static_cast<Eigen::Index>(i)) =
            GetFitModel(model_snapshot.node, atom_index).GetOffset();
    }
    if (previous_offset.size() == 0 || !previous_offset.allFinite())
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
            previous_offset,
            log_debug_diagnostics);
    }
    catch (const std::runtime_error &)
    {
        return JointOffsetSolveResult{
            JointOffsetSolveStatus::SystemBuildFailed,
            previous_offset
        };
    }
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
                kJointFittingResidualScaleMin)
        };
        for (Eigen::Index i = 0; i < residual.size(); i++)
        {
            weight(i) = algorithm::CalculateCauchyWeight(
                residual(i),
                residual_scale,
                kJointFittingRobustLossCutoffMultiplier);
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
            return JointOffsetSolveResult{
                JointOffsetSolveStatus::IrlsObjectiveDeteriorated,
                offset };
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
            return JointOffsetSolveResult{ JointOffsetSolveStatus::Converged, std::move(offset) };
        }
    }

    return JointOffsetSolveResult{ JointOffsetSolveStatus::IrlsMaximumIterationsReached, std::move(offset) };
}

std::optional<JointPolishParameterization> JointPolishParameterization::Build(
    const std::vector<GaussianModel3D> & base_model_list)
{
    return BuildActiveSet(
        base_model_list,
        std::vector<char>(base_model_list.size(), 1),
        std::vector<char>(base_model_list.size(), 1));
}

std::optional<JointPolishParameterization> JointPolishParameterization::BuildActiveSet(
    const std::vector<GaussianModel3D> & base_model_list,
    const std::vector<char> & shape_active_mask,
    const std::vector<char> & offset_active_mask)
{
    if (shape_active_mask.size() != base_model_list.size() ||
        offset_active_mask.size() != base_model_list.size())
    {
        return std::nullopt;
    }
    if (base_model_list.empty()) return std::nullopt;

    JointPolishParameterization parameterization;
    parameterization.m_shape_column_by_atom.assign(
        base_model_list.size(), kInactiveColumn);
    parameterization.m_offset_column_by_atom.assign(
        base_model_list.size(), kInactiveColumn);
    parameterization.m_base_shape_coordinate_by_atom.reserve(base_model_list.size());
    parameterization.m_base_offset_by_atom.reserve(base_model_list.size());
    Eigen::Index next_column{ 0 };
    for (std::size_t atom_position = 0; atom_position < base_model_list.size(); atom_position++)
    {
        if (shape_active_mask.at(atom_position) == 0) continue;
        parameterization.m_shape_column_by_atom.at(atom_position) = next_column;
        next_column += kShapeParameterSize;
    }
    for (std::size_t atom_position = 0; atom_position < base_model_list.size(); atom_position++)
    {
        if (offset_active_mask.at(atom_position) != 0)
        {
            parameterization.m_offset_column_by_atom.at(atom_position) = next_column++;
        }
    }
    parameterization.m_seed_parameter = Eigen::VectorXd::Zero(next_column);
    for (std::size_t atom_position = 0; atom_position < base_model_list.size(); atom_position++)
    {
        const auto & model{ base_model_list.at(atom_position) };
        const auto transformed{ model.ToTransformedCoordinates() };
        if (!transformed.has_value()) return std::nullopt;
        const Eigen::Vector2d shape{
            (*transformed)(static_cast<Eigen::Index>(GaussianModel3D::LogPeakHeightCoordinateIndex())),
            (*transformed)(static_cast<Eigen::Index>(GaussianModel3D::LogWidthCoordinateIndex()))
        };
        parameterization.m_base_shape_coordinate_by_atom.emplace_back(shape);
        parameterization.m_base_offset_by_atom.emplace_back(model.GetOffset());
        if (parameterization.HasShapeColumn(atom_position))
        {
            for (std::size_t coordinate = 0; coordinate < 2; coordinate++)
            {
                parameterization.m_seed_parameter(
                    parameterization.ShapeColumn(atom_position, coordinate)) =
                    shape(static_cast<Eigen::Index>(coordinate));
            }
        }
        if (parameterization.HasOffsetColumn(atom_position))
        {
            parameterization.m_seed_parameter(parameterization.OffsetColumn(atom_position)) =
                model.GetOffset();
        }
    }
    return parameterization;
}

static std::optional<Eigen::VectorXd> BuildJointPolishDirection(
    const SecondStageContext & context,
    const FitStateView & base_state,
    const ClusterKey & key,
    const std::vector<SampleRef> & sample_ref_list,
    const std::vector<double> & ridge_multiplier_list,
    const JointPolishParameterization & parameterization,
    const std::vector<GaussianModel3D> & seed_model_list,
    algorithm::WeightedRidgeSolver & reusable_solver)
{
    if (key.empty() || sample_ref_list.empty() ||
        parameterization.AtomCount() != key.size() ||
        seed_model_list.size() != key.size())
    {
        return std::nullopt;
    }

    const auto column_count{ parameterization.ParameterCount() };
    std::unordered_map<std::size_t, std::size_t> local_position_by_atom_index;
    local_position_by_atom_index.reserve(key.size());
    Eigen::VectorXd ridge_multiplier_by_column{ Eigen::VectorXd::Ones(column_count) };
    for (std::size_t local_position = 0; local_position < key.size(); local_position++)
    {
        const auto atom_index{ key.at(local_position) };
        local_position_by_atom_index.emplace(atom_index, local_position);
        const auto ridge_multiplier{ ridge_multiplier_list.at(atom_index) };
        if (parameterization.HasShapeColumn(local_position))
        {
            for (std::size_t parameter_index = 0;
                parameter_index < static_cast<std::size_t>(Eigen::Vector2d::SizeAtCompileTime);
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
            ridge_multiplier_by_column(offset_column) = std::max(
                ridge_multiplier_by_column(offset_column),
                ridge_multiplier);
        }
    }
    std::vector<Eigen::Triplet<double>> triplet_list;
    std::vector<double> residual_list;
    residual_list.reserve(sample_ref_list.size());
    for (const auto & sample_ref : sample_ref_list)
    {
        const auto & atom_context{ context.at(sample_ref.atom_index) };
        const auto & sample{
            atom_context.raw_sampling_entries.at(sample_ref.sample_index)
        };
        if (!std::isfinite(sample.response))
        {
            return std::nullopt;
        }

        const auto row_index{ static_cast<Eigen::Index>(residual_list.size()) };
        double predicted_response{ GetFrozenBackgroundResponse(context, sample_ref) };
        if (!std::isfinite(predicted_response)) return std::nullopt;
        const auto append_model = [&](std::size_t atom_index, double distance) -> bool
        {
            const auto local_position_iter{
                local_position_by_atom_index.find(atom_index)
            };
            const auto & model{
                local_position_iter != local_position_by_atom_index.end() ?
                    seed_model_list.at(local_position_iter->second) :
                    base_state.GetModel(atom_index)
            };
            const auto evaluation{ EvaluatePhysicalOffsetResponse(model, distance) };
            if (!evaluation.has_value()) return false;
            predicted_response += evaluation->response;

            if (local_position_iter == local_position_by_atom_index.end()) return true;
            const auto local_position{ local_position_iter->second };
            if (parameterization.HasShapeColumn(local_position))
            {
                for (std::size_t parameter_index = 0;
                    parameter_index < static_cast<std::size_t>(evaluation->shape_jacobian.size());
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
            if (parameterization.HasOffsetColumn(local_position))
            {
                triplet_list.emplace_back(
                    row_index,
                    parameterization.OffsetColumn(local_position),
                    evaluation->offset_jacobian);
            }
            return true;
        };
        if (!append_model(sample_ref.atom_index, sample.point.distance))
        {
            return std::nullopt;
        }
        for (const auto & neighbor_atom_sample : atom_context.Neighbors(sample_ref.sample_index))
        {
            const auto appended{ append_model(
                neighbor_atom_sample.atom_index,
                neighbor_atom_sample.distance) };
            if (!appended) return std::nullopt;
        }

        const auto residual{
            sample.response - predicted_response
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
            kJointFittingConditioningPivotRatioThreshold)
    };
    if (conditioning.guard_required)
    {
        ridge_multiplier_by_column.array() = ridge_multiplier_by_column.array().max(kJointFittingConditioningRidgeMultiplier);
    }
    for (Eigen::Index column_index = 0; column_index < column_count; column_index++)
    {
        system.ridge_diagonal(column_index) =
            CalculateJointFittingRidgeDiagonal(
                column_square_sum(column_index),
                kJointFittingRidgeRatio,
                ridge_multiplier_by_column(column_index));
    }

    const auto residual_scale{
        std::max(
            array_helper::ComputeMedianAbsoluteDeviationScale(residual_list),
            kJointFittingResidualScaleMin)
    };
    if (!std::isfinite(residual_scale)) return std::nullopt;
    Eigen::VectorXd weight{ Eigen::VectorXd::Ones(row_count) };
    for (Eigen::Index row_index = 0; row_index < row_count; row_index++)
    {
        weight(row_index) = algorithm::CalculateCauchyWeight(
            system.response(row_index),
            residual_scale,
            kJointFittingRobustLossCutoffMultiplier);
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
    algorithm::WeightedRidgeSolver & reusable_solver,
    double trust_region_radius)
{
    std::vector<GaussianModel3D> outer_previous_model_list;
    std::vector<GaussianModel3D> base_model_list;
    outer_previous_model_list.reserve(key.size());
    base_model_list.reserve(key.size());
    for (const auto atom_index : key)
    {
        outer_previous_model_list.emplace_back(base_state.GetBaseModel(atom_index));
        base_model_list.emplace_back(base_state.GetModel(atom_index));
    }
    const auto parameterization{
        JointPolishParameterization::Build(
            base_model_list)
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
            *seed_model_list,
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
                if (!HasMaterialJointFittingChange(
                        *candidate_model_list,
                        *seed_model_list))
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
    const FitStateView & endpoint_state,
    const std::vector<std::size_t> & shape_active_atom_index_list,
    const std::vector<std::size_t> & offset_active_atom_index_list,
    const std::vector<SampleRef> & sample_ref_list,
    const std::vector<double> & ridge_multiplier_list,
    const std::vector<BoundaryJointTrustRegion> & trust_region_list,
    algorithm::WeightedRidgeSolver & reusable_solver)
{
    BoundaryJointCorrectionResult result;
    if ((shape_active_atom_index_list.empty() && offset_active_atom_index_list.empty()) ||
        sample_ref_list.empty() ||
        endpoint_state.size() != context.size() ||
        ridge_multiplier_list.size() != context.size() ||
        trust_region_list.empty() ||
        !std::ranges::is_sorted(shape_active_atom_index_list) ||
        !std::ranges::is_sorted(offset_active_atom_index_list))
    {
        return result;
    }
    auto parameter_atom_index_list{ shape_active_atom_index_list };
    parameter_atom_index_list.insert(
        parameter_atom_index_list.end(),
        offset_active_atom_index_list.begin(), offset_active_atom_index_list.end());
    std::ranges::sort(parameter_atom_index_list);
    parameter_atom_index_list.erase(
        std::ranges::unique(parameter_atom_index_list).begin(), parameter_atom_index_list.end());
    for (const auto atom_index : parameter_atom_index_list)
    {
        if (atom_index >= context.size() ||
            !std::ranges::any_of(trust_region_list, [&](const auto & region)
            {
                return std::ranges::find(region.key, atom_index) != region.key.end();
            }))
        {
            return result;
        }
    }

    std::vector<GaussianModel3D> endpoint_model_list;
    std::vector<char> shape_active_mask;
    std::vector<char> offset_active_mask;
    endpoint_model_list.reserve(parameter_atom_index_list.size());
    shape_active_mask.reserve(parameter_atom_index_list.size());
    offset_active_mask.reserve(parameter_atom_index_list.size());
    for (const auto atom_index : parameter_atom_index_list)
    {
        endpoint_model_list.emplace_back(endpoint_state.GetModel(atom_index));
        shape_active_mask.emplace_back(
            std::ranges::binary_search(shape_active_atom_index_list, atom_index) ? 1 : 0);
        offset_active_mask.emplace_back(
            std::ranges::binary_search(offset_active_atom_index_list, atom_index) ? 1 : 0);
    }
    const auto parameterization{
        JointPolishParameterization::BuildActiveSet(
            endpoint_model_list,
            shape_active_mask,
            offset_active_mask)
    };
    if (!parameterization.has_value())
    {
        result.status = BoundaryJointCorrectionStatus::InvalidSeed;
        return result;
    }
    result.parameter_count = static_cast<std::size_t>(parameterization->ParameterCount());
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
            parameter_atom_index_list,
            sample_ref_list,
            ridge_multiplier_list,
            *parameterization,
            *seed_model_list,
            reusable_solver)
    };
    if (!direction.has_value())
    {
        result.status = BoundaryJointCorrectionStatus::SystemBuildFailed;
        return result;
    }

    std::vector<std::vector<GaussianModel3D>> previous_model_list_by_trust_region;
    previous_model_list_by_trust_region.reserve(trust_region_list.size());
    for (const auto & trust_region : trust_region_list)
    {
        if (!std::isfinite(trust_region.radius) ||
            trust_region.radius <= 0.0 || trust_region.key.empty())
        {
            return result;
        }
        auto & previous_model_list{
            previous_model_list_by_trust_region.emplace_back()
        };
        previous_model_list.reserve(trust_region.key.size());
        for (const auto atom_index : trust_region.key)
        {
            if (atom_index >= context.size()) return result;
            previous_model_list.emplace_back(endpoint_state.GetBaseModel(atom_index));
        }
    }

    double damping{ 1.0 };
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
        for (std::size_t trust_region_position = 0;
            trust_region_position < trust_region_list.size();
            trust_region_position++)
        {
            const auto & trust_region{ trust_region_list.at(trust_region_position) };
            const auto & previous_model_list{
                previous_model_list_by_trust_region.at(trust_region_position)
            };
            std::vector<GaussianModel3D> candidate_cluster_model_list;
            candidate_cluster_model_list.reserve(trust_region.key.size());
            for (const auto atom_index : trust_region.key)
            {
                const auto parameter_iter{
                    std::ranges::lower_bound(
                        parameter_atom_index_list,
                        atom_index)
                };
                if (parameter_iter != parameter_atom_index_list.end() && *parameter_iter == atom_index)
                {
                    const auto position{ static_cast<std::size_t>(
                        std::distance(
                            parameter_atom_index_list.begin(),
                            parameter_iter)) };
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
        if (!HasMaterialJointFittingChange(
                *candidate_model_list,
                endpoint_model_list))
        {
            result.status = BoundaryJointCorrectionStatus::NoMaterialChange;
            return result;
        }

        FitStatePatch patch;
        patch.atom_index_list = parameter_atom_index_list;
        patch.mdpde_list.reserve(parameter_atom_index_list.size());
        for (std::size_t position = 0; position < parameter_atom_index_list.size(); position++)
        {
            const auto atom_index{
                parameter_atom_index_list.at(position)
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
    result.status = BoundaryJointCorrectionStatus::TrustRegionUnavailable;
    return result;
}

} // namespace rhbm_gem::core::detail
