#pragma once

#include "core/detail/FittingModel.hpp"

#include <cstddef>
#include <map>
#include <optional>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <rhbm_gem/utils/algorithm/WeightedRidgeSolver.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>

namespace rhbm_gem::core::detail {

struct JointFittingConditioning
{
    bool guard_required{ false };
    double pivot_ratio{ 0.0 };
};

struct JointFittingGroupLayout
{
    std::vector<std::size_t> group_position_by_atom{};
    std::size_t group_count{ 0 };
};

std::optional<JointFittingGroupLayout> BuildJointFittingGroupLayout(
    const std::vector<std::size_t> & group_id_by_atom_position,
    std::size_t atom_count);

std::optional<double> CalculateJointFittingGroupMedian(
    std::vector<double> & value_list);

JointFittingConditioning EvaluateJointFittingConditioning(
    const Eigen::SparseMatrix<double> & design_matrix,
    double pivot_ratio_threshold);

double CalculateJointFittingRidgeDiagonal(
    double column_square_sum,
    double ridge_ratio,
    double multiplier);


enum class JointOffsetSolveStatus
{
    Converged,
    SystemBuildFailed,
    EmptySystem,
    InitialSolveFailed,
    IrlsSolveFailed,
    IrlsObjectiveDeteriorated,
    IrlsMaximumIterationsReached
};

bool IsJointOffsetSolveHardFailure(JointOffsetSolveStatus status);

const char * GetJointOffsetSolveStatusText(JointOffsetSolveStatus status);


class ReusableWeightedRidgeSolver
{
    algorithm::WeightedRidgeSolver m_solver{};
    Eigen::Index m_row_count{ -1 };
    Eigen::Index m_column_count{ -1 };
    std::vector<std::pair<Eigen::Index, Eigen::Index>> m_pattern{};

public:
    bool Solve(
        const algorithm::WeightedRidgeSystem & system,
        const Eigen::VectorXd & weight,
        Eigen::VectorXd & parameter);

    std::size_t GetSymbolicAnalysisCount() const
    {
        return m_solver.GetSymbolicAnalysisCount();
    }

private:
    static std::vector<std::pair<Eigen::Index, Eigen::Index>> BuildPattern(
        const Eigen::SparseMatrix<double> & matrix);
};

struct ClusterSolverWorkspace
{
    ReusableWeightedRidgeSolver joint_offset{};
    ReusableWeightedRidgeSolver joint_polish{};
};

using ClusterSolverWorkspaceMap = std::map<ClusterKey, ClusterSolverWorkspace>;

void ResetClusterSolverWorkspace(
    const std::vector<ClusterKey> & cluster_key_list,
    ClusterSolverWorkspaceMap & workspace_by_key);


struct ClusterHealth
{
    explicit ClusterHealth(JointOffsetSolveStatus status)
        : joint_offset_status(status)
    {
    }

    JointOffsetSolveStatus joint_offset_status;
    bool is_refit_stationarity_eligible{ true };

    bool IsStationarityEligible() const
    {
        return joint_offset_status == JointOffsetSolveStatus::Converged &&
            is_refit_stationarity_eligible;
    }
};

using ClusterHealthMap = std::map<ClusterKey, ClusterHealth>;

bool AreClustersStationarityEligible(const ClusterHealthMap & health_by_key);

bool IsLocalRefitStatusStationarityEligible(RHBMEstimationStatus status);


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

    Eigen::Index OffsetColumn(std::size_t atom_position) const
    {
        return static_cast<Eigen::Index>(group_position_by_atom.at(atom_position));
    }

    Eigen::VectorXd ExpandOffsets(const Eigen::VectorXd & group_offset) const;
};

std::optional<JointOffsetParameterization> BuildJointOffsetParameterization(
    const std::vector<std::size_t> & group_id_by_atom_position,
    const Eigen::VectorXd & atom_offset);

struct JointOffsetSolveResult
{
    JointOffsetSolveStatus status{ JointOffsetSolveStatus::SystemBuildFailed };
    Eigen::VectorXd offset{};
};

algorithm::WeightedRidgeSystem BuildJointOffsetSystem(
    const SecondStageContext & context,
    const std::vector<std::size_t> & active_index_list,
    const SecondStageModelSnapshot & model_snapshot,
    const std::vector<double> & ridge_multiplier_list,
    const JointOffsetParameterization & parameterization,
    bool log_debug_diagnostics);

double CalculateWeightedRidgeSurrogateObjective(
    const algorithm::WeightedRidgeSystem & system,
    const Eigen::VectorXd & weight,
    const Eigen::VectorXd & offset);

bool IsJointOffsetObjectiveDeteriorated(double updated_objective, double current_objective);

JointOffsetSolveResult EstimateJointOffsets(
    const SecondStageContext & context,
    const std::vector<std::size_t> & active_index_list,
    const SecondStageModelSnapshot & model_snapshot,
    const std::vector<double> & ridge_multiplier_list,
    ReusableWeightedRidgeSolver & reusable_solver,
    bool log_debug_diagnostics);


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
        const Eigen::VectorXd & parameter) const;

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
        double damping) const;

    std::optional<std::vector<GaussianModel3D>> DecodeSeedModels() const;
};

std::optional<JointPolishParameterization> BuildJointPolishParameterization(
    const std::vector<std::size_t> & group_id_by_atom_position,
    const std::vector<GaussianModel3D> & base_model_list);

std::optional<Eigen::VectorXd> BuildJointPolishDirection(
    const SecondStageContext & context,
    const FitStateView & base_state,
    const ClusterKey & key,
    const std::vector<SampleRef> & sample_ref_list,
    const std::vector<double> & ridge_multiplier_list,
    const JointPolishParameterization & parameterization,
    ReusableWeightedRidgeSolver & reusable_solver);

std::optional<FitStateProposal> BuildJointPolishProposal(
    const SecondStageContext & context,
    const FitStateView & base_state,
    const ClusterKey & key,
    const std::vector<SampleRef> & sample_ref_list,
    const std::vector<double> & ridge_multiplier_list,
    ReusableWeightedRidgeSolver & reusable_solver,
    double trust_region_radius);

} // namespace rhbm_gem::core::detail
