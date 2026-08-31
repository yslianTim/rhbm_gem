#pragma once

#include "core/detail/FittingModel.hpp"

#include <cstddef>
#include <map>
#include <optional>
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

JointFittingConditioning EvaluateJointFittingConditioning(
    const Eigen::SparseMatrix<double> & design_matrix,
    double pivot_ratio_threshold);

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

struct ClusterSolverWorkspace
{
    algorithm::WeightedRidgeSolver joint_offset{};
    algorithm::WeightedRidgeSolver joint_polish{};
};

using ClusterSolverWorkspaceMap = std::map<ClusterKey, ClusterSolverWorkspace>;

struct BoundaryJointCorrectionWorkspaceKey
{
    std::vector<std::size_t> shape_active_atom_index_list{};
    std::vector<std::size_t> offset_active_atom_index_list{};
    std::vector<std::size_t> offset_closure_atom_index_list{};
    std::vector<SampleRef> affected_sample_ref_list{};

    friend auto operator<=>(
        const BoundaryJointCorrectionWorkspaceKey &,
        const BoundaryJointCorrectionWorkspaceKey &) = default;
};

using BoundaryJointCorrectionWorkspaceMap =
    std::map<BoundaryJointCorrectionWorkspaceKey, algorithm::WeightedRidgeSolver>;

struct ClusterHealth
{
    explicit ClusterHealth(JointOffsetSolveStatus status)
        : joint_offset_status(status)
    {
    }

    JointOffsetSolveStatus joint_offset_status;
    bool all_local_refits_solver_qualified{ true };
    bool production_convergence_qualified{ true };

    bool IsSolverQualified() const
    {
        return joint_offset_status == JointOffsetSolveStatus::Converged &&
            all_local_refits_solver_qualified;
    }

};

using ClusterHealthMap = std::map<ClusterKey, ClusterHealth>;

bool IsLocalRefitStatusSolverQualified(RHBMEstimationStatus status);

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

JointOffsetSolveResult EstimateJointOffsets(
    const SecondStageContext & context,
    const std::vector<std::size_t> & active_index_list,
    const SecondStageModelSnapshot & model_snapshot,
    const std::vector<double> & ridge_multiplier_list,
    algorithm::WeightedRidgeSolver & reusable_solver,
    bool log_debug_diagnostics);

class JointPolishParameterization
{
    static constexpr std::size_t kShapeParameterSize{ 2 };
    std::vector<std::size_t> m_group_position_by_atom{};
    std::vector<std::size_t> m_shape_position_by_atom{};
    std::vector<std::size_t> m_offset_position_by_group{};
    std::vector<Eigen::Vector2d> m_base_shape_coordinate_by_atom{};
    std::vector<double> m_base_offset_by_group{};
    std::size_t m_shape_atom_count{ 0 };
    std::size_t m_offset_group_count{ 0 };
    Eigen::VectorXd m_seed_parameter{};

public:
    std::size_t AtomCount() const { return m_group_position_by_atom.size(); }
    Eigen::Index ParameterCount() const { return m_seed_parameter.size(); }

    Eigen::Index ShapeColumn(std::size_t atom_position, std::size_t shape_parameter_index) const
    {
        return static_cast<Eigen::Index>(
            m_shape_position_by_atom.at(atom_position) * kShapeParameterSize + shape_parameter_index);
    }

    bool HasShapeColumn(std::size_t atom_position) const
    {
        return m_shape_position_by_atom.at(atom_position) < m_shape_atom_count;
    }

    Eigen::Index OffsetColumn(std::size_t atom_position) const
    {
        return static_cast<Eigen::Index>(
            m_shape_atom_count * kShapeParameterSize +
            m_offset_position_by_group.at(m_group_position_by_atom.at(atom_position)));
    }

    bool HasOffsetColumn(std::size_t atom_position) const
    {
        return m_offset_position_by_group.at(m_group_position_by_atom.at(atom_position)) <
            m_offset_group_count;
    }

    std::optional<std::vector<GaussianModel3D>> DecodeModels(
        const Eigen::VectorXd & direction,
        double damping) const;

    std::optional<std::vector<GaussianModel3D>> DecodeSeedModels() const;

private:
    JointPolishParameterization() = default;
    std::optional<std::vector<GaussianModel3D>> DecodeParameter(const Eigen::VectorXd & parameter) const;
    friend std::optional<JointPolishParameterization>
    BuildActiveSetJointPolishParameterization(
        const std::vector<std::size_t> & group_id_by_atom_position,
        const std::vector<GaussianModel3D> & base_model_list,
        const std::vector<char> & shape_active_mask,
        const std::vector<char> & offset_active_mask);
};

std::optional<JointPolishParameterization> BuildJointPolishParameterization(
    const std::vector<std::size_t> & group_id_by_atom_position,
    const std::vector<GaussianModel3D> & base_model_list);

std::optional<JointPolishParameterization> BuildActiveSetJointPolishParameterization(
    const std::vector<std::size_t> & group_id_by_atom_position,
    const std::vector<GaussianModel3D> & base_model_list,
    const std::vector<char> & shape_active_mask,
    const std::vector<char> & offset_active_mask);

std::optional<FitStateProposal> BuildJointPolishProposal(
    const SecondStageContext & context,
    const FitStateView & base_state,
    const ClusterKey & key,
    const std::vector<SampleRef> & sample_ref_list,
    const std::vector<double> & ridge_multiplier_list,
    algorithm::WeightedRidgeSolver & reusable_solver,
    double trust_region_radius);

enum class BoundaryJointCorrectionStatus
{
    CandidateReady,
    InvalidInput,
    InvalidSeed,
    SystemBuildFailed,
    TrustRegionUnavailable,
    NoMaterialChange
};

const char * GetBoundaryJointCorrectionStatusText(BoundaryJointCorrectionStatus status);

struct BoundaryJointTrustRegion
{
    ClusterKey key{};
    double radius{ 0.0 };
};

struct BoundaryJointCorrectionResult
{
    BoundaryJointCorrectionStatus status{ BoundaryJointCorrectionStatus::InvalidInput };
    std::optional<FitStatePatch> patch{};
    double damping{ 0.0 };
    double maximum_normalized_trust_step{ 0.0 };
    std::size_t parameter_count{ 0 };
};

BoundaryJointCorrectionResult BuildBoundaryJointCorrection(
    const SecondStageContext & context,
    const FitStateView & endpoint_state,
    const std::vector<std::size_t> & shape_active_atom_index_list,
    const std::vector<std::size_t> & offset_active_atom_index_list,
    const std::vector<std::size_t> & offset_closure_atom_index_list,
    const std::vector<SampleRef> & sample_ref_list,
    const std::vector<double> & ridge_multiplier_list,
    const std::vector<BoundaryJointTrustRegion> & trust_region_list,
    algorithm::WeightedRidgeSolver & reusable_solver);

} // namespace rhbm_gem::core::detail
