#pragma once

#include <algorithm>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Sparse>

#include <rhbm_gem/utils/algorithm/WeightedRidgeSolver.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>

#include "core/detail/JointOffsetSolveStatus.hpp"
#include "core/detail/SecondStageContext.hpp"

namespace rhbm_gem::core::detail {

struct ClusterHealth
{
    JointOffsetSolveStatus joint_offset_status{ JointOffsetSolveStatus::SystemBuildFailed };
    bool is_refit_stationarity_eligible{ true };

    bool IsStationarityEligible() const
    {
        return IsJointOffsetSolveStationarityEligible(joint_offset_status) &&
            is_refit_stationarity_eligible;
    }
};

using ClusterHealthMap = std::map<ClusterKey, ClusterHealth>;

inline bool AreClustersStationarityEligible(const ClusterHealthMap & health_by_key)
{
    return std::all_of(
        health_by_key.begin(),
        health_by_key.end(),
        [](const auto & entry)
        {
            return entry.second.IsStationarityEligible();
        });
}

inline bool IsLocalRefitStatusStationarityEligible(RHBMEstimationStatus status)
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

    std::size_t GetSymbolicAnalysisCount() const
    {
        return m_solver.GetSymbolicAnalysisCount();
    }

private:
    static std::vector<std::pair<Eigen::Index, Eigen::Index>> BuildPattern(const Eigen::SparseMatrix<double> & matrix)
    {
        std::vector<std::pair<Eigen::Index, Eigen::Index>> pattern;
        pattern.reserve(static_cast<std::size_t>(matrix.nonZeros()));
        for (Eigen::Index outer = 0; outer < matrix.outerSize(); outer++)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator iter(matrix, outer); iter; ++iter)
            {
                pattern.emplace_back(iter.row(), iter.col());
            }
        }
        return pattern;
    }
};

struct ClusterSolverWorkspace
{
    ReusableWeightedRidgeSolver joint_offset{};
    ReusableWeightedRidgeSolver joint_polish{};
};

using ClusterSolverWorkspaceMap = std::map<ClusterKey, ClusterSolverWorkspace>;

inline void ResetClusterSolverWorkspace(
    const std::vector<ClusterKey> & cluster_key_list,
    ClusterSolverWorkspaceMap & workspace_by_key)
{
    workspace_by_key.clear();
    for (const auto & key : cluster_key_list)
    {
        workspace_by_key.try_emplace(key);
    }
}

} // namespace rhbm_gem::core::detail
